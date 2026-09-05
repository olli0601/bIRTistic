#!/usr/bin/env python3
"""
Focused-proposal training of the RAGGED deepset amortiser (§14.4.11).

Reframe (no super-population, no prior-predictive): anchor at interim 1
(weekly interim, n=57 ~ 15% of N=503). Its SVI posterior-predictive
future cohort (m=446 participants per draw s) IS a cohort simulated from
theta^(s); its endpoint pps_ratio_x IS rho^(s)=r(theta^(s)) (the exact
deployment target, already computed). We:

  - subsample cohorts of VARYING size n (x) and m (z) from the anchor's
    446-participant pool per draw (decoupled sizes -> breaks the §14.4.10
    n+m=N saturation);
  - label with rho^(s) (contracts, because the size-n cohort only
    estimates theta^(s));
  - feed the RAGGED encoder (segment-mean pool, no padding/mask);
  - train with the interval score (Gneiting-Raftery, scores width).

Then a proposal-probe: does predicted 90% width now contract with n on
held-out proposal draws (slope -> -0.5)? That is the checkpoint before
the deployment test.
"""
import os
import sys
import time
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'python'))
import numpy as np
import pandas as pd
import jax
import jax.numpy as jnp
import optax

warnings.filterwarnings('ignore')
from amortiser_common import save_trained_model
from amortiser_pps_features_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead import (
    Amortiser_PPS_features_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead as Net)
from model_pcm import PartialCreditModel

SB = "/Users/or105/sandbox/bIRTistic"
WK = f"{SB}/py-ukraine-interim-weekly-svi-260811"
dir_out = f"{SB}/py-ukraine-interim-amortise-deepsetXcompAtt-net-260812"
os.makedirs(dir_out, exist_ok=True)
file_prefix = "pcm_1_interim"
ANCHOR = 4                 # weekly interim, n=57 (~11%, nearest existing to 15%)
N_POOL = 446               # anchor future-cohort size (m at interim 4)
N_REF = 503                # deployment total, for aux n/N scaling
D = int(os.environ.get('RAG_D', 1500))       # proposal draws used
NET_STEPS = int(os.environ.get('RAG_STEPS', 6000))
B = int(os.environ.get('RAG_B', 32))         # batch elements per step
Q = 4                                        # queries per element
TAUS = np.array([0.05, 0.25, 0.5, 0.75, 0.95], np.float32)
CLIP = 20.0
SIZE_GRID = np.unique(np.round(np.geomspace(2, N_POOL, 10)).astype(int))  # n,m grid
CACHE = "/private/tmp/claude-501/-Users-or105-git-bIRTistic/49372e22-d11b-4882-a442-d1c60bcbdfb0/scratchpad/ragged_pool"
os.makedirs(CACHE, exist_ok=True)

# ---- item metadata (shared ordering) ----
dit = pd.read_csv(f"{WK}/{file_prefix}_1_data_dit.csv")
wa = pd.read_pickle(f"{WK}/{file_prefix}_i{ANCHOR}_regression_training.pkl")
items = (wa[['item_label', 'item_type', 'item_high_label']].drop_duplicates()
         .sort_values(['item_type', 'item_label']).reset_index(drop=True))
J = len(items); labels = items.item_label.tolist()
itype = items.item_type.map({'out-of-7': 0., 'categorical': 1.}).to_numpy(np.float32)
ihigh = items.item_high_label.map({'higher_is_better': 1., 'lower_is_better': 0.}).to_numpy(np.float32)
items = items.merge(dit[['item_label', 'cat_length']].drop_duplicates(), on='item_label', how='left')
kmax = items.cat_length.to_numpy(np.float32)
META = np.stack([itype, ihigh], -1).astype(np.float32)      # (J, 2)


def build_pool():
    """POOL (D, N_POOL, J, 2) rescaled + rho_target (D, J)."""
    xi = pd.read_csv(f"{WK}/{file_prefix}_{ANCHOR}_data_dp1.csv")
    model = PartialCreditModel(dit=dit, dcati=xi, x_formula="~ time - 1", seed=123)
    zi = model.get_interim_z_from_ypredi(f"{WK}/{file_prefix}_{ANCHOR}_draws.zarr",
                                         N_POOL, pps_z_total=D, seed=123, keep_order=True)
    cols = [f'ypred_{s}' for s in range(D)]
    pids = np.sort(zi.pid.unique())                 # canonical future-cohort pids
    NP = len(pids)
    POOL = np.zeros((D, NP, J, 2), np.float32)       # 0-fill per-item missing (as deploy _pivot)
    for j, l in enumerate(labels):
        km = kmax[j] - 1.0
        for t in (0, 1):
            blk = zi[(zi.item_label == l) & (zi.time == t)].set_index('pid')[cols].reindex(pids)
            v = blk.to_numpy(np.float32).T          # (D, NP), NaN for missing
            POOL[:, :, j, t] = np.nan_to_num((v - 1.0) / km)
    piv = wa.pivot_table(index='draw', columns='item_label', values='pps_ratio_x').reindex(columns=labels)
    rt = np.full((D, J), np.nan)
    idx = piv.index[piv.index < D].to_numpy()
    rt[idx] = np.clip(piv.loc[idx].to_numpy(np.float64), -CLIP, CLIP)
    return POOL, rt.astype(np.float32)


cf = f"{CACHE}/pool_D{D}.npz"
if os.path.exists(cf):
    z = np.load(cf); POOL, RT = z['pool'], z['rt']
else:
    print(f"building proposal pool (D={D}) ...")
    POOL, RT = build_pool(); np.savez(cf, pool=POOL, rt=RT)
ok_draw = np.isfinite(RT).all(1)
POOL, RT = POOL[ok_draw], RT[ok_draw]
Dn = POOL.shape[0]
N_POOL = POOL.shape[1]                              # actual future-cohort size
SIZE_GRID = np.unique(np.round(np.geomspace(2, N_POOL, 10)).astype(int))
sig = np.nanstd(RT, 0).clip(1e-3).astype(np.float32)        # per-item sd
np.save(f"{dir_out}/{file_prefix}_item_std.npy", sig)
print(f"pool {POOL.shape}, targets {RT.shape}, sigma {np.round(sig,3)[:5]}...")
D_train = Dn - 200
print(f"train draws {D_train}, held-out {Dn - D_train}")


def sample_batch(rng):
    """Homogeneous (n,m) per batch from the grid; ragged concat, no pad/mask."""
    n = int(rng.choice(SIZE_GRID)); m = int(rng.choice(SIZE_GRID))
    s = rng.integers(0, D_train, size=B)
    x_list, z_list = [], []
    for b in range(B):
        xi = rng.choice(N_POOL, size=n, replace=False)
        zi = rng.choice(N_POOL, size=m, replace=False)
        x_list.append(POOL[s[b], xi]); z_list.append(POOL[s[b], zi])
    x_flat = np.concatenate(x_list, 0)                       # (B*n, J, 2)
    z_flat = np.concatenate(z_list, 0)
    x_seg = np.repeat(np.arange(B), n).astype(np.int32)
    z_seg = np.repeat(np.arange(B), m).astype(np.int32)
    qidx = rng.integers(0, J, size=(B, Q)).astype(np.int32)
    aux = np.tile(np.array([1/np.sqrt(n), 1/np.sqrt(m), n/N_REF, m/N_REF], np.float32), (B, 1))
    y = RT[s[:, None], qidx] / sig[qidx]                     # (B, Q) standardised
    return (dict(x_flat=jnp.asarray(x_flat), x_seg=jnp.asarray(x_seg),
                 z_flat=jnp.asarray(z_flat), z_seg=jnp.asarray(z_seg),
                 item_metadata=jnp.asarray(np.tile(META[None], (B, 1, 1))),
                 query_idx=jnp.asarray(qidx), aux=jnp.asarray(aux)),
            jnp.asarray(y.astype(np.float32)))


net = Net(num_quantiles=5)
rng = np.random.default_rng(0)
b0, y0 = sample_batch(rng)
params = net.init(jax.random.PRNGKey(0), b0)
n_par = sum(int(np.prod(v.shape)) for v in jax.tree_util.tree_leaves(params))
print(f"ragged net params: {n_par}")


def interval_score(preds, y):
    pr = jnp.maximum.accumulate(preds, -1)                   # monotone quantiles
    l05, l25, med, u75, u95 = pr[..., 0], pr[..., 1], pr[..., 2], pr[..., 3], pr[..., 4]
    is10 = (u95 - l05) + 20.0 * jnp.maximum(l05 - y, 0) + 20.0 * jnp.maximum(y - u95, 0)
    is50 = (u75 - l25) + 4.0 * jnp.maximum(l25 - y, 0) + 4.0 * jnp.maximum(y - u75, 0)
    return jnp.mean(is10 + is50 + jnp.abs(y - med))


sched = optax.warmup_cosine_decay_schedule(1e-4, 1e-3, NET_STEPS // 10, NET_STEPS, 1e-5)
opt = optax.adam(sched); ost = opt.init(params)


@jax.jit
def step(params, ost, batch, y):
    loss, g = jax.value_and_grad(lambda p: interval_score(net.apply(p, batch), y))(params)
    up, ost = opt.update(g, ost, params)
    return optax.apply_updates(params, up), ost, loss


print(f"training ragged net ({NET_STEPS} steps, B={B}, Q={Q}) ...")
t0 = time.time()
for it in range(NET_STEPS):
    batch, y = sample_batch(rng)
    params, ost, loss = step(params, ost, batch, y)
    if (it + 1) % max(1, NET_STEPS // 20) == 0:
        print(f"  step {it+1}/{NET_STEPS}  IS={float(loss):.4f}")
mins = (time.time() - t0) / 60.0
fit = {'params': params, 'apply_fn': net.apply,
       'net_class_name': type(net).__name__,
       'net_class_module': type(net).__module__,
       'net_kwargs': {'num_quantiles': 5}, 'training_mins': mins,
       'pps_ProbH1_lwr_quantiles_mesh': tuple(TAUS.tolist())}
save_trained_model(fit, f"{dir_out}/{file_prefix}_amortised_pps_net.pkl")
print(f"saved ragged net ({mins:.1f} min) -> {dir_out}")


# ---- proposal-probe: width vs n on held-out draws (m fixed) ----
@jax.jit
def fwd(params, batch):
    return net.apply(params, batch)


ho = np.arange(D_train, Dn)
m_fix = 100
grid = np.unique(np.round(np.geomspace(6, N_POOL, 8)).astype(int))
rows = []
rng2 = np.random.default_rng(1)
for n in grid:
    ws = []
    for s in ho:
        xi = rng2.choice(N_POOL, size=int(n), replace=False)
        zi = rng2.choice(N_POOL, size=m_fix, replace=False)
        batch = dict(x_flat=jnp.asarray(POOL[s, xi]), x_seg=jnp.zeros(int(n), jnp.int32),
                     z_flat=jnp.asarray(POOL[s, zi]), z_seg=jnp.zeros(m_fix, jnp.int32),
                     item_metadata=jnp.asarray(META[None]),
                     query_idx=jnp.arange(J)[None], aux=jnp.asarray(
                         np.array([[1/np.sqrt(n), 1/np.sqrt(m_fix), n/N_REF, m_fix/N_REF]], np.float32)))
        pr = np.maximum.accumulate(np.asarray(fwd(params, batch))[0], 1)  # (J,5) in std units
        pr = pr * sig[:, None]
        ws.append(pr[:, 4] - pr[:, 0])                       # per-item 90% width
    ws = np.array(ws)                                         # (n_ho, J)
    rows.append(dict(n=int(n), width=float(np.median(ws))))
pw = pd.DataFrame(rows)
slope = float(np.polyfit(np.log(pw.n), np.log(pw.width), 1)[0])
pw.to_csv(f"{dir_out}/{file_prefix}_proposal_probe.csv", index=False)
print("\nPROPOSAL-PROBE width vs n (held-out proposal draws, m=100):")
print(pw.to_string(index=False))
print(f"\nwidth log-log slope: {slope:.3f}  (flat baseline -0.02; target ~ SVI -0.59)")
print(f"width n={grid[0]}: {pw.width.iloc[0]:.3f}  n={grid[-1]}: {pw.width.iloc[-1]:.3f}")
