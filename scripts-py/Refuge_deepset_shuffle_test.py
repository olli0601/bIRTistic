"""Shuffle test (§14.5.1) on REFUGE-ED (MSPSS), the highest-item-covariance application.

Same two marginal-preserving permutations as the Ukraine test
(scripts-py/Ukraine_deepset_shuffle_test.py):
  cross-item : permute participants per item, keep each (baseline,endline) pair -> preserves every item's
               marginal, destroys the cross-item joint. A plain deepsetXcompAtt is cross-item blind BY
               CONSTRUCTION, so its output must be invariant on ANY data -- the point here is to confirm
               the blindness holds on a scale (MSPSS) with ~2x Ukraine's item correlation (mean|r| 0.52
               vs 0.26), i.e. the deepset discards substantially MORE genuine covariance here.
  within-item: permute baseline and endline independently per item -> breaks the per-item joint (control
               that the net is not dead).
Uses the item-general J64 deepset (powerlaw head, scalar sigma) -- the cross-application asset deployed
on REFUGE -- applied to the REFUGE MSPSS cohorts. Run:
  .pixi/envs/default/bin/python scripts-py/Refuge_deepset_shuffle_test.py
"""
import os, sys, glob, numpy as np, pandas as pd, jax, jax.numpy as jnp
sys.path.insert(0, 'python')
from amortiser_common import load_fitted_model
from model_pcm import PartialCreditModel
SB = "/Users/or105/sandbox/bIRTistic"
WK = f"{SB}/py-refugee_interim_260831"                                   # REFUGE data (dp1 / draws / dit)
BASE = f"{SB}/py-ukraine-interim-amortise-deepsetXcompAtt-itemamortise-J64-260831"  # item-general net
PFX = "pcm_1_interim"; N_REF = 324; S = 100; KM = 10
dit = pd.read_csv(f"{WK}/{PFX}_1_data_dit.csv")
_wa = sorted(glob.glob(f"{WK}/{PFX}_i*_regression_training.pkl"),
             key=lambda p: int(p.split('_i')[-1].split('_')[0]))
items = (pd.read_pickle(_wa[0])[['item_label', 'item_type', 'item_high_label']].drop_duplicates()
         .sort_values(['item_type', 'item_label']).reset_index(drop=True))
J = len(items); labels = items.item_label.tolist()
itype = items.item_type.map({'out-of-7': 0., 'categorical': 1.}).to_numpy(np.float32)
ihigh = items.item_high_label.map({'higher_is_better': 1., 'lower_is_better': 0.}).to_numpy(np.float32)
items = items.merge(dit[['item_label', 'cat_length']].drop_duplicates(), on='item_label', how='left')
kmax = items.cat_length.to_numpy(np.float32)
META = np.stack([itype, ihigh, kmax / KM], -1).astype(np.float32)        # J64 item-general META (M=3)
fit = load_fitted_model(f"{BASE}/{PFX}_amortised_pps_net.pkl")
import importlib as _il
Net = getattr(_il.import_module(fit['net_class_module']), fit['net_class_name'])
net = Net(**dict(fit['net_kwargs']))
_sig0 = float(np.reshape(np.load(f"{BASE}/{PFX}_item_std.npy"), -1)[0])   # scalar global sigma
sig = np.full(J, _sig0, np.float32)


@jax.jit
def fwd(b):
    return net.apply(fit['params'], b)


def pivot(df, col, pids):
    piv = df.pivot_table(index='pid', columns=['item_label', 'group'], values=col).reindex(pids)
    out = np.zeros((len(pids), J, 2), np.float32)
    for j, l in enumerate(labels):
        for t in (0, 1):
            if (l, t) in piv.columns:
                out[:, j, t] = piv[(l, t)].to_numpy(np.float32)
    return np.nan_to_num((out - 1.) / (kmax[None, :, None] - 1.))


def sh_crossitem(a, rng):            # permute participants per item, keep (base,end) pair
    o = np.empty_like(a)
    for j in range(a.shape[1]):
        o[:, j, :] = a[rng.permutation(a.shape[0]), j, :]
    return o


def sh_withinitem(a, rng):           # permute base and end independently per item -> breaks within-item joint
    o = np.empty_like(a)
    for j in range(a.shape[1]):
        for t in (0, 1):
            o[:, j, t] = a[rng.permutation(a.shape[0]), j, t]
    return o


def run(xr, zl, metab, qidb, aux):
    qs = np.empty((S, J, 5))
    for s in range(S):
        b = dict(x_flat=jnp.asarray(xr), x_seg=jnp.zeros(len(xr), jnp.int32),
                 z_flat=jnp.asarray(zl[s]), z_seg=jnp.zeros(len(zl[s]), jnp.int32),
                 item_metadata=metab, query_idx=qidb, aux=aux)
        qs[s] = np.maximum.accumulate(np.asarray(fwd(b))[0], 1) * sig[:, None]
    return qs


def pps(q, eta0=0.5, etaH=0.89):
    taus = np.array([.05, .25, .5, .75, .95])
    PH = np.array([[1 - np.interp(eta0, q[s, j], taus) for j in range(J)] for s in range(S)])
    return (PH > etaH).mean(0)


metab = jnp.asarray(META[None]); qidb = jnp.arange(J)[None]
print(f"REFUGE-ED (MSPSS) shuffle test: J={J} items, item-general J64 net, S={S} draws")
for k in (2, 5, 8):
    xi = pd.read_csv(f"{WK}/{PFX}_{k}_data_dp1.csv"); xpids = np.sort(xi.pid.unique()); n = len(xpids)
    x_raw = pivot(xi, 'y_stan', xpids)
    model = PartialCreditModel(dit=dit, dcati=xi, x_formula="~ group - 1", seed=123)
    zi = model.get_interim_z_from_ypredi(f"{WK}/{PFX}_{k}_draws.zarr", max(N_REF - n, 1),
                                         pps_z_total=S, seed=123, keep_order=True)
    zpids = np.sort(zi.pid.unique()); mm = len(zpids)
    aux = jnp.asarray(np.array([[1 / np.sqrt(n), 1 / np.sqrt(mm), n / N_REF, mm / N_REF]], np.float32))
    zraws = [pivot(zi.rename(columns={f'ypred_{s}': '_y'}), '_y', zpids) for s in range(S)]
    rng = np.random.default_rng(0)
    base = run(x_raw, zraws, metab, qidb, aux)
    qA = run(sh_crossitem(x_raw, rng), [sh_crossitem(z, rng) for z in zraws], metab, qidb, aux)
    qB = run(sh_withinitem(x_raw, rng), [sh_withinitem(z, rng) for z in zraws], metab, qidb, aux)
    p0, pA, pB = pps(base), pps(qA), pps(qB)
    print(f"k={k:2d} n={n:3d} m={mm:3d}: cross-item  max|Δρ̂|={np.abs(qA - base).max():.2e} PPS max|Δ|={np.abs(pA - p0).max():.2e} | "
          f"within-item max|Δρ̂|={np.abs(qB - base).max():.3f} PPS max|Δ|={np.abs(pB - p0).max():.3f}")
