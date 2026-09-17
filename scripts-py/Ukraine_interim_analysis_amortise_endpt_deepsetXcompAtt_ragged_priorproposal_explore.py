#!/usr/bin/env python3
"""
EXPLORATORY (§14.4.12 corrected recipe). New script — does not touch the
focused-proposal build. Prior proposal + theta-level target + decoupled
(n,m) + ragged encoder + interval score.

  - proposal: item params (beta, tau, lambda) ~ PRIOR (wide theta range);
  - target rho^(s): theta-LEVEL endpoint = mean of the per-ability response
    functional over a large REFERENCE ability set (M_ref), i.e. the same
    quantity get_endpoints_per_draw computes over a design population, NOT
    the finite input cohort (that was the §14.4.10 saturation);
  - input: FRESH finite cohorts x_n, z_m simulated from the SAME params,
    decoupled n,m -> the size-n cohort only ESTIMATES rho^(s) -> contraction.

Prints an ORACLE slope first (plug-in endpoint error vs n on the theta-
level target): decisive + cheap. If it contracts (~ -0.5) the recipe has
the signal; then trains the ragged net + interval score + proposal-probe.
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

SB = "/Users/or105/sandbox/bIRTistic"
WK = f"{SB}/py-ukraine-interim-weekly-svi-260811"
OUT_TAG = os.environ.get('PP_TAG', '260812')
dir_out = f"{SB}/py-ukraine-interim-amortise-deepsetXcompAtt-{OUT_TAG}"
os.makedirs(dir_out, exist_ok=True)
file_prefix = "pcm_1_interim"
K_O7, K_CAT, C_CAT = 8, 4, 2
THRESH = 3.5
M_REF = 2000
N_REF = 503
NET_STEPS = int(os.environ.get('PP_STEPS', 6000))
B = int(os.environ.get('PP_B', 24)); Q = 4
CLIP = 20.0
# §14.4.14 capacity-fix knobs
EMBED_DIM = int(os.environ.get('PP_EMBED', 32))
HEAD_HIDDEN = tuple(int(x) for x in os.environ.get('PP_HEAD', '64,64').split(','))
HEAD_MODE = os.environ.get('PP_HEADMODE', 'plain')       # plain|factored|semiparam
PREC_POOL = os.environ.get('PP_PRECPOOL', '0') == '1'
SCHEDULE = os.environ.get('PP_SCHED', 'strat')           # strat|dirichlet
RUN_ORACLE = os.environ.get('PP_ORACLE', '0') == '1'
PREG = float(os.environ.get('PP_PREG', '0'))             # regularise powerlaw exponent p -> 1/2
DR = float(os.environ.get('PP_DR', '0'))                 # domain-randomise z-cohort (perturb time effect)
MEMPTY = float(os.environ.get('PP_MEMPTY', '0'))         # frac of batches with EMPTY future cohort (m=0)
#   §14.4.25 Option A (law-of-total-prob): train some batches with no z at all so the net
#   learns p(rho | x, empty) = p(rho|x) directly. At deploy the effect-size marginal is one
#   empty-z forward (no z^(s) mixing) -> clean n-contraction, no finite-m MC inflation.
ZCON = os.environ.get('PP_ZCON', '0') == '1'             # z-contrast token (pool_z-pool_x, pool_z*pool_x)
RAWPOOL = os.environ.get('PP_RAWPOOL', '0') == '1'       # explicit raw plug-in ratio channel (we/wb)
PERPART = os.environ.get('PP_PERPART', '0') == '1'       # per-participant rho-mapping channels
LINTAU = os.environ.get('PP_LINTAU', '0') == '1'         # linear (identity-capable) q_tau
HEADRAW = os.environ.get('PP_HEADRAW', '0') == '1'       # inject raw (wb,we,ratio) at the head
# arbitrary decision-threshold experiments:
#   PP_NQ    : number of quantile levels (Option 3: dense quantiles; default 5 = baseline)
#   PP_ETA0  : 1 => Option 4, amortise eta_0 with a success-prob head P(rho>eta0|x,z),
#              eta_0 (standardised) appended to aux, binary target, BCE loss.
NQ = int(os.environ.get('PP_NQ', '5'))
ETA0AM = os.environ.get('PP_ETA0', '0') == '1'
# §14.4.9 transfer fine-tune knobs: warm-start from a checkpoint, M=3 metadata (add K/K_max),
# single global sigma (so the deploy treats the net as item-amortised).
INIT = os.environ.get('PP_INIT', '')
META3 = os.environ.get('PP_META3', '0') == '1'
GLOBSIG = os.environ.get('PP_GLOBSIG', '0') == '1'
_num_q = 1 if ETA0AM else NQ
_head_mode = 'successprob' if ETA0AM else HEAD_MODE
NET_KW = dict(num_quantiles=_num_q, embed_dim=EMBED_DIM,
              q_tau_hidden=(EMBED_DIM, EMBED_DIM), q_tok_hidden=(EMBED_DIM, EMBED_DIM),
              q_query_hidden=(EMBED_DIM, EMBED_DIM), hidden_dims=HEAD_HIDDEN,
              head_mode=_head_mode, precision_pool=PREC_POOL, z_contrast=ZCON,
              raw_pool=RAWPOOL, perpart_map=PERPART, linear_tau=LINTAU, head_raw=HEADRAW)
print(f"CONFIG tag={OUT_TAG} sched={SCHEDULE} E={EMBED_DIM} head={HEAD_HIDDEN} "
      f"mode={HEAD_MODE} prec_pool={PREC_POOL} zcon={ZCON} dr={DR}")
TAUS = (np.linspace(0.05, 0.95, NQ).astype(np.float32) if (NQ != 5 and not ETA0AM)
        else np.array([0.05, 0.25, 0.5, 0.75, 0.95], np.float32))

# item metadata / directions (shared ordering)
dit = pd.read_csv(f"{WK}/{file_prefix}_1_data_dit.csv")
wa = pd.read_pickle(f"{WK}/{file_prefix}_i4_regression_training.pkl")
items = (wa[['item_label', 'item_type', 'item_high_label']].drop_duplicates()
         .sort_values(['item_type', 'item_label']).reset_index(drop=True))
J = len(items); labels = items.item_label.tolist()
itype = items.item_type.map({'out-of-7': 0., 'categorical': 1.}).to_numpy(np.float32)
ihigh = items.item_high_label.map({'higher_is_better': 1., 'lower_is_better': 0.}).to_numpy(np.float32)
if META3:                                        # M=3: add K/K_max (K=8 out-of-7, 4 categorical)
    _Karr = np.where(itype < .5, K_O7, K_CAT).astype(np.float32)
    META = np.stack([itype, ihigh, _Karr / 10.0], -1).astype(np.float32)
else:
    META = np.stack([itype, ihigh], -1).astype(np.float32)
o7 = np.flatnonzero(itype == 0.); ca = np.flatnonzero(itype == 1.)
GROUPS = [(o7, K_O7), (ca, K_CAT)]


def draw_params(rng):
    """One PCM parameter draw from the prior. Returns dict per group."""
    beta = rng.standard_normal(2)
    P = {}
    for js, K in GROUPS:
        tau = THRESH * rng.standard_normal((js.size, 2, K - 1))
        lam = np.abs(rng.standard_t(3.0, size=(js.size, 2))); lam[0, 0] = 1.0
        P[K] = (tau, lam)
    return beta, P


def _probs(theta, beta, tau, lam, K):
    """theta (Nab,), tau (Jt,2,K-1), lam (Jt,2) -> p (Nab, Jt, 2, K)."""
    latent = theta[:, None, None] + beta[None, None, :]           # (Nab, 1, 2)
    inc = lam[None, :, :, None] * (latent[:, :, :, None] - tau[None, :, :, :])  # (Nab, Jt, 2, K-1)
    eta = np.concatenate([np.zeros((*inc.shape[:3], 1)), np.cumsum(inc, -1)], -1)
    eta -= eta.max(-1, keepdims=True)
    p = np.exp(eta); p /= p.sum(-1, keepdims=True)
    return p                                                       # (Nab, Jt, 2, K)


def _mu(p, K):
    if K == K_O7:
        ks = np.arange(1, K + 1)
        return (p * ks[None, None, None, :]).sum(-1)               # (Nab, Jt, 2)
    return p[:, :, :, C_CAT - 1:].sum(-1)                          # caseness


def theta_level_rho(beta, P):
    """rho^(s)_j (J,) = theta-level endpoint over M_REF reference abilities."""
    th = np.random.default_rng(0).standard_normal(M_REF)  # fixed ref set per call is fine
    rho = np.zeros(J)
    for js, K in GROUPS:
        tau, lam = P[K]
        p = _probs(th, beta, tau, lam, K)
        w = _mu(p, K).mean(0)                                      # (Jt, 2)
        wb = np.maximum(w[:, 0], 1e-3); we = w[:, 1]
        hi = ihigh[js] > .5
        rho[js] = np.clip(np.where(hi, we/wb - 1, 1 - we/wb), -CLIP, CLIP)
    return rho


def sim_cohort(beta, P, nab, rng):
    """Simulate nab participants' rescaled responses -> (nab, J, 2)."""
    Y = np.zeros((nab, J, 2), np.float32)
    th = rng.standard_normal(nab)
    for js, K in GROUPS:
        tau, lam = P[K]
        p = _probs(th, beta, tau, lam, K)                         # (nab, Jt, 2, K)
        cdf = np.cumsum(p, -1); u = rng.random((nab, js.size, 2))
        k = (cdf < u[..., None]).sum(-1)                          # 0..K-1
        Y[:, js, :] = k.astype(np.float32) / (K - 1)
    return Y


# ---------- pilot: per-item sd for target standardisation (+ optional oracle) ----------
rng = np.random.default_rng(1)
params_list = [draw_params(rng) for _ in range(400)]
rhos = np.array([theta_level_rho(b, P) for b, P in params_list])   # (NP, J)
sig = np.maximum(np.std(rhos, 0), 1e-3).astype(np.float32)
if GLOBSIG:                                      # single global sigma (item-amortised-style)
    _sg = float(np.maximum(np.std(rhos), 1e-3)); sig = np.full(J, _sg, np.float32)
    np.save(f"{dir_out}/{file_prefix}_item_std.npy", np.array([_sg], np.float32))
else:
    np.save(f"{dir_out}/{file_prefix}_item_std.npy", sig)
os_slope = float('nan')
if RUN_ORACLE:
    print("ORACLE: plug-in endpoint error vs n (prior params, theta-level target)")
    grid = np.unique(np.round(np.geomspace(6, N_REF, 8)).astype(int))
    orows = []
    for n in grid:
        errs = []
        for (b, P), rho in zip(params_list, rhos):
            Y = sim_cohort(b, P, int(n), rng)
            wb = np.maximum(Y[:, :, 0].mean(0), 1e-3); we = Y[:, :, 1].mean(0)
            rhat = np.where(ihigh > .5, we/wb - 1, 1 - we/wb)
            errs.append(np.abs(np.clip(rhat, -CLIP, CLIP) - rho))
        orows.append((int(n), float(np.median(errs)))); print(f"  n={int(n):4d}  err={orows[-1][1]:.3f}")
    oa = np.array(orows); os_slope = float(np.polyfit(np.log(oa[:, 0]), np.log(oa[:, 1]), 1)[0])
    print(f"ORACLE error slope: {os_slope:.3f}")

# ---------- ragged training (prior proposal, theta-level target) ----------
# §14.4.19 Nband-noDR recipe knobs (defaults reproduce the §14.4.13 recipe):
#   PP_NHI : upper size for the per-sample size schedule (N-band cap)
#   PP_NX  : mean observed-cohort size  (T_X = B * PP_NX)
#   PP_MZ  : mean future-cohort size    (T_Z = B * PP_MZ)  -> big => clean z
# Nband big-m clean-z: PP_NHI=3000 PP_NX=250 PP_MZ=800 (and PP_DR=0, the default).
net = Net(**NET_KW)
N_HI = int(os.environ.get('PP_NHI', N_REF))                 # size-schedule upper bound
NX = int(os.environ.get('PP_NX', 120)); MZ = int(os.environ.get('PP_MZ', 120))
T_X = B * NX; T_Z = B * MZ            # fixed totals -> jit-stable ragged, no padding
_OCTS = np.geomspace(4, N_HI - 3, B)  # per-sample target sizes covering all octaves
print(f"SIZES N_HI={N_HI} T_X={T_X}(NX={NX}) T_Z={T_Z}(MZ={MZ}) DR={DR}")


def _fix_total(n, total, lo=2, hi=N_HI - 3, rng=None):
    n = np.clip(n, lo, hi).astype(int)
    idx = (rng.permutation(len(n)) if rng is not None else np.arange(len(n))); i = 0
    while n.sum() != total and i < 60 * len(n):
        j = idx[i % len(n)]; s = 1 if n.sum() < total else -1
        if lo <= n[j] + s <= hi:
            n[j] += s
        i += 1
    return n


def _compose(rng, total, lo=2, hi=N_HI - 3):
    """B per-sample sizes summing to `total`. SCHEDULE:
    'strat'  -> log-spread over all octaves + jitter (dense large-n gradient);
    'dirichlet' -> low-concentration Dirichlet (mass on few, favours small n).
    Real participants only (no padding)."""
    if SCHEDULE == 'strat':
        n = _OCTS * rng.uniform(0.7, 1.4, B)
        n = n * total / n.sum()
    else:
        n = rng.dirichlet(np.full(B, 0.3)) * total
    return _fix_total(np.round(n), total, lo, hi, rng)


def sample_batch(rng):
    nb = _compose(rng, T_X)                                        # (B,) observed sizes
    empty_z = MEMPTY > 0 and rng.random() < MEMPTY                 # whole-batch empty future cohort
    mb = np.zeros(B, dtype=int) if empty_z else _compose(rng, T_Z)
    xl, zl, ys = [], [], []
    for b in range(B):
        beta, P = draw_params(rng)
        ys.append(theta_level_rho(beta, P))
        xl.append(sim_cohort(beta, P, int(nb[b]), rng))
        # domain randomisation (Gap-2): draw z from a time-effect-perturbed
        # parameter so the future cohort is not a clean sample of theta^(s),
        # mimicking the deployment posterior-predictive z^(s) (target rho stays
        # the x-side truth). DR=0 recovers the plain §14.4.13 recipe.
        beta_z = beta + rng.normal(0.0, DR, size=2) if DR > 0 else beta
        if mb[b] > 0:
            zl.append(sim_cohort(beta_z, P, int(mb[b]), rng))
    rho_b = np.stack(ys)                                           # (B, J)
    x_flat = np.concatenate(xl, 0)
    z_flat = np.concatenate(zl, 0) if zl else np.zeros((0, J, 2), np.float32)
    x_seg = np.repeat(np.arange(B), nb).astype(np.int32)
    z_seg = np.repeat(np.arange(B), mb).astype(np.int32)           # empty when all mb=0
    qidx = rng.integers(0, J, size=(B, Q)).astype(np.int32)
    inv_m = np.where(mb > 0, 1 / np.sqrt(np.maximum(mb, 1)), 0.0)  # 0 flags "no future data"
    aux = np.stack([1/np.sqrt(nb), inv_m, nb/N_REF, mb/N_REF], -1).astype(np.float32)  # (B,4)
    y = rho_b[np.arange(B)[:, None], qidx] / sig[qidx]            # standardised rho (B,Q)
    if ETA0AM:                                                   # Option 4: amortise eta_0
        eta0 = rng.uniform(-2.0, 2.0, B).astype(np.float32)      # standardised threshold, per example
        aux = np.concatenate([aux, eta0[:, None]], -1)           # (B,5): eta_0 as last aux feature
        y = (y > eta0[:, None]).astype(np.float32)               # binary success target (B,Q)
    return (dict(x_flat=jnp.asarray(x_flat), x_seg=jnp.asarray(x_seg),
                 z_flat=jnp.asarray(z_flat), z_seg=jnp.asarray(z_seg),
                 item_metadata=jnp.asarray(np.tile(META[None], (B, 1, 1))),
                 query_idx=jnp.asarray(qidx), aux=jnp.asarray(aux)),
            jnp.asarray(y.astype(np.float32)))


b0, y0 = sample_batch(rng)
params = net.init(jax.random.PRNGKey(0), b0)
if INIT:                                         # warm-start (transfer fine-tune from a pretrained net)
    from amortiser_common import load_fitted_model as _lfm
    params = _lfm(INIT)['params']; print(f"warm-started params from {INIT}")


def interval_score(preds, y):
    pr = jnp.maximum.accumulate(preds, -1)
    l05, l25, med, u75, u95 = pr[..., 0], pr[..., 1], pr[..., 2], pr[..., 3], pr[..., 4]
    is10 = (u95 - l05) + 20.0*jnp.maximum(l05 - y, 0) + 20.0*jnp.maximum(y - u95, 0)
    is50 = (u75 - l25) + 4.0*jnp.maximum(l25 - y, 0) + 4.0*jnp.maximum(y - u75, 0)
    return jnp.mean(is10 + is50 + jnp.abs(y - med))


sched = optax.warmup_cosine_decay_schedule(1e-4, 1e-3, NET_STEPS // 10, NET_STEPS, 1e-5)
opt = optax.adam(sched); ost = opt.init(params)


_TAUS_J = jnp.asarray(TAUS)


def pinball(preds, y):                             # general quantile loss over the TAUS grid
    pr = jnp.maximum.accumulate(preds, -1)
    d = y[..., None] - pr
    return jnp.mean(jnp.maximum(_TAUS_J * d, (_TAUS_J - 1.0) * d))


def bce(preds, y):                                 # binary cross-entropy for the success-prob head
    p = jnp.clip(preds[..., 0], 1e-6, 1.0 - 1e-6)
    return -jnp.mean(y * jnp.log(p) + (1.0 - y) * jnp.log(1.0 - p))


def _loss(params, batch, y):
    if ETA0AM:                                     # Option 4: success-prob head
        return bce(net.apply(params, batch), y)
    if NQ != 5:                                    # Option 3: dense quantiles
        return pinball(net.apply(params, batch), y)
    if PREG > 0:                                   # powerlaw exponent -> 1/2 (BvM) prior
        out, st = net.apply(params, batch, mutable=['intermediates'])
        pe = st['intermediates']['p_exp'][0]
        return interval_score(out, y) + PREG * jnp.mean((pe - 0.5) ** 2)
    return interval_score(net.apply(params, batch), y)


@jax.jit
def step(params, ost, batch, y):
    loss, g = jax.value_and_grad(_loss)(params, batch, y)
    up, ost = opt.update(g, ost, params)
    return optax.apply_updates(params, up), ost, loss


print(f"\ntraining ragged net (prior proposal, {NET_STEPS} steps, B={B}) ...")
t0 = time.time()
for it in range(NET_STEPS):
    batch, y = sample_batch(rng)
    params, ost, loss = step(params, ost, batch, y)
    if (it + 1) % max(1, NET_STEPS // 20) == 0:
        print(f"  step {it+1}/{NET_STEPS}  IS={float(loss):.4f}")
mins = (time.time() - t0) / 60.0
fit = {'params': params, 'apply_fn': net.apply, 'net_class_name': type(net).__name__,
       'net_class_module': type(net).__module__, 'net_kwargs': dict(NET_KW),
       'training_mins': mins, 'pps_ProbH1_lwr_quantiles_mesh': tuple(TAUS.tolist())}
save_trained_model(fit, f"{dir_out}/{file_prefix}_amortised_pps_net.pkl")
print(f"saved ({mins:.1f} min) -> {dir_out}")

if ETA0AM:                     # success-prob head: the quantile-width probe below does not apply
    raise SystemExit(0)


# ---------- proposal-probe: width vs n (fresh prior draws, m fixed) ----------
@jax.jit
def fwd(params, batch):
    return net.apply(params, batch)


m_fix = 100
pg = np.unique(np.round(np.geomspace(6, N_REF - 3, 8)).astype(int))
pw = []
for n in pg:
    ws = []
    for _ in range(120):
        beta, P = draw_params(rng)
        xb = sim_cohort(beta, P, int(n), rng); zb = sim_cohort(beta, P, m_fix, rng)
        batch = dict(x_flat=jnp.asarray(xb), x_seg=jnp.zeros(int(n), jnp.int32),
                     z_flat=jnp.asarray(zb), z_seg=jnp.zeros(m_fix, jnp.int32),
                     item_metadata=jnp.asarray(META[None]), query_idx=jnp.arange(J)[None],
                     aux=jnp.asarray(np.array([[1/np.sqrt(n), 1/np.sqrt(m_fix), n/N_REF, m_fix/N_REF]], np.float32)))
        pr = np.maximum.accumulate(np.asarray(fwd(params, batch))[0], 1) * sig[:, None]
        ws.append(pr[:, 4] - pr[:, 0])
    pw.append((int(n), float(np.median(ws))))
pwa = np.array(pw); sl = np.polyfit(np.log(pwa[:, 0]), np.log(pwa[:, 1]), 1)[0]
pd.DataFrame(pw, columns=['n', 'width']).to_csv(f"{dir_out}/{file_prefix}_proposal_probe.csv", index=False)
print("\nPROPOSAL-PROBE width vs n:")
for n, w in pw:
    print(f"  n={n:4d}  width={w:.3f}")
print(f"width slope: {sl:.3f}  (focused-proposal build was -0.05; oracle {os_slope:.3f}; SVI -0.59)")
