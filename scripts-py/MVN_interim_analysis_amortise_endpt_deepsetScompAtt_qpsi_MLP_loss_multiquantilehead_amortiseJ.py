#!/usr/bin/env python3
"""
Amortise over the number of components J (§3.3.1 goal). ONE ragged
``deepsetXcompAtt`` net (§14.4.3) is trained with a per-step J-curriculum over
J in [J_MIN, J_MAX] = [2, 100], with the covariance structure held FIXED as
block-equicorrelation (§3.3.1 Cell B: rho_w=0.8, rho_b=0.1, block size 10, ragged
last block) and only J varying. The net's parameters are J-invariant by
construction (J is a runtime item axis; no positional J-embedding), so a single
trained net prices any J.

Evaluation sweeps a J-grid and, for each J, builds the interim schedule against
the CLOSED-FORM Gaussian reference (exact for the sigma-known MVN, so no HMC is
needed): per interim the posterior of rho_j = mu_j - 1 is
    mu_n,j = (mu_0/tau^2 + sum_i x_ij) / c_n,  c_n = 1/tau^2 + n/sigma^2,
    Var(mu_j|x) = K_jj / c_n = 1/c_n           (unit diagonal),
and the joint-draw target rho^(s) with matched future cohort z^(s) is drawn from
N(mu_n, K/c_n). The prior-trained net is deployed, calibrated with the expanding
head fine-tune + affine median-shift (§14.4.6), and scored against the closed-form
PPS (Phi tail) and the reference draws (PIT-KS, marg-KS). Metrics are reported as
a function of J to demonstrate J-invariance.

Outputs (dir_out):
  - ``mvn_amortiseJ_net_{head_mode}.pkl``          the J-amortised net
  - ``mvn_amortiseJ_metrics_by_J.csv``             per-J PPS-MSE / PIT-KS / marg-KS
  - ``mvn_amortiseJ_metrics_by_J.pdf``             metrics vs J

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/MVN_interim_analysis_amortise_endpt_deepsetXcompAtt_qpsi_MLP_loss_multiquantilehead_amortiseJ.py
Env: MVNJ_STEPS (8000), MVNJ_S (200 eval draws), MVNJ_JGRID ('2,5,10,20,35,50,75,100'),
     MVNJ_JMIN (2), MVNJ_JMAX (100), MVNJ_HEADMODE (powerlaw), MVNJ_SMOKE (1).
"""

# %%

import os
import sys
import time
import warnings
from functools import partial
from pathlib import Path

try:
    project_root = Path(__file__).resolve().parent.parent
except NameError:
    project_root = Path.cwd()
    if project_root.name == 'scripts-py':
        project_root = project_root.parent
sys.path.insert(0, str(project_root / 'python'))

import numpy as np
import pandas as pd
import jax
import jax.numpy as jnp
import optax
from scipy.stats import norm as _norm

warnings.filterwarnings('ignore')

from model_mvn import make_block_equicorr_chol_anyJ, sample_ragged_block_equicorr
from amortiser_common import save_trained_model, load_fitted_model
import amortiser_calibration as cal          # shared §14.4.6/§14.4.7 calibration
from amortiser_pps_features_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead import (
    Amortiser_PPS_features_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead as Net,
    _MLP,
)
from plotnine import (ggplot, aes, geom_line, geom_point, facet_wrap, theme_bw,
                      theme, labs, element_text, element_blank, scale_x_continuous)

print("Imports successful")

# %%
# =============================================================================
# Configuration (§3.3.1 Cell B, sigma^2 known Gaussian special case).
# =============================================================================

_sandbox = "/Users/or105/sandbox/bIRTistic"
dir_out = os.path.join(_sandbox, "py-mvn-interim-amortiseJ-260915")
os.makedirs(dir_out, exist_ok=True)

SMOKE = os.environ.get('MVNJ_SMOKE', '0') == '1'
STEPS = int(os.environ.get('MVNJ_STEPS', '30' if SMOKE else '8000'))
S_EVAL = int(os.environ.get('MVNJ_S', '20' if SMOKE else '200'))
J_MIN = int(os.environ.get('MVNJ_JMIN', '2'))
J_MAX = int(os.environ.get('MVNJ_JMAX', '100'))
J_GRID = [int(x) for x in os.environ.get(
    'MVNJ_JGRID', '2,5,10,20,35,50,75,100').split(',')]
HEAD_MODE = os.environ.get('MVNJ_HEADMODE', 'powerlaw')
if SMOKE:
    J_GRID = [2, 10, 50]

# §3.3.1 constants
N_TOTAL = 500                          # total cohort accrued over interims
INTERIM_N = list(range(50, N_TOTAL, 50))          # n = 50,100,...,450 (m>0)
SIGMA = 1.0; SIGMA2 = 1.0; TAU2 = 100.0; MU0 = 1.0
RHO_W, RHO_B, BLOCK = 0.8, 0.1, 10
DELTA = 0.3                            # true half-effect
ETA0 = 0.0; ETAH = 0.89
Z_ETA = float(_norm.ppf(ETAH))
TAUS = (0.05, 0.25, 0.5, 0.75, 0.95)
ZQ = np.array([-1.6449, -0.6745, 0.0, 0.6745, 1.6449])
NET_EMBED_DIM = 32; NET_HIDDEN = (64, 64); NET_LR = 1e-3; NET_BATCH = 32
NET_QUERIES = 4; HEADFT_STEPS = 40 if SMOKE else 600
seed = 123
if SMOKE:
    INTERIM_N = [50, 250, 450]
print(f"J-amortise: train J in [{J_MIN},{J_MAX}] head={HEAD_MODE} steps={STEPS}; "
      f"eval J_GRID={J_GRID} S={S_EVAL} interims n={INTERIM_N} smoke={SMOKE}")

# %%
# =============================================================================
# Train one J-invariant net (per-step J-curriculum, fixed block-equicorr K).
# =============================================================================


def train_amortiseJ():
    ckpt = os.path.join(dir_out, f'mvn_amortiseJ_net_{HEAD_MODE}.pkl')
    if os.path.exists(ckpt):
        print(f"load cached J-amortised net {ckpt}")
        fit = load_fitted_model(ckpt)
        return fit, float(fit.get('training_mins', 0.0))
    net_kwargs = dict(head_mode=HEAD_MODE, embed_dim=NET_EMBED_DIM,
                      hidden_dims=NET_HIDDEN, num_quantiles=len(TAUS))
    net = Net(**net_kwargs)
    taus = jnp.asarray(TAUS, jnp.float32)
    rng = np.random.default_rng(seed)

    def _sample(J, S):
        return sample_ragged_block_equicorr(
            rng, S, J, n_max=N_TOTAL, prior_tau=np.sqrt(TAU2), sigma=SIGMA,
            mu_0_baseline=MU0, rho_w=RHO_W, rho_b=RHO_B, block_size=BLOCK,
            queries_per_sample=NET_QUERIES, n_dist='log_uniform')

    peek, _ = _sample(10, 4)
    peek = jax.tree_util.tree_map(
        lambda a: jnp.asarray(a, jnp.float32) if np.asarray(a).dtype.kind == 'f'
        else jnp.asarray(a), peek)
    params = net.init(jax.random.PRNGKey(seed), peek)
    warmup = min(max(200, STEPS // 20), STEPS // 2) if STEPS >= 100 else 0
    sched = (optax.warmup_cosine_decay_schedule(0.0, NET_LR, warmup, STEPS, NET_LR * 0.05)
             if STEPS >= 100 else NET_LR)
    opt = optax.adam(sched); ost = opt.init(params)

    @jax.jit
    def step(params, ost, batch, y):
        def loss(p):
            pr = net.apply(p, batch)
            e = y[..., None] - pr
            return jnp.mean(jnp.maximum(taus * e, (taus - 1.0) * e))
        l, g = jax.value_and_grad(loss)(params)
        up, ost2 = opt.update(g, ost, params)
        return optax.apply_updates(params, up), ost2, l

    t0 = time.time(); hist = []
    for it in range(STEPS):
        J = int(rng.integers(J_MIN, J_MAX + 1))            # per-step J curriculum
        batch, rho = _sample(J, NET_BATCH)
        batch = jax.tree_util.tree_map(
            lambda a: jnp.asarray(a, jnp.float32) if np.asarray(a).dtype.kind == 'f'
            else jnp.asarray(a), batch)
        params, ost, l = step(params, ost, batch, jnp.asarray(rho, jnp.float32))
        hist.append(float(l))
        if (it + 1) % max(1, STEPS // 20) == 0 or it == 0:
            print(f"  step {it+1}/{STEPS}  J={J}  loss={float(l):.5f}")
    mins = (time.time() - t0) / 60.0
    fit = {'params': params, 'apply_fn': net.apply,
           'pps_ProbH1_lwr_quantiles_mesh': np.asarray(TAUS, np.float32),
           'net_class_module': Net.__module__, 'net_class_name': Net.__name__,
           'net_kwargs': net_kwargs, 'loss_history': hist, 'training_mins': mins}
    save_trained_model(fit, ckpt)
    print(f"trained J-amortised net in {mins:.1f} min -> {ckpt}")
    return fit, mins


# %%
# =============================================================================
# Closed-form reference per J + deploy + calibrate + metrics.
# =============================================================================


def eval_J(fit, J):
    """Closed-form interim schedule at this J, deploy + head-ft/affine, metrics."""
    net = Net(**dict(fit['net_kwargs']))

    @jax.jit
    def _fwd(params, batch):
        out, st = net.apply(params, batch, mutable=['intermediates'])
        return out, st['intermediates']['head_in'][0]

    rng = np.random.default_rng(seed + J)
    K_chol = make_block_equicorr_chol_anyJ(J, BLOCK, RHO_W, RHO_B)
    mu_true = np.where(np.arange(J) < J / 2, MU0 + DELTA, MU0 - DELTA)
    X = mu_true[None, :] + SIGMA * (rng.standard_normal((N_TOTAL, J)) @ K_chol.T)   # (N,J)

    QS, HD, TGT, PPS_CF, NOBS = {}, {}, {}, {}, {}
    for k, n in enumerate(INTERIM_N):
        m = N_TOTAL - n
        c_n = 1.0 / TAU2 + n / SIGMA2
        c_np = 1.0 / TAU2 + (n + m) / SIGMA2
        mu_n = (MU0 / TAU2 + X[:n].sum(0)) / c_n                     # (J,) per-component
        # closed-form PPS (H1: mu_j > MU0), Phi tail (fit_closed_form_pps)
        threshold = MU0 + Z_ETA * np.sqrt(1.0 / c_np)
        var_pred = 1.0 * m / (SIGMA2 * c_n * c_np)
        PPS_CF[k] = _norm.cdf((mu_n - threshold) / np.sqrt(var_pred))   # (J,)
        # joint draws: mu^(s) ~ N(mu_n, K/c_n); z^(s) from mu^(s)
        cov_chol = K_chol / np.sqrt(c_n)
        mu_s = mu_n[None, :] + rng.standard_normal((S_EVAL, J)) @ cov_chol.T   # (S,J)
        TGT[k] = (mu_s - MU0)                                          # rho^(s) (S,J)
        NOBS[k] = n
        # deploy: one cohort, Q=J queries; z^(s) changes across s
        x_flat = X[:n].astype(np.float32)[:, :, None]
        x_seg = np.zeros(n, np.int32)
        meta_b = np.ones((1, J, 1), np.float32)
        qidx = np.arange(J, dtype=np.int32)[None, :]
        aux_b = np.array([[1.0 / np.sqrt(n), 1.0 / np.sqrt(m), n / N_TOTAL, m / N_TOTAL]], np.float32)
        qs = np.empty((S_EVAL, J, 5)); hd = None
        for s in range(S_EVAL):
            Yz = (mu_s[s][None, :] + SIGMA * (rng.standard_normal((m, J)) @ K_chol.T))
            batch = dict(x_flat=x_flat, x_seg=x_seg, z_flat=Yz.astype(np.float32)[:, :, None],
                         z_seg=np.zeros(m, np.int32), item_metadata=meta_b,
                         query_idx=qidx, aux=aux_b)
            out, head = _fwd(fit['params'], batch)
            qs[s] = np.maximum.accumulate(np.asarray(out)[0], 1)
            if hd is None:
                hd = np.empty((S_EVAL, J, np.asarray(head).shape[-1]), np.float32)
            hd[s] = np.asarray(head)[0]
        QS[k], HD[k] = qs, hd.astype(np.float32)

    # ---- expanding head-ft + affine (§14.4.6) via the shared calibration module ----
    KS = list(range(len(INTERIM_N)))
    labels = np.array([f'mu_{j}' for j in range(J)])
    HIDDEN = tuple(dict(fit['net_kwargs']).get('hidden_dims', NET_HIDDEN))
    HK = cal.expanding_head_ft(HD, TGT, NOBS, KS, fit['params']['params']['q_psi'],
                               HIDDEN, HEAD_MODE, np.ones(J), taus=TAUS, steps=HEADFT_STEPS)
    HK = cal.affine_median_shift(HK, TGT, KS, taus=TAUS)
    summ = cal.calibration_summary(HK, TGT, NOBS, KS, labels, taus=TAUS)
    se = []
    for k in KS:
        for j in range(J):
            qs = HK[k][:, j, :]
            ph1 = 1.0 - np.array([np.interp(ETA0, qs[s], TAUS, 0., 1.) for s in range(qs.shape[0])])
            se.append((float(np.mean(ph1 > ETAH)) - PPS_CF[k][j]) ** 2)
    return dict(J=J, pit_ks=float(summ.pit_ks.mean()), marg_ks=float(summ.marg_ks.mean()),
                pps_mse=float(np.nanmean(se)), n_interims=len(KS))


# %%
# =============================================================================
# Main.
# =============================================================================

fit, mins_train = train_amortiseJ()
rows = []
for J in J_GRID:
    t0 = time.time()
    r = eval_J(fit, J); r['eval_min'] = round((time.time() - t0) / 60.0, 2)
    rows.append(r)
    print(f"  J={J:3d}: PPS-MSE={r['pps_mse']:.5f}  PIT-KS={r['pit_ks']:.3f}  "
          f"marg-KS={r['marg_ks']:.3f}  ({r['eval_min']}m)")
res = pd.DataFrame(rows)
res.to_csv(os.path.join(dir_out, 'mvn_amortiseJ_metrics_by_J.csv'), index=False)

long = res.melt(id_vars='J', value_vars=['pps_mse', 'pit_ks', 'marg_ks'],
                var_name='metric', value_name='value')
long['metric'] = pd.Categorical(long['metric'], ['pps_mse', 'pit_ks', 'marg_ks'], ordered=True)
p = (ggplot(long, aes('J', 'value')) + geom_line(colour='#1f77b4') + geom_point(size=1.5)
     + facet_wrap('~ metric', scales='free_y', ncol=3)
     + scale_x_continuous(breaks=J_GRID)
     + theme_bw() + theme(figure_size=(13, 4), strip_background=element_blank(),
                          strip_text=element_text(face='bold'),
                          axis_text_x=element_text(angle=45, hjust=1))
     + labs(x='J (number of components)', y='',
            title=f'J-amortised deepsetXcompAtt ({HEAD_MODE}) vs closed-form — metrics by J'))
p.save(os.path.join(dir_out, 'mvn_amortiseJ_metrics_by_J.pdf'), verbose=False)
print(f"\nSaved metrics + plot -> {dir_out}")
print(res.to_string(index=False))
