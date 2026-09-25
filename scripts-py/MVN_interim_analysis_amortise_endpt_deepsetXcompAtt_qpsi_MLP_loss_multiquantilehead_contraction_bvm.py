#!/usr/bin/env python3
"""
MVN sanity case for the §14.4.6/§14.4.7 deployment calibration stack:
expanding-window head fine-tune + per-item affine median-shift (§14.4.6),
the parametric contraction heads **A power-law** / **C floor** (§14.4.7), and
the Bernstein-von Mises self-consistency correction (§14.4.7). Reported in
§13.8 of ``dev/amortised_decision_making.md``.

Why MVN. Here the reference posterior of the endpoint ``rho_j = mu_j - 1`` is
**closed-form Gaussian**, so Bernstein-von Mises is EXACT: ``SD(rho_j|x) =
sqrt(K_jj / (1/tau^2 + n/sigma^2)) ~ n^{-1/2}``. The A power-law head should
recover the exponent ``p ~ 1/2``, the C floor head a plateau ``a ~ 0``, and the
BvM correction should collapse the deepset's finite-``m`` large-``n`` marginal
upturn exactly onto the closed-form contraction curve. This is the controlled
check that the Ukraine PCM machinery (§14.4) is correct where the truth is known.

The deepset encoder (``deepsetXcompAtt``) is reused unchanged; only the head is
parametric (``head_mode``). The reference marginal + joint-draw target both come
from the cached per-interim HMC ``mu_draws`` (aligned to ``zi.ypred_s``), i.e.
``rho^(s) = mu_draws[s] - mu_0_baseline`` is a valid draw from ``p(rho|x,z^s)``
(the §14.1.7 joint-draw construction). All calibration/BvM math is a verbatim
port of the ragged Ukraine deploy driver, with SVI replaced by the closed-form
reference and per-item standardiser ``sig = 1``.

Per config (SUF in {RGXP plain+affine, RGXA A-powerlaw+BvM, RGXC C-floor+BvM})
writes, in ``dir_out``:
  - ``mvn_J{J}_pps_{SUF}_p_h1_xz.pkl`` / ``.csv`` / ``_timing.csv`` / ``_perf.csv``
    (compare-methods schema; picked up by MVN_interim_analyses_compare_methods.py)
  - the four §14.4.5 diagnostic PDFs (pit box, marginal-quantile box,
    contraction-factor, contraction-law) via python/amortiser_diag_plots.py
  - ``mvn_J{J}_{SUF}_calibrated_qs.pkl`` (per-interim calibrated (S,J,5) arrays)
  - ``mvn_calibration_summary.csv`` (PIT-KS / marg-KS per config, for §13.8)

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/MVN_interim_analysis_amortise_endpt_deepsetXcompAtt_qpsi_MLP_loss_multiquantilehead_contraction_bvm.py
Env: MVN_STEPS (train steps, 4000), MVN_S (deploy draws, 200), MVN_J (20),
     MVN_CONFIGS ('plain,powerlaw,floor'), MVN_SMOKE (1 -> tiny steps/S/interims).
"""

# %%

import os
import sys
import time
import warnings
from functools import partial
from pathlib import Path

try:
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
except NameError:
    project_root = Path.cwd()
    if project_root.name == 'scripts-py':
        project_root = project_root.parent

python_path = str(project_root / 'python')
if python_path not in sys.path:
    sys.path.insert(0, python_path)

import numpy as np
import pandas as pd
import jax
import jax.numpy as jnp
import optax

warnings.filterwarnings('ignore')

from model_mvn import MVNModel
from amortiser_common import train, save_trained_model, load_fitted_model, _pinball_loss
import amortiser_diag_plots as adp
import amortiser_calibration as cal
import model_mvn_common as mc        # shared §14.4.6/§14.4.7 calibration (MVN + Ukraine)

# §14.4.3 ragged encoder switch. MVN_RAGGED=1 uses the ragged-participant-axis
# deepsetXcompAtt (no padding/mask, segmented mean-pool) instead of the padded
# variant. Same head modes; the K-row hand feature is dropped (learned instead).
RAGGED = os.environ.get('MVN_RAGGED', '0') == '1'
if RAGGED:
    from amortiser_pps_features_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead import (
        Amortiser_PPS_features_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead as Net,
        _MLP,
    )
else:
    from amortiser_pps_features_deepsetXcompAtt_qpsi_MLP_loss_multiquantilehead import (
        Amortiser_PPS_features_deepsetXcompAtt_qpsi_MLP_loss_multiquantilehead as Net,
        _MLP,
    )

print(f"Imports successful (RAGGED={RAGGED})")

# %%

# =============================================================================
# Configuration.
# =============================================================================

_sandbox = "/Users/or105/sandbox/bIRTistic"
DIR_SIM = os.path.join(_sandbox, "py-mvn-interim-simulations-260609")
DIR_PLAIN = os.path.join(                       # existing plain-head deepsetXcompAtt net
    _sandbox,
    "py-mvn-interim-amortise-endptx-on-wz-with-features-"
    "deepsetXcompAtt-qpsi-MLP-loss-multiquantilehead_260803",
)
dir_out = os.path.join(
    _sandbox,
    "py-mvn-interim-amortise-deepsetXcompAtt-ragged-bvm-260914" if RAGGED
    else "py-mvn-interim-amortise-deepsetXcompAtt-bvm-260914")
os.makedirs(dir_out, exist_ok=True)
print(f"Output dir: {dir_out}")

SMOKE = os.environ.get('MVN_SMOKE', '0') == '1'
NET_STEPS = int(os.environ.get('MVN_STEPS', '20' if SMOKE else '4000'))
PPS_Z_TOTAL = int(os.environ.get('MVN_S', '12' if SMOKE else '200'))
J_DEPLOY = int(os.environ.get('MVN_J', '20'))
CONFIGS = os.environ.get('MVN_CONFIGS', 'plain,powerlaw,floor').split(',')

TAUS = (0.05, 0.25, 0.5, 0.75, 0.95)
ZQ = np.array([-1.6449, -0.6745, 0.0, 0.6745, 1.6449])
NET_EMBED_DIM = 32
NET_HIDDEN = (64, 64)
NET_LR = 1e-3
NET_BATCH = 32
NET_QUERIES = 4
N_DIST = 'uniform'
K_FAMILIES = ('identity', 'ar1', 'block', 'factor')
HEADFT_STEPS = 40 if SMOKE else 800

# config -> (head_mode, use_bvm, SUF, human label)
CFG = {
    'plain':    ('plain',    False, 'RGXP', 'plain + affine'),
    'powerlaw': ('powerlaw', True,  'RGXA', 'A power-law + BvM'),
    'floor':    ('floor',    True,  'RGXC', 'C floor + BvM'),
}

# %%

# =============================================================================
# Load upstream simulation data (on-disk simu_params uses the OLD pps key names).
# =============================================================================

sim_data = pd.read_pickle(os.path.join(DIR_SIM, 'mvn_sim_data.pkl'))
simu_params = sim_data['simu_params']
seed = int(simu_params['seed'])
N_FULL = int(simu_params['N_full'])                     # 1050 == n_max == N_REF (n+m)
MU0 = float(simu_params['mu_0_baseline'])               # 1.0; rho_j = mu_j - MU0
ETA0 = float(simu_params.get('pps_H1_min_effect_size_thresh',
                             simu_params.get('pps_H1_def', 0.0)))          # H1: rho>0
ETAH = float(simu_params.get('pps_ProbH1_target_lwr_quantile',
                             simu_params.get('pps_ProbH1_thresh', 0.89)))
TRAIN_J = 20

interim_pkl = pd.read_pickle(os.path.join(DIR_SIM, f'mvn_J{J_DEPLOY}_interim_data.pkl'))
pps_cf = pd.read_pickle(os.path.join(DIR_SIM, 'mvn_pps_closed_form.pkl'))['pps_cf']

cell = sim_data['cells'][J_DEPLOY]
dit, R_chol = cell['dit'], cell['R_chol']
K = (R_chol @ R_chol.T).astype(np.float64)
K_diag = np.diag(K).astype(np.float64)
labels = dit['item_label'].to_numpy()
item_types = dit['item_type'].to_numpy()
item_highs = dit['item_high_label'].to_numpy()
sig = np.ones(J_DEPLOY, dtype=np.float64)              # MVN rho is unscaled

INTERIMS = sorted(int(k) for k in interim_pkl.keys())
if SMOKE:
    INTERIMS = INTERIMS[:2]
S = min(PPS_Z_TOTAL, min(interim_pkl[k]['mu_draws'].shape[0] for k in INTERIMS))
print(f"J={J_DEPLOY}  interims={INTERIMS}  S={S}  N_FULL={N_FULL}  "
      f"eta0={ETA0}  etaH={ETAH}  configs={CONFIGS}  smoke={SMOKE}")

# per-interim date labels for plots / outputs
DATEOF = {k: interim_pkl[k]['interim_date'] for k in INTERIMS}
MONTHYR = {k: interim_pkl[k]['interim_month_year'] for k in INTERIMS}
NOBS = {k: int(interim_pkl[k]['n_obs']) for k in INTERIMS}
MFUT = {k: int(interim_pkl[k]['interim_m']) for k in INTERIMS}


# %%

# =============================================================================
# Deploy-side head transform, mirroring the net's head_mode branches exactly.
# raw: (..., 5) q_psi output; prec: n^{-1/2} scalar or (...,1). Returns real-scale
# quantiles (sig applied by the caller).
# =============================================================================


_transform = cal.head_transform                        # shared deploy-side head map


# %%

# =============================================================================
# Build per-interim raw forward arrays for a given trained net:
#   qs  (S,J,5)  raw net output (real scale)
#   hd  (S,J,H)  head_in (frozen-encoder embedding, for the head fine-tune)
#   tgt (S,J)    joint-draw reference rho^(s) = mu_draws[s] - MU0
# =============================================================================


def build_forward(fit, head_mode):
    net = Net(**dict(fit['net_kwargs']))

    @jax.jit
    def _fwd(params, batch):
        out, st = net.apply(params, batch, mutable=['intermediates'])
        return out, st['intermediates']['head_in'][0]

    QS, HD, TGT = {}, {}, {}
    for k in INTERIMS:
        blk = interim_pkl[k]
        dpi = blk['dpi']; zi_full = blk['zi']; n_obs = NOBS[k]; m = MFUT[k]
        mu_draws = np.asarray(blk['mu_draws'][:S], dtype=np.float64)   # (S,J)
        TGT[k] = (mu_draws - MU0).astype(np.float64)                  # rho^(s)

        x_wide, ypred_arr = mc.xz_arrays(dpi, zi_full, J_DEPLOY, S, m)   # (n,J), (S,m,J)

        qs = np.empty((S, J_DEPLOY, 5), np.float64); hd = None
        if RAGGED:
            # one batch element (the observed cohort), Q=J queries; only z_flat
            # changes across the S posterior-predictive future cohorts.
            x_flat = x_wide.astype(np.float32)[:, :, None]             # (n,J,1)
            x_seg = np.zeros(n_obs, np.int32)
            meta_b = np.ones((1, J_DEPLOY, 1), np.float32)
            qidx = np.arange(J_DEPLOY, dtype=np.int32)[None, :]        # (1,J)
            aux_b = np.array([[1.0 / np.sqrt(n_obs), 1.0 / np.sqrt(max(m, 1)),
                               n_obs / N_FULL, m / N_FULL]], np.float32)
            for s in range(S):
                z_flat = ypred_arr[s].astype(np.float32)[:, :, None]  # (m,J,1)
                batch = dict(x_flat=x_flat, x_seg=x_seg, z_flat=z_flat,
                             z_seg=np.zeros(m, np.int32), item_metadata=meta_b,
                             query_idx=qidx, aux=aux_b)
                out, head = _fwd(fit['params'], batch)                # out (1,J,5); head (1,J,H)
                qs[s] = np.maximum.accumulate(np.asarray(out)[0], 1) * sig[:, None]
                if hd is None:
                    hd = np.empty((S, J_DEPLOY, np.asarray(head).shape[-1]), np.float32)
                hd[s] = np.asarray(head)[0]
        else:
            static = mc.deepset_static_batch(x_wide, K, K_diag, n_obs, m, N_FULL, J_DEPLOY)
            for s in range(S):
                batch = mc.with_z(static, ypred_arr[s])          # shared padded batch builder
                out, head = _fwd(fit['params'], batch)
                qs[s] = np.maximum.accumulate(np.asarray(out), 1) * sig[:, None]
                if hd is None:
                    hd = np.empty((S, J_DEPLOY, np.asarray(head).shape[-1]), np.float32)
                hd[s] = np.asarray(head)
        QS[k], HD[k] = qs, hd.astype(np.float32)
        print(f"    fwd interim {k} (n={n_obs}) done")
    return QS, HD, TGT


# %%

# =============================================================================
# Head fine-tune (expanding) + affine median-shift + BvM, per config.
# =============================================================================


def calibrate(fit, head_mode, use_bvm, QS, HD, TGT):
    """Head-ft expand + affine median-shift (+ BvM). Thin wrapper over the shared
    `amortiser_calibration` module so MVN and Ukraine run one implementation."""
    HK = cal.expanding_head_ft(
        HD, TGT, NOBS, INTERIMS, fit['params']['params']['q_psi'],
        tuple(dict(fit['net_kwargs']).get('hidden_dims', NET_HIDDEN)),
        head_mode, sig, taus=TAUS, steps=HEADFT_STEPS, verbose=True)
    HK = cal.affine_median_shift(HK, TGT, INTERIMS, taus=TAUS)
    if not use_bvm:
        return HK, HK
    HKV, _law = cal.bvm_correction(HK, TGT, NOBS, INTERIMS, N_FULL, verbose=True)
    return HK, HKV


# %%

# =============================================================================
# Scalar diagnostics: PIT-KS (conditional) + marg-KS (amortiser marginal vs
# reference), pooled over items/interims and per cohort-size band.
# =============================================================================


def scalar_summary(PPSSRC, TGT, label, suf):
    """Per (interim, item) PIT-KS + marg-KS via the shared amortiser_calibration layer
    (one implementation shared with Ukraine + by-architecture; adds config/suf columns)."""
    return cal.calibration_summary(PPSSRC, TGT, NOBS, INTERIMS, labels, taus=TAUS,
                                   extra=dict(config=label, suf=suf))


# %%

# =============================================================================
# PPS + p(H1|x,z) outputs in the compare-methods schema.
# =============================================================================


def emit_pps(PPSSRC, suf, mins_deploy, mins_train):
    p_rows, perf_rows, timing_rows = [], [], []
    for k in INTERIMS:
        for j in range(J_DEPLOY):
            q = PPSSRC[k][:, j, :]
            ph1 = 1.0 - np.array([np.interp(ETA0, q[s], TAUS, 0., 1.) for s in range(q.shape[0])])
            p_rows.append(pd.DataFrame({
                'item_label': labels[j], 'item_type': item_types[j],
                'item_high_label': item_highs[j], 'j': j, 's': np.arange(1, S + 1),
                'p_h1_xz': ph1.astype(np.float64), 'q_hat': q[:, 2].astype(np.float64),
                'J': J_DEPLOY, 'interim_id': k, 'interim_date': DATEOF[k],
                'interim_month_year': MONTHYR[k],
            }))
            perf_rows.append(dict(J=J_DEPLOY, interim_id=k, interim_date=DATEOF[k],
                                  interim_month_year=MONTHYR[k], j=j, item_label=labels[j],
                                  item_type=item_types[j], rho=float('nan'), r2=float('nan')))
        timing_rows.append(dict(J=J_DEPLOY, interim_id=k, interim_date=DATEOF[k],
                                interim_month_year=MONTHYR[k], n_obs=NOBS[k], m_future=MFUT[k],
                                mins_interim_id=round(mins_deploy / len(INTERIMS), 3)))
    dp = pd.concat(p_rows, ignore_index=True)
    dp['pps_H1_min_effect_size_thresh'] = ETA0
    dp['pps_ProbH1_target_lwr_quantile'] = ETAH
    dp['S'] = S
    dp.to_pickle(os.path.join(dir_out, f'mvn_J{J_DEPLOY}_pps_{suf}_p_h1_xz.pkl'))
    pd.DataFrame(perf_rows).to_csv(os.path.join(dir_out, f'mvn_J{J_DEPLOY}_pps_{suf}_perf.csv'), index=False)
    tdf = pd.DataFrame(timing_rows); tdf['mins_J_total'] = round(mins_deploy, 3)
    tdf.to_csv(os.path.join(dir_out, f'mvn_J{J_DEPLOY}_pps_{suf}_timing.csv'), index=False)
    pps_df = (dp.groupby(['J', 'interim_id', 'interim_date', 'interim_month_year',
                          'item_label', 'item_type', 'item_high_label', 'j'], observed=True)['p_h1_xz']
              .apply(lambda p: float((p > ETAH).mean())).reset_index(name='pps'))
    pps_df['eta'] = ETAH; pps_df['S'] = S
    pps_df.to_csv(os.path.join(dir_out, f'mvn_J{J_DEPLOY}_pps_{suf}.csv'), index=False)
    # accuracy vs closed-form PPS
    tab = pps_df.merge(pps_cf[pps_cf['J'] == J_DEPLOY][['interim_id', 'j', 'pps']]
                       .rename(columns={'pps': 'pps_analytic'}), on=['interim_id', 'j'], how='left')
    mse = float(((tab['pps'] - tab['pps_analytic']) ** 2).mean())
    return mse


# %%

# =============================================================================
# Main: per config -> net (load plain / train A,C) -> forward -> calibrate ->
# diagnostics -> outputs.
# =============================================================================


def _train_ragged(head_mode, steps):
    """Custom training loop for the ragged encoder: the sampler emits (B,Q)
    targets and the net returns (B,Q,K), so amortiser_common.train (which
    assumes (B,)/(B,K)) cannot be reused. Warmup-cosine Adam + multi-quantile
    pinball, mirroring amortiser_common.train."""
    cell_tr = sim_data['cells'][TRAIN_J]
    model = MVNModel(dit=cell_tr['dit'], dcati=MVNModel.get_interim_data_x(cell_tr['dp']),
                     seed=seed, J=TRAIN_J, K_chol=cell_tr['R_chol'], sigma=simu_params['sigma'],
                     prior_tau=simu_params['prior_tau'], mu_0_baseline=MU0)
    sample_fn = partial(model.make_training_data_ragged_random_K,
                        n_max=N_FULL, n_dist=N_DIST, queries_per_sample=NET_QUERIES,
                        K_families=K_FAMILIES)
    net_kwargs = dict(head_mode=head_mode, embed_dim=NET_EMBED_DIM,
                      hidden_dims=NET_HIDDEN, num_quantiles=len(TAUS))
    net = Net(**net_kwargs)
    taus = jnp.asarray(TAUS, jnp.float32)
    rng = np.random.default_rng(seed)
    peek, _ = sample_fn(rng, 4)
    peek = jax.tree_util.tree_map(lambda a: jnp.asarray(a, jnp.float32)
                                  if np.asarray(a).dtype.kind == 'f' else jnp.asarray(a), peek)
    params = net.init(jax.random.PRNGKey(seed), peek)
    warmup = min(max(200, steps // 20), steps // 2) if steps >= 100 else 0
    sched = (optax.warmup_cosine_decay_schedule(0.0, NET_LR, warmup, steps, NET_LR * 0.05)
             if steps >= 100 else NET_LR)
    opt = optax.adam(sched); ost = opt.init(params)

    @jax.jit
    def step(params, ost, batch, y):
        def loss(p):
            pr = net.apply(p, batch)                      # (B,Q,K)
            e = y[..., None] - pr
            return jnp.mean(jnp.maximum(taus * e, (taus - 1.0) * e))
        l, g = jax.value_and_grad(loss)(params)
        up, ost2 = opt.update(g, ost, params)
        return optax.apply_updates(params, up), ost2, l

    t0 = time.time(); hist = []
    for it in range(steps):
        batch, rho = sample_fn(rng, NET_BATCH)
        batch = jax.tree_util.tree_map(lambda a: jnp.asarray(a, jnp.float32)
                                       if np.asarray(a).dtype.kind == 'f' else jnp.asarray(a), batch)
        params, ost, l = step(params, ost, batch, jnp.asarray(rho, jnp.float32))
        hist.append(float(l))
        if (it + 1) % max(1, steps // 20) == 0 or it == 0:
            print(f"  step {it+1}/{steps}  loss={float(l):.5f}")
    return {
        'params': params, 'apply_fn': net.apply,
        'pps_ProbH1_lwr_quantiles_mesh': np.asarray(TAUS, np.float32),
        'net_class_module': Net.__module__, 'net_class_name': Net.__name__,
        'net_kwargs': net_kwargs, 'loss_history': hist,
        'training_mins': float((time.time() - t0) / 60.0),
    }


def get_net(head_mode):
    """Load a cached net or train fresh. Padded 'plain' reuses the 260803 net;
    ragged always trains (different architecture/inputs)."""
    if not RAGGED and head_mode == 'plain':
        ckpt = os.path.join(DIR_PLAIN, 'mvn_interim_amortised_pps_net.pkl')
        if os.path.exists(ckpt):
            print(f"  load plain net {ckpt}")
            return load_fitted_model(ckpt), 0.0
    ckpt = os.path.join(dir_out, f'mvn_interim_amortised_pps_net_{head_mode}.pkl')
    if os.path.exists(ckpt):
        print(f"  load cached {head_mode} net {ckpt}")
        fit = load_fitted_model(ckpt)
        return fit, float(fit.get('training_mins', 0.0))
    print(f"  train {head_mode} net ({NET_STEPS} steps, ragged={RAGGED}) ...")
    if RAGGED:
        fit = _train_ragged(head_mode, NET_STEPS)
    else:
        cell_tr = sim_data['cells'][TRAIN_J]
        model = MVNModel(dit=cell_tr['dit'], dcati=MVNModel.get_interim_data_x(cell_tr['dp']),
                         seed=seed, J=TRAIN_J, K_chol=cell_tr['R_chol'], sigma=simu_params['sigma'],
                         prior_tau=simu_params['prior_tau'], mu_0_baseline=MU0)
        sample_fn = partial(model.make_training_data_with_participant_tokens_random_K,
                            n_max=N_FULL, n_dist=N_DIST, queries_per_sample=NET_QUERIES,
                            K_families=K_FAMILIES)
        fit = train(sample_fn, Net,
                    net_kwargs=dict(head_mode=head_mode, n_max=float(N_FULL),
                                    embed_dim=NET_EMBED_DIM, hidden_dims=NET_HIDDEN),
                    pps_ProbH1_lwr_quantiles_mesh=TAUS, num_steps=NET_STEPS,
                    batch_size=NET_BATCH, lr=NET_LR, seed=seed, verbose=True)
    save_trained_model(fit, ckpt)
    return fit, float(fit.get('training_mins', 0.0))


summaries = []
qs_taus = np.asarray(TAUS)
for cfg in CONFIGS:
    head_mode, use_bvm, suf, label = CFG[cfg]
    print(f"\n{'#' * 72}\n# CONFIG {cfg}: head_mode={head_mode} bvm={use_bvm} suf={suf}\n{'#' * 72}")
    fit, mins_train = get_net(head_mode)
    t0 = time.time()
    QS, HD, TGT = build_forward(fit, head_mode)
    HK, PPSSRC = calibrate(fit, head_mode, use_bvm, QS, HD, TGT)
    mins_deploy = (time.time() - t0) / 60.0
    # calibrated arrays + PPS/p(H1) outputs
    pd.to_pickle({k: PPSSRC[k] for k in INTERIMS},
                 os.path.join(dir_out, f'mvn_J{J_DEPLOY}_{suf}_calibrated_qs.pkl'))
    mse = emit_pps(PPSSRC, suf, mins_deploy, mins_train)
    # diagnostics (§14.4.5): calibrated (deployed) + raw-network base
    adp.all_plots(lambda k: PPSSRC[k], lambda k: TGT[k], lambda k: NOBS[k],
                  labels, INTERIMS, qs_taus, dir_out, f'mvn_J{J_DEPLOY}', suf, None)
    adp.all_plots(lambda k: QS[k], lambda k: TGT[k], lambda k: NOBS[k],
                  labels, INTERIMS, qs_taus, dir_out, f'mvn_J{J_DEPLOY}', f'{suf}base', None)
    # scalar summary for §13.8
    ss = scalar_summary(PPSSRC, TGT, label, suf)
    summaries.append(ss)
    print(f"  [{label}] PIT-KS={ss.pit_ks.mean():.3f}  marg-KS={ss.marg_ks.mean():.3f}  "
          f"PPS-MSE(vs closed-form)={mse:.5f}  train={mins_train:.1f}m deploy={mins_deploy:.1f}m")

# %%

# =============================================================================
# Combined calibration summary table for §13.8.
# =============================================================================

allss = pd.concat(summaries, ignore_index=True)
allss.to_csv(os.path.join(dir_out, 'mvn_calibration_detail.csv'), index=False)


def _band(n):
    return 'small (n<=300)' if n <= 300 else ('mid (300<n<=700)' if n <= 700 else 'large (n>700)')


allss['band'] = allss['n'].map(_band)
tbl = (allss.groupby(['suf', 'config'], observed=True)
       .agg(pit_ks=('pit_ks', 'mean'), marg_ks=('marg_ks', 'mean')).reset_index())
band = (allss.groupby(['suf', 'band'], observed=True)['marg_ks'].mean().unstack().reset_index())
summary = tbl.merge(band, on='suf', how='left')
summary.to_csv(os.path.join(dir_out, 'mvn_calibration_summary.csv'), index=False)
print("\n===== §13.8 calibration summary =====")
print(summary.to_string(index=False))
print(f"\nAll configs done -> {dir_out}")
