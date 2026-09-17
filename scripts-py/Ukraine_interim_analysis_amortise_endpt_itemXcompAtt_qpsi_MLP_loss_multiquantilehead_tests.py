#!/usr/bin/env python3
"""
Calibration tests for the Ukraine `itemXcompAtt` amortiser (§14.2.1 of
``dev/amortised_decision_making.md``): does the predictive distribution
of rho | x, z^(s), approximated through the Q = 5 quantiles, match the
SVI posterior?

Since (z^(s), rho^(s)) are JOINT draws given x (theta^(s) ~ p(theta|x)
generates both), a calibrated multi-quantile predictive must satisfy:

1. **Quantile coverage.** For each level eta^inpol_q,
       Cov_q^(j*) = S^-1 sum_s 1{ rho_true^(j*,s) <= rho_hat_q^(j*,s) }
   should equal eta^inpol_q.
2. **PIT uniformity.** u^(j*,s) = F_hat(rho_true^(j*,s) | x, z^(s)),
   obtained by CDF interpolation of the 5 predicted quantiles, should
   be Uniform(0, 1) over s. Deviations are summarised by the
   Kolmogorov-Smirnov distance max_u |ecdf(u) - u|.

This deliberately does NOT compare per-draw medians point-by-point
(noisy, and only tests location); it tests the whole conditional
distribution.

Reuses the trained net + deployment-token construction of the main
script (loaded from its output dir; run the main script first).

Outputs (in the main script's output dir):

  - ``pcm_1_interim_pps_RGEG_tests_coverage.csv``
  - ``pcm_1_interim_pps_RGEG_tests_coverage.pdf``
  - ``pcm_1_interim_pps_RGEG_tests_pit.pdf``
  - ``pcm_1_interim_pps_RGEG_tests_pit_ks.csv``

Usage:
    cd /Users/or105/git/bIRTistic
    PYTHONUNBUFFERED=1 pixi run python -u \
      scripts-py/Ukraine_interim_analysis_amortise_endptx_on_wz_with_features_itemXcompAtt_qpsi_MLP_loss_multiquantilehead_tests.py
"""

# %%

import os
import sys
import warnings
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
from plotnine import (
    ggplot, aes, geom_abline, geom_histogram, geom_hline, geom_line,
    geom_point,
    facet_wrap,
    theme_bw, theme, element_blank, element_text, labs,
)

warnings.filterwarnings('ignore')

from amortiser_common import load_fitted_model, _jax_pytree
from model_pcm import PartialCreditModel

print("Imports successful")

# %%

# =============================================================================
# Configuration (mirrors the main script).
# =============================================================================

DIR_RGE = os.path.join(
    "/Users/or105/sandbox/bIRTistic",
    "py-ukraine-interim-with-regression-on-endptx-wz-260601",
)
dir_out = os.environ.get('UKR_DIR', os.path.join(
    "/Users/or105/sandbox/bIRTistic",
    "py-ukraine-interim-amortise-itemXcompAtt-net-260716",
))
file_prefix = "pcm_1_interim"
SUFFIX = "RGEG"

N_FULL = 503
KHAT_SHRINK_N0 = 50.0
F_TOK = 8
TARGET_CLIP = 20.0
NET_TAUS = np.array([0.05, 0.25, 0.5, 0.75, 0.95])   # eta^inpol_q
INTERIM_IDS = list(range(1, 9))
PPS_Z_TOTAL = int(os.environ.get('UKR_TESTS_S', 4000))

# %%

# =============================================================================
# Item metadata + trained net + item_std (from the main script's run).
# =============================================================================

interim_wa = {}
for i in INTERIM_IDS:
    p = os.path.join(DIR_RGE, f"{file_prefix}_i{i}_regression_training.pkl")
    if os.path.exists(p):
        interim_wa[i] = pd.read_pickle(p)
INTERIM_IDS = sorted(interim_wa.keys())

items_df = (
    interim_wa[INTERIM_IDS[0]][['item_label', 'item_type', 'item_high_label']]
    .drop_duplicates()
    .sort_values(['item_type', 'item_label'])
    .reset_index(drop=True)
)
items_df['item_pos'] = np.arange(len(items_df), dtype=int)
J = len(items_df)
item_pos = dict(zip(items_df['item_label'], items_df['item_pos']))
item_type_idx = items_df['item_type'].map(
    {'out-of-7': 0.0, 'categorical': 1.0}
).to_numpy(dtype=np.float32)
item_high_idx = items_df['item_high_label'].map(
    {'higher_is_better': 1.0, 'lower_is_better': 0.0}
).to_numpy(dtype=np.float32)

dit_df = pd.read_csv(os.path.join(DIR_RGE, f"{file_prefix}_1_data_dit.csv"))
items_df = items_df.merge(
    dit_df[['item_label', 'cat_length']].drop_duplicates(),
    on='item_label', how='left',
)
item_klevels = items_df['cat_length'].to_numpy(dtype=np.float32)

fit = load_fitted_model(
    os.path.join(dir_out, f"{file_prefix}_amortised_pps_net.pkl"),
)
item_std = np.load(
    os.path.join(dir_out, f"{file_prefix}_item_std.npy"),
).astype(np.float32)
print(f"Loaded net + item_std from {dir_out}")

# %%

# =============================================================================
# Deployment-token construction (identical to the main script).
# =============================================================================


def _interim_x_features(interim_id):
    xi = pd.read_csv(
        os.path.join(DIR_RGE, f"{file_prefix}_{interim_id}_data_dp1.csv"),
    )
    piv = xi.pivot_table(index='pid', columns=['item_label', 'time'],
                         values='y')
    n_pid = piv.index.size
    w_x = np.zeros((J, 2), dtype=np.float32)
    chg = np.zeros((n_pid, J), dtype=np.float64)
    for j, lbl in enumerate(items_df['item_label']):
        y0 = piv[(lbl, 0)].to_numpy(dtype=np.float64)
        y1 = piv[(lbl, 1)].to_numpy(dtype=np.float64)
        kmax = item_klevels[j] - 1.0
        w_x[j, 0] = np.nanmean(y0) / kmax
        w_x[j, 1] = np.nanmean(y1) / kmax
        chg[:, j] = y1 - y0
    rho = pd.DataFrame(chg).corr(method='spearman').to_numpy()
    rho = np.nan_to_num(rho, nan=0.0)
    np.fill_diagonal(rho, 1.0)
    lam = n_pid / (n_pid + KHAT_SHRINK_N0)
    khat = (lam * rho + (1.0 - lam) * np.eye(J)).astype(np.float32)
    return w_x, khat, n_pid


def _interim_z_means(interim_id, D):
    xi = pd.read_csv(
        os.path.join(DIR_RGE, f"{file_prefix}_{interim_id}_data_dp1.csv"),
    )
    n_pid = xi['pid'].nunique()
    interim_m = N_FULL - n_pid
    model = PartialCreditModel(dit=dit_df, dcati=xi, x_formula="~ time - 1",
                               seed=123)
    zi = model.get_interim_z_from_ypredi(
        os.path.join(DIR_RGE, f"{file_prefix}_{interim_id}_draws.zarr"),
        interim_m, pps_z_total=D, seed=123, keep_order=True,
    )
    ypred_cols = [f'ypred_{s}' for s in range(D)]
    out = np.zeros((D, J, 2), dtype=np.float32)
    for j, lbl in enumerate(items_df['item_label']):
        kmax = item_klevels[j] - 1.0
        for tt in (0, 1):
            sub = zi[(zi['item_label'] == lbl) & (zi['time'] == tt)]
            vals = sub[ypred_cols].to_numpy(dtype=np.float32)
            out[:, j, tt] = (vals - 1.0).mean(axis=0) / kmax
    return out


def _tensorise_interim(wa, interim_id):
    w_x, khat, n_pid = _interim_x_features(interim_id)
    wa = wa.copy()
    wa['item_pos'] = wa['item_label'].map(item_pos)
    piv = wa.pivot_table(index='draw', columns='item_pos',
                         values='pps_ratio_x')
    D = min(piv.index.size, PPS_Z_TOTAL)
    target = piv.reindex(columns=range(J)).to_numpy(dtype=np.float32)[:D]
    target = np.nan_to_num(target, nan=0.0, posinf=TARGET_CLIP,
                           neginf=-TARGET_CLIP)
    target = np.clip(target, -TARGET_CLIP, TARGET_CLIP)
    zm = _interim_z_means(interim_id, D)
    wb = np.maximum(zm[:, :, 0], 1e-3)
    we = zm[:, :, 1]
    high = (item_high_idx > 0.5)[None, :]
    w_rat = np.clip(np.where(high, we / wb - 1.0, 1.0 - we / wb),
                    -TARGET_CLIP, TARGET_CLIP)
    w_x_b = np.broadcast_to(w_x[None, :, 0], (D, J))
    w_x_e = np.broadcast_to(w_x[None, :, 1], (D, J))
    feats = np.stack([w_x_b, w_x_e, zm[:, :, 0], zm[:, :, 1], w_rat],
                     axis=-1)
    feats = np.nan_to_num(feats, nan=0.0, posinf=0.0, neginf=0.0)
    metadata_b = np.broadcast_to(
        np.stack([item_type_idx, item_high_idx], axis=-1)[None], (D, J, 2),
    ).astype(np.float32)
    tokens = np.concatenate([feats, metadata_b], axis=-1)
    aux_row = np.array([n_pid / N_FULL, (N_FULL - n_pid) / N_FULL],
                       dtype=np.float32)
    return tokens.astype(np.float32), target, khat, aux_row


def _fanout_with_khat(tok, khat, aux_row):
    D = tok.shape[0]
    B = D * J
    tokens_rep = np.empty((D, J, J, F_TOK), dtype=np.float32)
    tokens_rep[..., :7] = tok[:, None, :, :]
    tokens_rep[..., 7] = np.broadcast_to(khat[None, :, :], (D, J, J))
    tokens_rep = tokens_rep.reshape(B, J, F_TOK)
    mask_rep = np.ones((B, J), dtype=np.float32)
    query_idx = np.broadcast_to(
        np.arange(J)[None, :], (D, J),
    ).reshape(B).astype(np.int32)
    aux_rep = np.broadcast_to(aux_row[None, :], (B, 2)).astype(np.float32)
    return tokens_rep, mask_rep, query_idx, aux_rep

# %%

# =============================================================================
# Per-interim: full quantile grid + coverage + PIT.
# =============================================================================

cov_rows = []
pit_rows = []
marg_ks_rows = []
marg_cdf_rows = []
for interim_id in INTERIM_IDS:
    tokens, target, khat_i, aux_i = _tensorise_interim(
        interim_wa[interim_id], interim_id,
    )
    D = tokens.shape[0]
    tokens_rep, mask_rep, query_idx, aux_rep = _fanout_with_khat(
        tokens, khat_i, aux_i,
    )
    batch = {'tokens': tokens_rep, 'mask': mask_rep,
             'query_idx': query_idx, 'aux': aux_rep}
    preds = np.asarray(fit['apply_fn'](fit['params'], _jax_pytree(batch)))
    preds = preds * item_std[query_idx][:, None]            # unstandardise
    preds = np.maximum.accumulate(preds, axis=1)            # monotone
    preds = preds.reshape(D, J, len(NET_TAUS))              # (D, J, Q)

    for j in range(J):
        lbl = items_df['item_label'].iloc[j]
        y = target[:, j]                                    # (D,) SVI rho
        qs = preds[:, j, :]                                 # (D, Q)
        # 1. conditional-test coverage per quantile level (§14.1.7)
        for qi, eta in enumerate(NET_TAUS):
            cov_rows.append({
                'interim_id': interim_id, 'item_label': lbl,
                'eta_inpol': float(eta),
                'coverage': float((y <= qs[:, qi]).mean()),
            })
        # 2. conditional-test PIT via CDF interpolation (§14.1.7)
        u = np.array([
            np.interp(y[s], qs[s], NET_TAUS, left=0.0, right=1.0)
            for s in range(D)
        ])
        uu = np.sort(u)
        ecdf = np.arange(1, D + 1) / D
        ks = float(np.max(np.abs(ecdf - uu)))
        pit_rows.append({
            'interim_id': interim_id, 'item_label': lbl, 'ks': ks,
            'pit': u,
        })
        # 3. MARGINAL diagnostic (§14.1.8): does the amortiser's
        #    p_hat(rho|x) = E_{z|x}[p_hat(rho|x,z)] -- the item-variant
        #    analogue of "empty z" -- match the SVI p(rho|x) = {rho_SVI-x}?
        #    Amortiser marginal CDF on a grid = mean over s of the per-s
        #    piecewise-linear quantile CDFs; SVI CDF = ecdf of the targets.
        lo = float(min(y.min(), qs.min()))
        hi = float(max(y.max(), qs.max()))
        r_grid = np.linspace(lo, hi, 200)
        F_marg = np.zeros_like(r_grid)
        for s in range(D):
            F_marg += np.interp(r_grid, qs[s], NET_TAUS, left=0.0, right=1.0)
        F_marg /= D                                         # amortiser p_hat(rho|x)
        F_svi = (y[:, None] <= r_grid[None, :]).mean(axis=0)  # SVI p(rho|x)
        ks_marg = float(np.max(np.abs(F_marg - F_svi)))
        marg_ks_rows.append({'interim_id': interim_id, 'item_label': lbl,
                             'ks_marginal': ks_marg})
        # store a subsampled curve pair for plotting
        idx = np.linspace(0, len(r_grid) - 1, 60).astype(int)
        for source, F in (('SVI  p(rho|x)', F_svi),
                          ('amortiser  p_hat(rho|x)', F_marg)):
            for r_v, c_v in zip(r_grid[idx], F[idx]):
                marg_cdf_rows.append({
                    'interim_id': interim_id, 'item_label': lbl,
                    'source': source, 'rho': float(r_v), 'cdf': float(c_v),
                })
    print(f"  interim {interim_id}: coverage + PIT + marginal done (D={D})")

# %%

# =============================================================================
# Save + plots.
# =============================================================================

cov = pd.DataFrame(cov_rows)
cov_path = os.path.join(dir_out, f"{file_prefix}_pps_{SUFFIX}_tests_coverage.csv")
cov.to_csv(cov_path, index=False)
print(f"Saved coverage CSV: {cov_path}")

ks_df = pd.DataFrame([{k: r[k] for k in ('interim_id', 'item_label', 'ks')}
                      for r in pit_rows])
ks_path = os.path.join(dir_out, f"{file_prefix}_pps_{SUFFIX}_tests_pit_ks.csv")
ks_df.to_csv(ks_path, index=False)
print(f"Saved PIT-KS CSV: {ks_path}")
print("\nKS by item (mean over interims):")
print(ks_df.groupby('item_label')['ks'].mean().round(3).to_string())

# Calibration curves: empirical coverage vs nominal, one line per interim.
p = (
    ggplot(cov, aes(x='eta_inpol', y='coverage',
                    colour='factor(interim_id)', group='interim_id'))
    + geom_abline(intercept=0, slope=1, linetype='dashed', colour='black')
    + geom_line(size=0.4, alpha=0.8)
    + geom_point(size=0.8)
    + facet_wrap('~ item_label', ncol=4)
    + theme_bw()
    + theme(figure_size=(12, 12.5), legend_position='top',
            strip_background=element_blank(),
            strip_text=element_text(face='bold'))
    + labs(x='nominal quantile level eta^inpol_q',
           y='empirical coverage  P(rho_true <= rho_hat_q)',
           colour='interim',
           title='Quantile-coverage calibration of q_psi(rho | x, z^s)')
)
cov_pdf = os.path.join(dir_out, f"{file_prefix}_pps_{SUFFIX}_tests_coverage.pdf")
p.save(cov_pdf, verbose=False, limitsize=False)
print(f"Saved coverage plot: {cov_pdf}")

# PIT histograms pooled over interims, per item.
pit_long = pd.concat([
    pd.DataFrame({'item_label': r['item_label'], 'u': r['pit']})
    for r in pit_rows
], ignore_index=True)
p = (
    ggplot(pit_long, aes(x='u'))
    + geom_histogram(aes(y='..density..'), bins=20,
                     fill='#1f77b4', colour='white', size=0.2)
    + geom_hline(yintercept=1.0, linetype='dashed', colour='black')
    + facet_wrap('~ item_label', ncol=4)
    + theme_bw()
    + theme(figure_size=(12, 12.5),
            strip_background=element_blank(),
            strip_text=element_text(face='bold'))
    + labs(x='PIT  u = F_hat(rho_true | x, z^s)', y='density',
           title='PIT uniformity of q_psi(rho | x, z^s) '
                 '(pooled over interims; dashed = Uniform(0,1))')
)
pit_pdf = os.path.join(dir_out, f"{file_prefix}_pps_{SUFFIX}_tests_pit.pdf")
p.save(pit_pdf, verbose=False, limitsize=False)
print(f"Saved PIT plot: {pit_pdf}")

# --- Marginal diagnostic (§14.1.8): SVI p(rho|x) vs amortiser p_hat(rho|x) ---
marg_ks = pd.DataFrame(marg_ks_rows)
marg_ks_path = os.path.join(
    dir_out, f"{file_prefix}_pps_{SUFFIX}_tests_marginal_ks.csv")
marg_ks.to_csv(marg_ks_path, index=False)
print(f"Saved marginal-KS CSV: {marg_ks_path}")
print("\nMarginal KS (SVI vs amortiser p(rho|x)) by item (mean over interims):")
print(marg_ks.groupby('item_label')['ks_marginal'].mean().round(3).to_string())
print("Marginal KS by interim (mean over items):")
print(marg_ks.groupby('interim_id')['ks_marginal'].mean().round(3).to_string())

marg_cdf = pd.DataFrame(marg_cdf_rows)
# One CDF-overlay plot per interim (facet per item, two lines).
for interim_id in INTERIM_IDS:
    sub = marg_cdf[marg_cdf['interim_id'] == interim_id]
    if sub.empty:
        continue
    p = (
        ggplot(sub, aes(x='rho', y='cdf', colour='source', group='source'))
        + geom_line(size=0.6)
        + facet_wrap('~ item_label', ncol=4, scales='free_x')
        + theme_bw()
        + theme(figure_size=(15, 14), legend_position='top',
                strip_background=element_blank(),
                strip_text=element_text(face='bold'))
        + labs(x=f'rho at interim {interim_id}', y='CDF',
               colour='',
               title=f'Marginal p(rho | x): SVI vs amortiser at interim '
                     f'{interim_id}')
    )
    marg_pdf = os.path.join(
        dir_out,
        f"{file_prefix}_i{interim_id}_svi_vs_amortiser_marginal_cdf.pdf")
    p.save(marg_pdf, verbose=False, limitsize=False)
    print(f"  saved marginal-CDF plot: {marg_pdf}")

print("\nCalibration tests complete.")
