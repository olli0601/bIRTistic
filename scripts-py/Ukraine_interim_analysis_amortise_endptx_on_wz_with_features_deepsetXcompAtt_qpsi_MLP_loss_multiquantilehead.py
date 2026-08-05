#!/usr/bin/env python3
"""
Ukraine interim analysis: **fully amortised** PPS per §8-10 of
``dev/amortised_decision_making.md`` using the nested-DeepSets +
cross-attention amortiser (`deepsetXcompAtt`, §13.7) on the partial
credit model.

Unlike the §14.1 item-level variant (regression-style: tokens built
from ``z^(s)`` summaries only), this script implements the §10
deployment faithfully: the network consumes the RAW observed cohort
``x_obs`` (n participants x J items x 2 time-points) AND the raw
posterior-predictive future cohort ``z^(s)`` (m shadow participants x
J x 2) per posterior draw ``s``. Cohort axis (x vs z) and time axis
(baseline vs endline; R = 2 response features) are kept separate by
construction.

Training data are prior-predictive draws from the PCM itself (§8):
``PartialCreditModel.make_training_data_with_participant_tokens_prior``
draws theta from the model prior, simulates both cohorts jointly and
labels each item with the population endpoint ratio rho_j(theta).
No SVI fits are consumed at training time.

Deployment (§10) reuses the cached per-interim artifacts of the RGE
script (xi csv, dit csv, draws.zarr, wa pkl for evaluation targets).

Outputs (``_RGDX_`` suffix):

  - ``pcm_1_interim_pps_RGDX_p_h1_xz.pkl``
  - ``pcm_1_interim_pps_RGDX.csv``
  - ``pcm_1_interim_pps_RGDX_timing.csv``
  - ``pcm_1_interim_pps_RGDX_perf.csv``
  - ``pcm_1_interim_pps_RGDX_svi_vs_amortiser_rho_all.pdf``

Usage:
    cd /Users/or105/git/bIRTistic
    PYTHONUNBUFFERED=1 pixi run python -u \
      scripts-py/Ukraine_interim_analysis_amortise_endptx_on_wz_with_features_deepsetXcompAtt_qpsi_MLP_loss_multiquantilehead.py
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
from plotnine import (
    ggplot, aes, geom_boxplot,
    facet_wrap, position_dodge,
    scale_fill_manual,
    theme_bw, theme, element_blank, element_text, labs,
)

warnings.filterwarnings('ignore')

from amortiser_common import (
    train as amortiser_train,
    save_trained_model,
    load_fitted_model,
    predict_amortised_p_h1_for_one_xz,
)
from amortiser_pps_features_deepsetXcompAtt_qpsi_MLP_loss_multiquantilehead import (
    Amortiser_PPS_features_deepsetXcompAtt_qpsi_MLP_loss_multiquantilehead,
)
from model_pcm import PartialCreditModel

print("Imports successful")

# %%

# =============================================================================
# Configuration
# =============================================================================

DIR_RGE = os.path.join(
    "/Users/or105/sandbox/bIRTistic",
    "py-ukraine-interim-with-regression-on-endptx-wz-260601",
)
dir_out = os.path.join(
    "/Users/or105/sandbox/bIRTistic",
    "py-ukraine-interim-amortise-endptx-on-wz-with-features-"
    "deepsetXcompAtt-qpsi-MLP-loss-multiquantilehead-260805",
)
os.makedirs(dir_out, exist_ok=True)
print(f"Output dir: {dir_out}")

file_prefix = "pcm_1_interim"
SUFFIX = "RGDX"

pps_H1_min_effect_size_thresh = 0.5
pps_ProbH1_target_lwr_quantile = 0.89
categorical_threshold = 2        # on 1-indexed categories (matches RGE get_w)
seed = 123
x_formula = "~ time - 1"

N_FULL = 503                     # total trial participants (n + m per interim)
# Response-level counts are read from dit['cat_length'] below (Ukraine:
# 8 for out-of-7, 4 for categorical) -- study-specific, not hardcoded.
PPS_Z_TOTAL = int(os.environ.get('UKR_DEEPSET_S', 200))

NET_STEPS = int(os.environ.get('UKR_DEEPSET_STEPS', 8000))
NET_BATCH = 32                   # S=32 prior samples * Q=4 queries -> B=128
NET_QUERIES = 4
NET_LR = 1e-3
NET_TAUS = (0.05, 0.25, 0.5, 0.75, 0.95)
NET_SEED = seed
TARGET_CLIP = 20.0

INTERIM_IDS = list(range(1, 9))

# %%

# =============================================================================
# Item metadata: stable ordering shared by sampler + deployment pivots.
# =============================================================================

dit = pd.read_csv(os.path.join(DIR_RGE, f"{file_prefix}_1_data_dit.csv"))
# Items actually fitted (20 of the 26 dit rows) -- derive from the wa
# frame like the §14.1 item scripts so ordering matches downstream.
_wa1 = pd.read_pickle(
    os.path.join(DIR_RGE, f"{file_prefix}_i1_regression_training.pkl"),
)
items_df = (
    _wa1[['item_label', 'item_type', 'item_high_label']]
    .drop_duplicates()
    .sort_values(['item_type', 'item_label'])
    .reset_index(drop=True)
)
J = len(items_df)
item_labels = items_df['item_label'].tolist()
item_type_idx = items_df['item_type'].map(
    {'out-of-7': 0.0, 'categorical': 1.0}
).to_numpy(dtype=np.float32)
item_high_idx = items_df['item_high_label'].map(
    {'higher_is_better': 1.0, 'lower_is_better': 0.0}
).to_numpy(dtype=np.float32)
items_df = items_df.merge(
    dit[['item_label', 'cat_length']].drop_duplicates(),
    on='item_label', how='left',
)
item_kmax = items_df['cat_length'].to_numpy(dtype=np.float32)
assert not np.isnan(item_kmax).any(), "cat_length missing for some items"
K_O7 = int(items_df.loc[items_df['item_type'] == 'out-of-7',
                        'cat_length'].iloc[0])
K_CAT = int(items_df.loc[items_df['item_type'] == 'categorical',
                         'cat_length'].iloc[0])
print(f"J = {J} items:\n{items_df}")


def _pivot_responses(df, value_col):
    """Long-form (pid, item_label, time, value) -> (n_pid, J, 2) float32,
    rescaled to [0, 1] via (y_stan - 1) / (K - 1) per item type."""
    piv = df.pivot_table(index='pid', columns=['item_label', 'time'],
                         values=value_col)
    pids = piv.index.to_numpy()
    out = np.zeros((len(pids), J, 2), dtype=np.float32)
    for j, lbl in enumerate(item_labels):
        for t in (0, 1):
            out[:, j, t] = piv[(lbl, t)].to_numpy(dtype=np.float32)
    out = (out - 1.0) / (item_kmax[None, :, None] - 1.0)
    return np.nan_to_num(out, nan=0.0), pids

# %%

# =============================================================================
# Train the amortised net on prior-predictive draws (§8 + §9), or load.
# =============================================================================

NET_CKPT = os.path.join(dir_out, f"{file_prefix}_amortised_pps_net.pkl")

if os.path.exists(NET_CKPT):
    print(f"\nLoading cached amortised net from {NET_CKPT}")
    fit = load_fitted_model(NET_CKPT)
else:
    print(f"\nTraining deepsetXcompAtt amortised net on PCM prior draws "
          f"({NET_STEPS} steps, S={NET_BATCH}, Q={NET_QUERIES}, "
          f"N_full={N_FULL}) ...")
    sample_fn = partial(
        PartialCreditModel.make_training_data_with_participant_tokens_prior,
        item_type_idx=item_type_idx,
        item_high_idx=item_high_idx,
        n_max=N_FULL,
        K_out_of_7=K_O7,
        K_categorical=K_CAT,
        categorical_threshold=categorical_threshold,
        queries_per_sample=NET_QUERIES,
        target_clip=TARGET_CLIP,
    )
    fit = amortiser_train(
        sample_fn,
        Amortiser_PPS_features_deepsetXcompAtt_qpsi_MLP_loss_multiquantilehead,
        net_kwargs={},
        pps_ProbH1_lwr_quantiles_mesh=NET_TAUS,
        num_steps=NET_STEPS,
        batch_size=NET_BATCH,
        lr=NET_LR,
        seed=NET_SEED,
        verbose=True,
    )
    save_trained_model(fit, NET_CKPT)
    print(f"Saved amortised net to {NET_CKPT}")
mins_train = float(fit.get('training_mins', float('nan')))
print(f"Amortised net trained in {mins_train:.2f} min.")

# %%

# =============================================================================
# Per-interim deployment (§10). For each interim: raw x_obs from the
# cached xi csv, raw z^(s) rebuilt from the cached draws.zarr via
# get_interim_z_from_ypredi, forward-pass (x_obs, z^(s)) per draw.
# =============================================================================

p_h1_xz_rows = []
pps_timing_rows = []
t_all0 = time.time()
for interim_id in INTERIM_IDS:
    t0 = time.time()
    xi_path = os.path.join(DIR_RGE, f"{file_prefix}_{interim_id}_data_dp1.csv")
    zarr_path = os.path.join(DIR_RGE, f"{file_prefix}_{interim_id}_draws.zarr")
    wa_path = os.path.join(
        DIR_RGE, f"{file_prefix}_i{interim_id}_regression_training.pkl",
    )
    if not (os.path.exists(xi_path) and os.path.exists(zarr_path)):
        print(f"[interim {interim_id}] skipped: missing cached artifacts")
        continue
    xi = pd.read_csv(xi_path)
    x_resp_raw, x_pids = _pivot_responses(xi, 'y_stan')      # (n, J, 2)
    n_obs = x_resp_raw.shape[0]
    interim_m = N_FULL - n_obs
    print(f"\n{'=' * 70}\n{SUFFIX} interim {interim_id}: n={n_obs} "
          f"m={interim_m} S={PPS_Z_TOTAL}\n{'=' * 70}")

    # Rebuild the posterior-predictive future block z^(s) from the cached
    # SVI draws (keep_order=True so draw s aligns with the wa targets).
    model = PartialCreditModel(dit=dit, dcati=xi, x_formula=x_formula,
                               seed=seed)
    zi = model.get_interim_z_from_ypredi(
        zarr_path, interim_m, pps_z_total=PPS_Z_TOTAL, seed=seed,
        keep_order=True,
    )

    # x cohort tensors: shared across draws + queries.
    x_pad = np.zeros((N_FULL, J, 2), dtype=np.float32)
    x_pad[:n_obs] = x_resp_raw
    mask_x_row = np.zeros((N_FULL,), dtype=np.float32)
    mask_x_row[:n_obs] = 1.0
    mask_z_row = np.zeros((N_FULL,), dtype=np.float32)
    mask_z_row[:interim_m] = 1.0

    x_resp_b = np.broadcast_to(
        x_pad[None], (J, N_FULL, J, 2),
    ).astype(np.float32)
    mask_x_b = np.broadcast_to(mask_x_row[None], (J, N_FULL)).astype(np.float32)
    mask_z_b = np.broadcast_to(mask_z_row[None], (J, N_FULL)).astype(np.float32)
    meta_b = np.broadcast_to(
        np.stack([item_type_idx, item_high_idx], axis=-1)[None], (J, J, 2),
    ).astype(np.float32)
    aux_b = np.broadcast_to(
        np.array([n_obs / N_FULL, interim_m / N_FULL],
                 dtype=np.float32)[None], (J, 2),
    ).astype(np.float32)
    qidx_b = np.arange(J, dtype=np.int32)

    preds_p = np.empty((PPS_Z_TOTAL, J), dtype=np.float64)
    preds_q = np.empty((PPS_Z_TOTAL, J), dtype=np.float64)
    med_idx = len(NET_TAUS) // 2
    for s in range(PPS_Z_TOTAL):
        z_resp_raw, _ = _pivot_responses(
            zi.rename(columns={f'ypred_{s}': '_y'}), '_y',
        )                                                    # (m, J, 2)
        z_pad = np.zeros((N_FULL, J, 2), dtype=np.float32)
        z_pad[:interim_m] = z_resp_raw
        z_resp_b = np.broadcast_to(
            z_pad[None], (J, N_FULL, J, 2),
        ).astype(np.float32)
        batch = {
            'x_responses':   x_resp_b,
            'mask_x':        mask_x_b,
            'z_responses':   z_resp_b,
            'mask_z':        mask_z_b,
            'item_metadata': meta_b,
            'query_idx':     qidx_b,
            'aux':           aux_b,
        }
        p_h1_xz_s, _q_s, preds_s = predict_amortised_p_h1_for_one_xz(
            fit, batch, pps_H1_min_effect_size_thresh,
        )
        preds_p[s] = p_h1_xz_s
        preds_q[s] = preds_s[:, med_idx]

    # Evaluation targets from the cached wa frame (draw-aligned).
    target_per_pos = np.full((PPS_Z_TOTAL, J), np.nan)
    if os.path.exists(wa_path):
        wa = pd.read_pickle(wa_path)
        wa_piv = wa.pivot_table(index='draw', columns='item_label',
                                values='pps_ratio_x')
        wa_piv = wa_piv.reindex(columns=item_labels)
        common = wa_piv.index[wa_piv.index < PPS_Z_TOTAL]
        target_per_pos[common.to_numpy()] = np.clip(
            wa_piv.loc[common].to_numpy(dtype=np.float64),
            -TARGET_CLIP, TARGET_CLIP,
        )

    rows = []
    for j in range(J):
        rows.append(pd.DataFrame({
            'item_label':      [items_df.iloc[j]['item_label']] * PPS_Z_TOTAL,
            'item_type':       [items_df.iloc[j]['item_type']] * PPS_Z_TOTAL,
            'item_high_label': [items_df.iloc[j]['item_high_label']] * PPS_Z_TOTAL,
            's':               np.arange(1, PPS_Z_TOTAL + 1, dtype=int),
            'p_h1_xz':         preds_p[:, j],
            'q_hat':           preds_q[:, j],
            'pps_ratio_x':     target_per_pos[:, j],
        }))
    df = pd.concat(rows, ignore_index=True)
    df['interim_id'] = interim_id
    p_h1_xz_rows.append(df)

    mins_interim = (time.time() - t0) / 60.0
    pps_timing_rows.append({
        'interim_id':      interim_id,
        'mins_interim_id': round(mins_interim, 3),
    })
    print(f"  interim {interim_id} {SUFFIX} amortised done in "
          f"{mins_interim:.3f} min")

# %%

# =============================================================================
# Aggregate: p_h1_xz pkl, PPS table, perf (rho / R^2 vs SVI target),
# timing, combined boxplot.
# =============================================================================

dp_h1_xz = pd.concat(p_h1_xz_rows, ignore_index=True)
dp_h1_xz['pps_H1_min_effect_size_thresh'] = pps_H1_min_effect_size_thresh
dp_h1_xz['pps_ProbH1_target_lwr_quantile'] = pps_ProbH1_target_lwr_quantile
dp_h1_xz['S'] = PPS_Z_TOTAL
pkl_path = os.path.join(dir_out, f"{file_prefix}_pps_{SUFFIX}_p_h1_xz.pkl")
dp_h1_xz.to_pickle(pkl_path)
print(f"\nSaved {SUFFIX} p(H_1 | x, z) samples to: {pkl_path}")

pps = (
    dp_h1_xz
    .assign(h1=(dp_h1_xz['p_h1_xz'] > pps_ProbH1_target_lwr_quantile)
            .astype(float))
    .groupby(['interim_id', 'item_label', 'item_type', 'item_high_label'])
    .agg(pps=('h1', 'mean'))
    .reset_index()
)
pps['eta'] = pps_ProbH1_target_lwr_quantile
pps['S'] = PPS_Z_TOTAL
csv_path = os.path.join(dir_out, f"{file_prefix}_pps_{SUFFIX}.csv")
pps.to_csv(csv_path, index=False)
print(f"Saved {SUFFIX} PPS table to: {csv_path}")

pps_timing = pd.DataFrame(pps_timing_rows)
pps_timing['mins_total'] = round((time.time() - t_all0) / 60.0, 3)
timing_path = os.path.join(dir_out, f"{file_prefix}_pps_{SUFFIX}_timing.csv")
pps_timing.to_csv(timing_path, index=False)
print(f"Saved {SUFFIX} timing to: {timing_path}")

perf_rows = []
for (interim_id, item), g in dp_h1_xz.groupby(['interim_id', 'item_label']):
    m_ok = np.isfinite(g['pps_ratio_x']) & np.isfinite(g['q_hat'])
    if m_ok.sum() >= 2:
        y = g.loc[m_ok, 'pps_ratio_x'].to_numpy()
        yh = g.loc[m_ok, 'q_hat'].to_numpy()
        rho_v = float(np.corrcoef(yh, y)[0, 1])
        ss_res = float(np.sum((y - yh) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        r2_v = 1.0 - ss_res / ss_tot if ss_tot > 0 else float('nan')
    else:
        rho_v, r2_v = float('nan'), float('nan')
    perf_rows.append({'interim_id': interim_id, 'item_label': item,
                      'rho': rho_v, 'r2': r2_v})
perf = pd.DataFrame(perf_rows).merge(
    items_df[['item_label', 'item_type', 'item_high_label']],
    on='item_label', how='left',
)
perf_path = os.path.join(dir_out, f"{file_prefix}_pps_{SUFFIX}_perf.csv")
perf.to_csv(perf_path, index=False)
print(f"Saved {SUFFIX} perf CSV to: {perf_path}")
print(perf[perf['interim_id'] == 1][
    ['item_label', 'item_type', 'rho', 'r2']].to_string(index=False))

# Combined boxplot: SVI rho|x vs amortiser q_psi(rho|x) per (interim, item).
_source_svi = 'SVI  p(rho|x)'
_source_amo = 'amortiser  q_psi(rho|x)'
svi = dp_h1_xz[['interim_id', 'item_label', 'pps_ratio_x']].rename(
    columns={'pps_ratio_x': 'rho'})
svi['source'] = _source_svi
amo = dp_h1_xz[['interim_id', 'item_label', 'q_hat']].rename(
    columns={'q_hat': 'rho'})
amo['source'] = _source_amo
box_all = pd.concat([svi, amo], ignore_index=True).dropna(subset=['rho'])
box_all['source'] = pd.Categorical(
    box_all['source'], categories=[_source_svi, _source_amo], ordered=True,
)
box_stats = (
    box_all.groupby(['interim_id', 'item_label', 'source'], observed=True)['rho']
    .quantile([0.025, 0.25, 0.5, 0.75, 0.975])
    .unstack().reset_index()
    .rename(columns={0.025: 'q025', 0.25: 'q25', 0.5: 'q50',
                     0.75: 'q75', 0.975: 'q975'})
)
box_stats['grp'] = (
    box_stats['interim_id'].astype(str) + '|' + box_stats['source'].astype(str)
)
_n_col = 4
_n_row = int(np.ceil(J / _n_col))
p = (
    ggplot(box_stats, aes(x='factor(interim_id)',
                          ymin='q025', lower='q25', middle='q50',
                          upper='q75', ymax='q975',
                          fill='source', group='grp'))
    + geom_boxplot(stat='identity',
                   position=position_dodge(width=0.8), width=0.7,
                   colour='#404040', size=0.3)
    + scale_fill_manual(values={_source_svi: '#1f77b4',
                                 _source_amo: '#d62728'})
    + facet_wrap('~ item_label', ncol=_n_col, scales='free_y')
    + theme_bw()
    + theme(axis_text_x=element_text(angle=45, vjust=1, hjust=1),
            figure_size=(3.0 * _n_col, 2.5 * _n_row),
            legend_position='top',
            strip_background=element_blank(),
            strip_text=element_text(face='bold'))
    + labs(x='Interim', y='rho  (box 25-75%, whiskers 2.5-97.5%)',
           title='SVI posterior of rho|x vs deepsetXcompAtt q_psi(rho|x)')
)
pdf_path = os.path.join(
    dir_out, f"{file_prefix}_pps_{SUFFIX}_svi_vs_amortiser_rho_all.pdf",
)
p.save(pdf_path, verbose=False, limitsize=False)
print(f"Saved combined boxplot: {pdf_path}")

print(f"\nAll interims done -- Ukraine {SUFFIX} deepset amortised PPS "
      f"complete in {(time.time() - t_all0) / 60.0:.1f} min.")
