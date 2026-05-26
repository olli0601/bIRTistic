#!/usr/bin/env python3
"""
Ukraine interim analyses, Python implementation.

Mirrors ``scripts-py/Colombia_interim_analyses.py`` but on the Ukraine data:
- ncats partial credit model (``src/numpyro/partial_credit_model_ncats_v260413.pyro``);
- NumPyro SVI with ``AutoLowRankMultivariateNormal``;
- ``get_endpoints`` for primary endpoint summaries (``categorical_threshold = 2`` for Ukraine);
- cross-interim plots in plotnine.

The monthly interim grid is auto-detected from the endline submission dates in
the data (first to last endline month).

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/Ukraine_interim_analyses.py
"""

# %%

import os
import sys
import time
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
    ggplot, aes, geom_col, geom_line, geom_linerange, geom_hline,
    geom_rect, geom_vline, facet_wrap, scale_x_datetime, scale_y_continuous,
    scale_color_manual, scale_fill_manual, theme_bw, theme, element_text, labs,
    position_dodge2, position_dodge, guides, guide_legend, ggsave,
)
import warnings
warnings.filterwarnings('ignore')

from data_loading import read_data_ukraine
from fit_partial_credit_model_ncats_pyrosvi import (
    fit_partial_credit_model_ncats_pyrosvi,
)
from get_endpoints import get_endpoints
from utils import _futurama_palette, _build_interim_dcati, _per_draw_ratio

print("✓ Imports successful")

# %%

# =============================================================================
# Configuration
# =============================================================================

np.random.seed(42)
seed = 123

dir_data = "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/2025_project_Hope_Groups/data"
file_data = os.path.join(dir_data, "Ukraine_Hope_Groups_Baseline_Endline_Wide_Aug6.csv")

dir_out = "/Users/or105/sandbox/bIRTistic/py-ukraine-interim-260526"
dir_logs = os.path.join(dir_out, "logs")
os.makedirs(dir_out, exist_ok=True)
os.makedirs(dir_logs, exist_ok=True)

file_prefix = "pcm_1_interim"
svi_algorithm = 'AutoLowRankMultivariateNormal'

timing_file = os.path.join(dir_out, f"{file_prefix}_timing.csv")

# Monthly interim grid: auto-detected from endline submission dates after the
# data is loaded (see `di = ...` below).
di = None

print(f"Output dir: {dir_out}")

# %%

# =============================================================================
# Load + preprocess full Ukraine data
# =============================================================================

print("\nLoading Ukraine data...")
raw = read_data_ukraine(file_data)
dp = raw['dp'].copy()
dit = raw['dit'].copy()
dmeta = raw['dmeta'].copy()
print(f"  dp: {len(dp):,} rows | dit: {len(dit):,} rows | dmeta: {len(dmeta):,} rows")

# drop participants with unknown displacement status
tmp = (
    dmeta[['pid', 'time_label', 'displacement_status']]
    .drop_duplicates()
    .rename(columns={'pid': 'pid_label'})
    .dropna(subset=["displacement_status"], how='all')
)
dp = dp.merge(tmp, on=['pid_label', 'time_label'], how='inner', validate='many_to_one')

dp1 = dp[~dp['item_label'].str.contains('agg')].copy()
dp1['y_stan'] = dp1['y'] + 1
dp1 = dp1.merge(dit[['item_label', 'item_type']], on='item_label', how='left')

# item_time_id: sequential id of (item_label, time) within each item_type.
item_time_df = (
    dp1[['item_type', 'item_label', 'time']]
    .drop_duplicates()
    .sort_values(['item_type', 'time', 'item_label'])
    .reset_index(drop=True)
)
item_time_df['item_time_id'] = item_time_df.groupby('item_type').cumcount() + 1
dp1 = dp1.merge(item_time_df, on=['item_label', 'time', 'item_type'], how='left')

dp1 = dp1.merge(
    dit[['item_type', 'item_type_id']].drop_duplicates(),
    on='item_type', how='left',
)
dp1 = dp1.sort_values(['item_type_id', 'pid', 'time', 'item_label']).reset_index(drop=True)
dp1['oid']  = range(1, len(dp1) + 1)
dp1['oidt'] = dp1.groupby('item_type').cumcount() + 1

print(f"  Pre-processed dp1: {len(dp1):,} observations | participants: {dp1['pid'].nunique()}")

# Auto-detect monthly interim grid from the endline submission dates.
# Each interim row spans one full calendar month: ``month_start`` (inclusive)
# to ``interim_date`` (month-end, the cumulative cutoff used by the fits).
_endline_dates = pd.to_datetime(
    dp1.loc[dp1['time_label'] == 'Endline', 'submission_date']
).dropna()
_start = _endline_dates.min().replace(day=1)
_end = (_endline_dates.max() + pd.offsets.MonthEnd(0)).normalize()
month_starts = pd.date_range(_start, _end, freq='MS')
di = pd.DataFrame({
    'interim_id': range(1, len(month_starts) + 1),
    'month_start': month_starts,
    'interim_date': month_starts + pd.offsets.MonthEnd(0),
})
di['interim_month_year'] = di['month_start'].dt.strftime('%Y-%b')
print(f"Interim grid: {len(di)} months from {di['month_start'].min().date()} to {di['interim_date'].max().date()}")

# %%

# =============================================================================
# Data viz: accruing participants over time
# =============================================================================

print("Plotting accruing participants...")
acc = (
    dp1[dp1['time_label'] == 'Endline']
    .groupby('submission_date')['pid'].nunique().reset_index(name='n_pid')
    .sort_values('submission_date')
)
acc['cn_pid'] = acc['n_pid'].cumsum()
acc['submission_date'] = pd.to_datetime(acc['submission_date'])

# Build per-interim background bands. One band per calendar month, spanning
# ``month_start`` (inclusive) to ``interim_date`` (month-end).
di_bands = di.copy()
di_bands['month_start'] = pd.to_datetime(di_bands['month_start'])
di_bands['interim_date'] = pd.to_datetime(di_bands['interim_date'])
di_bands = di_bands.sort_values('month_start').reset_index(drop=True)
di_bands['xmin'] = di_bands['month_start']
di_bands['xmax'] = di_bands['interim_date']
y_max = float(acc['cn_pid'].max() * 1.05)
di_bands['ymin'] = 0.0
di_bands['ymax'] = y_max
di_bands['interim_month_year'] = pd.Categorical(
    di_bands['interim_month_year'],
    categories=di_bands['interim_month_year'].tolist(),
    ordered=True,
)
interim_pal = dict(zip(
    di_bands['interim_month_year'].cat.categories,
    _futurama_palette(len(di_bands)),
))

p = (
    ggplot()
    + geom_rect(
        aes(xmin='xmin', xmax='xmax', ymin='ymin', ymax='ymax', fill='interim_month_year'),
        data=di_bands,
        alpha=0.25,
        inherit_aes=False,
    )
    + geom_vline(
        aes(xintercept='interim_date'),
        data=di_bands,
        colour='#666666',
        linetype='dotted',
        size=0.3,
    )
    + geom_line(aes(x='submission_date', y='cn_pid'), data=acc, size=0.8)
    + scale_fill_manual(values=interim_pal)
    + scale_y_continuous(expand=[0, 0])
    + scale_x_datetime(expand=[0, 0])
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        legend_position='right',
        figure_size=(9, 7),
    )
    + guides(fill=guide_legend(ncol=1))
    + labs(x='Endline assessment completed', 
           y='total number participants', 
           fill = 'interim analysis period')
)
ggsave(p, os.path.join(dir_out, 'ukraine_interim_accrueing_number_participants.pdf'),
       width=9, height=7)

# %%

# =============================================================================
# Data viz: per-interim displacement-status mix
# =============================================================================

print("Plotting displacement-status mix per interim...")
_disp_rows = []
for _, _row in di.iterrows():
    _interim_id = int(_row['interim_id'])
    _interim_date = _row['interim_date']
    _sub = dp1[dp1['submission_date'] <= _interim_date]
    if _sub.empty:
        continue
    # Keep only participants with both baseline + endline (same cohort used by the fits).
    _n_per = (
        _sub.groupby(['pid', 'item_label'])['time']
        .nunique().reset_index(name='n_times')
    )
    _complete = (
        _n_per.groupby('pid')['n_times']
        .apply(lambda x: bool((x == 2).all())).reset_index(name='complete')
    )
    _keep_pids = _complete.loc[_complete['complete'], 'pid'].tolist()
    _sub = _sub[_sub['pid'].isin(_keep_pids)]
    if _sub.empty:
        continue
    # One displacement_status per pid (drop NaN); count + normalise.
    _pid_disp = (
        _sub.dropna(subset=['displacement_status'])
            .drop_duplicates('pid')[['pid', 'displacement_status']]
    )
    _counts = _pid_disp.groupby('displacement_status').size().reset_index(name='n')
    _counts['interim_id'] = _interim_id
    _counts['interim_date'] = _interim_date
    _counts['proportion'] = _counts['n'] / _counts['n'].sum()
    _disp_rows.append(_counts)

disp_df = pd.concat(_disp_rows, ignore_index=True)
disp_df['interim_date'] = pd.to_datetime(disp_df['interim_date'])
_disp_levels = sorted(disp_df['displacement_status'].unique())
disp_pal = dict(zip(_disp_levels, _futurama_palette(len(_disp_levels))))

p = (
    ggplot(disp_df, aes(x='interim_date', y='proportion', fill='displacement_status'))
    + geom_col(position='stack', width=20)
    + scale_y_continuous(labels=lambda l: [f'{v:.0%}' for v in l], expand=[0, 0])
    + scale_x_datetime(expand=[0, 0])
    + scale_fill_manual(values=disp_pal)
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        legend_position='right',
        figure_size=(9, 5),
    )
    + guides(fill=guide_legend(ncol=1))
    + labs(
        x='interim analysis end',
        y='proportion of participants',
        fill='displacement status',
    )
)
ggsave(p, os.path.join(dir_out, 'ukraine_interim_displacement_status_proportions.pdf'),
       width=9, height=5)

# %%

# =============================================================================
# Fit PCM at each interim via AutoLowRankMVN SVI
# =============================================================================

print(f"\n{'='*70}\nFitting PCM (AutoLowRankMVN) per interim\n{'='*70}")

interim_dcati = {}      # interim_id -> dcati
interim_zarr = {}       # interim_id -> draws_file path
interim_endpoints = []  # list of get_endpoints DataFrames with interim_id

t0 = time.time() 
for i, row in di.iterrows():
    interim_id = int(row['interim_id'])
    interim_date = row['interim_date']
    print(f"\n--- Interim {interim_id}: {interim_date.date()} ---")

    dcati = _build_interim_dcati(dp1, interim_date)
    if dcati.empty or dcati['pid'].nunique() < 2:
        print(f"  Skipping (insufficient complete participants)")
        continue
    print(f"  n_obs={len(dcati):,} | n_pid={dcati['pid'].nunique()} | n_items={dcati['item_label'].nunique()}")

    interim_prefix = os.path.join(dir_out, f"{file_prefix}_{interim_id}")
    result = fit_partial_credit_model_ncats_pyrosvi(
        dit,
        dcati,
        output_file_prefix=interim_prefix,
        algorithm=svi_algorithm,
        lr=0.01,
        num_steps=10000,
        output_samples=4000,
        seed=seed,
        x_formula="~ time - 1",
        resume=True,
        with_core_analyses=True,
        with_additional_analyses=False,
    )

    zarr_path = f"{interim_prefix}_draws.zarr"
    interim_dcati[interim_id] = dcati
    interim_zarr[interim_id] = zarr_path

    pos = get_endpoints(
        dp1=dcati,
        dit=dit,
        draws_file=zarr_path,
        categorical_threshold=2,
        endpoint_type='items',
        param_name='ordered_prob_by_cat_qu_fit',
    )
    pos['interim_id'] = interim_id
    interim_endpoints.append(pos)

if not interim_endpoints:
    raise RuntimeError("No interim endpoints computed.")

t1 = time.time() 
post_minutes = (t1 - t0) / 60.0
print(f"  Posterior sampling completed in {post_minutes:.2f} minutes")

print("Extracting timing information...")
timing_data = pd.DataFrame({
    'mins_total': [round(post_minutes, 3)],
    })
timing_data.to_csv(timing_file, index=False)
print(f"Saved timing information to: {timing_file}")
# %%

ipos = pd.concat(interim_endpoints, ignore_index=True)
ipos = ipos.merge(di, on='interim_id')
# get_endpoints output omits endpoint_measure; pull it from dit.
ipos = ipos.merge(
    dit[['item_label', 'endpoint_measure']].drop_duplicates('item_label'),
    on='item_label', how='left',
)
ipos['item_label_long'] = ipos['group_label_long'] + np.where(
    ipos['item_label_short'].notna(), ' --- ' + ipos['item_label_short'].fillna(''), ''
)
ipos.to_csv(os.path.join(dir_out, f"{file_prefix}_primary_endpoints_b.csv"), index=False)
print(f"\nSaved primary endpoints to: {file_prefix}_primary_endpoints_b.csv")

# %%

# =============================================================================
# Plot: ratio per item over interim time
# =============================================================================

print("\nPlotting ratio over interim time...")
ratio_df = ipos[ipos['variable'] == 'ratio'].copy()
ratio_df['facet'] = ratio_df['group_label_long'] + '\n( ' + ratio_df['endpoint_measure'].fillna('') + ' )'
ratio_pal = dict(zip(sorted(ratio_df['item_label_long'].unique()),
                     _futurama_palette(ratio_df['item_label_long'].nunique())))

p = (
    ggplot(ratio_df, aes(x='interim_month_year', y='median',
                         colour='item_label_long', group='item_label_long'))
    + geom_hline(yintercept=0, colour='black')
    + geom_line()
    + geom_linerange(aes(ymin='q_lower', ymax='q_upper'), colour='#7f7f7f',
                     position=position_dodge(width=0.5))
    + scale_color_manual(values=ratio_pal)
    + scale_y_continuous(labels=lambda l: [f'{v:.0%}' for v in l], limits=[-2, 2])
    + facet_wrap('~ facet', ncol=4)
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        legend_position='top',
        figure_size=(14, 16),
    )
    + guides(colour=guide_legend(ncol=3))
    + labs(
        x='\ninterim analysis date',
        y='improvement in outcomes\n(1-Endline/Baseline or Endline/Baseline-1)',
        colour='survey items',
    )
)
ggsave(p, os.path.join(dir_out, f"{file_prefix}_percentchanges_over_time.pdf"),
       width=14, height=16, limitsize=False)

# %%

# =============================================================================
# Plot: diff per item over interim time
# =============================================================================

print("Plotting diff over interim time...")
diff_df = ipos[ipos['variable'] == 'diff'].copy()
diff_df['facet'] = diff_df['group_label_long'] + '\n( ' + diff_df['endpoint_measure'].fillna('') + ' )'
diff_pal = dict(zip(sorted(diff_df['item_label_long'].unique()),
                    _futurama_palette(diff_df['item_label_long'].nunique())))

p = (
    ggplot(diff_df, aes(x='interim_month_year', y='median',
                        colour='item_label_long', group='item_label_long'))
    + geom_hline(yintercept=0, colour='black')
    + geom_line()
    + geom_linerange(aes(ymin='q_lower', ymax='q_upper'), colour='#7f7f7f',
                     position=position_dodge(width=0.5))
    + scale_color_manual(values=diff_pal)
    + facet_wrap('~ facet', ncol=4)
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        legend_position='top',
        figure_size=(14, 16),
    )
    + guides(colour=guide_legend(ncol=3))
    + labs(
        x='\ninterim analysis date',
        y='improvement in outcomes\n(Baseline - Endline or Endline - Baseline)',
        colour='survey items',
    )
)
ggsave(p, os.path.join(dir_out, f"{file_prefix}_differences_over_time.pdf"),
       width=14, height=16, limitsize=False)

# %%

# =============================================================================
# Per-draw relative improvement (ratio / per-draw average ratio across items)
# =============================================================================

print("\nComputing per-draw relative improvement per interim...")
rel_rows = []
for interim_id, zarr_path in interim_zarr.items():
    dcati = interim_dcati[interim_id]
    per_draw = _per_draw_ratio(zarr_path, dcati, dit, categorical_threshold=2)
    # Normalise to per-draw mean ratio across items.
    mean_per_draw = per_draw.groupby('draw')['ratio'].mean().rename('ratio_avg').reset_index()
    per_draw = per_draw.merge(mean_per_draw, on='draw')
    per_draw['rel'] = per_draw['ratio'] / per_draw['ratio_avg']

    quantiles = (
        per_draw.groupby(['item_type', 'item_label', 'item_high_label'])['rel']
        .quantile([0.025, 0.25, 0.5, 0.75, 0.975])
        .unstack()
        .reset_index()
    )
    quantiles.columns = ['item_type', 'item_label', 'item_high_label',
                         'q_lower', 'iqr_lower', 'median', 'iqr_upper', 'q_upper']
    quantiles['interim_id'] = interim_id
    rel_rows.append(quantiles)

rel_df = pd.concat(rel_rows, ignore_index=True)
rel_df = rel_df.merge(di, on='interim_id')
rel_df = rel_df.merge(
    dit[['item_label', 'group_label_long', 'item_label_short', 'endpoint_measure']]
       .drop_duplicates('item_label'),
    on='item_label', how='left',
)
rel_df['item_label_long'] = rel_df['group_label_long'] + np.where(
    rel_df['item_label_short'].notna(), ' --- ' + rel_df['item_label_short'].fillna(''), ''
)
rel_df['facet'] = rel_df['group_label_long'] + '\n( ' + rel_df['endpoint_measure'].fillna('') + ' )'

rel_pal = dict(zip(sorted(rel_df['item_label_long'].unique()),
                   _futurama_palette(rel_df['item_label_long'].nunique())))

rel_df.to_csv(os.path.join(dir_out, f"{file_prefix}_primary_endpoints_c.csv"), index=False)
print(f"Saved relative endpoints to: {file_prefix}_primary_endpoints_c.csv")

p = (
    ggplot(rel_df, aes(x='interim_month_year', y='median',
                       colour='item_label_long', group='item_label_long'))
    + geom_hline(yintercept=1, colour='black')
    + geom_line()
    + geom_linerange(aes(ymin='q_lower', ymax='q_upper'), colour='#7f7f7f',
                     position=position_dodge(width=0.5))
    + scale_color_manual(values=rel_pal)
    + facet_wrap('~ facet', ncol=4)
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        legend_position='top',
        figure_size=(14, 12),
    )
    + guides(colour=guide_legend(ncol=3))
    + labs(
        x='\ninterim analysis date',
        y='improvement in outcomes relative to average improvement',
        colour='survey items',
    )
)
ggsave(p, os.path.join(dir_out, f"{file_prefix}_relativeimprovement_over_time.pdf"),
       width=14, height=12, limitsize=False)

# %%

# =============================================================================
# Predictive probability of success (PPS) at interim 1
# =============================================================================
#
# Numerical approximation of PPS per item_label, evaluated at the first interim.
# Strategy follows dev/amortised_decision_making.md §2 Target 2.
#
# H_1 (per item_label): the directional ratio endpoint > 0 (i.e. an improvement
# in outcomes vs Baseline, where direction respects ``item_high_label``).
#
# Procedure:
#   1. Get ypred draws from the interim-1 fit; sample S = 10 hypothetical
#      future data sets z by selecting S draws from ``ypred``.
#   2. For each s, augment the interim-1 dcati with z and refit the PCM via
#      ``fit_partial_credit_model_ncats_pyrosvi`` with
#      ``algorithm='AutoLowRankMultivariateNormal'``, no extra analyses.
#   3. Compute p(H_1 | x, z) per item from per-draw ratios on the refit.
#   4. Save raw p(H_1 | x, z) samples to a pickle.
#   5. Approximate PPS as the fraction of s with p(H_1 | x, z) > eta and save
#      to a CSV/DataFrame.

import arviz as az  # local import to keep top of file lean

S = 10
eta_pps = 0.89  # decision threshold from the doc
pps_interim_id = 1

pps_zarr = interim_zarr.get(pps_interim_id)
pps_dcati = interim_dcati.get(pps_interim_id)
if pps_zarr is None or pps_dcati is None or pps_dcati.empty:
    print(f"\n[PPS] Skipping: no fit available for interim {pps_interim_id}.")
else:
    print(f"\n{'='*70}\nPPS for interim {pps_interim_id}\n{'='*70}")

    # Current P(H_1 | x) per item using the interim-1 posterior.
    ratio_x = _per_draw_ratio(pps_zarr, pps_dcati, dit, categorical_threshold=2)
    p_h1_x = (
        ratio_x.groupby(['item_label', 'item_type', 'item_high_label'])['ratio']
        .apply(lambda r: float((r > 0).mean()))
        .reset_index(name='p_h1_x')
    )
    print(f"  P(H_1 | x) computed for {len(p_h1_x)} items")

    # Pull ypred draws (1-indexed integers) from the interim-1 fit and pick
    # S equally-spaced draws to use as hypothetical z.
    _idata = az.from_zarr(pps_zarr)
    if 'ypred' not in _idata.posterior.data_vars:
        raise RuntimeError("ypred not found in interim-1 posterior; cannot draw z.")
    _ypred = _idata.posterior['ypred'].values  # (chain, draw, N_total)
    _n_draw = _ypred.shape[1]
    _draw_idx = np.linspace(0, _n_draw - 1, S, dtype=int)

    # For each s, augment the dataset with one ypred draw (one extra "shadow"
    # participant per existing pid using the ypred outcomes) and refit.
    _pid_offset = int(pps_dcati['pid'].max())
    p_h1_xz_rows = []

    for s_idx, draw_i in enumerate(_draw_idx):
        s_label = s_idx + 1
        print(f"\n--- PPS sample {s_label}/{S} (ypred draw {draw_i}) ---")
        yz = np.asarray(_ypred[0, draw_i, :]).astype(int)  # (N_total,)

        z_dcati = pps_dcati.copy()
        z_dcati['y_stan'] = yz
        z_dcati['y'] = yz - 1
        z_dcati['pid'] = z_dcati['pid'] + _pid_offset * s_label  # fresh pid block per s

        aug = pd.concat([pps_dcati, z_dcati], ignore_index=True)
        aug = aug.sort_values(['item_type_id', 'pid', 'time', 'item_label']).reset_index(drop=True)
        aug['oid'] = range(1, len(aug) + 1)
        aug['oidt'] = aug.groupby('item_type').cumcount() + 1

        aug_prefix = os.path.join(dir_out, f"{file_prefix}_pps_i{pps_interim_id}_s{s_label}")
        fit_partial_credit_model_ncats_pyrosvi(
            dit,
            aug,
            output_file_prefix=aug_prefix,
            algorithm=svi_algorithm,
            lr=0.01,
            num_steps=10000,
            output_samples=4000,
            seed=seed + s_label,
            x_formula="~ time - 1",
            resume=True,
            with_core_analyses=False,
            with_additional_analyses=False,
        )

        aug_zarr = f"{aug_prefix}_draws.zarr"
        ratio_xz = _per_draw_ratio(aug_zarr, aug, dit, categorical_threshold=2)
        sample_p = (
            ratio_xz.groupby(['item_label', 'item_type', 'item_high_label'])['ratio']
            .apply(lambda r: float((r > 0).mean()))
            .reset_index(name='p_h1_xz')
        )
        sample_p['s'] = s_label
        p_h1_xz_rows.append(sample_p)

    p_h1_xz = pd.concat(p_h1_xz_rows, ignore_index=True)
    pkl_path = os.path.join(dir_out, f"{file_prefix}_pps_i{pps_interim_id}_p_h1_xz.pkl")
    p_h1_xz.to_pickle(pkl_path)
    print(f"\nSaved P(H_1 | x, z) samples to: {pkl_path}")

    pps_df = (
        p_h1_xz.groupby(['item_label', 'item_type', 'item_high_label'])['p_h1_xz']
        .apply(lambda p: float((p > eta_pps).mean()))
        .reset_index(name='pps')
    )
    pps_df = pps_df.merge(p_h1_x, on=['item_label', 'item_type', 'item_high_label'], how='left')
    pps_df['eta'] = eta_pps
    pps_df['S'] = S
    pps_df['interim_id'] = pps_interim_id

    csv_path = os.path.join(dir_out, f"{file_prefix}_pps_i{pps_interim_id}.csv")
    pps_df.to_csv(csv_path, index=False)
    print(f"Saved PPS table to: {csv_path}")
    print(pps_df.head(10).to_string(index=False))

print("\n✓ Interim analyses complete")
