#!/usr/bin/env python3
"""
Colombia interim analyses, Python refactor.

Mirrors ``scripts-R/Colombia_interim_analyses_v251007.Rmd``, but:
- replaces the 2cats Stan credit model with the ncats partial credit model
  (``src/numpyro/partial_credit_model_ncats_v260413.pyro``);
- replaces the Stan HMC fits with NumPyro SVI using
  ``AutoLowRankMultivariateNormal`` via
  ``fit_partial_credit_model_ncats_pyrosvi``;
- reuses ``get_endpoints`` for primary endpoint summaries;
- implements the cross-interim plots in plotnine.

Iteration cost: AutoLowRankMVN fits Colombia data in ~10-15s each, so the full
monthly interim sweep (~11 fits) finishes in a few minutes.

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/Colombia_interim_analyses.py
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

from data_loading import read_data_colombia
from model_pcm import PartialCreditModel
from utils import _futurama_palette, _build_interim_dcati

print("✓ Imports successful")

# %%

# =============================================================================
# Configuration
# =============================================================================

np.random.seed(42)
seed = 123

dir_data = "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/2025_project_Hope_Groups/data"
file_data = os.path.join(dir_data, "Colombia_data_baseline_endline_itemised_250927.csv")

dir_out = "/Users/or105/sandbox/bIRTistic/py-colombia-interim-260526"
dir_logs = os.path.join(dir_out, "logs")
os.makedirs(dir_out, exist_ok=True)
os.makedirs(dir_logs, exist_ok=True)

file_prefix = "pcm_1_interim"
svi_algorithm = 'AutoLowRankMultivariateNormal'

timing_file = os.path.join(dir_out, f"{file_prefix}_timing.csv")

# Monthly interim grid: 2024-07-01 .. 2025-05-01 (matches the Rmd).
interim_dates = pd.date_range("2024-07-01", "2025-05-01", freq="MS")
di = pd.DataFrame({
    'interim_id': range(1, len(interim_dates) + 1),
    'interim_date': interim_dates,
})
di['interim_month_year'] = di['interim_date'].dt.strftime('%Y-%b')

print(f"Output dir: {dir_out}")
print(f"Interim grid: {len(di)} dates from {di['interim_date'].min().date()} to {di['interim_date'].max().date()}")

# %%

# =============================================================================
# Load + preprocess full Colombia data
# =============================================================================

print("\nLoading Colombia data...")
raw = read_data_colombia(file_data)
dp = raw['dp'].copy()
dit = raw['dit'].copy()
dmeta = raw['dmeta'].copy()
print(f"  dp: {len(dp):,} rows | dit: {len(dit):,} rows | dmeta: {len(dmeta):,} rows")

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

# Build per-interim background bands. Each band runs from the previous interim
# date (or earliest endline submission) to the current interim date.
di_bands = di.copy()
di_bands['interim_date'] = pd.to_datetime(di_bands['interim_date'])
di_bands = di_bands.sort_values('interim_date').reset_index(drop=True)
acc_min = acc['submission_date'].min()
di_bands['xmin'] = di_bands['interim_date'].shift(1).fillna(acc_min)
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
           fill = 'interim analysisperiod')
)
ggsave(p, os.path.join(dir_out, 'colombia_interim_accrueing_number_participants.pdf'),
       width=9, height=7)

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
    _pcm = PartialCreditModel(
        dit=dit, dcati=dcati, x_formula="~ time - 1", seed=seed,
    )
    result = _pcm.fit_pyro_svi(
        output_file_prefix=interim_prefix,
        algorithm=svi_algorithm,
        lr=0.01,
        num_steps=10000,
        output_samples=4000,
        resume=True,
        with_core_analyses=True,
        with_additional_analyses=False,
    )

    zarr_path = f"{interim_prefix}_draws.zarr"
    interim_dcati[interim_id] = dcati
    interim_zarr[interim_id] = zarr_path

    pos = _pcm.get_endpoints(
        draws_file=zarr_path,
        categorical_threshold=3,
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
    + scale_y_continuous(labels=lambda l: [f'{v:.0%}' for v in l])
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
    _pcm = PartialCreditModel(dit=dit, dcati=dcati)
    per_draw = _pcm.get_endpoints_per_draw(
        draws_file=zarr_path,
        categorical_threshold=3,
        endpoint_type='items',
        param_name='ordered_prob_by_cat_qu_fit',
    )
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

print("\n✓ Interim analyses complete")
