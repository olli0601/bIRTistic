#!/usr/bin/env python3
"""
Ukraine interim analysis: BATCHED SMC sampler with resample-move (Case A).

Prototype of the vectorised SMC: instead of S independent SMC runs per interim
(``Ukraine_interim_analysis_with_SMC_resample.py``), all S future samples share
ONE tempering schedule and are moved together as a single ``(S, K, D)`` particle
tensor under a vmapped MALA kernel that compiles once (no per-sample recompile,
no process pool). The shared schedule advances by the most conservative sample
(``min_s`` tempering ESS hits the target). Via
``fit_interim_SMC_batched_PPS_of_posterior_xz_from_x``.

Outputs mirror the other interim scripts (``_SMCbatched`` suffix): p_h1_xz pkl,
PPS csv, box-stats pkl, p_h1_xz boxplot, perf long-form pkl, timing csv.

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/Ukraine_interim_analysis_SMC_batched.py
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
    ggplot, aes, geom_boxplot, geom_hline, facet_wrap,
    scale_y_continuous, scale_fill_manual,
    theme_bw, theme, element_text, labs,
)
import warnings
warnings.filterwarnings('ignore')

from data_loading import read_data_ukraine
from fit_partial_credit_model import fit_partial_credit_model_ncats_pyrosvi
from fit_interim import (
    get_interim_x,
    get_interim_z_from_ypredi,
    fit_interim_SMC_batched_PPS_of_posterior_xz_from_x,
)
from utils import _futurama_palette

print("✓ Imports successful")

# %%

np.random.seed(42)
seed = 123

dir_data = "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/2025_project_Hope_Groups/data"
file_data = os.path.join(dir_data, "Ukraine_Hope_Groups_Baseline_Endline_Wide_Aug6.csv")
dir_out = "/Users/or105/sandbox/bIRTistic/py-ukraine-interim-with-SMC-batched-260526"
os.makedirs(dir_out, exist_ok=True)

file_prefix = "pcm_1_interim"
svi_algorithm = 'AutoLowRankMultivariateNormal'
x_formula = "~ time - 1"

pps_z_total = 12           # PROTOTYPE: bump to 200 once validated
pps_H1_def = 0.5
pps_ProbH1_thresh = 0.89
categorical_threshold = 2

# Batched-SMC fitting-method arguments (one shared schedule across all S).
fitting_method_args = dict(
    n_particles=128,
    ess_frac_target=0.5,
    n_move_steps=20,
    init_step_size=0.02,
    max_temps=300,
    x_formula=x_formula,
    seed=seed,
)

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

item_time_df = (
    dp1[['item_type', 'item_label', 'time']]
    .drop_duplicates()
    .sort_values(['item_type', 'time', 'item_label'])
    .reset_index(drop=True)
)
item_time_df['item_time_id'] = item_time_df.groupby('item_type').cumcount() + 1
dp1 = dp1.merge(item_time_df, on=['item_label', 'time', 'item_type'], how='left')
dp1 = dp1.merge(
    dit[['item_type', 'item_type_id']].drop_duplicates(), on='item_type', how='left',
)
dp1 = dp1.sort_values(['item_type_id', 'pid', 'time', 'item_label']).reset_index(drop=True)
dp1['oid'] = range(1, len(dp1) + 1)
dp1['oidt'] = dp1.groupby('item_type').cumcount() + 1
print(f"  Pre-processed dp1: {len(dp1):,} observations | participants: {dp1['pid'].nunique()}")

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
# Batched SMC PPS at every interim
# =============================================================================

n_full = dp1['pid'].nunique()
p_h1_xz_rows = []
sched_rows = []
pps_timing_rows = []
t_all0 = time.time()
for interim_id in di['interim_id']:
    interim_id = int(interim_id)
    interim_date = pd.to_datetime(di.loc[di['interim_id'] == interim_id, 'interim_date'].iloc[0])

    xi = get_interim_x(dp1, interim_date)
    if xi.empty or xi['pid'].nunique() < 2:
        continue
    interim_m = n_full - xi['pid'].nunique()
    if interim_m <= 0:
        continue
    print(f"\n{'='*70}\nbatched-SMC-PPS for interim {interim_id} ({interim_date.date()})\n{'='*70}")
    print(f"  n_obs={len(xi):,} | n_pid={xi['pid'].nunique()} | n_items={xi['item_label'].nunique()}")

    t0 = time.time()
    interim_prefix = os.path.join(dir_out, f"{file_prefix}_{interim_id}")
    fit_partial_credit_model_ncats_pyrosvi(
        dit, xi,
        output_file_prefix=interim_prefix,
        algorithm=svi_algorithm,
        lr=0.01, num_steps=10000, output_samples=4000, seed=seed,
        x_formula=x_formula, resume=True,
        with_core_analyses=True, with_additional_analyses=False, verbose=False,
    )
    zi = get_interim_z_from_ypredi(
        xi, f"{interim_prefix}_draws.zarr", interim_m, pps_z_total=pps_z_total, seed=seed,
    )
    p, schedule = fit_interim_SMC_batched_PPS_of_posterior_xz_from_x(
        xi=xi, zi=zi, dit=dit,
        fitting_method_args=fitting_method_args,
        draws=None, draws_file=f"{interim_prefix}_draws.zarr",
        pps_z_total=pps_z_total,
        pps_H1_def=pps_H1_def,
        pps_ProbH1_thresh=pps_ProbH1_thresh,
        categorical_threshold=categorical_threshold,
        save_to_file=False,
        verbose=True,
    )
    mins_interim_id = (time.time() - t0) / 60.0
    p['interim_id'] = interim_id
    p['interim_date'] = interim_date
    schedule['interim_id'] = interim_id
    p_h1_xz_rows.append(p)
    sched_rows.append(schedule)
    pps_timing_rows.append({'interim_id': interim_id, 'mins_interim_id': round(mins_interim_id, 3),
                            'n_temps': int(len(schedule))})
    print(f"  interim {interim_id} batched-SMC-PPS done in {mins_interim_id:.2f} min ({len(schedule)} temps)")

if not p_h1_xz_rows:
    raise RuntimeError("[batched-SMC-PPS] no interims with missing participants to evaluate.")

p_h1_xz = pd.concat(p_h1_xz_rows, ignore_index=True)
sched_all = pd.concat(sched_rows, ignore_index=True)
mins_total = (time.time() - t_all0) / 60.0
pps_timing = pd.DataFrame(pps_timing_rows)
pps_timing['mins_total'] = round(mins_total, 3)
print(f"\nAll interims batched-SMC-PPS done in {mins_total:.2f} min")

sched_all.to_csv(os.path.join(dir_out, f"{file_prefix}_pps_SMCbatched_schedule.csv"), index=False)
pkl_path = os.path.join(dir_out, f"{file_prefix}_pps_SMCbatched_p_h1_xz.pkl")
p_h1_xz.to_pickle(pkl_path)
print(f"Saved batched-SMC P(H_1 | x, z) samples to: {pkl_path}")
pps_timing.to_csv(os.path.join(dir_out, f"{file_prefix}_pps_SMCbatched_timing.csv"), index=False)

# %%

pps_df = (
    p_h1_xz.groupby(['interim_id', 'interim_date', 'item_label', 'item_type', 'item_high_label'])['p_h1_xz']
    .apply(lambda p: float((p > pps_ProbH1_thresh).mean()))
    .reset_index(name='pps')
)
pps_df['eta'] = pps_ProbH1_thresh
pps_df['S'] = pps_z_total
csv_path = os.path.join(dir_out, f"{file_prefix}_pps_SMCbatched.csv")
pps_df.to_csv(csv_path, index=False)
print(f"Saved batched-SMC PPS table to: {csv_path}")
print(pps_df.head(10).to_string(index=False))

# %%

# =============================================================================
# Plot: distribution of p(H_1 | x, z) per item across the S samples
# =============================================================================

print("\nPlotting batched-SMC p(H_1 | x, z) distribution...")
tmp = p_h1_xz.merge(
    dit[['item_label', 'group_label_long', 'item_label_short']].drop_duplicates(),
    on='item_label', how='left',
)
tmp['item_label_long'] = tmp['group_label_long'] + np.where(
    tmp['item_label_short'].notna(), '\n' + tmp['item_label_short'], ''
)
all_items = list(pd.unique(tmp['item_label_long']))
color_dict = dict(zip(all_items, _futurama_palette(len(all_items))))

tmp = tmp.merge(di[['interim_id', 'interim_month_year']], on='interim_id', how='left')
order = (
    di[di['interim_id'].isin(tmp['interim_id'].unique())]
    .sort_values('interim_date')['interim_month_year'].tolist()
)
tmp['interim_month_year'] = pd.Categorical(tmp['interim_month_year'], categories=order, ordered=True)

box_stats = (
    tmp.groupby(['item_label_long', 'interim_month_year'], observed=True)['p_h1_xz']
    .agg(
        min='min',
        q025=lambda s: s.quantile(0.025),
        q25=lambda s: s.quantile(0.25),
        q50=lambda s: s.quantile(0.50),
        q75=lambda s: s.quantile(0.75),
        q975=lambda s: s.quantile(0.975),
        max='max',
    )
    .reset_index()
)
pkl_path = os.path.join(dir_out, f"{file_prefix}_pps_SMCbatched_p_h1_xz_boxplot.pkl")
box_stats.to_pickle(pkl_path)
print(f"Saved batched-SMC p(H_1 | x, z) box-stats to: {pkl_path}")

p = (
    ggplot(box_stats, aes(x='interim_month_year', fill='item_label_long'))
    + geom_boxplot(
        aes(ymin='q025', lower='q25', middle='q50', upper='q75', ymax='q975',
            group='interim_month_year'),
        stat='identity',
    )
    + geom_hline(yintercept=pps_ProbH1_thresh, colour='black', size=1.5)
    + facet_wrap('~ item_label_long', ncol=4)
    + scale_fill_manual(values=color_dict)
    + scale_y_continuous(
        limits=[0, 1],
        breaks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        labels=['0%', '20%', '40%', '60%', '80%', '100%'],
    )
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        legend_position='none',
        figure_size=(15, 15),
        strip_text_y=element_text(angle=0),
    )
    + labs(x='Interim', y='p(H_1 | x, z)  [SMC batched]')
)
box_pdf = os.path.join(dir_out, f"{file_prefix}_pps_SMCbatched_p_h1_xz_boxplot.pdf")
p.save(box_pdf, verbose=False, limitsize=False)
print(f"Saved batched-SMC p(H_1 | x, z) boxplot to: {box_pdf}")

# %%

# =============================================================================
# Performance long-form pkl (post-move particles uniformly weighted:
# ESS/particle = 1, ESS = K, E(w^2) = 1/K^2; time = per-interim stage).
# =============================================================================

print("\nBuilding batched-SMC performance long-form table...")
K = int(fitting_method_args['n_particles'])
perf = pps_timing.merge(di[['interim_id', 'interim_month_year']], on='interim_id', how='left')
perf_order = (
    di[di['interim_id'].isin(perf['interim_id'].unique())]
    .sort_values('interim_date')['interim_month_year'].tolist()
)
perf['interim_month_year'] = pd.Categorical(perf['interim_month_year'], categories=perf_order, ordered=True)
perf['ess'] = float(K)
perf['ess_per_particle'] = 1.0
perf['ew2'] = 1.0 / (K ** 2)
perf['time_min'] = perf['mins_interim_id']
perf['s'] = 1  # one shared schedule per interim

metric_labels = {
    'ess': 'ESS',
    'ess_per_particle': 'ESS / particle',
    'ew2': 'E(w^2)',
    'time_min': 'time (min)',
}
perf_long = perf.melt(
    id_vars=['interim_id', 'interim_month_year', 's'],
    value_vars=list(metric_labels.keys()),
    var_name='metric', value_name='value',
)
perf_long['metric'] = pd.Categorical(
    perf_long['metric'].map(metric_labels),
    categories=list(metric_labels.values()), ordered=True,
)
perf_long['method'] = 'SMC (batched)'
pkl_path = os.path.join(dir_out, f"{file_prefix}_pps_SMCbatched_perf_long.pkl")
perf_long.to_pickle(pkl_path)
print(f"Saved batched-SMC performance long-form to: {pkl_path}")

print("\n✓ batched-SMC interim analysis complete")
