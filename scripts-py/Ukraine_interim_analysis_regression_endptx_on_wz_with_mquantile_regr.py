#!/usr/bin/env python3
"""
Ukraine interim analysis: regression-based EVSI labels (Strong & Oakley, 2014),
continuous-target variant.

This is the more discriminative sibling of
``Ukraine_interim_analysis_regression_H1x_on_wz.py``. Instead of binarising the
endpoint ratio to ``pps_H1_x = 1{pps_ratio_x > pps_H1_def}`` and fitting a
binomial GLM, here we fit a Gaussian GLM on the continuous ``pps_ratio_x``:

  - draw  theta^(s) ~ p(theta | x)               (x-fit posterior draws)
  - draw  z^(s) ~ p(z | theta^(s))               (posterior-predictive ypred)
  - y^(s) := ratio_item(theta^(s))               (continuous endpoint ratio from x)
  - W^(s) := per-item summary of z^(s) at baseline and endline (see
             :func:`get_interim_endpt_and_w_from_poi` for details).

Per item we fit linear multi-quantile regressions at 11 levels (0.05..0.95)
via ``statsmodels.QuantReg``, monotone-correct the predicted quantiles, and
interpolate the conditional CDF at ``pps_H1_def`` to produce a continuous
``p_h1_xz = 1 - F_hat(pps_H1_def | w) in [0, 1]`` (Rao-Blackwellised). The PPS is the Monte-Carlo average
``S^-1 sum_s 1{p_h1_xz(W(z^(s))) > pps_ProbH1_thresh}``.

Outputs (``_RGEM_`` suffix to distinguish from the H1x ``_RG_`` variant):
p_h1_xz pkl, PPS csv, box-stats pkl, p_h1_xz boxplot, perf long-form pkl
(rho / R^2 / time(min)), timing csv, plus a per-interim diagnostic scatter
(w_ratio vs pps_ratio_x with Gaussian-GLM smoother + rho / R^2 annotations).

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/Ukraine_interim_analysis_regression_endptx_on_wz_with_mquantile_regr.py
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
    ggplot, aes, geom_boxplot, geom_hline, geom_point,
    geom_smooth, geom_text, facet_wrap,
    scale_colour_manual, scale_fill_manual, scale_y_continuous,
    theme_bw, theme, element_text, element_rect, labs,
)
import warnings
warnings.filterwarnings('ignore')

from data_loading import read_data_ukraine
from model_pcm import PartialCreditModel
from fit_interim import (
    fit_interim_regress_endptx_on_wz_with_mquantile_regr,
    fit_interim_regress_H1x_on_wz_per_item_summary,
)
from utils import _futurama_palette

print("✓ Imports successful")

# %%

# =============================================================================
# Configuration
# =============================================================================

np.random.seed(42)
seed = 123

dir_data = "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/2025_project_Hope_Groups/data"
file_data = os.path.join(dir_data, "Ukraine_Hope_Groups_Baseline_Endline_Wide_Aug6.csv")
dir_out = "/Users/or105/sandbox/bIRTistic/py-ukraine-interim-with-regression-on-endptx-wz-mquantile-regr-260702"
os.makedirs(dir_out, exist_ok=True)

file_prefix = "pcm_1_interim"
svi_algorithm = 'AutoLowRankMultivariateNormal'
x_formula = "~ time - 1"

pps_z_total = 4000                   # one z-sample per x-posterior draw
pps_H1_def = 0.5
pps_ProbH1_thresh = 0.89
categorical_threshold = 2


print(f"Output dir: {dir_out}")

# %%

# =============================================================================
# Load + preprocess data
# =============================================================================

print("\nLoading Ukraine data...")
raw = read_data_ukraine(file_data)
dp = raw['dp'].copy()
dit = raw['dit'].copy()
dmeta = raw['dmeta'].copy()

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

n_full = dp1['pid'].nunique()

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
print(f"  n_full={n_full} | looping over {len(di)} interim rounds")
# %%

# =============================================================================
# Per-interim loop: fit (resume) -> wa via get_interim_endpt_and_w_from_poi ->
# per-item Gaussian GLM on pps_ratio_x ~ w_ratio with predictive-Normal label
# conversion. Diagnostic scatter plot is built in a separate loop below.
# =============================================================================

p_h1_xz_rows = []
perf_rows = []
pps_timing_rows = []
t_all0 = time.time()
for interim_id in di['interim_id']:
    interim_id = int(interim_id)
    interim_date = pd.to_datetime(di.loc[di['interim_id'] == interim_id, 'interim_date'].iloc[0])

    xi = PartialCreditModel.get_interim_data_x(dp1, interim_date)
    if xi.empty or xi['pid'].nunique() < 2:
        print(f"\n[interim {interim_id}] skipped: empty xi"); continue
    interim_m = n_full - xi['pid'].nunique()
    if interim_m <= 0:
        print(f"\n[interim {interim_id}] skipped: no missing participants"); continue
    print(f"\n{'='*70}\nRegression-endptx for interim {interim_id} ({interim_date.date()})\n{'='*70}")
    print(f"  n_obs={len(xi):,} | n_pid={xi['pid'].nunique()} | m={interim_m}")

    t0 = time.time()
    
    # fit model on today's data x
    interim_prefix = os.path.join(dir_out, f"{file_prefix}_{interim_id}")
    model = PartialCreditModel(dit=dit, dcati=xi,
                                    x_formula=x_formula, seed=seed)
    fit = model.fit_pyro_svi(
        output_file_prefix=interim_prefix,
        algorithm=svi_algorithm,
        lr=0.01, num_steps=10000, output_samples=pps_z_total,
        resume=True,
        with_core_analyses=True, with_additional_analyses=False, verbose=False,
    )

    # expand to future data keeping draw index linked to posterior draws, and get w(z^s)
    zi = model.get_interim_z_from_ypredi(
        f"{interim_prefix}_draws.zarr", interim_m,
        pps_z_total=pps_z_total, seed=seed, keep_order=True,
    )
    wa = model.get_w(zi, categorical_threshold=categorical_threshold)

    # calculate endpoint on theta | x
    x_ratio = model.get_endpoints_per_draw(
        draws=fit['draws'], categorical_threshold=categorical_threshold,
        endpoint_type='items',
    )

    # calculate Prob(H1|x) ie mean over samples theta^s
    tmp = (
        model.get_p_h1(x_ratio, pps_H1_def=pps_H1_def)
        .rename(columns={'p_h1': 'pps_ProbH1_x'})
    )
    tmp['pps_H1_yes'] = (tmp['pps_ProbH1_x'] > pps_ProbH1_thresh).astype(int)

    # calculate 1{theta^s in H1} per draw
    x_ratio = x_ratio.rename(columns={'ratio': 'pps_ratio_x'})
    x_ratio['pps_H1_x'] = (x_ratio['pps_ratio_x'] > pps_H1_def).astype(int)
    x_ratio = x_ratio.merge(
        tmp[['item_label', 'pps_ProbH1_x', 'pps_H1_yes']],
        on='item_label', how='left',
    )
    wa = wa.merge(
        x_ratio[['draw', 'item_label', 'item_type', 'pps_ratio_x', 'pps_H1_x']],
        on=['draw', 'item_label', 'item_type'], how='inner',
    )

    # save data for regression
    pkl_path = os.path.join(dir_out, f"{file_prefix}_i{interim_id}_regression_training.pkl")
    wa.to_pickle(pkl_path)
    print(f"  saved training pkl: {pkl_path}")

    # do regression of endpoint on w(x)
    p_h1_xz_interim, perf_interim = fit_interim_regress_endptx_on_wz_with_mquantile_regr(
        wa, pps_H1_def=pps_H1_def, pps_ProbH1_thresh=pps_ProbH1_thresh,
    )
    p_h1_xz_interim['interim_id'] = interim_id
    p_h1_xz_interim['interim_date'] = interim_date
    p_h1_xz_rows.append(p_h1_xz_interim)

    perf_interim['interim_id'] = interim_id
    perf_interim['interim_date'] = interim_date
    perf_rows.append(perf_interim)

    mins_interim_id = (time.time() - t0) / 60.0
    pps_timing_rows.append({'interim_id': interim_id,
                            'mins_interim_id': round(mins_interim_id, 3)})
    print(f"  interim {interim_id} regression-endptx done in {mins_interim_id:.2f} min")

if not p_h1_xz_rows:
    raise RuntimeError("[RGE-PPS] no interims with missing participants to evaluate.")

# %%

# =============================================================================
# Aggregate across interims; save artifacts mirroring the other interim scripts.
# =============================================================================

p_h1_xz = pd.concat(p_h1_xz_rows, ignore_index=True)
p_h1_xz['pps_H1_def'] = pps_H1_def
p_h1_xz['pps_ProbH1_thresh'] = pps_ProbH1_thresh
p_h1_xz['S'] = pps_z_total
perf_all = pd.concat(perf_rows, ignore_index=True)
mins_total = (time.time() - t_all0) / 60.0
pps_timing = pd.DataFrame(pps_timing_rows)
pps_timing['mins_total'] = round(mins_total, 3)
print(f"\nAll interims RGE-PPS done in {mins_total:.2f} min")

perf_all.to_csv(os.path.join(dir_out, f"{file_prefix}_pps_RGEM_perf.csv"), index=False)
pkl_path = os.path.join(dir_out, f"{file_prefix}_pps_RGEM_p_h1_xz.pkl")
p_h1_xz.to_pickle(pkl_path)
print(f"Saved RGE P(H_1 | x, z) samples to: {pkl_path}")
pps_timing.to_csv(os.path.join(dir_out, f"{file_prefix}_pps_RGEM_timing.csv"), index=False)

pps_df = (
    p_h1_xz.groupby(['interim_id', 'interim_date', 'item_label', 'item_type',
                     'item_high_label'])['p_h1_xz']
    .apply(lambda p: float((p > pps_ProbH1_thresh).mean()))
    .reset_index(name='pps')
)
pps_df['eta'] = pps_ProbH1_thresh
pps_df['S'] = pps_z_total
pps_df.to_csv(os.path.join(dir_out, f"{file_prefix}_pps_RGEM.csv"), index=False)
print(pps_df.head(10).to_string(index=False))

perf = perf_all.merge(di[['interim_id', 'interim_month_year']], on='interim_id', how='left')
perf = perf.merge(pps_timing[['interim_id', 'mins_interim_id']], on='interim_id', how='left')
perf_order = (
    di[di['interim_id'].isin(perf['interim_id'].unique())]
    .sort_values('interim_date')['interim_month_year'].tolist()
)
perf['interim_month_year'] = pd.Categorical(perf['interim_month_year'],
                                            categories=perf_order, ordered=True)
perf['time_min'] = perf['mins_interim_id']
metric_labels = {'rho': 'rho', 'r2': 'R^2', 'time_min': 'time (min)'}
perf_long = perf.melt(
    id_vars=['interim_id', 'interim_month_year', 'item_label'],
    value_vars=list(metric_labels.keys()),
    var_name='metric', value_name='value',
)
perf_long['metric'] = pd.Categorical(
    perf_long['metric'].map(metric_labels),
    categories=list(metric_labels.values()), ordered=True,
)
perf_long['method'] = 'Regression endpt mquantile (Strong-Oakley)'
perf_long.to_pickle(os.path.join(dir_out, f"{file_prefix}_pps_RGEM_perf_long.pkl"))
print(f"Saved RGE performance long-form pkl.")

# %%

# =============================================================================
# Plot: distribution of p(H_1 | x, z) per item across the S samples.
# =============================================================================

tmp = p_h1_xz.merge(
    dit[['item_label', 'group_label_long', 'item_label_short']].drop_duplicates(),
    on='item_label', how='left',
)
tmp['item_label_long'] = tmp['group_label_long'] + np.where(
    tmp['item_label_short'].notna(), '\n' + tmp['item_label_short'], ''
)
tmp = tmp.merge(di[['interim_id', 'interim_month_year']], on='interim_id', how='left')
order = (
    di[di['interim_id'].isin(tmp['interim_id'].unique())]
    .sort_values('interim_date')['interim_month_year'].tolist()
)
tmp['interim_month_year'] = pd.Categorical(tmp['interim_month_year'],
                                           categories=order, ordered=True)

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
pkl_path = os.path.join(dir_out, f"{file_prefix}_pps_RGEM_p_h1_xz_boxplot.pkl")
box_stats.to_pickle(pkl_path)
print(f"Saved RGE p(H_1 | x, z) box-stats to: {pkl_path}")

all_items = list(pd.unique(tmp['item_label_long']))
color_dict = dict(zip(all_items, _futurama_palette(len(all_items))))

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
    + labs(x='Interim', y='p(H_1 | x, z)  [Regression endpt (Strong-Oakley)]')
)
box_pdf = os.path.join(dir_out, f"{file_prefix}_pps_RGEM_p_h1_xz_boxplot.pdf")
p.save(box_pdf, verbose=False, limitsize=False)
print(f"Saved RGE p(H_1 | x, z) boxplot to: {box_pdf}")

# %%

# =============================================================================
# Diagnostic-plot loop (moved from the H1x script): per-interim
# w_ratio vs pps_ratio_x scatter with the Gaussian-GLM smoother + per-item
# rho / R^2 annotations. This is the regression target for fit_interim_regress
# _endptx_on_wz, so the diagnostic naturally belongs to this script.
# =============================================================================

print("\nBuilding per-interim diagnostic plots (Gaussian fit on pps_ratio_x)...")
for interim_id in di['interim_id']:
    interim_id = int(interim_id)
    wa_pkl = os.path.join(dir_out, f"{file_prefix}_i{interim_id}_regression_training.pkl")
    if not os.path.exists(wa_pkl):
        continue
    wa = pd.read_pickle(wa_pkl)

    _ann = wa.groupby('item_label').apply(fit_interim_regress_H1x_on_wz_per_item_summary).reset_index()
    _ann['rho_label'] = _ann['rho'].apply(lambda r: f"ρ={r:.2f}")
    _ann['r2_label'] = _ann['r2'].apply(lambda r: f"R2={100 * r:.1f}%")

    _items = sorted(wa['item_label'].unique())
    _item_colors = dict(zip(_items, _futurama_palette(len(_items))))

    p = (
        ggplot(wa, aes(x='w_ratio', y='pps_ratio_x', colour='item_label'))
        + geom_point(alpha=0.3, size=0.5)
        + geom_smooth(method='glm', method_args={'family': 'gaussian'},
                      se=True, colour='black', size=0.8)
        + geom_text(_ann, aes(x='x_left', y='y_top', label='rho_label'),
                    ha='left', va='top', size=8, inherit_aes=False)
        + geom_text(_ann, aes(x='x_right', y='y_top', label='r2_label'),
                    ha='right', va='top', size=8, inherit_aes=False)
        + facet_wrap('~ item_label', ncol=4, scales='free')
        + scale_colour_manual(values=_item_colors)
        + theme_bw()
        + theme(
            figure_size=(15, 14),
            strip_background=element_rect(fill='none', colour='none'),
            panel_background=element_rect(fill='none', colour='none'),
            legend_position='none',
        )
        + labs(x=f'empirical endpoint summary w(z) at interim id {interim_id}',
               y=f'actual endpoint from posterior p(theta | x) draws at interim id {interim_id}')
    )
    pdf_path = os.path.join(dir_out, f"{file_prefix}_i{interim_id}_endpointratio_vs_wratio.pdf")
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"  saved diagnostic plot: {pdf_path}")

print("\n✓ regression-endptx-mquantile PPS + diagnostic plots complete for all interims")
