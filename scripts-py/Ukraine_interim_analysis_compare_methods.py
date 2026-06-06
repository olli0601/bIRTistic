#!/usr/bin/env python3
"""
Ukraine interim analysis: cross-method comparison of weight quality + timing.

Loads the per-method performance artifacts written by the individual analysis
scripts and plots, per interim round, the importance-sampling diagnostics
ESS, ESS / particle and E(w^2) (the second-order weight moment) plus the
per-stage inference time (minutes), with one Futurama fill colour per method:

  - IS (reweight)            : ``..._pps_IS_perf_long.pkl``   (Ukraine_interim_analysis_with_IS_from_x.py)
  - IS (moment-match)        : ``..._pps_MM_perf_long.pkl``   (Ukraine_interim_analysis_with_IS_moment_matching_from_x.py)
  - SMC (resample-move)      : ``..._pps_SMC_perf_long.pkl``  (Ukraine_interim_analysis_with_SMC_resample_from_x.py)
  - Regression (Strong-Oakley): ``..._pps_RG_perf_long.pkl``   (Ukraine_interim_analysis_regression.py)

The first three are IS-style methods with weight metrics (ESS, ESS/particle,
E(w^2)); the regression-based method has no weights and reports per-item rho /
R^2 instead. Only the per-stage inference time (minutes) is common across all
four, so the bar plot shows time for all methods and weight metrics for the IS
family only.

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/Ukraine_interim_analysis_dc_methods.py
"""

# %%

import os
import sys
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
    ggplot, aes, geom_boxplot, geom_col, geom_hline, geom_text, facet_wrap,
    scale_fill_manual, scale_y_continuous, scale_y_sqrt, position_dodge,
    theme_bw, theme, element_text, labs,
)
import warnings
warnings.filterwarnings('ignore')

from utils import _futurama_palette

print("✓ Imports successful")

# %%

_sandbox = "/Users/or105/sandbox/bIRTistic"
dir_hmc = os.path.join(_sandbox, "py-ukraine-interim-260526")
dir_is = os.path.join(_sandbox, "py-ukraine-interim-with-IS-260526")
dir_mm = os.path.join(_sandbox, "py-ukraine-interim-with-IS-moment-matching-260526")
dir_smc = os.path.join(_sandbox, "py-ukraine-interim-with-SMC-resample-260526")
dir_rg = os.path.join(_sandbox, "py-ukraine-interim-with-regression-on-H1x-wz-260601")
dir_rge = os.path.join(_sandbox, "py-ukraine-interim-with-regression-on-endptx-wz-260601")
dir_out = os.path.join(_sandbox, "py-ukraine-interim-compare-methods-260526")
os.makedirs(dir_out, exist_ok=True)
file_prefix = "pcm_1_interim"

# %%

# =============================================================================
# Load each method's long-form perf pkl (all interims) and concatenate.
# Each pkl has columns [interim_id, interim_month_year, s, metric, value, method].
# =============================================================================

is_long = pd.read_pickle(os.path.join(dir_is, f"{file_prefix}_pps_perf_long.pkl"))
mm_long = pd.read_pickle(os.path.join(dir_mm, f"{file_prefix}_pps_perf_long.pkl"))
smc_long = pd.read_pickle(os.path.join(dir_smc, f"{file_prefix}_pps_perf_long.pkl"))
rg_long = pd.read_pickle(os.path.join(dir_rg, f"{file_prefix}_pps_RG_perf_long.pkl"))
rge_long = pd.read_pickle(os.path.join(dir_rge, f"{file_prefix}_pps_RGE_perf_long.pkl"))

cols = ['interim_month_year', 'metric', 'value', 'method']
dc = pd.concat(
    [is_long[cols], mm_long[cols], smc_long[cols], rg_long[cols], rge_long[cols]],
    ignore_index=True,
)

interim_order = list(pd.Categorical(is_long['interim_month_year']).categories)
method_order = ['IS (reweight)', 'IS (moment-match)', 'SMC (resample-move)',
                'Regression (Strong-Oakley)', 'Regression endpt (Strong-Oakley)']
dc['interim_month_year'] = pd.Categorical(
    dc['interim_month_year'], categories=interim_order, ordered=True)
dc['method'] = pd.Categorical(dc['method'], categories=method_order, ordered=True)
# Only 'time (min)' is shared across all 4 methods; the section below recasts
# the metric column to a single-level Categorical for the timing-only bar plot.

csv_path = os.path.join(dir_out, f"{file_prefix}_dc_methods.csv")
dc.to_csv(csv_path, index=False)
print(f"Saved combined comparison table to: {csv_path}")

# %%

# =============================================================================
# Plot: time per interim
# =============================================================================

dc2 = (
    dc[dc['metric'] == 'time (min)']
    .groupby(['method', 'interim_month_year', 'metric'], observed=True)['value']
    .median()
    .reset_index()
)
# HMC: known wall-clock of 10 min per interim; no IS-style weight metrics.
hmc_time = pd.DataFrame({
    'method': 'HMC (refit)',
    'interim_month_year': interim_order,
    'metric': 'time (min)',
    'value': 10.0,
})
dc2 = pd.concat([dc2, hmc_time], ignore_index=True)

# Rename methods to match the box_paths labels used in the boxplot section.
_method_rename = {
    'HMC (refit)': 'HMC',
    'Regression (Strong-Oakley)': 'Regression H1(x) on w(z)',
    'Regression endpt (Strong-Oakley)': 'Regression endpt(x) on w(z)',
}
dc2['method'] = dc2['method'].astype(str).map(lambda m: _method_rename.get(m, m))

perf_method_order = ['HMC', 'IS (reweight)', 'IS (moment-match)',
                     'SMC (resample-move)', 'Regression H1(x) on w(z)',
                     'Regression endpt(x) on w(z)']
dc2['method'] = pd.Categorical(dc2['method'], categories=perf_method_order, ordered=True)
dc2['interim_month_year'] = pd.Categorical(
    dc2['interim_month_year'], categories=interim_order, ordered=True)
dc2['value_label'] = dc2['value'].map(lambda v: f"{v:.2f} mins")
perf_method_colors = dict(zip(perf_method_order, _futurama_palette(len(perf_method_order))))

p = (
    ggplot(dc2, aes(x='interim_month_year', y='value', fill='method'))
    + geom_col(position=position_dodge(width=0.8), width=0.7)
    + geom_text(
        aes(label='value_label'),
        position=position_dodge(width=0.8),
        size=6, angle=30, ha='left', va='bottom',
    )
    + scale_fill_manual(values=perf_method_colors)
    + scale_y_sqrt(breaks=[0.1, 0.5, 1, 5, 10, 30, 60],
                   labels=['0.1', '0.5', '1', '5', '10', '30', '60'],
                   expand=(0, 0, 0.1, 0))
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        legend_position='top',
        figure_size=(13, 6),
    )
    + labs(x='Interim date', y='time (mins)', fill='method')
)
pdf_path = os.path.join(dir_out, f"{file_prefix}_compare_methods_timing.pdf")
p.save(pdf_path, verbose=False, limitsize=False)
print(f"\nSaved cross-method comparison plot to: {pdf_path}")


# %%

# =============================================================================
# Plot: ESS / particle + E(w^2) per interim for the IS-style methods only
# (regression methods have no weights; HMC has no IS diagnostics).
# =============================================================================

is_methods = ['IS (reweight)', 'IS (moment-match)', 'SMC (resample-move)']
ess_metrics = ['ESS / particle', 'E(w^2)']
dc3 = (
    dc[dc['method'].isin(is_methods) & dc['metric'].isin(ess_metrics)]
    .groupby(['method', 'interim_month_year', 'metric'], observed=True)['value']
    .median()
    .reset_index()
)
dc3['method'] = pd.Categorical(dc3['method'], categories=is_methods, ordered=True)
dc3['interim_month_year'] = pd.Categorical(
    dc3['interim_month_year'], categories=interim_order, ordered=True)
dc3['metric'] = pd.Categorical(
    dc3['metric'].astype(str), categories=ess_metrics, ordered=True)
is_method_colors = {m: perf_method_colors[m] for m in is_methods}

p = (
    ggplot(dc3, aes(x='interim_month_year', y='value', fill='method'))
    + geom_col(position=position_dodge(width=0.8), width=0.7)
    + facet_wrap('~ metric', ncol=1, scales='free_y')
    + scale_fill_manual(values=is_method_colors)
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        legend_position='top',
        figure_size=(11, 7),
    )
    + labs(x='Interim date', y='', fill='method')
)
pdf_path = os.path.join(dir_out, f"{file_prefix}_compare_methods_ess.pdf")
p.save(pdf_path, verbose=False, limitsize=False)
print(f"\nSaved IS-method ESS/particle + E(w^2) plot to: {pdf_path}")

# %%

# =============================================================================
# Per-item p(H_1 | x, z) by method (same layout as the IS-from-x boxplot,
# but methods dodged in Futurama fill colours -- includes the HMC refit).
# =============================================================================

box_paths = {
    'HMC': os.path.join(dir_hmc, f"{file_prefix}_pps_p_h1_xz_boxplot.pkl"),
    'IS (reweight)': os.path.join(dir_is, f"{file_prefix}_pps_p_h1_xz_boxplot.pkl"),
    'IS (moment-match)': os.path.join(dir_mm, f"{file_prefix}_pps_p_h1_xz_boxplot.pkl"),
    'SMC (resample-move)': os.path.join(dir_smc, f"{file_prefix}_pps_p_h1_xz_boxplot.pkl"),
    'Regression H1(x) on w(z)': os.path.join(dir_rg, f"{file_prefix}_pps_RG_p_h1_xz_boxplot.pkl"),
    'Regression endpt(x) on w(z)': os.path.join(dir_rge, f"{file_prefix}_pps_RGE_p_h1_xz_boxplot.pkl"),
}
tmp = []
for m, bp in box_paths.items():
    bs = pd.read_pickle(bp).copy()
    bs['method'] = m
    tmp.append(bs)
box_all = pd.concat(tmp, ignore_index=True)

box_method_order = list(box_paths.keys())
box_all['method'] = pd.Categorical(box_all['method'], categories=box_method_order, ordered=True)
box_all['interim_month_year'] = pd.Categorical(
    box_all['interim_month_year'], categories=interim_order, ordered=True)
box_all['grp'] = box_all['interim_month_year'].astype(str) + '|' + box_all['method'].astype(str)
box_method_colors = dict(zip(box_method_order, _futurama_palette(len(box_method_order))))
pps_ProbH1_thresh = 0.89

def _boxplot_subset(methods, suffix):
    tmp = box_all[box_all['method'].isin(methods)].copy()
    tmp['method'] = pd.Categorical(
        tmp['method'].astype(str), categories=methods, ordered=True,
    )
    sub_colors = {m: box_method_colors[m] for m in methods}
    p = (
        ggplot(tmp, aes(x='interim_month_year', fill='method'))
        + geom_boxplot(
            aes(ymin='q025', lower='q25', middle='q50', upper='q75', ymax='q975', group='grp'),
            stat='identity', position=position_dodge(width=0.8),
        )
        + geom_hline(yintercept=pps_ProbH1_thresh, colour='black', size=1.5)
        + facet_wrap('~ item_label_long', ncol=4)
        + scale_fill_manual(values=sub_colors)
        + scale_y_continuous(
            limits=[0, 1],
            breaks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            labels=['0%', '20%', '40%', '60%', '80%', '100%'],
        )
        + theme_bw()
        + theme(
            axis_text_x=element_text(angle=45, vjust=1, hjust=1),
            legend_position='top',
            figure_size=(15, 15),
            strip_text_y=element_text(angle=0),
        )
        + labs(x='Interim', y='p(H_1 | x, z)', fill='method')
    )
    box_pdf = os.path.join(dir_out, f"{file_prefix}_dc_methods_p_h1_xz_boxplot_{suffix}.pdf")
    p.save(box_pdf, verbose=False, limitsize=False)
    print(f"Saved cross-method p(H_1 | x, z) boxplot to: {box_pdf}")


_boxplot_subset(
    ['HMC', 'IS (reweight)', 'IS (moment-match)', 'SMC (resample-move)'],
    'IS_methods',
)
_boxplot_subset(
    ['HMC', 'Regression H1(x) on w(z)', 'Regression endpt(x) on w(z)'],
    'regression_methods',
)

print("\n✓ method comparison complete")
