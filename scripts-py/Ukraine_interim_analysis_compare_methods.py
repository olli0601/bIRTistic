#!/usr/bin/env python3
"""
Ukraine interim analysis: cross-method comparison of weight quality + timing.

Loads the per-method performance artifacts written by the individual analysis
scripts and plots, per interim round, the importance-sampling diagnostics
ESS, ESS / particle and E(w^2) (the second-order weight moment) plus the
per-stage inference time (minutes), with one Futurama fill colour per method:

  - IS (reweight)       : ``..._pps_IS_perf_long.pkl``   (Ukraine_interim_analysis_with_IS_from_x.py)
  - IS (moment-match)   : ``..._pps_MM_perf_long.pkl``   (Ukraine_interim_analysis_with_IS_moment_matching_from_x.py)
  - SMC (resample-move) : ``..._pps_SMC_perf_long.pkl``  (Ukraine_interim_analysis_with_SMC_resample_from_x.py)

All three are full analyses over every interim. ESS and E(w^2) depend on the
particle budget N (IS/MM use N=4000 x-posterior draws; SMC uses K=128 moved
particles), so ESS / particle is the scale-free metric to compare across methods.

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/Ukraine_interim_analysis_compare_methods.py
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
    ggplot, aes, geom_boxplot, geom_col, geom_hline, facet_wrap,
    scale_fill_manual, scale_y_continuous, position_dodge,
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
dir_out = os.path.join(_sandbox, "py-ukraine-interim-compare-methods-260526")
os.makedirs(dir_out, exist_ok=True)
file_prefix = "pcm_1_interim"

METRIC_ORDER = ['time (min)', 'ESS / particle', 'E(w^2)']

# %%

# =============================================================================
# Load each method's long-form perf pkl (all interims) and concatenate.
# Each pkl has columns [interim_id, interim_month_year, s, metric, value, method].
# =============================================================================

is_long = pd.read_pickle(os.path.join(dir_is, f"{file_prefix}_pps_perf_long.pkl"))
mm_long = pd.read_pickle(os.path.join(dir_mm, f"{file_prefix}_pps_perf_long.pkl"))
smc_long = pd.read_pickle(os.path.join(dir_smc, f"{file_prefix}_pps_perf_long.pkl"))

cols = ['interim_month_year', 'metric', 'value', 'method']
compare = pd.concat([is_long[cols], mm_long[cols], smc_long[cols]], ignore_index=True)

interim_order = list(pd.Categorical(is_long['interim_month_year']).categories)
method_order = ['IS (reweight)', 'IS (moment-match)', 'SMC (resample-move)']
compare['interim_month_year'] = pd.Categorical(
    compare['interim_month_year'], categories=interim_order, ordered=True)
compare['method'] = pd.Categorical(compare['method'], categories=method_order, ordered=True)
compare['metric'] = pd.Categorical(
    compare['metric'].astype(str), categories=METRIC_ORDER, ordered=True)

csv_path = os.path.join(dir_out, f"{file_prefix}_compare_methods.csv")
compare.to_csv(csv_path, index=False)
print(f"Saved combined comparison table to: {csv_path}")
print(compare.groupby(['method', 'metric'], observed=True)['value'].mean().to_string())

# %%

# =============================================================================
# Plot: ESS / ESS-per-particle / E(w^2) / time per interim, methods dodged in
# Futurama colours (whisker outlier dots removed).
# =============================================================================

median = (
    compare.groupby(['method', 'interim_month_year', 'metric'], observed=True)['value']
    .median().reset_index()
)
# HMC: known wall-clock of 10 min per interim; no IS-style weight metrics.
hmc_time = pd.DataFrame({
    'method': 'HMC (refit)',
    'interim_month_year': interim_order,
    'metric': 'time (min)',
    'value': 10.0,
})
median = pd.concat([median, hmc_time], ignore_index=True)

perf_method_order = ['HMC (refit)', 'IS (reweight)', 'IS (moment-match)', 'SMC (resample-move)']
median['method'] = pd.Categorical(median['method'], categories=perf_method_order, ordered=True)
median['interim_month_year'] = pd.Categorical(
    median['interim_month_year'], categories=interim_order, ordered=True)
median['metric'] = pd.Categorical(
    median['metric'].astype(str), categories=METRIC_ORDER, ordered=True)
perf_method_colors = dict(zip(perf_method_order, _futurama_palette(len(perf_method_order))))

p = (
    ggplot(median, aes(x='interim_month_year', y='value', fill='method'))
    + geom_col(position=position_dodge(width=0.8), width=0.7)
    + facet_wrap('~ metric', ncol=1, scales='free_y')
    + scale_fill_manual(values=perf_method_colors)
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        legend_position='top',
        figure_size=(11, 8),
    )
    + labs(x='Interim date', y='', fill='method')
)
pdf_path = os.path.join(dir_out, f"{file_prefix}_compare_methods.pdf")
p.save(pdf_path, verbose=False, limitsize=False)
print(f"\nSaved cross-method comparison plot to: {pdf_path}")

# %%

# =============================================================================
# Per-item p(H_1 | x, z) by method (same layout as the IS-from-x boxplot,
# but methods dodged in Futurama fill colours -- includes the HMC refit).
# =============================================================================

box_paths = {
    'HMC (refit)': os.path.join(dir_hmc, f"{file_prefix}_pps_p_h1_xz_boxplot.pkl"),
    'IS (reweight)': os.path.join(dir_is, f"{file_prefix}_pps_p_h1_xz_boxplot.pkl"),
    'IS (moment-match)': os.path.join(dir_mm, f"{file_prefix}_pps_p_h1_xz_boxplot.pkl"),
    'SMC (resample-move)': os.path.join(dir_smc, f"{file_prefix}_pps_p_h1_xz_boxplot.pkl"),
}
box_frames = []
for m, bp in box_paths.items():
    bs = pd.read_pickle(bp).copy()
    bs['method'] = m
    box_frames.append(bs)
box_all = pd.concat(box_frames, ignore_index=True)

box_method_order = list(box_paths.keys())
box_all['method'] = pd.Categorical(box_all['method'], categories=box_method_order, ordered=True)
box_all['interim_month_year'] = pd.Categorical(
    box_all['interim_month_year'], categories=interim_order, ordered=True)
box_all['grp'] = box_all['interim_month_year'].astype(str) + '|' + box_all['method'].astype(str)
box_method_colors = dict(zip(box_method_order, _futurama_palette(len(box_method_order))))
pps_ProbH1_thresh = 0.89

pbox = (
    ggplot(box_all, aes(x='interim_month_year', fill='method'))
    + geom_boxplot(
        aes(ymin='q025', lower='q25', middle='q50', upper='q75', ymax='q975', group='grp'),
        stat='identity', position=position_dodge(width=0.8),
    )
    + geom_hline(yintercept=pps_ProbH1_thresh, colour='black', size=1.5)
    + facet_wrap('~ item_label_long', ncol=4)
    + scale_fill_manual(values=box_method_colors)
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
box_pdf = os.path.join(dir_out, f"{file_prefix}_compare_methods_p_h1_xz_boxplot.pdf")
pbox.save(box_pdf, verbose=False, limitsize=False)
print(f"Saved cross-method p(H_1 | x, z) boxplot to: {box_pdf}")

print("\n✓ method comparison complete")
