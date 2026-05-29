#!/usr/bin/env python3
"""
Ukraine interim analysis: cross-method comparison of weight quality + timing.

Loads the per-method performance artifacts written by the individual analysis
scripts and plots, per interim round, the importance-sampling diagnostics
ESS, ESS / particle and E(w^2) (the second-order weight moment) plus the
per-stage inference time (minutes), with one Futurama fill colour per method:

  - IS (reweight)       : ``..._pps_IS_perf_long.pkl``   (Ukraine_interim_analysis_with_IS_from_x.py)
  - IS (moment-match)   : ``..._pps_MM_perf_long.pkl``   (Ukraine_interim_analysis_with_IS_moment_matching.py)
  - SMC (resample-move) : ``..._pps_SMC_perf_long.pkl``  (Ukraine_interim_analysis_with_SMC_resample.py)

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
    ggplot, aes, geom_boxplot, facet_wrap, scale_fill_manual, position_dodge,
    theme_bw, theme, element_text, labs,
)
import warnings
warnings.filterwarnings('ignore')

from utils import _futurama_palette

print("✓ Imports successful")

# %%

_sandbox = "/Users/or105/sandbox/bIRTistic"
dir_is = os.path.join(_sandbox, "py-ukraine-interim-with-IS-260526")
dir_mm = os.path.join(_sandbox, "py-ukraine-interim-with-IS-moment-matching-260526")
dir_smc = os.path.join(_sandbox, "py-ukraine-interim-with-SMC-resample-260526")
dir_out = os.path.join(_sandbox, "py-ukraine-interim-compare-methods-260526")
os.makedirs(dir_out, exist_ok=True)
file_prefix = "pcm_1_interim"

METRIC_ORDER = ['ESS', 'ESS / particle', 'E(w^2)', 'time (min)']

# %%

# =============================================================================
# Load each method's long-form perf pkl (all interims) and concatenate.
# Each pkl has columns [interim_id, interim_month_year, s, metric, value, method].
# =============================================================================

is_long = pd.read_pickle(os.path.join(dir_is, f"{file_prefix}_pps_IS_perf_long.pkl"))
mm_long = pd.read_pickle(os.path.join(dir_mm, f"{file_prefix}_pps_MM_perf_long.pkl"))
smc_long = pd.read_pickle(os.path.join(dir_smc, f"{file_prefix}_pps_SMC_perf_long.pkl"))

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

method_colors = dict(zip(method_order, _futurama_palette(len(method_order))))

p = (
    ggplot(compare, aes(x='interim_month_year', y='value', fill='method'))
    + geom_boxplot(outlier_alpha=0, position=position_dodge(width=0.8))
    + facet_wrap('~ metric', ncol=1, scales='free_y')
    + scale_fill_manual(values=method_colors)
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        legend_position='bottom',
        figure_size=(11, 13),
    )
    + labs(x='Interim', y='value (boxes across S importance-sampling draws)',
           fill='method',
           title='Weight quality + timing by method, per interim')
)
pdf_path = os.path.join(dir_out, f"{file_prefix}_compare_methods.pdf")
p.save(pdf_path, verbose=False, limitsize=False)
print(f"\nSaved cross-method comparison plot to: {pdf_path}")

print("\n✓ method comparison complete")
