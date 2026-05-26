#!/usr/bin/env python3
"""
Investigate ADVI vs HMC for partial credit model on Ukraine data

Date: 2026-04-30

This script applies the partial credit model (ncats v260413) to the Ukraine data
for comparison with other IRT models. It mirrors the Colombia analysis script
(``scripts-py/Colombia_analysis_partial_credit_ADVI-HMC-pyro.py``).

Ported from: ``scripts-R/Colombia_analysis_for_HGpaper_pcm_v251224.Rmd``
            (Ukraine sections starting around line 498).

Usage:
    # Option 1: Run as script
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/Ukraine_analysis_partial_credit_ADVI-HMC-pyro.py

    # Option 2: Run interactively line-by-line
    cd /Users/or105/git/bIRTistic
    pixi run python
    >>> exec(open('scripts-py/__init__.py').read())  # Setup paths
    >>> exec(open('scripts-py/Ukraine_analysis_partial_credit_ADVI-HMC-pyro.py').read())
"""

# %%

# =============================================================================
# Setup and Imports
# =============================================================================

import sys
import os
import subprocess
from pathlib import Path

# Add python module to path (handles both script and interactive mode)
try:
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
except NameError:
    # Interactive mode: assume we're in project root
    project_root = Path.cwd()
    if project_root.name == 'scripts-py':
        project_root = project_root.parent

python_path = str(project_root / 'python')
if python_path not in sys.path:
    sys.path.insert(0, python_path)

import pandas as pd
import numpy as np
import arviz as az
import xarray as xr
from plotnine import (
    ggplot, aes, geom_col, geom_errorbar, geom_density, geom_boxplot, geom_text, facet_wrap, facet_grid,
    scale_fill_manual, scale_color_manual, theme_minimal, theme_bw, theme, element_text,
    element_blank, labs, position_dodge, coord_flip, ggsave, scale_y_continuous
)
from ggsci import scale_fill_futurama
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from PIL import Image
import warnings
warnings.filterwarnings('ignore')
# Import bIRTistic functions
from data_loading import read_data_ukraine
from fit_partial_credit_model_ncats_stanadvi import fit_partial_credit_model_ncats_stanadvi
from fit_partial_credit_model_ncats_stanhmc import fit_partial_credit_model_ncats_stanhmc
from fit_partial_credit_model_ncats_pyrosvi import fit_partial_credit_model_ncats_pyrosvi
from get_endpoints import get_endpoints
from utils import _summarize_ordered_prob_quantiles

print("✓ Imports successful")

# %%

# =============================================================================
# Configuration
# =============================================================================

# Set random seed
np.random.seed(42)
seed = 123

# Define directories
dir_data = "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/2025_project_Hope_Groups/data"
file_data = os.path.join(dir_data, "Ukraine_Hope_Groups_Baseline_Endline_Wide_Aug6.csv")

# Output directories
dir_out_pcm = "/Users/or105/sandbox/bIRTistic/py-ukraine-partial_credit-260430-vanilla-advi-hmc-pyro"
dir_logs_pcm = os.path.join(dir_out_pcm, "logs")

# Create output directories
os.makedirs(dir_out_pcm, exist_ok=True)
os.makedirs(dir_logs_pcm, exist_ok=True)

output_file_prefix = os.path.join(dir_out_pcm, "pcm_1")
output_file_prefix_hmc = os.path.join(dir_out_pcm, "pcm_1_hmc")
output_file_prefix_advi = os.path.join(dir_out_pcm, "pcm_1_advi")
svi_algorithm_autodiagnormal = "AutoDiagonalNormal"
svi_algorithm_autolaplaceapproximation = 'AutoLaplaceApproximation'
svi_algorithm_automultivariatenormal = 'AutoMultivariateNormal'
svi_algorithm_autoiafnormal = 'AutoIAFNormal'
output_file_svi_autodiagnormal = os.path.join(dir_out_pcm, "pcm_1_svi_autodiagnormal")
output_file_svi_autolaplaceapproximation = os.path.join(dir_out_pcm, "pcm_1_svi_autolaplaceapproximation")
output_file_svi_automultivariatenormal = os.path.join(dir_out_pcm, "pcm_1_svi_automultivariatenormal")
output_file_svi_autoiafnormal = os.path.join(dir_out_pcm, "pcm_1_svi_autoiafnormal")

# Endpoint summarization: Ukraine uses categorical_threshold=2 (Colombia uses 3)
categorical_threshold = 2

print(f"Data file: {file_data}")
print(f"Output directory: {dir_out_pcm}")

# %%

# =============================================================================
# Load and Preprocess Data
# =============================================================================

# Read Ukraine data
print("\nLoading Ukraine data...")
result = read_data_ukraine(file_data)
dp_ukr = result['dp'].copy()
dit_ukr = result['dit'].copy()
dmeta_ukr = result['dmeta'].copy()

print(f"\nLoaded data:")
print(f"  dp: {len(dp_ukr):,} rows")
print(f"  dit: {len(dit_ukr):,} rows")
print(f"  dmeta: {len(dmeta_ukr):,} rows")

# Preprocess data
print("\nPreprocessing data...")
dp1_ukr = dp_ukr[~dp_ukr['item_label'].str.contains('agg')].copy()
dp1_ukr['y_stan'] = dp1_ukr['y'] + 1

# Merge item_type to avoid indexing issues
dp1_ukr = dp1_ukr.merge(dit_ukr[['item_label', 'item_type']], on='item_label', how='left')

# Create item_time mapping
item_time_df = (
    dp1_ukr[['item_type', 'item_label', 'time']]
    .drop_duplicates()
    .sort_values(['item_type', 'time', 'item_label'])
    .reset_index(drop=True)
)
item_time_df['item_time_id'] = item_time_df.groupby('item_type').cumcount() + 1


# Merge item_time_id
dp1_ukr = dp1_ukr.merge(item_time_df, on=['item_label', 'time', 'item_type'], how='left')

# Add item_type_id
dp1_ukr = dp1_ukr.merge(
    dit_ukr[['item_type', 'item_type_id']].drop_duplicates(),
    on='item_type',
    how='left'
)

# Sort and create oid
dp1_ukr = dp1_ukr.sort_values(['item_type_id', 'pid', 'time', 'item_label']).reset_index(drop=True)
dp1_ukr['oid'] = range(1, len(dp1_ukr) + 1)

# Create oidt (observation id within item type)
dp1_ukr['oidt'] = dp1_ukr.groupby('item_type').cumcount() + 1

print(f"\nPreprocessed data:")
print(f"  Observations: {len(dp1_ukr):,}")
print(f"  Participants: {dp1_ukr['pid'].nunique()}")
print(f"  Items: {dp1_ukr['item_label'].nunique()}")
print(f"  Item types: {sorted(dp1_ukr['item_type'].unique())}")

# Save data for model fitting
tmp = os.path.join(dir_out_pcm, "pcm_1_data.pkl")
print(f"\nSaved preprocessed data to: {tmp}")
pd.to_pickle({'dp1': dp1_ukr, 'dit': dit_ukr, 'dmeta': dmeta_ukr}, tmp)

# %%

# =============================================================================
# Model Fitting - HMC
# =============================================================================

result_hmc = fit_partial_credit_model_ncats_stanhmc(
    dit_ukr,
    dp1_ukr,
    output_file_prefix=output_file_prefix_hmc,
    stan_file="/Users/or105/git/bIRTistic/src/stan/partial_credit_model_ncats_v260413.stan",
    chains=4,
    parallel_chains=4,
    threads_per_chain=1,
    iter_warmup=500,
    iter_sampling=1000,
    seed=seed,
    x_formula="~ time - 1",
    resume=True,
    with_core_analyses=True,
    with_additional_analyses=True,
    show_messages=True
)

print(f"\n✓ HMC fitting complete")
print(f"  Good chains: {result_hmc['good_chains']}")
print(f"  Draws file: {output_file_prefix_hmc}_draws.zarr")

# %%

# =============================================================================
# Model Fitting - ADVI
# =============================================================================

result_advi = fit_partial_credit_model_ncats_stanadvi(
    dit_ukr,
    dp1_ukr,
    output_file_prefix=output_file_prefix_advi,
    stan_file="/Users/or105/git/bIRTistic/src/stan/partial_credit_model_ncats_v260413.stan",
    iter=10000,
    grad_samples=1,
    elbo_samples=100,
    output_samples=4000,
    seed=seed,
    x_formula="~ time - 1",
    resume=True,
    with_core_analyses=True,
    with_additional_analyses=True,
    show_messages=True
)

print(f"\n✓ ADVI fitting complete")
print(f"  Draws file: {output_file_prefix_advi}_draws.zarr")

# %%

# =============================================================================
# Model Fitting - SVI AutoDiagonalNormal
# =============================================================================

result_svi_autodiagnormal = fit_partial_credit_model_ncats_pyrosvi(
    dit_ukr,
    dp1_ukr,
    output_file_prefix=output_file_svi_autodiagnormal,
    algorithm=svi_algorithm_autodiagnormal,
    lr=0.01,
    num_steps=10000,
    output_samples=4000,
    seed=seed,
    x_formula="~ time - 1",
    resume=True,
    with_core_analyses=True,
    with_additional_analyses=True,
)

print(f"\n✓ SVI fitting complete")
print(f"  Algorithm: {result_svi_autodiagnormal['algorithm']}")
print(f"  Draws file: {output_file_svi_autodiagnormal}_draws.zarr")

# %%

# =============================================================================
# Model Fitting - SVI AutoLaplaceApproximation
# =============================================================================

result_svi_autolaplaceapproximation = fit_partial_credit_model_ncats_pyrosvi(
    dit_ukr,
    dp1_ukr,
    output_file_prefix=output_file_svi_autolaplaceapproximation,
    algorithm=svi_algorithm_autolaplaceapproximation,
    lr=0.01,
    num_steps=10000,
    output_samples=4000,
    seed=seed,
    x_formula="~ time - 1",
    resume=True,
    with_core_analyses=True,
    with_additional_analyses=True,
)

print(f"\n✓ SVI fitting complete")
print(f"  Algorithm: {result_svi_autolaplaceapproximation['algorithm']}")
print(f"  Draws file: {output_file_svi_autolaplaceapproximation}_draws.zarr")

# %%

# =============================================================================
# Model Fitting - SVI AutoMultivariateNormal
# =============================================================================

result_svi_automultivariatenormal = fit_partial_credit_model_ncats_pyrosvi(
    dit_ukr,
    dp1_ukr,
    output_file_prefix=output_file_svi_automultivariatenormal,
    algorithm=svi_algorithm_automultivariatenormal,
    lr=0.01,
    num_steps=10000,
    output_samples=4000,
    seed=seed,
    x_formula="~ time - 1",
    resume=True,
    with_core_analyses=True,
    with_additional_analyses=True,
)

print(f"\n✓ SVI fitting complete")
print(f"  Algorithm: {result_svi_automultivariatenormal['algorithm']}")
print(f"  Draws file: {output_file_svi_automultivariatenormal}_draws.zarr")

# %%

# =============================================================================
# Model Fitting - SVI Auto IAF Normal
# =============================================================================

result_svi_autoiafnormal = fit_partial_credit_model_ncats_pyrosvi(
    dit_ukr,
    dp1_ukr,
    output_file_prefix=output_file_svi_autoiafnormal,
    algorithm=svi_algorithm_autoiafnormal,
    lr=0.01,
    num_steps=10000,
    output_samples=4000,
    seed=seed,
    x_formula="~ time - 1",
    resume=True,
    with_core_analyses=True,
    with_additional_analyses=True,
)

print(f"\n✓ SVI fitting complete")
print(f"  Algorithm: {result_svi_autoiafnormal['algorithm']}")
print(f"  Draws file: {output_file_svi_autoiafnormal}_draws.zarr")


# %%

# =============================================================================
# Compare HMC vs ADVI Results - Compute Endpoints
# =============================================================================

# Load preprocessed data
tmp = pd.read_pickle(f"{output_file_prefix}_data.pkl")
dp1_ukr = tmp['dp1']
dit_ukr = tmp['dit']
dmeta_ukr = tmp['dmeta']

# Collect all fitted methods in one place and compute endpoints uniformly.
method_cfg = {
    'stan_hmc': {
        'draws_file': f"{output_file_prefix_hmc}_draws.zarr",
        'param_name': None,
    },
    'stan_advi': {
        'draws_file': f"{output_file_prefix_advi}_draws.zarr",
        'param_name': 'ordered_prob_by_cat_qu_pr',
    },
    'numpyro_svi_autodiagnormal': {
        'draws_file': f"{output_file_svi_autodiagnormal}_draws.zarr",
        'param_name': 'ordered_prob_by_cat_qu_pr',
    },
    'numpyro_svi_autolaplaceapproximation': {
        'draws_file': f"{output_file_svi_autolaplaceapproximation}_draws.zarr",
        'param_name': 'ordered_prob_by_cat_qu_pr',
    },
    'numpyro_svi_automultivariatenormal': {
        'draws_file': f"{output_file_svi_automultivariatenormal}_draws.zarr",
        'param_name': 'ordered_prob_by_cat_qu_pr',
    },
    'numpyro_svi_autoiafnormal': {
        'draws_file': f"{output_file_svi_autoiafnormal}_draws.zarr",
        'param_name': 'ordered_prob_by_cat_qu_pr',
    },
}

endpoints_by_method = {}
for method, cfg in method_cfg.items():
    print(f"\nComputing {method} endpoints...")
    ep_kwargs = {
        'dp1': dp1_ukr,
        'dit': dit_ukr,
        'draws_file': cfg['draws_file'],
        'categorical_threshold': categorical_threshold,
        'endpoint_type': 'items',
    }
    if cfg['param_name'] is not None:
        ep_kwargs['param_name'] = cfg['param_name']

    pos = get_endpoints(**ep_kwargs)
    pos['method'] = method
    endpoints_by_method[method] = pos

# Keep legacy variable names used later in the script.
endpoints_hmc = endpoints_by_method['stan_hmc']
endpoints_advi = endpoints_by_method['stan_advi']
endpoints_svi_autodiagnormal = endpoints_by_method['numpyro_svi_autodiagnormal']

print(f"\n✓ Computed endpoints")
for method, pos in endpoints_by_method.items():
    print(f"  {method}: {len(pos)} rows")

# Reuse posterior draws for all downstream comparisons.
method_order = list(method_cfg.keys())
posterior_by_method = {
    'stan_hmc': result_hmc['draws'].posterior,
    'stan_advi': result_advi['draws'].posterior,
    'numpyro_svi_autodiagnormal': result_svi_autodiagnormal['draws'].posterior,
    'numpyro_svi_autolaplaceapproximation': result_svi_autolaplaceapproximation['draws'].posterior,
    'numpyro_svi_automultivariatenormal': result_svi_automultivariatenormal['draws'].posterior,
    'numpyro_svi_autoiafnormal': result_svi_autoiafnormal['draws'].posterior,
}

# %%

# =============================================================================
# Create Comparison Table
# =============================================================================

# Combine all methods
pos = pd.concat(endpoints_by_method.values(), ignore_index=True)

# Pivot to wide format for comparison
pos = pos.pivot_table(
    index=['item_type_id','item_type','item_label','item_label_short','group_label','group_label_long','item_high_label','variable'],
    columns='method',
    values=['median', 'q_lower', 'iqr_lower', 'iqr_upper', 'q_upper']
).reset_index()

# Flatten column names from multi-index
pos.columns = ['_'.join(col).strip('_') if isinstance(col, tuple) else col
                          for col in pos.columns]

# Difference in medians vs stan_hmc reference, sort by largest absolute gap.
ref_col = 'median_stan_hmc'
diff_cols = []
if ref_col in pos.columns:
    for m in method_order:
        if m == 'stan_hmc':
            continue
        col = f'median_{m}'
        if col in pos.columns:
            dcol = f'median_diff_vs_hmc_{m}'
            pos[dcol] = pos[col] - pos[ref_col]
            diff_cols.append(dcol)
    if diff_cols:
        pos['max_abs_median_diff_vs_hmc'] = pos[diff_cols].abs().max(axis=1)
        pos = pos.sort_values('max_abs_median_diff_vs_hmc', ascending=False)

# Save comparison
tmp = os.path.join(dir_out_pcm, "comparison_endpoints_all_methods.csv")
pos.to_csv(tmp, index=False)
print(f"\nSaved comparison to: {tmp}")

# %%

# =============================================================================
# Visualization: Endpoint Comparison by Item
# =============================================================================

# Filter for difference variable and add item_label_short
# Merge with dit to get item_label_short for better labels
pos = pd.concat(endpoints_by_method.values(), ignore_index=True)
endpoints_plot = pos[pos['variable'] == 'diff'].copy()

# Shift estimates relative to stan_hmc median per unique (item, group, high-label, variable).
tmp = ['item_type_id', 'item_type', 'item_label', 'item_label_short',
       'group_label', 'group_label_long', 'item_high_label', 'variable']
# Fill NaN to allow merge equality on item_label_short.
for k in tmp:
    if endpoints_plot[k].isna().any():
        endpoints_plot[k] = endpoints_plot[k].fillna('__NA__')
hmc_ref = (
    endpoints_plot.loc[endpoints_plot['method'] == 'stan_hmc', tmp + ['median']]
    .rename(columns={'median': 'median_hmc'})
)
assert not hmc_ref.duplicated(subset=tmp).any(), \
    "stan_hmc rows not unique on merge keys — check get_endpoints output"
endpoints_plot = endpoints_plot.merge(hmc_ref, on=tmp, how='left', validate='many_to_one')
endpoints_plot['median'] = endpoints_plot['median'] - endpoints_plot['median_hmc']
endpoints_plot['iqr_lower'] = endpoints_plot['iqr_lower'] - endpoints_plot['median_hmc']
endpoints_plot['iqr_upper'] = endpoints_plot['iqr_upper'] - endpoints_plot['median_hmc']

# Drop stan_hmc itself — it's the reference, always zero.
endpoints_plot = endpoints_plot[endpoints_plot['method'] != 'stan_hmc'].copy()

# Create composite label for facets
endpoints_plot['facet_label'] = (
    endpoints_plot['group_label_long'] +
    np.where(endpoints_plot['item_label_short'] != '__NA__',
             '---' + endpoints_plot['item_label_short'] ,
             ''
             ) + '\n' +
    np.where(endpoints_plot['item_type'] == 'categorical',
             'Diff in probability per week vs stan_hmc\n(Baseline - Endline)',
             'Diff in mean days per week vs stan_hmc\n(Baseline - Endline)'
             )
)

# Create plot with individual items as facets (matching R version)
print("\nCreating endpoint difference plot...")

p_diff = (
    ggplot(endpoints_plot, aes(x='method', y='median', fill='method')) +
    geom_col(alpha=0.8, width=0.7) +
    geom_errorbar(aes(ymin='iqr_lower', ymax='iqr_upper'), width=0.3) +
    facet_wrap('~facet_label', scales='free', ncol=3) +
    scale_fill_futurama() +
    theme_bw() +
    theme(
        axis_text_x=element_text(angle=45, hjust=1, size=8),
        axis_title_x=element_blank(),
        strip_text=element_text(size=7, weight='bold'),
        strip_background=element_blank(),
        legend_position='top',
        figure_size=(14, 18)
    ) +
    labs(
        title='Endpoint Comparison: Difference vs stan_hmc',
        y='Estimate − stan_hmc median',
        fill='Method'
    )
)

# Save plot
tmp = os.path.join(dir_out_pcm, 'comparison_effect_size_difference.pdf')
ggsave(p_diff, filename=tmp, width=14, height=30, limitsize=False)
print(f"Saved plot to: {tmp}")

# %%

# Create ratio plot
print("Creating endpoint ratio plot...")
endpoints_plot = pos[pos['variable'] == 'ratio'].copy()

# Shift estimates relative to stan_hmc median per unique (item, group, high-label, variable).
tmp = ['item_type_id', 'item_type', 'item_label', 'item_label_short',
       'group_label', 'group_label_long', 'item_high_label', 'variable']
for k in tmp:
    if endpoints_plot[k].isna().any():
        endpoints_plot[k] = endpoints_plot[k].fillna('__NA__')
hmc_ref = (
    endpoints_plot.loc[endpoints_plot['method'] == 'stan_hmc', tmp + ['median']]
    .rename(columns={'median': 'median_hmc'})
)
assert not hmc_ref.duplicated(subset=tmp).any(), \
    "stan_hmc rows not unique on merge keys — check get_endpoints output"
endpoints_plot = endpoints_plot.merge(hmc_ref, on=tmp, how='left', validate='many_to_one')
endpoints_plot['median'] = endpoints_plot['median'] - endpoints_plot['median_hmc']
endpoints_plot['iqr_lower'] = endpoints_plot['iqr_lower'] - endpoints_plot['median_hmc']
endpoints_plot['iqr_upper'] = endpoints_plot['iqr_upper'] - endpoints_plot['median_hmc']

# Drop stan_hmc itself — it's the reference, always zero.
endpoints_plot = endpoints_plot[endpoints_plot['method'] != 'stan_hmc'].copy()

endpoints_plot['facet_label'] = (
    endpoints_plot['group_label_long'] +
    np.where(endpoints_plot['item_label_short'] != '__NA__',
             '---' + endpoints_plot['item_label_short'] ,
             ''
             ) + '\n' +
    np.where(endpoints_plot['item_type'] == 'categorical',
             'Ratio in probability per week vs stan_hmc\n(1- Baseline/Endline)',
             'Ratio in mean days per week vs stan_hmc\n(1- Baseline/Endline)'
             )
)

# Create plot with individual items as facets (matching R version)
print("\nCreating endpoint ratio plot...")

p_ratio = (
    ggplot(endpoints_plot, aes(x='method', y='median', fill='method')) +
    geom_col(alpha=0.8, width=0.7) +
    geom_errorbar(aes(ymin='iqr_lower', ymax='iqr_upper'), width=0.3) +
    facet_wrap('~facet_label', scales='free', ncol=3) +
    scale_fill_futurama() +
    theme_bw() +
    theme(
        axis_text_x=element_text(angle=45, hjust=1, size=8),
        axis_title_x=element_blank(),
        strip_text=element_text(size=7, weight='bold'),
        strip_background=element_blank(),
        legend_position='top',
        figure_size=(14, 18)
    ) +
    labs(
        title='Endpoint Comparison: Ratio vs stan_hmc',
        y='Estimate − stan_hmc median',
        fill='Method'
    )
)

tmp = os.path.join(dir_out_pcm, 'comparison_effect_size_ratio.pdf')
ggsave(p_ratio, filename=tmp, width=14, height=30, limitsize=False)
print(f"Saved plot to: {tmp}")

# %%

# =============================================================================
# Compare ordered_prob_by_cat_qu_fit parameters
# =============================================================================

print("\n" + "="*70)
print("COMPARING PARTIAL CREDIT PROBABILITIES")
print("="*70)

pos = []
for method in method_order:
    po = posterior_by_method[method]
    tmp = 'ordered_prob_by_cat_qu_fit'
    if tmp not in po:
        tmp = 'ordered_prob_by_cat_qu_pr'
    tmp2 = _summarize_ordered_prob_quantiles(po[tmp].values, dp1_ukr, dit_ukr)
    tmp2['method'] = method
    pos.append(tmp2)
pos = pd.concat(pos, ignore_index=True)

pos = pos.pivot_table(
    index=['cq_id', 'item_type_id', 'item_time_id', 'y', 'item_label', 'time_label',
           'group_label_long', 'item_label_short', 'endpoint_measure'],
    columns='method',
    values='median',
    aggfunc='first'
).reset_index()

tmp = [m for m in method_order if m in pos.columns]
pos['abs_median_range'] = pos[tmp].max(axis=1) - pos[tmp].min(axis=1)
pos = pos.sort_values('abs_median_range', ascending=False)

print(f"\nTop 10 items with largest median probability differences across all methods:")
print(pos.head(10)[['item_label', 'time_label', 'y'] + tmp + ['abs_median_range']])

# Save combined results
tmp = os.path.join(dir_out_pcm, "comparison_pcm_prob_all_methods.csv")
pos.to_csv(tmp, index=False)
print(f"\nSaved ordered_prob comparison to: {tmp}")

# %%

# =============================================================================
# Generate Ordered Probability Comparison Plots (Integrated)
# =============================================================================

print("\n" + "="*70)
print("GENERATING ORDERED PROBABILITY COMPARISON PLOTS")
print("="*70)

pos = []
for method in method_order:
    po = posterior_by_method[method]
    tmp = 'ordered_prob_by_cat_qu_fit'
    if tmp not in po:
        tmp = 'ordered_prob_by_cat_qu_pr'
    tmp2 = _summarize_ordered_prob_quantiles(po[tmp].values, dp1_ukr, dit_ukr)
    tmp2['method'] = method
    pos.append(tmp2)
pos = pd.concat(pos, ignore_index=True)

tmp_emp = dp1_ukr.groupby(
    ['time_label', 'item_label', 'y_label', 'item_type_id', 'item_time_id', 'y']
).size().reset_index(name='n')
tmp_emp_totals = dp1_ukr.groupby(
    ['time_label', 'item_label', 'item_type_id', 'item_time_id']
).size().reset_index(name='total')
tmp_emp = tmp_emp.merge(
    tmp_emp_totals,
    on=['time_label', 'item_label', 'item_type_id', 'item_time_id'],
)
tmp_emp['median'] = tmp_emp['n'] / tmp_emp['total']
tmp_emp['iqr_lower'] = np.nan
tmp_emp['iqr_upper'] = np.nan
tmp_emp['method'] = 'Empirical'
tmp_emp = tmp_emp.merge(
    dit_ukr[['item_type_id', 'item_label', 'group_label_long', 'item_label_short']],
    on=['item_type_id', 'item_label'],
    how='left',
)

plot_cols = [
    'item_type_id', 'item_time_id', 'item_label', 'time_label', 'y', 'y_label',
    'group_label_long', 'item_label_short', 'median', 'iqr_lower', 'iqr_upper', 'method'
]
pos_plot = pd.concat([pos[plot_cols], tmp_emp[plot_cols]], ignore_index=True)
pos_plot['item_label_long'] = pos_plot['group_label_long'] + np.where(
    pos_plot['item_label_short'].notna(), '\n' + pos_plot['item_label_short'], ''
)
pos_plot['method'] = pd.Categorical(
    pos_plot['method'],
    categories=['Empirical'] + method_order,
    ordered=True,
)
pos_plot = pos_plot.dropna(subset=['item_label_long', 'time_label'])

p = (
    ggplot(pos_plot, aes(x='y_label', y='median', fill='method'))
    + geom_col(position=position_dodge(width=0.9), alpha=0.75, width=0.8)
    + geom_errorbar(
        aes(ymin='iqr_lower', ymax='iqr_upper'),
        position=position_dodge(width=0.9),
        width=0.25,
        na_rm=True,
    )
    + facet_wrap('~ item_label_long +time_label', ncol = 3, scales='free')
    + scale_y_continuous(labels=lambda l: [f'{v:.0%}' for v in l], limits=[0, None])
    + scale_fill_futurama()
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, va='top', ha='right', size=7),
        strip_text=element_text(weight='bold'),
        strip_background=element_blank(),
        legend_position='top',
        figure_size=(14, max(8, 2.5 * pos_plot['item_label_long'].nunique())),
    )
    + labs(x='', y='probability', fill='Method')
)

tmp = os.path.join(dir_out_pcm, "comparison_pcm_prob_by_item.pdf")
print(f"Printed combined plot: {tmp}")
p.save(tmp, width=14, height=max(8, 2.5 * pos_plot['item_label_long'].nunique()), units='in', limitsize=False)


# %%

# =============================================================================
# Compare Model Parameters
# =============================================================================

print("\n" + "="*70)
print("COMPARING MODEL PARAMETERS")
print("="*70)

# Define parameter groups to extract
model_pars = [
    'latent_factor_unit',
    'skill_thresholds',
    'loadings_questions_m1',
]

# Use posterior draws from fit results
po_hmc = posterior_by_method['stan_hmc']

# compute quantiles
q_levels = np.array([0.025, 0.25, 0.5, 0.75, 0.975], dtype=float)
tmp2 = {
    m: np.quantile(
        np.concatenate([posterior_by_method[m][p].values for p in model_pars], axis=2),
        q_levels,
        axis=(0, 1)
    )
    for m in method_order
}

# make data frame
tmp = [
    f"{p}[{i+1}]"
    for p in model_pars
    for i in range(po_hmc[p].shape[2])
]
pos = pd.DataFrame({
    'model_par': tmp,
    })
for method in method_order:
    pos[f'q_lower_{method}'] = tmp2[method][0]
    pos[f'iqr_lower_{method}'] = tmp2[method][1]
    pos[f'median_{method}'] = tmp2[method][2]
    pos[f'iqr_upper_{method}'] = tmp2[method][3]
    pos[f'q_upper_{method}'] = tmp2[method][4]

tmp2 = [f'median_{m}' for m in method_order]
pos['abs_median_range'] = pos[tmp2].max(axis=1) - pos[tmp2].min(axis=1)
pos = pos.sort_values('abs_median_range', ascending=True)

tmp = os.path.join(dir_out_pcm, "comparison_model_params_all_methods.csv")
print(f"\nSaved model params comparison to: {tmp}")
pos.to_csv(tmp, index=False)

# =============================================================================
# Plot Comparisong of Model Parameters for top 8 worst parameters
# =============================================================================

print("\nCreating model parameter boxplots for worst 8 parameters...")
pos = pos.tail(8)
tmp = '|'.join(method_order)
pos = pd.wide_to_long(
    pos.drop('abs_median_range', axis=1),
    stubnames=['q_lower', 'iqr_lower', 'median', 'iqr_upper', 'q_upper'],
    i='model_par',
    j='method',
    sep='_',
    suffix=f'({tmp})'
).reset_index()
pos['method'] = pd.Categorical(pos['method'], categories=method_order, ordered=True)

p = (
    ggplot(
        pos,
        aes(x='model_par', fill='method',
            ymin='q_lower', lower='iqr_lower',
            middle='median', upper='iqr_upper',
            ymax='q_upper'
        ),
    )
    + geom_boxplot(stat='identity', position=position_dodge(width=0.8), width=0.7)
    + scale_fill_futurama()
    + coord_flip()
    + theme_bw()
    + theme(
        legend_position='top',
        axis_text_y=element_text(size=7),
        axis_text_x=element_text(size=8),
        axis_title=element_text(size=9),
        figure_size=(12, 8),
    )
    + labs(
        x='Worst 8 model parameters',
        y='Posterior quantiles',
        fill='Method',
    )
)

tmp = os.path.join(dir_out_pcm, 'comparison_model_params_boxplots.pdf')
ggsave(p, filename=tmp, width=12, height=6)
print(f"Saved model parameter boxplot to: {tmp}")

# %%

# =============================================================================
# Compare ypred (Posterior Predictions of Outcomes) - using efficient xarray
# =============================================================================

print("\n" + "="*70)
print("COMPARING YPRED (POSTERIOR PREDICTIONS)")
print("="*70)

print("Computing quantiles for ypred parameters...")

# Use posterior draws from fit results
po_hmc = posterior_by_method['stan_hmc']['ypred']

# Compute quantiles
q_levels = np.array([0.025, 0.25, 0.5, 0.75, 0.975], dtype=float)
tmp2 = {
    m: np.quantile(posterior_by_method[m]['ypred'], q_levels, axis=(0, 1))
    for m in method_order
}

# make data frame
tmp = [ f"{'ypred'}[{i+1}]" for i in range(po_hmc.shape[2]) ]
pos = pd.DataFrame({
    'parameter': tmp,
})
for method in method_order:
    pos[f'q_lower_{method}'] = tmp2[method][0]
    pos[f'iqr_lower_{method}'] = tmp2[method][1]
    pos[f'median_{method}'] = tmp2[method][2]
    pos[f'iqr_upper_{method}'] = tmp2[method][3]
    pos[f'q_upper_{method}'] = tmp2[method][4]

tmp2 = [f'median_{m}' for m in method_order]
pos['abs_median_range'] = pos[tmp2].max(axis=1) - pos[tmp2].min(axis=1)
pos = pos.sort_values('abs_median_range', ascending=False)

# Save comparison
tmp = os.path.join(dir_out_pcm, "comparison_ypred_all_methods.csv")
pos.to_csv(tmp, index=False)
print(f"Saved ypred comparison to: {tmp}")

# =============================================================================
# Plot Comparison of ypred for top 8 worst y predictions
# =============================================================================

print("\nCreating ypred comparison plot...")

pos = pos.head(8)
tmp = '|'.join(method_order)
pos = pd.wide_to_long(
    pos.drop('abs_median_range', axis=1),
    stubnames=['q_lower', 'iqr_lower', 'median', 'iqr_upper', 'q_upper'],
    i='parameter',
    j='method',
    sep='_',
    suffix=f'({tmp})'
).reset_index()
pos['method'] = pd.Categorical(pos['method'], categories=method_order, ordered=True)

p = (
    ggplot(
        pos,
        aes(x='parameter', fill='method',
            ymin='q_lower', lower='iqr_lower',
            middle='median', upper='iqr_upper',
            ymax='q_upper'
        ),
    )
    + geom_boxplot(stat='identity', position=position_dodge(width=0.8), width=0.7)
    + scale_fill_futurama()
    + coord_flip()
    + theme_bw()
    + theme(
        legend_position='top',
        axis_text_y=element_text(size=6),
        axis_text_x=element_text(size=8),
        axis_title=element_text(size=9),
        figure_size=(12, 6)
        )
    + labs(
        x='Worst 8 y predictions',
        y='Posterior quantiles',
        fill='Method'
    )
)

tmp = os.path.join(dir_out_pcm, 'comparison_ypred_boxplots.pdf')
ggsave(p, filename=tmp, width=12, height=6)
print(f"Saved ypred comparison plot to: {tmp}")

# %%

# =============================================================================
# Compare Computation Time Across Methods
# =============================================================================

print("\n" + "="*70)
print("COMPARING COMPUTATION TIME")
print("="*70)

timing_files = {
    'stan_hmc':                          f"{output_file_prefix_hmc}_timing.csv",
    'stan_advi':                         f"{output_file_prefix_advi}_timing.csv",
    'numpyro_svi_autodiagnormal':        f"{output_file_svi_autodiagnormal}_timing.csv",
    'numpyro_svi_autolaplaceapproximation': f"{output_file_svi_autolaplaceapproximation}_timing.csv",
    'numpyro_svi_automultivariatenormal':f"{output_file_svi_automultivariatenormal}_timing.csv",
    'numpyro_svi_autoiafnormal':         f"{output_file_svi_autoiafnormal}_timing.csv",
}

# use mean mins_total across chains for each method,
# since computations are done in parallel for all chains.
pos = []
for method, fpath in timing_files.items():
    if os.path.exists(fpath):
        tmp2 = pd.read_csv(fpath)
        pos.append({'method': method, 'mins_total': tmp2['mins_total'].mean()})
    else:
        print(f"  WARNING: timing file not found for {method}: {fpath}")
pos = pd.DataFrame(pos)
pos['method'] = pd.Categorical(pos['method'], categories=method_order, ordered=True)
pos = pos.sort_values('method')

print("\nTotal mins_total summed across chains (n_sample = 4000):")
print(pos.to_string(index=False))

tmp = os.path.join(dir_out_pcm, "comparison_timing_all_methods.csv")
pos.to_csv(tmp, index=False)
print(f"\nSaved timing comparison to: {tmp}")

p = (
    ggplot(pos, aes(x='method', y='mins_total', fill='method')) +
    geom_col(alpha=0.9, width=0.7) +
    geom_text(aes(label='mins_total.round(1).astype(str)'), va='bottom', size=9, nudge_y=0.02 * pos['mins_total'].max()) +
    scale_fill_futurama() +
    theme_bw() +
    theme(
        axis_text_x=element_text(angle=45, hjust=1, size=9),
        axis_title_x=element_blank(),
        strip_background=element_blank(),
        legend_position='top',
        figure_size=(10, 6),
    ) +
    labs(
        title='Total computation time by method (n_sample = 4000)',
        y='Total minutes (summed across chains)',
        fill='Method',
    )
)

tmp = os.path.join(dir_out_pcm, 'comparison_timing_all_methods.pdf')
ggsave(p, filename=tmp, width=10, height=6)
print(f"Saved timing plot to: {tmp}")
