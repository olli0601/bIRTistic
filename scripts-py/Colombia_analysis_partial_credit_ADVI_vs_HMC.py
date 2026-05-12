#!/usr/bin/env python3
"""
Investigate ADVI vs HMC for partial credit model on Colombia data

Date: 2026-04-30

This script applies the partial credit model (ncats v260413) to the Colombia data 
for comparison with other IRT models. It generates plots and tables for the Hope 
Groups paper using the partial credit model v260413.

Ported from: Colombia_analysis_partial_credit_ADVI_vs_HMC.Rmd

Usage:
    # Option 1: Run as script
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/Colombia_analysis_partial_credit_ADVI_vs_HMC.py
    
    # Option 2: Run interactively line-by-line
    cd /Users/or105/git/bIRTistic
    pixi run python
    >>> exec(open('scripts-py/__init__.py').read())  # Setup paths
    >>> exec(open('scripts-py/Colombia_analysis_partial_credit_ADVI_vs_HMC.py').read())
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
    ggplot, aes, geom_col, geom_errorbar, geom_density, geom_boxplot, facet_wrap, facet_grid,
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
from data_loading import read_data_colombia
from fit_partial_credit_model_ncats_advi import fit_partial_credit_model_ncats_advi
from fit_partial_credit_model_ncats import fit_partial_credit_model_ncats
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
file_data = os.path.join(dir_data, "Colombia_data_baseline_endline_itemised_250927.csv")

# Output directories
dir_out_pcm = "/Users/or105/sandbox/bIRTistic/py-colombia-partial_credit-260430-vanilla-advi-vs-hmc"
dir_logs_pcm = os.path.join(dir_out_pcm, "logs")

# Create output directories
os.makedirs(dir_out_pcm, exist_ok=True)
os.makedirs(dir_logs_pcm, exist_ok=True)

output_file_prefix = os.path.join(dir_out_pcm, "pcm_1")
output_file_prefix_hmc = os.path.join(dir_out_pcm, "pcm_1_hmc")
output_file_prefix_advi = os.path.join(dir_out_pcm, "pcm_1_advi")

print(f"Data file: {file_data}")
print(f"Output directory: {dir_out_pcm}")

# %%

# =============================================================================
# Load and Preprocess Data
# =============================================================================

# Read Colombia data
print("\nLoading Colombia data...")
result = read_data_colombia(file_data)
dp_col = result['dp'].copy()
dit_col = result['dit'].copy()
dmeta_col = result['dmeta'].copy()

print(f"\nLoaded data:")
print(f"  dp: {len(dp_col):,} rows")
print(f"  dit: {len(dit_col):,} rows")
print(f"  dmeta: {len(dmeta_col):,} rows")

# Preprocess data 
print("\nPreprocessing data...")
dp1_col = dp_col[~dp_col['item_label'].str.contains('agg')].copy()
dp1_col['y_stan'] = dp1_col['y'] + 1

# Merge item_type to avoid indexing issues
dp1_col = dp1_col.merge(dit_col[['item_label', 'item_type']], on='item_label', how='left')

# Create item_time mapping
item_time_df = (
    dp1_col[['item_type', 'item_label', 'time']]
    .drop_duplicates()
    .sort_values(['item_type', 'time', 'item_label'])
    .reset_index(drop=True)
)
item_time_df['item_time_id'] = item_time_df.groupby('item_type').cumcount() + 1
    

# Merge item_time_id
dp1_col = dp1_col.merge(item_time_df, on=['item_label', 'time', 'item_type'], how='left')

# Add item_type_id
dp1_col = dp1_col.merge(
    dit_col[['item_type', 'item_type_id']].drop_duplicates(),
    on='item_type',
    how='left'
)

# Sort and create oid
dp1_col = dp1_col.sort_values(['item_type_id', 'pid', 'time', 'item_label']).reset_index(drop=True)
dp1_col['oid'] = range(1, len(dp1_col) + 1)

# Create oidt (observation id within item type)
dp1_col['oidt'] = dp1_col.groupby('item_type').cumcount() + 1

print(f"\nPreprocessed data:")
print(f"  Observations: {len(dp1_col):,}")
print(f"  Participants: {dp1_col['pid'].nunique()}")
print(f"  Items: {dp1_col['item_label'].nunique()}")
print(f"  Item types: {sorted(dp1_col['item_type'].unique())}")

# Save data for model fitting
tmp = os.path.join(dir_out_pcm, "pcm_1_data.pkl")
print(f"\nSaved preprocessed data to: {tmp}")
pd.to_pickle( {'dp1': dp1_col, 'dit': dit_col, 'dmeta': dmeta_col}, tmp)

# %%

# =============================================================================
# Model Fitting - HMC
# =============================================================================

result_hmc = fit_partial_credit_model_ncats(
    dit_col,
    dp1_col,
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
    with_core_analyses=False,
    with_additional_analyses=False,
    show_messages=False
)

print(f"\n✓ HMC fitting complete")
print(f"  Good chains: {result_hmc['good_chains']}")
print(f"  Draws file: {output_file_prefix_hmc}_draws.zarr")

# %%

# =============================================================================
# Model Fitting - ADVI
# =============================================================================

result_advi = fit_partial_credit_model_ncats_advi(
    dit_col,
    dp1_col,
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
# Compare HMC vs ADVI Results - Compute Endpoints
# =============================================================================

# Load preprocessed data
tmp = pd.read_pickle(f"{output_file_prefix}_data.pkl")
dp1_col = tmp['dp1']
dit_col = tmp['dit']
dmeta_col = tmp['dmeta']

# Get endpoints for HMC
print("\nComputing HMC endpoints...")
endpoints_hmc = get_endpoints(
    dp1=dp1_col,
    dit=dit_col,
    draws_file=f"{output_file_prefix_hmc}_draws.zarr",
    categorical_threshold=3,
    endpoint_type="items"
)
endpoints_hmc['method'] = 'HMC'

# Get endpoints for ADVI
print("Computing ADVI endpoints...")
endpoints_advi = get_endpoints(
    dp1=dp1_col,
    dit=dit_col,
    draws_file=f"{output_file_prefix_advi}_draws.zarr",
    param_name="ordered_prob_by_cat_qu_pr",
    categorical_threshold=3,
    endpoint_type="items"
)
endpoints_advi['method'] = 'ADVI'

print(f"\n✓ Computed endpoints")
print(f"  HMC: {len(endpoints_hmc)} rows")
print(f"  ADVI: {len(endpoints_advi)} rows")

# %%

# =============================================================================
# Create Comparison Table
# =============================================================================

# Combine HMC and ADVI results
pos = pd.concat([endpoints_hmc, endpoints_advi], ignore_index=True)

# Pivot to wide format for comparison
pos = pos.pivot_table(
    index=['item_type_id','item_type','item_label','item_label_short','group_label','group_label_long','item_high_label','variable'],
    columns='method',
    values=['median', 'q_lower', 'iqr_lower', 'iqr_upper', 'q_upper']
).reset_index()

# Flatten column names from multi-index
pos.columns = ['_'.join(col).strip('_') if isinstance(col, tuple) else col 
                          for col in pos.columns]

# Calculate median difference
pos['median_diff'] = pos['median_HMC'] - pos['median_ADVI']
pos['abs_median_diff'] = pos['median_diff'].abs()

# Save comparison
tmp = os.path.join(dir_out_pcm, "comparison_endpoints_hmc_vs_advi.csv")
pos.to_csv(tmp, index=False)
print(f"\nSaved comparison to: {tmp}")

# %%

# =============================================================================
# Visualization: Endpoint Comparison by Item
# =============================================================================

# Filter for difference variable and add item_label_short
# Merge with dit to get item_label_short for better labels
pos = pd.concat([endpoints_hmc, endpoints_advi], ignore_index=True)
endpoints_plot = pos[pos['variable'] == 'diff'].copy()

# Create composite label for facets
endpoints_plot['facet_label'] = (
    endpoints_plot['group_label_long'] + 
    np.where(endpoints_plot['item_label_short'].notna(), 
             '---' + endpoints_plot['item_label_short'] , 
             ''
             ) + '\n' +
    np.where(endpoints_plot['item_type'] == 'categorical', 
             'Difference in probability per week\n(Baseline - Endline)', 
             'Difference in mean days per week\n(Baseline - Endline)'
             )
)

# Create plot with individual items as facets (matching R version)
print("\nCreating endpoint difference plot...")

p_diff = (
    ggplot(endpoints_plot, aes(x='method', y='median', fill='method')) +
    geom_col(alpha=0.8, width=0.7) +
    geom_errorbar(aes(ymin='iqr_lower', ymax='iqr_upper'), width=0.3) +
    facet_wrap('~facet_label', scales='free', ncol=3) +
    scale_fill_manual(values={'HMC': '#008080', 'ADVI': '#CA562C'}) +
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
        title='Endpoint Comparison: Difference',
        y='Estimate',
        fill='Method'
    )
)

# Save plot
tmp = os.path.join(dir_out_pcm, 'comparison_effect_size_difference.pdf')
ggsave(p_diff, filename=tmp, width=14, height=18)
print(f"Saved plot to: {tmp}")

# %%

# Create ratio plot
print("Creating endpoint ratio plot...")
endpoints_plot = pos[pos['variable'] == 'ratio'].copy()

endpoints_plot['facet_label'] = (
    endpoints_plot['group_label_long'] + 
    np.where(endpoints_plot['item_label_short'].notna(), 
             '---' + endpoints_plot['item_label_short'] , 
             ''
             ) + '\n' +
    np.where(endpoints_plot['item_type'] == 'categorical', 
             'Ratio in probability per week\n(1- Baseline/Endline)', 
             'Ratio in mean days per week\n(1- Baseline/Endline)'
             )
)

# Create plot with individual items as facets (matching R version)
print("\nCreating endpoint ratio plot...")

p_ratio = (
    ggplot(endpoints_plot, aes(x='method', y='median', fill='method')) +
    geom_col(alpha=0.8, width=0.7) +
    geom_errorbar(aes(ymin='iqr_lower', ymax='iqr_upper'), width=0.3) +
    facet_wrap('~facet_label', scales='free', ncol=3) +
    scale_fill_manual(values={'HMC': '#008080', 'ADVI': '#CA562C'}) +
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
        title='Endpoint Comparison: Ratio',
        y='Estimate',
        fill='Method'
    )
)

tmp = os.path.join(dir_out_pcm, 'comparison_effect_size_ratio.pdf')
ggsave(p_ratio, filename=tmp, width=14, height=18)
print(f"Saved plot to: {tmp}")

# %%

# =============================================================================
# Compare ordered_prob_by_cat_qu_fit parameters
# =============================================================================

print("\n" + "="*70)
print("COMPARING PARTIAL CREDIT PROBABILITIES")
print("="*70)

po_hmc = result_hmc['draws'].posterior['ordered_prob_by_cat_qu_fit'].values
po_advi = result_advi['draws'].posterior['ordered_prob_by_cat_qu_fit'].values

# Summarize both methods
pos_hmc = _summarize_ordered_prob_quantiles(po_hmc, dp1_col, dit_col)
pos_hmc['method'] = 'HMC'
pos_advi = _summarize_ordered_prob_quantiles(po_advi, dp1_col, dit_col)
pos_advi['method'] = 'ADVI'
pos = pd.concat([pos_hmc, pos_advi], ignore_index=True)

pos = pos.pivot_table(
    index=['cq_id', 'item_type_id', 'item_time_id', 'y', 'item_label', 'time_label',
           'group_label_long', 'item_label_short', 'endpoint_measure'],
    columns='method',
    values='median',
    aggfunc='first'
).reset_index()

pos['median_diff'] = pos['HMC'] - pos['ADVI']
pos['abs_median_diff'] = pos['median_diff'].abs()
pos = pos.sort_values('abs_median_diff', ascending=False)
        
print(f"\nTop 10 items with largest median probability differences (HMC vs ADVI):")
print(pos.head(10)[['item_label', 'time_label', 'y', 'HMC', 'ADVI', 'median_diff']])
    
# Save combined results
tmp = os.path.join(dir_out_pcm, "comparison_pcm_prob_hmc_vs_advi.csv")
pos.to_csv(tmp, index=False)
print(f"\nSaved ordered_prob comparison to: {tmp}")

# %%

# =============================================================================
# Generate Ordered Probability Comparison Plots (Integrated)
# =============================================================================

print("\n" + "="*70)
print("GENERATING ORDERED PROBABILITY COMPARISON PLOTS")
print("="*70)

pos = pd.concat([pos_hmc, pos_advi], ignore_index=True)

tmp_emp = dp1_col.groupby(
    ['time_label', 'item_label', 'y_label', 'item_type_id', 'item_time_id', 'y']
).size().reset_index(name='n')
tmp_emp_totals = dp1_col.groupby(
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
    dit_col[['item_type_id', 'item_label', 'group_label_long', 'item_label_short']],
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
    categories=['Empirical', 'HMC', 'ADVI'],
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
po_hmc = result_hmc['draws'].posterior
po_advi = result_advi['draws'].posterior

# compute quantiles
q_levels = np.array([0.025, 0.25, 0.5, 0.75, 0.975], dtype=float)
q_hmc = np.quantile(np.concatenate([po_hmc[p].values for p in model_pars], axis=2), 
                    q_levels, 
                    axis=(0, 1)
                    )
q_advi = np.quantile(np.concatenate([po_advi[p].values for p in model_pars], axis=2), 
                     q_levels, 
                     axis=(0, 1)
                     )

# make data frame
tmp = [
    f"{p}[{i+1}]"
    for p in model_pars
    for i in range(po_hmc[p].shape[2])
]
pos = pd.DataFrame({
    'model_par': tmp,
    'q_lower_hmc': q_hmc[0],
    'iqr_lower_hmc': q_hmc[1],
    'median_hmc': q_hmc[2],
    'iqr_upper_hmc': q_hmc[3],
    'q_upper_hmc': q_hmc[4],
    'q_lower_advi': q_advi[0],
    'iqr_lower_advi': q_advi[1],
    'median_advi': q_advi[2],
    'iqr_upper_advi': q_advi[3],
    'q_upper_advi': q_advi[4],
    })
pos['abs_median_diff'] = abs(pos['median_hmc'] - pos['median_advi'])
pos = pos.sort_values('abs_median_diff', ascending=True)

tmp = os.path.join(dir_out_pcm, "comparison_model_params_hmc_vs_advi.csv")
print(f"\nSaved model params comparison to: {tmp}")
pos.to_csv(tmp, index=False)
    
# =============================================================================
# Plot Comparisong of Model Parameters for top 8 worst parameters
# =============================================================================

print("\nCreating model parameter boxplots for worst 8 parameters...")
pos = pos.tail(8)
pos = pd.wide_to_long(
    pos.drop('abs_median_diff', axis=1),
    stubnames=['q_lower', 'iqr_lower', 'median', 'iqr_upper', 'q_upper'],
    i='model_par',
    j='method',
    sep='_',
    suffix='(hmc|advi)'
).reset_index()
pos['method'] = pos['method'].str.upper()

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
po_hmc = result_hmc['draws'].posterior['ypred']
po_advi = result_advi['draws'].posterior['ypred']

# Compute quantiles 
q_levels = np.array([0.025, 0.25, 0.5, 0.75, 0.975], dtype=float)
q_hmc = np.quantile(po_hmc, q_levels, axis=(0, 1))
q_advi = np.quantile(po_advi, q_levels, axis=(0, 1))
    
# make data frame
tmp = [ f"{'ypred'}[{i+1}]" for i in range(po_hmc.shape[2]) ]
pos = pd.DataFrame({
    'parameter': tmp,
    'q_lower_hmc': q_hmc[0],
    'iqr_lower_hmc': q_hmc[1],
    'median_hmc': q_hmc[2],
    'iqr_upper_hmc': q_hmc[3],
    'q_upper_hmc': q_hmc[4],
    'q_lower_advi': q_advi[0],
    'iqr_lower_advi': q_advi[1],
    'median_advi': q_advi[2],
    'iqr_upper_advi': q_advi[3],
    'q_upper_advi': q_advi[4],
})
pos['abs_median_diff'] = abs(pos['median_hmc'] - pos['median_advi'])
pos = pos.sort_values('abs_median_diff', ascending=False)
    
# Save comparison
tmp = os.path.join(dir_out_pcm, "comparison_ypred_hmc_vs_advi.csv")
pos.to_csv(tmp, index=False)
print(f"Saved ypred comparison to: {tmp}")    
    
# =============================================================================
# Plot Comparison of ypred for top 8 worst y predictions
# =============================================================================

print("\nCreating ypred comparison plot...")

pos = pos.head(8)
pos = pd.wide_to_long(
    pos.drop('abs_median_diff', axis=1),
    stubnames=['q_lower', 'iqr_lower', 'median', 'iqr_upper', 'q_upper'],
    i='parameter',
    j='method',
    sep='_',
    suffix='(hmc|advi)'
).reset_index()
pos['method'] = pos['method'].str.upper()
    
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
