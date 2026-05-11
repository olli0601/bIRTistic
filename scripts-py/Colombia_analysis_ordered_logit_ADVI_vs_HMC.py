#!/usr/bin/env python3
"""
Investigate ADVI vs HMC for ordered logit model on Colombia data

Date: 2026-04-30

This script applies the ordered logit model (ncats v260413) to the Colombia data 
for comparison with other IRT models. It generates plots and tables for the Hope 
Groups paper using the ordered logit model v260413.

Ported from: Colombia_analysis_ordered_logit_ADVI_vs_HMC.Rmd

Usage:
    # Option 1: Run as script
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/Colombia_analysis_ordered_logit_ADVI_vs_HMC.py
    
    # Option 2: Run interactively line-by-line
    cd /Users/or105/git/bIRTistic
    pixi run python
    >>> exec(open('scripts-py/__init__.py').read())  # Setup paths
    >>> exec(open('scripts-py/Colombia_analysis_ordered_logit_ADVI_vs_HMC.py').read())
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
    ggplot, aes, geom_col, geom_errorbar, geom_density, facet_wrap, facet_grid,
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
from fit_ordered_logit_model_ncats_advi import fit_ordered_logit_model_ncats_advi
from fit_ordered_logit_model_ncats import fit_ordered_logit_model_ncats
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
dir_out_ol = "/Users/or105/sandbox/bIRTistic/py-colombia-ordered_logit-260430-vanilla-advi-vs-hmc"
dir_logs_ol = os.path.join(dir_out_ol, "logs")

# Create output directories
os.makedirs(dir_out_ol, exist_ok=True)
os.makedirs(dir_logs_ol, exist_ok=True)

output_file_prefix = os.path.join(dir_out_ol, "ol_1")
output_file_prefix_hmc = os.path.join(dir_out_ol, "ol_1_hmc")
output_file_prefix_advi = os.path.join(dir_out_ol, "ol_1_advi")

print(f"Data file: {file_data}")
print(f"Output directory: {dir_out_ol}")

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
tmp = os.path.join(dir_out_ol, "ol_1_data.pkl")
print(f"\nSaved preprocessed data to: {tmp}")
pd.to_pickle( {'dp1': dp1_col, 'dit': dit_col, 'dmeta': dmeta_col}, tmp)

# %%

# =============================================================================
# Model Fitting - HMC
# =============================================================================

result_hmc = fit_ordered_logit_model_ncats(
    dit_col,
    dp1_col,
    output_file_prefix=output_file_prefix_hmc,
    stan_file="/Users/or105/git/bIRTistic/src/stan/ordered_logit_ncats_v260413.stan",
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

result_advi = fit_ordered_logit_model_ncats_advi(
    dit_col,
    dp1_col,
    output_file_prefix=output_file_prefix_advi,
    stan_file="/Users/or105/git/bIRTistic/src/stan/ordered_logit_ncats_v260413.stan",
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
tmp = os.path.join(dir_out_ol, "endpoints_comparison_hmc_vs_advi.csv")
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
tmp = os.path.join(dir_out_ol, 'effect_size_difference_comparison.pdf')
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

tmp = os.path.join(dir_out_ol, 'effect_size_ratio_comparison.pdf')
ggsave(p_ratio, filename=tmp, width=14, height=18)
print(f"Saved plot to: {tmp}")

# %%

# =============================================================================
# Compare ordered_prob_by_cat_qu_fit parameters
# =============================================================================

print("\n" + "="*70)
print("COMPARING ORDERED PROBABILITIES")
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
tmp = os.path.join(dir_out_ol, "ordered_prob_comparison_hmc_vs_advi.csv")
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

tmp = os.path.join(dir_out_ol, "ordered_prob_comparison_by_item.pdf")
print(f"Printed combined plot: {tmp}")
p.save(tmp, width=14, height=max(8, 2.5 * pos_plot['item_label_long'].nunique()), units='in', limitsize=False)

    
# %%

# =============================================================================
# Compare Model Parameters (using efficient xarray slicing)
# =============================================================================

print("\n" + "="*70)
print("COMPARING MODEL PARAMETERS")
print("="*70)

# Define parameter groups to extract
param_patterns = {
    'latent_factor_unit': 'latent_factor_unit',
    'skill_thresholds_1': 'skill_thresholds_1',
    'skill_thresholds_incs': 'skill_thresholds_incs',
    'loadings_questions_m1': 'loadings_questions_m1'
}

print("\nParameter counts:")
param_groups = {}
for group, pattern in param_patterns.items():
    vars_in_group = [v for v in idata_hmc.posterior.data_vars if v.startswith(pattern)]
    param_groups[group] = vars_in_group
    print(f"  {group}: {len(vars_in_group)}")

# Compute quantiles for all model parameters
all_model_params = []
for params in param_groups.values():
    all_model_params.extend(params)

if len(all_model_params) > 0:
    print(f"\nComputing quantiles for {len(all_model_params)} model parameters...")
    
    # Select parameters and compute quantiles
    model_hmc = idata_hmc.posterior[all_model_params]
    model_advi = idata_advi.posterior[all_model_params]
    
    model_hmc_q = model_hmc.quantile(quantiles, dim=['chain', 'draw'])
    model_advi_q = model_advi.quantile(quantiles, dim=['chain', 'draw'])
    
    # Convert to DataFrame
    results = []
    for var in all_model_params:
        results.append({
            'parameter': var,
            'median_hmc': float(model_hmc_q[var].sel(quantile=0.5).values),
            'median_advi': float(model_advi_q[var].sel(quantile=0.5).values),
            'q_05_hmc': float(model_hmc_q[var].sel(quantile=0.05).values),
            'q_05_advi': float(model_advi_q[var].sel(quantile=0.05).values),
            'q_25_hmc': float(model_hmc_q[var].sel(quantile=0.25).values),
            'q_25_advi': float(model_advi_q[var].sel(quantile=0.25).values),
            'q_75_hmc': float(model_hmc_q[var].sel(quantile=0.75).values),
            'q_75_advi': float(model_advi_q[var].sel(quantile=0.75).values),
            'q_975_hmc': float(model_hmc_q[var].sel(quantile=0.975).values),
            'q_975_advi': float(model_advi_q[var].sel(quantile=0.975).values),
        })
    
    model_params_combined = pd.DataFrame(results)
    model_params_combined['median_diff'] = model_params_combined['median_hmc'] - model_params_combined['median_advi']
    model_params_combined['abs_median_diff'] = model_params_combined['median_diff'].abs()
    
    # Add parameter group
    model_params_combined['param_group'] = ''
    for group, params in param_groups.items():
        model_params_combined.loc[model_params_combined['parameter'].isin(params), 'param_group'] = group
    
    # Sort by group and absolute difference
    model_params_combined = model_params_combined.sort_values(['param_group', 'abs_median_diff'], ascending=[True, False])
    
    print("\nTop 5 parameters per group with largest median differences:")
    for group in param_groups.keys():
        if len(param_groups[group]) > 0:
            group_data = model_params_combined[model_params_combined['param_group'] == group]
            print(f"\n{group}:")
            print(group_data.head(5)[['parameter', 'median_hmc', 'median_advi', 'median_diff']])
    
    # Save results
    model_params_file = os.path.join(dir_out_ol, "model_params_comparison_hmc_vs_advi.csv")
    model_params_combined.to_csv(model_params_file, index=False)
    print(f"\nSaved model params comparison to: {model_params_file}")
    
    # Create density plot for top parameters from each group (sample efficiently)
    print("\nCreating model parameter density plots...")
    top_per_group = 8
    top_params = []
    for group in param_groups.keys():
        if len(param_groups[group]) > 0:
            group_data = model_params_combined[model_params_combined['param_group'] == group]
            top_params.extend(group_data.head(top_per_group)['parameter'].tolist())
    
    if len(top_params) > 0:
        # Sample draws for plotting (convert to DataFrame only for selected params)
        plot_data_model = []
        for param in top_params:
            # Stack chains and draws, then convert to numpy
            hmc_vals = idata_hmc.posterior[param].stack(sample=('chain', 'draw')).values
            advi_vals = idata_advi.posterior[param].stack(sample=('chain', 'draw')).values
            plot_data_model.extend([
                {'parameter': param, 'value': v, 'method': 'HMC'} for v in hmc_vals
            ])
            plot_data_model.extend([
                {'parameter': param, 'value': v, 'method': 'ADVI'} for v in advi_vals
            ])
        
        plot_data_model = pd.DataFrame(plot_data_model)
        
        # Add parameter group
        plot_data_model['param_group'] = ''
        for group, params in param_groups.items():
            plot_data_model.loc[plot_data_model['parameter'].isin(params), 'param_group'] = group
        
        p_model = (
            ggplot(plot_data_model, aes(x='value', color='method')) +
            geom_density(size=1.0, alpha=0.8) +
            facet_wrap('~parameter', scales='free', ncol=4) +
            scale_color_manual(values={'HMC': '#008080', 'ADVI': '#CA562C'}) +
            theme_bw() +
            theme(
                legend_position='top',
                strip_text=element_text(size=7),
                axis_text=element_text(size=6),
                axis_title=element_text(size=9),
                figure_size=(16, 12)
            ) +
            labs(
                title='Comparison of Posterior Densities: Model Parameters',
                subtitle=f'Top {top_per_group} parameters per group ordered by largest difference',
                x='Parameter value',
                y='Density',
                color='Method'
            )
        )
        
        model_plot_file = os.path.join(dir_out_ol, 'model_params_comparison_densities.pdf')
        ggsave(p_model, filename=model_plot_file, width=16, height=12)
        print(f"Saved model parameter density plot to: {model_plot_file}")
else:
    print("No model parameters found")

# %%

# =============================================================================
# Compare ypred (Posterior Predictions of Outcomes) - using efficient xarray
# =============================================================================

print("\n" + "="*70)
print("COMPARING YPRED (POSTERIOR PREDICTIONS)")
print("="*70)

# Extract ypred parameters
ypred_vars = [v for v in idata_hmc.posterior.data_vars if v.startswith('ypred')]
print(f"Found {len(ypred_vars)} ypred parameters")

if len(ypred_vars) > 0:
    print("Computing quantiles for ypred parameters...")
    
    # Select ypred parameters
    ypred_hmc = idata_hmc.posterior[ypred_vars]
    ypred_advi = idata_advi.posterior[ypred_vars]
    
    # Compute mean, median, and quantiles
    ypred_hmc_mean = ypred_hmc.mean(dim=['chain', 'draw'])
    ypred_advi_mean = ypred_advi.mean(dim=['chain', 'draw'])
    
    ypred_hmc_q = ypred_hmc.quantile([0.25, 0.5, 0.75], dim=['chain', 'draw'])
    ypred_advi_q = ypred_advi.quantile([0.25, 0.5, 0.75], dim=['chain', 'draw'])
    
    # Convert to DataFrame
    results = []
    for var in ypred_vars:
        results.append({
            'parameter': var,
            'mean_hmc': float(ypred_hmc_mean[var].values),
            'mean_advi': float(ypred_advi_mean[var].values),
            'median_hmc': float(ypred_hmc_q[var].sel(quantile=0.5).values),
            'median_advi': float(ypred_advi_q[var].sel(quantile=0.5).values),
            'q25_hmc': float(ypred_hmc_q[var].sel(quantile=0.25).values),
            'q25_advi': float(ypred_advi_q[var].sel(quantile=0.25).values),
            'q75_hmc': float(ypred_hmc_q[var].sel(quantile=0.75).values),
            'q75_advi': float(ypred_advi_q[var].sel(quantile=0.75).values),
        })
    
    ypred_combined = pd.DataFrame(results)
    ypred_combined['median_diff'] = ypred_combined['median_hmc'] - ypred_combined['median_advi']
    ypred_combined['mean_diff'] = ypred_combined['mean_hmc'] - ypred_combined['mean_advi']
    ypred_combined['abs_median_diff'] = ypred_combined['median_diff'].abs()
    ypred_combined = ypred_combined.sort_values('abs_median_diff', ascending=False)
    
    print(f"\nTop 10 ypred with largest median differences:")
    print(ypred_combined.head(10)[['parameter', 'median_hmc', 'median_advi', 'median_diff', 'mean_hmc', 'mean_advi']])
    
    # Find disagreements (different rounded medians)
    disagreements = ypred_combined[ypred_combined['median_hmc'].round() != ypred_combined['median_advi'].round()]
    print(f"\nYpred parameters with different rounded medians: {len(disagreements)} out of {len(ypred_combined)}")
    print(f"Disagreement rate: {100 * len(disagreements) / len(ypred_combined):.2f}%")
    
    # Create bar plot for top 30 ypred
    print("\nCreating ypred comparison plot...")
    top_n_ypred = 30
    top_ypred = ypred_combined.head(top_n_ypred)
    
    plot_data_ypred = pd.concat([
        pd.DataFrame({
            'parameter': top_ypred['parameter'],
            'method': 'HMC',
            'mean': top_ypred['mean_hmc'],
            'q25': top_ypred['q25_hmc'],
            'q75': top_ypred['q75_hmc']
        }),
        pd.DataFrame({
            'parameter': top_ypred['parameter'],
            'method': 'ADVI',
            'mean': top_ypred['mean_advi'],
            'q25': top_ypred['q25_advi'],
            'q75': top_ypred['q75_advi']
        })
    ])
    
    # Reverse order for better visualization
    plot_data_ypred['parameter'] = pd.Categorical(
        plot_data_ypred['parameter'],
        categories=top_ypred['parameter'].tolist()[::-1],
        ordered=True
    )
    
    p_ypred = (
        ggplot(plot_data_ypred, aes(x='parameter', y='mean', fill='method')) +
        geom_col(position=position_dodge(width=0.9), alpha=0.7, width=0.8) +
        geom_errorbar(
            aes(ymin='q25', ymax='q75'),
            position=position_dodge(width=0.9),
            width=0.3
        ) +
        scale_fill_manual(values={'HMC': '#008080', 'ADVI': '#CA562C'}) +
        coord_flip() +
        theme_bw() +
        theme(
            legend_position='top',
            axis_text_y=element_text(size=6),
            axis_text_x=element_text(size=8),
            axis_title=element_text(size=10),
            figure_size=(14, 10)
        ) +
        labs(
            title='Comparison of Posterior Means: ypred',
            subtitle=f'Top {top_n_ypred} parameters ordered by largest median difference\nError bars show 25%-75% quantiles',
            x='Parameter',
            y='Posterior Mean',
            fill='Method'
        )
    )
    
    ypred_plot_file = os.path.join(dir_out_ol, 'ypred_comparison_geom_col.pdf')
    ggsave(p_ypred, filename=ypred_plot_file, width=14, height=10)
    print(f"Saved ypred comparison plot to: {ypred_plot_file}")
else:
    print("No ypred parameters found")

# %%

# =============================================================================
# Summary
# =============================================================================

print("\n" + "="*70)
print("=== SUMMARY OF HMC VS ADVI COMPARISON ===")
print("="*70)

if len(prob_vars) > 0:
    print("\n1. ordered_prob_by_cat_qu_fit:")
    print(f"   - Total parameters: {len(probs_combined)}")
    print(f"   - Max absolute median difference: {probs_combined['abs_median_diff'].max():.6f}")
    print(f"   - Mean absolute median difference: {probs_combined['abs_median_diff'].mean():.6f}")
    print(f"   - Results saved to: ordered_prob_comparison_hmc_vs_advi.csv")

if len(all_model_params) > 0:
    print("\n2. Model parameters:")
    for group in param_groups.keys():
        if len(param_groups[group]) > 0:
            group_data = model_params_combined[model_params_combined['param_group'] == group]
            print(f"   {group}:")
            print(f"     - Parameters: {len(group_data)}")
            print(f"     - Max abs diff: {group_data['abs_median_diff'].max():.6f}")
            print(f"     - Mean abs diff: {group_data['abs_median_diff'].mean():.6f}")
    print(f"   - Results saved to: model_params_comparison_hmc_vs_advi.csv")

if len(ypred_vars) > 0:
    print("\n3. ypred:")
    print(f"   - Total parameters: {len(ypred_combined)}")
    print(f"   - Disagreements (different rounded medians): {len(disagreements)}")
    print(f"   - Disagreement rate: {100 * len(disagreements) / len(ypred_combined):.2f}%")
    print(f"   - Max absolute median difference: {ypred_combined['median_diff'].abs().max():.4f}")



# %%

# =============================================================================
# Generate Endpoint Comparison Plots
# =============================================================================

print("\n[2/2] Generating endpoint comparison bar plots...")
try:
    result = subprocess.run(
        [sys.executable, str(project_root / 'scripts-py' / 'Colombia_endpoints_plots.py')],
        capture_output=True,
        text=True,
        timeout=120
    )
    if result.returncode == 0:
        print("  ✓ Endpoint comparison plots generated successfully")
    else:
        print(f"  ⚠ Warning: Endpoint plots failed with code {result.returncode}")
        if result.stderr:
            print(f"    {result.stderr[:500]}")
except Exception as e:
    print(f"  ⚠ Warning: Could not generate endpoint plots: {e}")

print("\n" + "="*70)
print("ANALYSIS COMPLETE")
print("="*70)
print(f"\n✓ All outputs saved to: {dir_out_ol}")
print("\n✨ Done!")
