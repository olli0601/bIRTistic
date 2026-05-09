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
    ggplot, aes, geom_col, geom_errorbar, geom_density, facet_wrap,
    scale_fill_manual, scale_color_manual, theme_minimal, theme_bw, theme, element_text,
    element_blank, labs, position_dodge, coord_flip, ggsave, scale_y_continuous
)
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

# %%

# Save data for model fitting
tmp = os.path.join(dir_out_ol, "ol_1_data_dp.pkl")
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
    show_messages=True
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

print("\n" + "="*70)
print("CREATING COMPARISON TABLE")
print("="*70)

# Combine HMC and ADVI results
endpoints_combined = pd.concat([endpoints_hmc, endpoints_advi], ignore_index=True)

# Create merge ID for comparison
if 'item_time_id' in endpoints_combined.columns:
    endpoints_combined['merge_id'] = (
        endpoints_combined['item_type'] + '___' + 
        endpoints_combined['item_label'] + '___' + 
        endpoints_combined['item_time_id'].fillna('NA').astype(str) + '___' + 
        endpoints_combined['variable']
    )
else:
    endpoints_combined['merge_id'] = (
        endpoints_combined['item_type'] + '___' + 
        endpoints_combined['item_label'] + '___' + 
        endpoints_combined['variable']
    )

# Pivot to wide format for comparison
endpoints_wide = endpoints_combined.pivot_table(
    index='merge_id',
    columns='method',
    values=['median', 'q_lower', 'iqr_lower', 'iqr_upper', 'q_upper']
).reset_index()

# Flatten column names from multi-index
endpoints_wide.columns = ['_'.join(col).strip('_') if isinstance(col, tuple) else col 
                          for col in endpoints_wide.columns]

# Calculate median difference
endpoints_wide['median_diff'] = endpoints_wide['median_HMC'] - endpoints_wide['median_ADVI']
endpoints_wide['abs_median_diff'] = endpoints_wide['median_diff'].abs()

# Save comparison
comparison_file = os.path.join(dir_out_ol, "endpoints_comparison_hmc_vs_advi.csv")
endpoints_wide.to_csv(comparison_file, index=False)
print(f"\nSaved comparison to: {comparison_file}")

# Show top differences
print("\nTop 10 endpoints with largest median differences:")
display_cols = ['merge_id', 'median_HMC', 'median_ADVI', 'median_diff']
print(endpoints_wide.nlargest(10, 'abs_median_diff')[display_cols])

# %%

# =============================================================================
# Visualization: Endpoint Comparison by Item
# =============================================================================

print("\n" + "="*70)
print("CREATING ENDPOINT VISUALIZATIONS")
print("="*70)

# Filter for difference variable and add item_label_short
# Merge with dit to get item_label_short for better labels
endpoints_plot = endpoints_combined[endpoints_combined['variable'] == 'diff'].copy()
endpoints_plot = endpoints_plot.merge(
    dit[['item_label', 'item_label_short']].drop_duplicates(),
    on='item_label',
    how='left'
)

# Create composite label for facets
endpoints_plot['facet_label'] = (
    endpoints_plot['item_label_short'] + '\n' +
    'Difference in mean days per week\n(Baseline - Endline)'
)

# Create plot with individual items as facets (matching R version)
print("\nCreating endpoint difference plot...")
from plotnine import facet_wrap, theme_bw

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
        strip_text=element_text(size=7),
        legend_position='bottom',
        figure_size=(14, 18)
    ) +
    labs(
        title='Endpoint Comparison: Difference',
        y='Estimate',
        fill='Method'
    )
)

# Save plot
plot_file = os.path.join(dir_out_ol, 'effect_size_difference_comparison.pdf')
ggsave(p_diff, filename=plot_file, width=14, height=18)
print(f"Saved plot to: {plot_file}")

# Create ratio plot
print("Creating endpoint ratio plot...")
endpoints_ratio = endpoints_combined[endpoints_combined['variable'] == 'ratio'].copy()
endpoints_ratio = endpoints_ratio.merge(
    dit[['item_label', 'item_label_short']].drop_duplicates(),
    on='item_label',
    how='left'
)
endpoints_ratio['facet_label'] = (
    endpoints_ratio['item_label_short'] + '\n' +
    'Ratio: 1 - (Endline/Baseline)'
)

p_ratio = (
    ggplot(endpoints_ratio, aes(x='method', y='median', fill='method')) +
    geom_col(alpha=0.8, width=0.7) +
    geom_errorbar(aes(ymin='iqr_lower', ymax='iqr_upper'), width=0.3) +
    facet_wrap('~facet_label', scales='free', ncol=3) +
    scale_fill_manual(values={'HMC': '#008080', 'ADVI': '#CA562C'}) +
    theme_bw() +
    theme(
        axis_text_x=element_text(angle=45, hjust=1, size=8),
        axis_title_x=element_blank(),
        strip_text=element_text(size=7),
        legend_position='bottom',
        figure_size=(14, 18)
    ) +
    labs(
        title='Endpoint Comparison: Ratio',
        y='Estimate',
        fill='Method'
    )
)

plot_file_ratio = os.path.join(dir_out_ol, 'effect_size_ratio_comparison.pdf')
ggsave(p_ratio, filename=plot_file_ratio, width=14, height=18)
print(f"Saved plot to: {plot_file_ratio}")

# %%

# =============================================================================
# Compare HMC vs ADVI (Load from pickle files for speed)
# =============================================================================

print("\n" + "="*70)
print("COMPARING HMC VS ADVI")
print("="*70)

# Load draws from Stan pickle files (much faster than NetCDF!)
# Following R script approach: load .rds -> convert to DataFrame -> compute
print("\nLoading HMC and ADVI fit objects from pickle files...")
import time
t0 = time.time()

with open(f"{output_file_prefix_hmc}_stan.pkl", 'rb') as f:
    fit_hmc = pickle.load(f)
with open(f"{output_file_prefix_advi}_stan.pkl", 'rb') as f:
    fit_advi = pickle.load(f)

print(f"Loaded pickle files in {time.time() - t0:.2f} seconds")

# Convert to DataFrames (fast: ~0.15s each)
print("Converting to DataFrames...")
t0 = time.time()
df_hmc = fit_hmc.draws_pd()
df_advi = fit_advi.variational_sample if hasattr(fit_advi, 'variational_sample') else fit_advi.draws_pd()
print(f"Converted in {time.time() - t0:.2f} seconds")
print(f"HMC draws: {df_hmc.shape[0]} samples x {df_hmc.shape[1]} parameters")
print(f"ADVI draws: {df_advi.shape[0]} samples x {df_advi.shape[1]} parameters")

# %%

# =============================================================================
# Compare ordered_prob_by_cat_qu_fit parameters
# =============================================================================

print("\n" + "="*70)
print("COMPARING ORDERED PROBABILITIES")
print("="*70)

# Get list of ordered_prob parameters
prob_params = [c for c in df_hmc.columns if c.startswith('ordered_prob_by_cat_qu_fit')]
print(f"Found {len(prob_params)} ordered_prob_by_cat_qu_fit parameters")

if len(prob_vars) > 0:
    # Compute quantiles directly on xarray (much faster than pandas)
    print("Computing quantiles using xarray (efficient)...")
    
    # Select only ordered_prob parameters
    prob_hmc = idata_hmc.posterior[prob_vars]
    prob_advi = idata_advi.posterior[prob_vars]
    
    # Compute quantiles across chains and draws
    quantiles = [0.05, 0.25, 0.5, 0.75, 0.975]
    prob_hmc_q = prob_hmc.quantile(quantiles, dim=['chain', 'draw'])
    prob_advi_q = prob_advi.quantile(quantiles, dim=['chain', 'draw'])
    
    # Convert to DataFrame for comparison (only summary stats, not full draws)
    results = []
    for var in prob_vars:
        results.append({
            'parameter': var,
            'median_hmc': float(prob_hmc_q[var].sel(quantile=0.5).values),
            'median_advi': float(prob_advi_q[var].sel(quantile=0.5).values),
            'q_05_hmc': float(prob_hmc_q[var].sel(quantile=0.05).values),
            'q_05_advi': float(prob_advi_q[var].sel(quantile=0.05).values),
            'q_25_hmc': float(prob_hmc_q[var].sel(quantile=0.25).values),
            'q_25_advi': float(prob_advi_q[var].sel(quantile=0.25).values),
            'q_75_hmc': float(prob_hmc_q[var].sel(quantile=0.75).values),
            'q_75_advi': float(prob_advi_q[var].sel(quantile=0.75).values),
            'q_975_hmc': float(prob_hmc_q[var].sel(quantile=0.975).values),
            'q_975_advi': float(prob_advi_q[var].sel(quantile=0.975).values),
        })
    
    probs_combined = pd.DataFrame(results)
    probs_combined['median_diff'] = probs_combined['median_hmc'] - probs_combined['median_advi']
    probs_combined['abs_median_diff'] = probs_combined['median_diff'].abs()
    probs_combined = probs_combined.sort_values('abs_median_diff', ascending=False)
    
    print(f"\nTop 10 ordered_prob with largest median differences:")
    print(probs_combined.head(10)[['parameter', 'median_hmc', 'median_advi', 'median_diff']])
    
    # Save results
    probs_file = os.path.join(dir_out_ol, "ordered_prob_comparison_hmc_vs_advi.csv")
    probs_combined.to_csv(probs_file, index=False)
    print(f"\nSaved ordered_prob comparison to: {probs_file}")
else:
    print("No ordered_prob_by_cat_qu_fit parameters found")

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
# Generate Ordered Probability Comparison Plots (Integrated)
# =============================================================================

print("\n" + "="*70)
print("GENERATING ORDERED PROBABILITY COMPARISON PLOTS")
print("="*70)

# Load InferenceData objects for ordered_prob extraction
print("\n[1/7] Loading InferenceData files...")
hmc_draws_file = os.path.join(dir_out_ol, "ol_1_hmc_draws.nc")
advi_draws_file = os.path.join(dir_out_ol, "ol_1_advi_draws.nc")

if not os.path.exists(hmc_draws_file):
    print(f"  ⚠ Warning: HMC draws file not found: {hmc_draws_file}")
    print("  Skipping ordered_prob plots")
elif not os.path.exists(advi_draws_file):
    print(f"  ⚠ Warning: ADVI draws file not found: {advi_draws_file}")
    print("  Skipping ordered_prob plots")
else:
    idata_hmc = az.from_zarr(hmc_draws_file)
    idata_advi = az.from_zarr(advi_draws_file)
    print(f"  ✓ Loaded HMC and ADVI draws")

    # ========================================================================
    # Extract ordered_prob parameters
    # ========================================================================

    print("\n[2/7] Extracting ordered_prob parameters from HMC...")
    ordered_prob_hmc = idata_hmc.posterior['ordered_prob_by_cat_qu_fit']
    print(f"  Shape: {ordered_prob_hmc.shape}")
    
    # Compute quantiles
    q25_hmc = ordered_prob_hmc.quantile(0.25, dim=['chain', 'draw']).values
    median_hmc = ordered_prob_hmc.quantile(0.50, dim=['chain', 'draw']).values
    q75_hmc = ordered_prob_hmc.quantile(0.75, dim=['chain', 'draw']).values
    print(f"  ✓ Extracted {len(median_hmc)} parameters")

    print("\n[3/7] Extracting ordered_prob parameters from ADVI...")
    try:
        ordered_prob_advi = idata_advi.posterior['ordered_prob_by_cat_qu_fit']
        print(f"  Shape: {ordered_prob_advi.shape}")
    except KeyError:
        # ADVI format (individual variables) - reconstruct array
        print("  Detected ADVI format (individual variables)")
        ordered_prob_vars = sorted(
            [v for v in idata_advi.posterior.data_vars if 'ordered_prob_by_cat_qu_fit[' in v],
            key=lambda x: int(x.split('[')[1].split(']')[0])
        )
        print(f"  Found {len(ordered_prob_vars)} variables")
        arrays = [idata_advi.posterior[v] for v in ordered_prob_vars]
        ordered_prob_advi = xr.concat(arrays, dim='ordered_prob_by_cat_qu_fit_dim_0')
        ordered_prob_advi = ordered_prob_advi.assign_coords(
            ordered_prob_by_cat_qu_fit_dim_0=np.arange(len(ordered_prob_vars))
        )
        print(f"  Reconstructed shape: {ordered_prob_advi.shape}")
    
    q25_advi = ordered_prob_advi.quantile(0.25, dim=['chain', 'draw']).values
    median_advi = ordered_prob_advi.quantile(0.50, dim=['chain', 'draw']).values
    q75_advi = ordered_prob_advi.quantile(0.75, dim=['chain', 'draw']).values
    print(f"  ✓ Extracted {len(median_advi)} parameters")

    # ========================================================================
    # Map cq_id to items and response categories
    # ========================================================================

    print("\n[4/7] Mapping cq_id to items and response categories...")
    
    tmp_map = dp1_col[['item_type_id', 'item_label', 'item_time_id']].drop_duplicates()
    tmp_map = tmp_map.merge(
        dit_col[['item_type_id', 'item_label', 'cat_length']],
        on=['item_type_id', 'item_label']
    )
    tmp_map = tmp_map.sort_values(['item_type_id', 'item_time_id'])
    tmp_map['cq_id'] = tmp_map['cat_length'].cumsum() - tmp_map['cat_length']
    
    map_rows = []
    for _, row in tmp_map.iterrows():
        for y in range(row['cat_length']):
            map_rows.append({
                'cq_id': row['cq_id'] + y,
                'item_type_id': row['item_type_id'],
                'item_time_id': row['item_time_id'],
                'y': y
            })
    tmp_map_expanded = pd.DataFrame(map_rows)
    
    tmp_map_expanded = tmp_map_expanded.merge(
        dp1_col[['item_type_id', 'item_time_id', 'y', 'y_label', 'item_label', 'time_label']].drop_duplicates(),
        on=['item_type_id', 'item_time_id', 'y'],
        how='left'
    )
    
    print(f"  ✓ Created mapping for {len(tmp_map_expanded)} combinations")
    
    # Create pos DataFrames
    pos_hmc = tmp_map_expanded.copy()
    pos_hmc['median'] = pos_hmc['cq_id'].map(lambda x: median_hmc[x] if x < len(median_hmc) else np.nan)
    pos_hmc['q25'] = pos_hmc['cq_id'].map(lambda x: q25_hmc[x] if x < len(q25_hmc) else np.nan)
    pos_hmc['q75'] = pos_hmc['cq_id'].map(lambda x: q75_hmc[x] if x < len(q75_hmc) else np.nan)
    pos_hmc['method'] = 'HMC'
    
    pos_advi = tmp_map_expanded.copy()
    pos_advi['median'] = pos_advi['cq_id'].map(lambda x: median_advi[x] if x < len(median_advi) else np.nan)
    pos_advi['q25'] = pos_advi['cq_id'].map(lambda x: q25_advi[x] if x < len(q25_advi) else np.nan)
    pos_advi['q75'] = pos_advi['cq_id'].map(lambda x: q75_advi[x] if x < len(q75_advi) else np.nan)
    pos_advi['method'] = 'ADVI'
    
    pos_df = pd.concat([pos_hmc, pos_advi], ignore_index=True)
    pos_df = pos_df.merge(dit_col[['item_type_id', 'item_label', 'group_label_long', 'item_label_short']], 
                          on=['item_type_id', 'item_label'], how='left')
    pos_df['item_label_long'] = pos_df['group_label_long'] + ' --- ' + pos_df['item_label_short']
    
    print(f"  ✓ Merged with item metadata")

    # ========================================================================
    # Compute empirical probabilities
    # ========================================================================

    print("\n[5/7] Computing empirical probabilities...")
    
    tmp_emp = dp1_col.groupby(['time_label', 'item_label', 'y_label', 'item_type_id', 'item_time_id', 'y']).size().reset_index(name='n')
    tmp2 = dp1_col.groupby(['time_label', 'item_label', 'item_type_id', 'item_time_id']).size().reset_index(name='total')
    tmp_emp = tmp_emp.merge(tmp2, on=['time_label', 'item_label', 'item_type_id', 'item_time_id'])
    tmp_emp['median'] = tmp_emp['n'] / tmp_emp['total']
    tmp_emp['q25'] = np.nan
    tmp_emp['q75'] = np.nan
    tmp_emp['method'] = 'Empirical'
    tmp_emp = tmp_emp.drop(columns=['n', 'total'])
    
    tmp_emp = tmp_emp.merge(dit_col[['item_type_id', 'item_label', 'group_label_long', 'item_label_short']], 
                            on=['item_type_id', 'item_label'], how='left')
    tmp_emp['item_label_long'] = tmp_emp['group_label_long'] + ' --- ' + tmp_emp['item_label_short']
    
    print(f"  ✓ Computed empirical probabilities for {len(tmp_emp)} combinations")

    # ========================================================================
    # Combine and create plots
    # ========================================================================

    common_cols = ['item_label_long', 'time_label', 'y_label', 'median', 'q25', 'q75', 'method']
    pos_combined = pd.concat([
        pos_df[common_cols],
        tmp_emp[common_cols]
    ], ignore_index=True)
    
    pos_combined['method'] = pd.Categorical(
        pos_combined['method'], 
        categories=['Empirical', 'HMC', 'ADVI'], 
        ordered=True
    )
    
    print(f"  ✓ Combined all methods: {len(pos_combined)} rows")

    print("\n[6/7] Creating individual subplots...")
    
    pos_combined = pos_combined.dropna(subset=['item_label_long', 'time_label'])
    items = sorted(pos_combined['item_label_long'].unique())
    times = sorted(pos_combined['time_label'].unique())
    
    print(f"  Items: {len(items)}, Times: {len(times)}")
    
    plot_files = []
    for item in items:
        for time in times:
            plot_data = pos_combined[
                (pos_combined['item_label_long'] == item) & 
                (pos_combined['time_label'] == time)
            ].copy()
            
            if len(plot_data) == 0:
                continue
                
            p = (
                ggplot(plot_data, aes(x='y_label', y='median', fill='method', group='method'))
                + geom_col(position=position_dodge(width=0.9), alpha=0.7, width=0.8)
                + geom_errorbar(
                    aes(ymin='q25', ymax='q75'),
                    position=position_dodge(width=0.9),
                    width=0.25
                )
                + scale_y_continuous(labels=lambda l: [f'{v:.0%}' for v in l], limits=[0, None])
                + scale_fill_manual(values={'Empirical': '#B4C8A8', 'HMC': '#008080', 'ADVI': '#CA562C'})
                + ggtitle(f"{item} - {time}")
                + theme_bw()
                + theme(
                    axis_text_x=element_text(angle=45, va='top', ha='right', size=7),
                    legend_position='none',
                    plot_title=element_text(size=8),
                    figure_size=(6, 4)
                )
                + labs(x='', y='probability', fill='Method')
            )
            
            subplot_file = os.path.join(dir_out_ol, f"_subplot_{item.replace('/', '_').replace(' ', '_')}_{time}.png")
            p.save(subplot_file, dpi=150, width=6, height=4)
            plot_files.append(subplot_file)
    
    print(f"  ✓ Created {len(plot_files)} individual subplots")

    print("\n[7/7] Combining subplots into multi-panel grid...")
    
    n_items = len(items)
    n_cols = 2
    n_rows = n_items
    
    fig = plt.figure(figsize=(12, 4 * n_items))
    gs = gridspec.GridSpec(n_rows + 1, n_cols, figure=fig, height_ratios=[0.3] + [1]*n_rows, hspace=0.3, wspace=0.2)
    
    ax_legend = fig.add_subplot(gs[0, :])
    ax_legend.axis('off')
    
    from matplotlib.patches import Rectangle
    legend_colors = {'Empirical': '#B4C8A8', 'HMC': '#008080', 'ADVI': '#CA562C'}
    legend_handles = [Rectangle((0, 0), 1, 1, fc=color, alpha=0.7) for color in legend_colors.values()]
    ax_legend.legend(legend_handles, legend_colors.keys(), loc='center', ncol=3, frameon=False, fontsize=12)
    
    for i, item in enumerate(items):
        for j, time in enumerate(times):
            subplot_file = os.path.join(dir_out_ol, f"_subplot_{item.replace('/', '_').replace(' ', '_')}_{time}.png")
            if os.path.exists(subplot_file):
                ax = fig.add_subplot(gs[i+1, j])
                img = Image.open(subplot_file)
                ax.imshow(img)
                ax.axis('off')
    
    combined_file = os.path.join(dir_out_ol, "ordered_prob_comparison_by_item.pdf")
    plt.savefig(combined_file, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved combined plot: {combined_file}")
    
    # Clean up temporary files
    print("\n  Cleaning up temporary subplot files...")
    for f in plot_files:
        if os.path.exists(f):
            os.remove(f)
    print(f"  ✓ Removed {len(plot_files)} temporary files")
    
    print("\n" + "="*70)
    print("ORDERED PROBABILITY COMPARISON PLOT COMPLETED")
    print("="*70)

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
