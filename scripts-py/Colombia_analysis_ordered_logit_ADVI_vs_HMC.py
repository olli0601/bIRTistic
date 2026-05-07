#!/usr/bin/env python3
"""
Investigate ADVI vs HMC for ordered logit model on Colombia data

Date: 2026-04-30

This script applies the ordered logit model (ncats v260413) to the Colombia data 
for comparison with other IRT models. It generates plots and tables for the Hope 
Groups paper using the ordered logit model v260413.

Ported from: Colombia_analysis_ordered_logit_ADVI_vs_HMC.Rmd

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python vignettes/Colombia_analysis_ordered_logit_ADVI_vs_HMC.py
"""

# =============================================================================
# Setup and Imports
# =============================================================================

import sys
import os
from pathlib import Path

# Add python module to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'python'))

import pandas as pd
import numpy as np
import arviz as az
from plotnine import *
import warnings
warnings.filterwarnings('ignore')

# Import bIRTistic functions
from data_loading import read_data_colombia
from model_fitting import (
    fit_ordered_logit_model_ncats,
    fit_ordered_logit_model_ncats_advi
)
from analysis import get_endpoints

print("✓ Imports successful")

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

output_file_prefix_hmc = os.path.join(dir_out_ol, "ol_1_hmc")
output_file_prefix_advi = os.path.join(dir_out_ol, "ol_1_advi")

print(f"Data file: {file_data}")
print(f"Output directory: {dir_out_ol}")

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

# Preprocess data (following R logic)
print("\nPreprocessing data...")
dp1_col = dp_col[~dp_col['item_label'].str.contains('agg')].copy()
dp1_col['y_stan'] = dp1_col['y'] + 1

# Merge item_type first to avoid indexing issues
dp1_col = dp1_col.merge(dit_col[['item_label', 'item_type']], on='item_label', how='left')

# Create item_time mapping
item_time_list = []
for item_type in sorted(dp1_col['item_type'].unique()):
    item_type_data = dp1_col[dp1_col['item_type'] == item_type]
    
    for time in sorted(item_type_data['time'].unique()):
        time_data = item_type_data[item_type_data['time'] == time]
        for item_label in sorted(time_data['item_label'].unique()):
            item_time_list.append({
                'item_type': item_type,
                'item_label': item_label,
                'time': time
            })

item_time_df = pd.DataFrame(item_time_list)
item_time_df['item_time_id'] = range(1, len(item_time_df) + 1)

# Merge item_time_id
dp1_col = dp1_col.merge(item_time_df, on=['item_label', 'time', 'item_type'], how='left')

# Add item_type_id (item_type already merged above)
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

# =============================================================================
# Save Preprocessed Data
# =============================================================================

# Save data for model fitting
data_file = os.path.join(dir_out_ol, "ol_1_hmc_data.pkl")
pd.to_pickle({'dp1': dp1_col, 'dit': dit_col}, data_file)
print(f"\nSaved preprocessed data to: {data_file}")

# =============================================================================
# Model Fitting - HMC
# =============================================================================

print("\n" + "="*70)
print("FITTING HMC MODEL")
print("="*70)

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
print(f"  Draws file: {output_file_prefix_hmc}_draws.pkl")

# =============================================================================
# Model Fitting - ADVI
# =============================================================================

print("\n" + "="*70)
print("FITTING ADVI MODEL")
print("="*70)

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
print(f"  Draws file: {output_file_prefix_advi}_draws.pkl")

# =============================================================================
# Compare HMC vs ADVI Results - Compute Endpoints
# =============================================================================

print("\n" + "="*70)
print("COMPUTING ENDPOINTS")
print("="*70)

# Load preprocessed data
data = pd.read_pickle(f"{output_file_prefix_hmc}_data.pkl")
dp1 = data['dp1']
dit = data['dit']

# Get endpoints for HMC
print("\nComputing HMC endpoints...")
endpoints_hmc = get_endpoints(
    dp1=dp1,
    dit=dit,
    draws_file=f"{output_file_prefix_hmc}_draws.pkl",
    categorical_threshold=3,
    endpoint_type="items"
)
endpoints_hmc['method'] = 'HMC'

# Get endpoints for ADVI
print("Computing ADVI endpoints...")
endpoints_advi = get_endpoints(
    dp1=dp1,
    dit=dit,
    draws_file=f"{output_file_prefix_advi}_draws.pkl",
    categorical_threshold=3,
    endpoint_type="items"
)
endpoints_advi['method'] = 'ADVI'

print(f"\n✓ Computed endpoints")
print(f"  HMC: {len(endpoints_hmc)} rows")
print(f"  ADVI: {len(endpoints_advi)} rows")

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

# =============================================================================
# Visualization: Effect Size Comparison
# =============================================================================

print("\n" + "="*70)
print("CREATING VISUALIZATIONS")
print("="*70)

# Filter for key variables
endpoints_plot = endpoints_combined[endpoints_combined['variable'].isin(['diff', 'ratio'])].copy()

# Create plot for difference
print("\nCreating difference plot...")
p_diff = (
    ggplot(endpoints_plot[endpoints_plot['variable'] == 'diff'], 
           aes(x='item_label', y='median', fill='method')) +
    geom_col(position=position_dodge(width=0.9), alpha=0.8) +
    geom_errorbar(aes(ymin='iqr_lower', ymax='iqr_upper'),
                  position=position_dodge(width=0.9),
                  width=0.3) +
    facet_wrap('~item_type', scales='free') +
    scale_fill_manual(values={'HMC': '#008080', 'ADVI': '#CA562C'}) +
    theme_minimal() +
    theme(axis_text_x=element_text(angle=45, hjust=1),
          figure_size=(12, 8)) +
    labs(title='Effect Size Comparison: Difference (Baseline - Endline)',
         x='Item',
         y='Median Difference',
         fill='Method')
)

# Save plot
plot_file = os.path.join(dir_out_ol, 'effect_size_difference_comparison.png')
ggsave(p_diff, filename=plot_file, width=12, height=8, dpi=300)
print(f"Saved plot to: {plot_file}")

# Create plot for ratio
print("Creating ratio plot...")
p_ratio = (
    ggplot(endpoints_plot[endpoints_plot['variable'] == 'ratio'], 
           aes(x='item_label', y='median', fill='method')) +
    geom_col(position=position_dodge(width=0.9), alpha=0.8) +
    geom_errorbar(aes(ymin='iqr_lower', ymax='iqr_upper'),
                  position=position_dodge(width=0.9),
                  width=0.3) +
    facet_wrap('~item_type', scales='free') +
    scale_fill_manual(values={'HMC': '#008080', 'ADVI': '#CA562C'}) +
    theme_minimal() +
    theme(axis_text_x=element_text(angle=45, hjust=1),
          figure_size=(12, 8)) +
    labs(title='Effect Size Comparison: Ratio (1 - Endline/Baseline)',
         x='Item',
         y='Median Ratio',
         fill='Method')
)

# Save plot
plot_file_ratio = os.path.join(dir_out_ol, 'effect_size_ratio_comparison.png')
ggsave(p_ratio, filename=plot_file_ratio, width=12, height=8, dpi=300)
print(f"Saved plot to: {plot_file_ratio}")

# =============================================================================
# Summary
# =============================================================================

print("\n" + "="*70)
print("ANALYSIS COMPLETE")
print("="*70)
print(f"\n✓ All outputs saved to: {dir_out_ol}")
print("\nGenerated files:")
print(f"  - {output_file_prefix_hmc}_draws.pkl")
print(f"  - {output_file_prefix_advi}_draws.pkl")
print(f"  - {comparison_file}")
print(f"  - {plot_file}")
print(f"  - {plot_file_ratio}")
print("\n✨ Done!")
