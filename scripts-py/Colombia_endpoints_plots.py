"""
Generate endpoint comparison plots for Colombia ADVI vs HMC analysis
Produces faceted bar plots comparing HMC and ADVI endpoint estimates with error bars
Matches the R markdown output format
"""

import pandas as pd
import numpy as np
from plotnine import *
import sys
import os

# Output directory
output_dir = "/Users/or105/sandbox/bIRTistic/py-colombia-ordered_logit-260430-vanilla-advi-vs-hmc"

print("="*70)
print("GENERATING ENDPOINT COMPARISON PLOTS")
print("="*70)

# Load endpoints comparison data
df = pd.read_csv(os.path.join(output_dir, "endpoints_comparison_hmc_vs_advi.csv"))

# Parse merge_id to extract metadata
df[['item_type', 'item_label', 'item_time_id', 'variable']] = df['merge_id'].str.split('___', expand=True)

# Filter to only diff and ratio variables
df_plot = df[df['variable'].isin(['diff', 'ratio'])].copy()

# Load dit (data item types) to get group_label_long and item_label_short
dit = pd.read_csv(os.path.join(output_dir, "ol_1_advi_data_dit.csv"))

# Merge with dit to get group_label_long and item_label_short
df_plot = df_plot.merge(
    dit[['item_type_id', 'item_label', 'group_label_long', 'item_label_short']],
    left_on='item_label',
    right_on='item_label',
    how='left'
)

# Restructuring data for faceted plotting

# Extract from comparison file
endpoints_combined = []
for _, row in df_plot.iterrows():
    # HMC row
    endpoints_combined.append({
        'item_type': row['item_type'],
        'item_label': row['item_label'],
        'item_label_short': row.get('item_label_short', row['item_label']),
        'group_label_long': row.get('group_label_long', row['item_type']),
        'item_time_id': row['item_time_id'],
        'variable': row['variable'],
        'median': row['median_HMC'],
        'iqr_lower': row['iqr_lower_HMC'],
        'iqr_upper': row['iqr_upper_HMC'],
        'q_lower': row['q_lower_HMC'],
        'q_upper': row['q_upper_HMC'],
        'method': 'HMC'
    })
    # ADVI row
    endpoints_combined.append({
        'item_type': row['item_type'],
        'item_label': row['item_label'],
        'item_label_short': row.get('item_label_short', row['item_label']),
        'group_label_long': row.get('group_label_long', row['item_type']),
        'item_time_id': row['item_time_id'],
        'variable': row['variable'],
        'median': row['median_ADVI'],
        'iqr_lower': row['iqr_lower_ADVI'],
        'iqr_upper': row['iqr_upper_ADVI'],
        'q_lower': row['q_lower_ADVI'],
        'q_upper': row['q_upper_ADVI'],
        'method': 'ADVI'
    })
endpoints_plot = pd.DataFrame(endpoints_combined)

print(f"  ✓ Combined data: {len(endpoints_plot)} rows")

print(f"  ✓ Combined data: {len(endpoints_plot)} rows")

# Create item identifier and facet labels
print("\n[4/5] Creating facet labels...")
if 'item_label_short' in endpoints_plot.columns:
    endpoints_plot['item_id'] = endpoints_plot['item_label_short']
else:
    endpoints_plot['item_id'] = endpoints_plot['item_label']

# Set method factor for consistent ordering
endpoints_plot['method'] = pd.Categorical(
    endpoints_plot['method'], 
    categories=['HMC', 'ADVI'], 
    ordered=True
)

# Create y-axis labels based on item_type and variable (matching R markdown logic)
def create_y_label(row):
    if row['variable'] == 'diff':
        if row['item_type'] == 'out-of-7':
            return "Difference in mean days per week\n(Baseline - Endline)"
        else:
            return "Difference in frequency\nof categories 'most of the time' and 'all of the time'"
    else:  # ratio
        if row['item_type'] == 'out-of-7':
            return "% change in mean days per week\n(1 - Endline/Baseline)"
        else:
            return "% change in frequency\nof categories 'most of the time' and 'all of the time'\n(1 - Endline/Baseline)"

endpoints_plot['y_label'] = endpoints_plot.apply(create_y_label, axis=1)

# Create facet label combining group and y-axis description
endpoints_plot['facet_label'] = endpoints_plot['group_label_long'] + '\n' + endpoints_plot['y_label']

print(f"  ✓ Created {endpoints_plot['facet_label'].nunique()} unique facet labels")

# ============================================================================
# Plot 1: Diff variable (Baseline - Endline) with facet_wrap
# ============================================================================

print("\n[5/5] Creating faceted plots...")
print("  Creating diff comparison plot...")

p_diff = (
    ggplot(endpoints_plot[endpoints_plot['variable'] == 'diff'], 
           aes(x='item_id', y='median', fill='method'))
    + geom_col(position=position_dodge(width=0.9), alpha=0.8, width=0.8)
    + geom_errorbar(
        aes(ymin='iqr_lower', ymax='iqr_upper'),
        position=position_dodge(width=0.9),
        width=0.3,
        size=0.4
    )
    + facet_wrap('~ facet_label', scales='free', ncol=2)
    + scale_fill_manual(values={'HMC': '#008080', 'ADVI': '#CA562C'})
    + theme_minimal()
    + theme(
        axis_text_x=element_text(angle=45, va='top', ha='right', size=9),
        axis_text_y=element_text(size=9),
        legend_position='bottom',
        legend_title=element_text(size=10, weight='bold'),
        legend_text=element_text(size=9),
        strip_text=element_text(size=9, weight='bold', linespacing=1.1),
        strip_background=element_blank(),
        panel_spacing=1,
        plot_title=element_text(size=12, weight='bold', ha=0.5),
        plot_margin=10,
        figure_size=(12, 14)
    )
    + labs(
        title='Endpoint Comparison: Difference',
        x='Item',
        y='Estimate',
        fill='Method'
    )
)

# Save plot
diff_plot_file = os.path.join(output_dir, "endpoints_comparison_barplot_diff.pdf")
p_diff.save(diff_plot_file, width=12, height=14, limitsize=False)
print(f"  ✓ Saved: {diff_plot_file}")

# ============================================================================
# Plot 2: Ratio variable (1 - Endline/Baseline) with facet_wrap
# ============================================================================

print("  Creating ratio comparison plot...")

p_ratio = (
    ggplot(endpoints_plot[endpoints_plot['variable'] == 'ratio'], 
           aes(x='item_id', y='median', fill='method'))
    + geom_col(position=position_dodge(width=0.9), alpha=0.8, width=0.8)
    + geom_errorbar(
        aes(ymin='iqr_lower', ymax='iqr_upper'),
        position=position_dodge(width=0.9),
        width=0.3,
        size=0.4
    )
    + facet_wrap('~ facet_label', scales='free_y', ncol=1)
    + scale_fill_manual(values={'HMC': '#008080', 'ADVI': '#CA562C'})
    + theme_minimal()
    + theme(
        axis_text_x=element_text(angle=45, va='top', ha='right', size=9),
        axis_text_y=element_text(size=9),
        legend_position='bottom',
        legend_title=element_text(size=10, weight='bold'),
        legend_text=element_text(size=9),
        strip_text=element_text(size=9, weight='bold', linespacing=1.1),
        strip_background=element_blank(),
        panel_spacing=1,
        plot_title=element_text(size=12, weight='bold', ha=0.5),
        plot_margin=10,
        figure_size=(12, 14)
    )
    + labs(
        title='Endpoint Comparison: Ratio',
        x='Item',
        y='Estimate',
        fill='Method'
    )
)

# Save plot
ratio_plot_file = os.path.join(output_dir, "endpoints_comparison_barplot_ratio.pdf")
p_ratio.save(ratio_plot_file, width=12, height=14, limitsize=False)
print(f"  ✓ Saved: {ratio_plot_file}")

print("\n" + "="*70)
print("ENDPOINT COMPARISON PLOTS COMPLETED")
print("="*70)
print(f"Generated 2 faceted plots:")
print(f"  1. {diff_plot_file}")
print(f"  2. {ratio_plot_file}")
print("="*70)
