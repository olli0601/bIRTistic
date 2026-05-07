"""
Analysis functions for Bayesian IRT models.

This module provides functions for analyzing posterior draws from IRT models,
including endpoint extraction and effect size computation.
"""

from typing import Optional, Literal
import os
import re

import pandas as pd
import numpy as np


def get_endpoints(
    dp1: pd.DataFrame,
    dit: pd.DataFrame,
    draws_file: str,
    categorical_threshold: int = 3,
    endpoint_type: Literal["items", "item_groups"] = "items",
    param_name: str = "ordered_prob_by_cat_qu_pr",
    fit_file: Optional[str] = None,
) -> pd.DataFrame:
    """
    Get endpoints from ordered logit model draws.
    
    This function processes posterior draws from an ordered logit model to compute
    endpoint statistics for both categorical and out-of-7 response items.
    For categorical items, it aggregates probabilities above a threshold. For out-of-7
    items, it computes means. Results include baseline, endline, differences, and ratios
    with credible intervals.
    
    Parameters
    ----------
    dp1 : pd.DataFrame
        Input data with columns: item_type_id, item_label, item_time_id, item_type, time_label, y, y_label.
    dit : pd.DataFrame
        Item definitions with columns: item_type_id, item_label, cat_length, item_type, group_label,
        item_label_short, group_label_long, item_high_label.
    draws_file : str
        Path to pickle file containing posterior draws (array format with ordered_prob_by_cat_qu_pr
        or ordered_prob_by_cat_qu_fit parameters).
    fit_file : str, optional
        Path to pickle file containing the cmdstanpy fit object. If provided, column names will be
        extracted from the fit object. If not provided, will attempt to infer from draws_file path.
    categorical_threshold : int, default 3
        Threshold for aggregating categorical responses. Categories >= this value are summed.
        Default is 3, which aggregates "most" and "all of the time" for 5-category items.
    endpoint_type : {"items", "item_groups"}, default "items"
        Aggregation level:
        - "items": Keeps individual items (item_type_id, item_time_id preserved)
        - "item_groups": Aggregates across items within groups (group_label level)
    param_name : str, default "ordered_prob_by_cat_qu_pr"
        Which ordered probability parameters to use:
        - "ordered_prob_by_cat_qu_pr": Posterior predictions
        - "ordered_prob_by_cat_qu_fit": Fitted values
    
    Returns
    -------
    pd.DataFrame
        Endpoint summaries with columns:
        - item_type, group_label, group_label_long: Item type and groupings
        - item_label, item_label_short (if endpoint_type="items"): Item identifiers
        - item_time_id (if endpoint_type="items"): Time point identifier (categorical items only)
        - item_high_label: Direction indicator (higher_is_better/lower_is_better)
        - variable: Endpoint type (Baseline, Endline, diff, ratio, diff2, ratio2)
        - q_lower, iqr_lower, median, iqr_upper, q_upper: Quantile summaries (2.5%, 25%, 50%, 75%, 97.5%)
    
    Notes
    -----
    The function processes two types of items:
    
    **Categorical items**: Aggregates probabilities for categories >= categorical_threshold.
    - endpoint_type="items": Computes for each individual item
    - endpoint_type="item_groups": Averages probabilities across items within each group
    
    **Out-of-7 items**: Computes expected means (weighted sum of category probabilities).
    - endpoint_type="items": Computes for each individual item  
    - endpoint_type="item_groups": Averages means across items within each group
    
    Computed measures:
    - diff: Baseline - Endline (reduction, for categorical or lower_is_better items)
    - ratio: 1 - Endline/Baseline (proportional reduction)
    - diff2: Endline - Baseline (increase)
    - ratio2: Endline/Baseline - 1 (proportional increase)
    
    For out-of-7 items, direction of diff/ratio respects item_high_label when endpoint_type="items".
    """
    # Validate inputs
    if not os.path.exists(draws_file):
        raise FileNotFoundError(f"Draws file not found: {draws_file}")
    
    if endpoint_type not in ["items", "item_groups"]:
        raise ValueError("endpoint_type must be either 'items' or 'item_groups'")
    
    if param_name not in ["ordered_prob_by_cat_qu_pr", "ordered_prob_by_cat_qu_fit"]:
        raise ValueError("param_name must be either 'ordered_prob_by_cat_qu_pr' or 'ordered_prob_by_cat_qu_fit'")
    
    # Load draws
    print(f"Loading draws from: {draws_file}")
    poa = pd.read_pickle(draws_file)
    
    # Handle both HMC (3D: iterations × chains × parameters) and ADVI (2D: samples × parameters)
    if poa.ndim == 2:
        # ADVI format: reshape to 3D with single chain
        print("Detected ADVI format (2D array), reshaping to 3D...")
        n_samples, n_params = poa.shape
        poa = poa.reshape(n_samples, 1, n_params)
    elif poa.ndim != 3:
        raise ValueError(f"Unexpected draws array shape: {poa.shape}. Expected 2D (ADVI) or 3D (HMC)")
    
    # Get column names from fit object
    if fit_file is None:
        # Infer fit file path from draws file path
        fit_file = draws_file.replace('_draws.pkl', '_stan.pkl')
    
    if not os.path.exists(fit_file):
        raise FileNotFoundError(f"Fit file not found: {fit_file}. Please provide fit_file parameter.")
    
    print(f"Loading fit object from: {fit_file}")
    fit = pd.read_pickle(fit_file)
    param_names = fit.column_names
    
    # =========================================================================
    # Part 1: Categorical items - aggregate probabilities for high categories
    # =========================================================================
    
    print("Processing categorical items...")
    
    # Extract ordered probability parameters from draws array
    # poa is shape (iterations, chains, parameters)
    # We need to convert to long format
    n_iter, n_chains, n_params = poa.shape
    
    # Find indices of ordered_prob parameters using actual column names
    ordered_prob_indices = [i for i, name in enumerate(param_names) if param_name in name]
    
    if len(ordered_prob_indices) == 0:
        raise ValueError(f"No parameters found matching '{param_name}' in the draws")
    
    # Create long-format dataframe
    po_data = []
    for chain_idx in range(n_chains):
        for iter_idx in range(n_iter):
            for param_idx in ordered_prob_indices:
                po_data.append({
                    '.draw': iter_idx + chain_idx * n_iter + 1,
                    '.chain': chain_idx + 1,
                    '.iteration': iter_idx + 1,
                    'variable': param_names[param_idx],
                    'prob': poa[iter_idx, chain_idx, param_idx]
                })
    
    po = pd.DataFrame(po_data)
    
    # Extract cq_id from variable name
    po['cq_id'] = po['variable'].str.extract(r'\[(\d+)\]')[0].astype(int)
    po = po.drop(columns=['variable'])
    
    # Map cq_id back to item structure
    tmp = dp1[['item_type_id', 'item_label', 'item_time_id']].drop_duplicates()
    tmp = tmp.merge(dit[['item_type_id', 'item_label', 'cat_length']], 
                   on=['item_type_id', 'item_label'])
    tmp = tmp.sort_values(['item_type_id', 'item_time_id']).reset_index(drop=True)
    tmp['cq_id'] = tmp.groupby(['item_type_id'])['cat_length'].cumsum() - tmp['cat_length'] + 1
    
    # Expand to get all category levels
    cq_mapping = []
    for _, row in tmp.iterrows():
        for y in range(row['cat_length']):
            cq_mapping.append({
                'item_type_id': row['item_type_id'],
                'item_time_id': row['item_time_id'],
                'cq_id': row['cq_id'] + y,
                'y': y
            })
    cq_mapping = pd.DataFrame(cq_mapping)
    
    po = po.merge(cq_mapping, on='cq_id')
    
    # Filter for categorical items and aggregate high categories
    tmp = dp1[dp1['item_type'] == 'categorical'][
        ['item_type_id', 'item_label', 'item_time_id', 'item_type', 'time_label']
    ].drop_duplicates()
    
    po = po.merge(tmp, on=['item_type_id', 'item_time_id'])
    
    # Aggregate categories >= threshold
    po = po[po['y'] >= categorical_threshold]
    
    po = po.groupby(['.draw', 'item_type_id', 'item_time_id', 'item_label', 'item_type', 'time_label']).agg({
        'prob': 'sum'
    }).reset_index()
    po['y_stan'] = 45  # Marker value
    
    po = po.merge(dit[['item_type', 'item_label', 'group_label']], 
                 on=['item_type', 'item_label'])
    
    # Aggregate across items if endpoint_type="item_groups"
    if endpoint_type == "item_groups":
        po = po.groupby(['.draw', 'item_type', 'time_label', 'group_label']).agg({
            'prob': 'mean'
        }).reset_index()
        
        # Reshape to wide format (Baseline vs Endline)
        po_wide = po.pivot_table(
            index=['.draw', 'item_type', 'group_label'],
            columns='time_label',
            values='prob'
        ).reset_index()
    else:
        # Save item_time_id from Baseline measurements
        item_time_mapping = po[po['time_label'] == 'Baseline'][
            ['item_type_id', 'item_label', 'item_time_id']
        ].drop_duplicates()
        
        # Reshape to wide format
        po_wide = po.pivot_table(
            index=['.draw', 'item_type_id', 'item_label', 'item_type', 'group_label'],
            columns='time_label',
            values='prob'
        ).reset_index()
        
        # Add back item_time_id from Baseline
        po_wide = po_wide.merge(item_time_mapping, on=['item_type_id', 'item_label'], how='left')
    
    po_wide = po_wide.dropna(subset=['Baseline', 'Endline'])
    
    # Compute differences and ratios
    po_wide['diff'] = po_wide['Baseline'] - po_wide['Endline']
    po_wide['ratio'] = 1 - po_wide['Endline'] / po_wide['Baseline']
    po_wide['diff2'] = po_wide['Endline'] - po_wide['Baseline']
    po_wide['ratio2'] = po_wide['Endline'] / po_wide['Baseline'] - 1
    
    # Reshape to long format for summarization
    if endpoint_type == "item_groups":
        id_vars = ['.draw', 'item_type', 'group_label']
    else:
        id_vars = ['.draw', 'item_type_id', 'item_time_id', 'item_label', 'item_type', 'group_label']
    
    po_long = po_wide.melt(
        id_vars=id_vars,
        value_vars=['diff', 'ratio', 'diff2', 'ratio2', 'Baseline', 'Endline'],
        var_name='variable',
        value_name='value'
    )
    
    # Compute quantile summaries
    quantiles = [0.025, 0.25, 0.5, 0.75, 0.975]
    quantile_names = ['q_lower', 'iqr_lower', 'median', 'iqr_upper', 'q_upper']
    
    if endpoint_type == "item_groups":
        group_cols = ['item_type', 'group_label', 'variable']
    else:
        group_cols = ['item_type_id', 'item_time_id', 'item_label', 'item_type', 'group_label', 'variable']
    
    pos = po_long.groupby(group_cols)['value'].quantile(quantiles).unstack()
    pos.columns = quantile_names
    pos = pos.reset_index()
    
    # Merge with item metadata
    if endpoint_type == "item_groups":
        tmp = dit[['item_type', 'group_label', 'group_label_long', 'item_high_label']].drop_duplicates()
        pos = pos.merge(tmp, on=['item_type', 'group_label'])
    else:
        tmp = dit[['item_type', 'item_label', 'item_label_short', 'group_label', 
                  'group_label_long', 'item_high_label']].drop_duplicates()
        pos = pos.merge(tmp, on=['item_type', 'item_label'])
    
    pos_cat1 = pos.copy()
    
    # =========================================================================
    # Part 2: Out-of-7 items - compute means
    # =========================================================================
    
    print("Processing out-of-7 items...")
    
    # Re-extract ordered probability parameters
    po = pd.DataFrame(po_data)
    po['cq_id'] = po['variable'].str.extract(r'\[(\d+)\]')[0].astype(int)
    po = po.drop(columns=['variable'])
    
    # Map cq_id back to item structure
    po = po.merge(cq_mapping, on='cq_id')
    
    # Filter for out-of-7 items and compute means
    tmp = dp1[dp1['item_type'] == 'out-of-7'][
        ['item_type_id', 'item_label', 'item_time_id', 'item_type', 'time_label']
    ].drop_duplicates()
    
    po = po.merge(tmp, on=['item_type_id', 'item_time_id'])
    
    # Compute weighted mean
    po['weighted'] = po['y'] * po['prob']
    po = po.groupby(['.draw', 'item_label', 'item_time_id', 'item_type', 'time_label']).agg({
        'weighted': 'sum'
    }).reset_index()
    po = po.rename(columns={'weighted': 'mean'})
    
    po = po.merge(dit[['item_type', 'item_label', 'group_label']], 
                 on=['item_type', 'item_label'])
    
    # Aggregate across items if endpoint_type="item_groups"
    if endpoint_type == "item_groups":
        po = po.groupby(['.draw', 'item_type', 'time_label', 'group_label']).agg({
            'mean': 'mean'
        }).reset_index()
        
        # Reshape to wide format
        po_wide = po.pivot_table(
            index=['.draw', 'item_type', 'group_label'],
            columns='time_label',
            values='mean'
        ).reset_index()
        
        po_wide = po_wide.dropna(subset=['Baseline', 'Endline'])
        po_wide = po_wide.merge(
            dit[['item_type', 'group_label', 'item_high_label']].drop_duplicates(),
            on=['item_type', 'group_label']
        )
        
        # Compute differences and ratios (direction depends on item_high_label)
        po_wide['diff'] = np.nan
        po_wide['ratio'] = np.nan
        
        lower_mask = po_wide['item_high_label'] == 'lower_is_better'
        po_wide.loc[lower_mask, 'diff'] = po_wide.loc[lower_mask, 'Baseline'] - po_wide.loc[lower_mask, 'Endline']
        po_wide.loc[lower_mask, 'ratio'] = 1 - po_wide.loc[lower_mask, 'Endline'] / po_wide.loc[lower_mask, 'Baseline']
        
        higher_mask = po_wide['item_high_label'] == 'higher_is_better'
        po_wide.loc[higher_mask, 'diff'] = po_wide.loc[higher_mask, 'Endline'] - po_wide.loc[higher_mask, 'Baseline']
        po_wide.loc[higher_mask, 'ratio'] = po_wide.loc[higher_mask, 'Endline'] / po_wide.loc[higher_mask, 'Baseline'] - 1
        
        # Reshape to long format
        id_vars = ['.draw', 'item_type', 'group_label', 'item_high_label']
        po_long = po_wide.melt(
            id_vars=id_vars,
            value_vars=['diff', 'ratio', 'Baseline', 'Endline'],
            var_name='variable',
            value_name='value'
        )
        
        # Compute quantile summaries
        pos2 = po_long.groupby(['item_type', 'group_label', 'variable'])['value'].quantile(quantiles).unstack()
        
        if len(pos2) > 0:
            pos2.columns = quantile_names
            pos2 = pos2.reset_index()
            
            # Merge with item metadata
            tmp = dit[['item_type', 'group_label', 'group_label_long', 'item_high_label']].drop_duplicates()
            pos2 = pos2.merge(tmp, on=['item_type', 'group_label'])
        else:
            print("Warning: No out-of-7 items found for endpoint computation (item_groups mode)")
            pos2 = pd.DataFrame()
    else:
        # Item-level processing
        # Reshape to wide format
        po_wide = po.pivot_table(
            index=['.draw', 'item_label', 'item_time_id', 'item_type', 'group_label'],
            columns='time_label',
            values='mean'
        ).reset_index()
        
        po_wide = po_wide.dropna(subset=['Baseline', 'Endline'])
        po_wide = po_wide.merge(
            dit[['item_type', 'item_label', 'item_high_label']].drop_duplicates(),
            on=['item_type', 'item_label']
        )
        
        # Compute differences and ratios
        po_wide['diff'] = np.nan
        po_wide['ratio'] = np.nan
        
        lower_mask = po_wide['item_high_label'] == 'lower_is_better'
        po_wide.loc[lower_mask, 'diff'] = po_wide.loc[lower_mask, 'Baseline'] - po_wide.loc[lower_mask, 'Endline']
        po_wide.loc[lower_mask, 'ratio'] = 1 - po_wide.loc[lower_mask, 'Endline'] / po_wide.loc[lower_mask, 'Baseline']
        
        higher_mask = po_wide['item_high_label'] == 'higher_is_better'
        po_wide.loc[higher_mask, 'diff'] = po_wide.loc[higher_mask, 'Endline'] - po_wide.loc[higher_mask, 'Baseline']
        po_wide.loc[higher_mask, 'ratio'] = po_wide.loc[higher_mask, 'Endline'] / po_wide.loc[higher_mask, 'Baseline'] - 1
        
        # Reshape to long format
        id_vars = ['.draw', 'item_label', 'item_time_id', 'item_type', 'group_label', 'item_high_label']
        po_long = po_wide.melt(
            id_vars=id_vars,
            value_vars=['diff', 'ratio', 'Baseline', 'Endline'],
            var_name='variable',
            value_name='value'
        )
        
        # Compute quantile summaries
        pos2 = po_long.groupby(['item_label', 'item_time_id', 'item_type', 'group_label', 'variable'])['value'].quantile(quantiles).unstack()
        
        if len(pos2) > 0:
            pos2.columns = quantile_names
            pos2 = pos2.reset_index()
            
            # Merge with item metadata
            tmp = dit[['item_type', 'item_label', 'item_label_short', 'group_label', 
                      'group_label_long', 'item_high_label']].drop_duplicates()
            pos2 = pos2.merge(tmp, on=['item_type', 'item_label'])
        else:
            print("Warning: No out-of-7 items found for endpoint computation")
            pos2 = pd.DataFrame()
    
    # Combine categorical and out-of-7 results
    if len(pos2) > 0:
        pos = pd.concat([pos_cat1, pos2], ignore_index=True)
    else:
        pos = pos_cat1
    
    print(f"Computed endpoints for {len(pos)} item-variable combinations")
    
    return pos
