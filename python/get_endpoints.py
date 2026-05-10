"""
Analysis functions for Bayesian IRT models.

This module provides functions for analyzing posterior draws from IRT models,
specifically endpoint extraction and effect size computation.
"""

from typing import Optional, Literal
import os
import re

import pandas as pd
import numpy as np


def get_endpoints(
    dp1: pd.DataFrame,
    dit: pd.DataFrame,
    fit_file: str,
    draws_file: str,
    categorical_threshold: int = 3,
    endpoint_type: Literal["items", "item_groups"] = "items",
    param_name: str = "ordered_prob_by_cat_qu_pr"
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

    import arviz as az
    
    # Read samples using ArviZ
    idata = az.from_zarr(draws_file)    
    tmp = [v for v in idata.posterior.data_vars if v == param_name]
    if not tmp:
        raise ValueError(f"No data_vars '{param_name}'. Available: {list(idata.posterior.data_vars)}")
                                
    # %%                                
    # Reshape to dataframe with columns: .draw, cq_id, prob
    po = idata.posterior[param_name].values
    n_chains, n_draws, n_cq = po.shape
    po = po.reshape(-1, n_cq)
    po = pd.DataFrame.from_records(po)
    po = po.melt(var_name='cq_id', value_name='prob', ignore_index=False).reset_index().rename(columns={"index":".draw"})	

    # Map cq_id back to item structure
    tmp = dp1[['item_type_id', 'item_label', 'item_time_id']].drop_duplicates()
    tmp = tmp.merge(dit[['item_type_id', 'item_label', 'cat_length']], 
                   on=['item_type_id', 'item_label'])
    tmp = tmp.sort_values(['item_type_id', 'item_time_id']).reset_index(drop=True)
    tmp['cq_id'] = tmp['cat_length'].cumsum() - tmp['cat_length']
    tmp = tmp.assign(y=tmp['cat_length'].map(lambda n: np.arange(n))).explode('y', ignore_index=True)
    tmp['y'] = tmp['y'].astype(int)
    tmp['cq_id'] = tmp['cq_id'] + tmp['y']
    tmp.drop(columns=['cat_length','item_label'], inplace=True)
    po = po.merge(tmp, on='cq_id')

    # =========================================================================
    # Part 1: Categorical items - aggregate probabilities for high categories
    # =========================================================================
    
    print("Processing categorical items...")
    
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
        id_vars = ['item_type_id', 'item_type', 'group_label']
        po = po.groupby(['.draw', 'item_type', 'time_label', 'group_label']).agg({
            'prob': 'mean'
        }).reset_index()
    else:                
        id_vars = ['item_type_id', 'item_type', 'item_label', 'group_label']
    
    # Reshape to wide format (Baseline vs Endline)    
    po = po.pivot_table(
            index=['.draw'] + id_vars,
            columns='time_label',
            values='prob'
        ).reset_index()            
    po = po.dropna(subset=['Baseline', 'Endline'])
    
    # Compute differences and ratios
    po['diff'] = po['Baseline'] - po['Endline']
    po['ratio'] = 1 - po['Endline'] / po['Baseline']
    po['diff2'] = po['Endline'] - po['Baseline']
    po['ratio2'] = po['Endline'] / po['Baseline'] - 1
    
    # Reshape to long format for summarization    
    po = po.melt(
        id_vars=['.draw'] + id_vars,
        value_vars=['diff', 'ratio', 'diff2', 'ratio2', 'Baseline', 'Endline'],
        var_name='variable',
        value_name='value'
    )
    
    # Compute quantile summaries
    quantiles = [0.025, 0.25, 0.5, 0.75, 0.975]
    quantile_names = ['q_lower', 'iqr_lower', 'median', 'iqr_upper', 'q_upper']    
    tmp = id_vars + ['variable']
    pos = po.groupby(tmp)['value'].quantile(quantiles).unstack()
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
    # %%                                

    # =========================================================================
    # Part 2: Out-of-7 items - compute means
    # =========================================================================
    
    # Reshape to dataframe with columns: .draw, cq_id, prob
    po = idata.posterior[param_name].values
    n_chains, n_draws, n_cq = po.shape
    po = po.reshape(-1, n_cq)
    po = pd.DataFrame.from_records(po)
    po = po.melt(var_name='cq_id', value_name='prob', ignore_index=False).reset_index().rename(columns={"index":".draw"})	

    # Map cq_id back to item structure
    tmp = dp1[['item_type_id', 'item_label', 'item_time_id']].drop_duplicates()
    tmp = tmp.merge(dit[['item_type_id', 'item_label', 'cat_length']], 
                   on=['item_type_id', 'item_label'])
    tmp = tmp.sort_values(['item_type_id', 'item_time_id']).reset_index(drop=True)
    tmp['cq_id'] = tmp['cat_length'].cumsum() - tmp['cat_length']
    tmp = tmp.assign(y=tmp['cat_length'].map(lambda n: np.arange(n))).explode('y', ignore_index=True)
    tmp['y'] = tmp['y'].astype(int)
    tmp['cq_id'] = tmp['cq_id'] + tmp['y']
    tmp.drop(columns=['cat_length','item_label'], inplace=True)
    po = po.merge(tmp, on='cq_id')

    print("Processing out-of-7 items...")
    
    # Filter for out-of-7 items
    tmp = dp1[dp1['item_type'] == 'out-of-7'][
        ['item_type_id', 'item_label', 'item_time_id', 'item_type', 'time_label']
    ].drop_duplicates()
    po = po.merge(tmp, on=['item_type_id', 'item_time_id'])
    
    # Compute weighted mean
    po['mean'] = po['y'] * po['prob']
    po = po.groupby(['.draw', 'item_type_id', 'item_label', 'item_time_id', 'item_type', 'time_label']).agg({
        'mean': 'sum'
    }).reset_index()
        
    po = po.merge(dit[['item_type', 'item_label', 'group_label']], on=['item_type', 'item_label'])
    
    # Aggregate across items if endpoint_type="item_groups"
    if endpoint_type == "item_groups":
        id_vars = ['item_type_id', 'item_type', 'group_label']
        po = po.groupby(['.draw', 'item_type', 'time_label', 'group_label']).agg({
            'prob': 'mean'
        }).reset_index()
    else:                
        id_vars = ['item_type_id', 'item_type', 'item_label', 'group_label']
    
    # Reshape to wide format
    po = po.pivot_table(
        index=['.draw'] + id_vars,
        columns='time_label',
        values='mean'
    ).reset_index()
        
    po = po.dropna(subset=['Baseline', 'Endline'])
    po = po.merge(
            dit[['item_type', 'group_label', 'item_high_label']].drop_duplicates(),
            on=['item_type', 'group_label']
    )
    
    # Compute differences and ratios (direction depends on item_high_label)
    po['diff'] = np.nan
    po['ratio'] = np.nan        
    tmp = po['item_high_label'] == 'lower_is_better'
    po.loc[tmp, 'diff'] = po.loc[tmp, 'Baseline'] - po.loc[tmp, 'Endline']
    po.loc[tmp, 'ratio'] = 1 - po.loc[tmp, 'Endline'] / po.loc[tmp, 'Baseline']        
    tmp = po['item_high_label'] == 'higher_is_better'
    po.loc[tmp, 'diff'] = po.loc[tmp, 'Endline'] - po.loc[tmp, 'Baseline']
    po.loc[tmp, 'ratio'] = po.loc[tmp, 'Endline'] / po.loc[tmp, 'Baseline'] - 1

    # Reshape to long format
    po = po.melt(
        id_vars=['.draw'] + id_vars,
        value_vars=['diff', 'ratio', 'Baseline', 'Endline'],
        var_name='variable',
        value_name='value'
    )
        
    # Compute quantile summaries
    tmp = id_vars + ['variable']
    pos = po.groupby(tmp)['value'].quantile(quantiles).unstack()
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
        
    # %%                                
    # Combine categorical and out-of-7 results
    pos = pd.concat([pos_cat1, pos], ignore_index=True)
    
    print(f"Computed endpoints for {len(pos)} item-variable combinations")
    
    return pos
