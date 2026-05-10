"""
Model fitting functions for Bayesian IRT analysis.

This module provides functions for fitting ordered logit models
using ADVI Stan via cmdstanpy.
"""

from typing import Dict, Optional, Union
import os
from pathlib import Path
import warnings
import re
import time

import pandas as pd
import numpy as np
import patsy
import matplotlib.pyplot as plt
from cmdstanpy import CmdStanModel
from plotnine import *
import arviz as az

def fit_ordered_logit_model_ncats_advi(
    dit: pd.DataFrame,
    dcati: pd.DataFrame,
    output_file_prefix: str,
    stan_file: Optional[str] = None,
    x_formula: str = "~ time - 1",
    x_formula_ignore_regex: Optional[str] = None,
    iter: int = 10000,
    grad_samples: int = 1,
    elbo_samples: int = 100,
    output_samples: int = 4000,
    seed: int = 123,
    show_messages: bool = False,
    show_exceptions: bool = False,
    resume: bool = False,
    with_core_analyses: bool = True,
    with_additional_analyses: bool = False,
) -> Dict:
    """
    Run ordered logit model analysis (ncats version) using ADVI.
    
    This function performs Bayesian IRT analysis using the flexible ncats ordered logit model
    with ADVI (Automatic Differentiation Variational Inference) instead of HMC.
    
    Parameters
    ----------
    dit : pd.DataFrame
        Item metadata table with item types and labels.
    dcati : pd.DataFrame
        Pre-processed data with observations for analysis. Must include
        item_type, y_stan, item_time_id, pid, and oidt columns.
    output_file_prefix : str
        Full path prefix for output files (without extension).
    stan_file : str, optional
        Path to Stan model file (.stan). If None, uses default ordered_logit_ncats_v260413.stan
    x_formula : str, default "~ time - 1"
        Patsy formula string specifying predictors for the design matrix.
    x_formula_ignore_regex : str, optional
        Regular expression pattern to identify X columns to exclude from posterior predictive.
    iter : int, default 10000
        Number of ADVI iterations.
    grad_samples : int, default 1
        Number of samples for ELBO gradient estimation.
    elbo_samples : int, default 100
        Number of samples for ELBO evaluation.
    output_samples : int, default 4000
        Number of samples to draw from the approximate posterior.
    seed : int, default 123
        Random seed for reproducibility.
    show_messages : bool, default False
        If True, show Stan informational messages.
    show_exceptions : bool, default False
        If True, show Stan exception messages when errors occur.
    resume : bool, default False
        If True and all output files exist, skip ADVI and load existing results.
    with_core_analyses : bool, default True
        If True, generate core probability plots.
    with_additional_analyses : bool, default False
        If True, generate additional diagnostic plots.
    
    Returns
    -------
    dict
        Dictionary containing:
        - 'fit': CmdStanVB object (if not resumed)
        - 'draws': Draws array from approximate posterior
        - 'timing': DataFrame with timing information
    """
    # Set default stan_file if not provided
    if stan_file is None:
        # Assume we're running from the bIRTistic root directory
        stan_file = str(Path(__file__).parents[2] / "src" / "stan" / "ordered_logit_ncats_v260413.stan")
    
    # Print configuration
    print("\n" + "=" * 40)
    print("Ordered Logit Model (ncats) ADVI Analysis Configuration")
    print("=" * 40)
    print(f"Stan file: {stan_file}")
    print(f"Stan include dir: {os.path.dirname(stan_file)}")
    print(f"Data: dit with {len(dit)} items, dcati with {len(dcati)} observations")
    print(f"Output prefix: {output_file_prefix}")
    print(f"ADVI iterations: {iter}")
    print(f"Gradient samples: {grad_samples}")
    print(f"ELBO samples: {elbo_samples}")
    print(f"Output samples: {output_samples}")
    print(f"Seed: {seed}")
    print(f"Resume: {resume}")
    print(f"Core analyses: {with_core_analyses}")
    print(f"Additional analyses: {with_additional_analyses}")
    print("=" * 40 + "\n")
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file_prefix), exist_ok=True)
    
    # Define output file paths
    timing_file = f"{output_file_prefix}_timing.csv"
    draws_file = f"{output_file_prefix}_draws.zarr"
    output_file = f"{output_file_prefix}_stan.pkl"
    data_file = f"{output_file_prefix}_data.csv"
    
    # Save the input data for later use
    dcati.to_csv(data_file.replace('.csv', '_dp1.csv'), index=False)
    dit.to_csv(data_file.replace('.csv', '_dit.csv'), index=False)
    print(f"Saved preprocessed data to: {data_file.replace('.csv', '_dp1.csv')} and _dit.csv")
    
    # Check if we should resume from existing outputs
    can_resume = (resume and 
                  os.path.exists(timing_file) and 
                  os.path.exists(draws_file) and 
                  os.path.exists(data_file.replace('.csv', '_dp1.csv')) and
                  os.path.exists(data_file.replace('.csv', '_dit.csv')) and
                  os.path.exists(output_file))
    
    # Prepare data in ncats format (same as HMC version)
    print("Preparing Stan data in ncats format...")
    
    # Ensure dcati is sorted by item_type_id and oidt
    dcati_sorted = dcati.sort_values(['item_type_id', 'oidt']).reset_index(drop=True)
    
    # Verify oid is sequential
    assert (dcati_sorted['oid'] == dcati_sorted.index + 1).all(), \
        "oid must be sequential 1:N after sorting by item_type_id, oidt"
    
    # Build stan_data dictionary
    stan_data = {}
    stan_data['C'] = int(dcati_sorted['item_type_id'].max())
    stan_data['U'] = int(dcati_sorted['pid'].max())
    
    # Arrays for each category type
    category_stats = dcati.groupby('item_type_id').agg({
        'oid': 'count',
        'item_time_id': 'max',
        'y_stan': lambda x: len(x.unique())
    }).rename(columns={'oid': 'N', 'item_time_id': 'Q', 'y_stan': 'K'})
    
    stan_data['N'] = category_stats['N'].astype(int).tolist()
    stan_data['Q'] = category_stats['Q'].astype(int).tolist()
    stan_data['K'] = category_stats['K'].astype(int).tolist()
    
    # Total counts
    stan_data['N_total'] = int(sum(stan_data['N']))
    stan_data['Q_total'] = int(sum(stan_data['Q']))
    
    # Concatenated observation data
    stan_data['y'] = dcati['y_stan'].astype(int).tolist()
    stan_data['unit_of_obs'] = dcati['pid'].astype(int).tolist()
    stan_data['question_of_obs'] = dcati['item_time_id'].astype(int).tolist()
    stan_data['cat_type'] = dcati['item_type_id'].astype(int).tolist()
    
    # Create design matrix using patsy
    design_matrix = patsy.dmatrix(x_formula, data=dcati, return_type='dataframe')
    
    # Determine which columns to use for predictive probabilities
    if x_formula_ignore_regex is not None:
        xpr_id = [i + 1 for i, col in enumerate(design_matrix.columns) 
                  if not re.search(x_formula_ignore_regex, col)]
    else:
        xpr_id = list(range(1, len(design_matrix.columns) + 1))
    
    stan_data['X'] = design_matrix.values.tolist()
    stan_data['Xpr_id'] = xpr_id
    stan_data['P'] = len(design_matrix.columns)
    stan_data['P_pr'] = len(xpr_id)
    
    print(f"Design matrix number of predictors: P = {stan_data['P']}")
    print(f"Design matrix column names (for fitting): {', '.join(design_matrix.columns)}")
    print(f"Design matrix column names (for prediction): {', '.join([design_matrix.columns[i-1] for i in xpr_id])}")
    print(f"Number of category types: C = {stan_data['C']}")
    print(f"Category type labels: {', '.join(dit['item_type'].unique())}")
    print(f"N per category: {', '.join(map(str, stan_data['N']))}")
    print(f"Q per category: {', '.join(map(str, stan_data['Q']))}")
    print(f"K per category: {', '.join(map(str, stan_data['K']))}")
    
    # Check if resuming from existing outputs
    if can_resume:
        print("\n" + "=" * 40)
        print("RESUMING from existing outputs")
        print("=" * 40)
        
        print(f"Loading draws from: {draws_file}")
        idata = az.from_zarr(draws_file)
        
        print(f"Loading timing data from: {timing_file}")
        timing_data = pd.read_csv(timing_file)
        print("=" * 40 + "\n")
        
        return {
            'draws': idata,
            'timing': timing_data
        }
    else:
        if resume:
            print("\nNote: Resume requested but not all output files exist. Running full analysis.\n")
        
        # Compile Stan model
        print("Compiling Stan model...")
        model = CmdStanModel(
            stan_file=stan_file,
            cpp_options={'STAN_THREADS': 'TRUE'}
        )
        
        # Run ADVI
        print("Running ADVI...")
        start_time = time.time()
        
        fit = model.variational(
            data=stan_data,
            algorithm='meanfield',
            iter=iter,
            grad_samples=grad_samples,
            elbo_samples=elbo_samples,
            output_samples=output_samples,
            seed=seed,
            show_console=show_messages
        )
        
        elapsed_time = time.time() - start_time
        
        # Extract timing information
        print("\nExtracting timing information...")
        timing_data = pd.DataFrame({
            'total_minutes': [elapsed_time / 60]
        })
        
        timing_data.to_csv(timing_file, index=False)
        print(f"Saved timing information to: {timing_file}")
        
        # Save model fit
        print(f"Saving model fit to: {output_file}")
        pd.to_pickle(fit, output_file)
        
        # Convert to ArviZ InferenceData and save
        print("Converting to ArviZ format...")
        # For ADVI, we need to manually extract draws since az.from_cmdstanpy() 
        # doesn't support variational fits (no chain structure)
        # Get the variational sample as a DataFrame
        try:
            draws_df = fit.variational_sample
            if isinstance(draws_df, np.ndarray):
                # If it's a numpy array, convert to DataFrame
                column_names = fit.column_names
                draws_df = pd.DataFrame(draws_df, columns=column_names)
        except Exception as e:
            # Fallback: extract as DataFrame directly
            draws_df = fit.draws_pd()
        
        # Create InferenceData from ADVI draws while preserving Stan parameter structure.
        # CmdStan variational draws flatten array parameters into columns like var[1], var[2], ...
        # Reconstruct those into a single variable with explicit dimensions for ArviZ.
        n_draws = len(draws_df)
        col_pattern = re.compile(r"^(?P<name>[^\[]+)(?:\[(?P<idx>[0-9,]+)\])?$")

        grouped_cols = {}
        for col in draws_df.columns:
            match = col_pattern.match(col)
            if match is None:
                raise ValueError(f"Unrecognized draw column format: {col}")

            base_name = match.group("name")
            idx_str = match.group("idx")
            idx_tuple = None if idx_str is None else tuple(int(i) for i in idx_str.split(","))
            grouped_cols.setdefault(base_name, []).append((idx_tuple, col))

        posterior_dict = {}
        dims = {}

        for base_name, entries in grouped_cols.items():
            first_idx = entries[0][0]

            if first_idx is None:
                # Scalar parameter, shape -> (chain, draw)
                posterior_dict[base_name] = draws_df[entries[0][1]].to_numpy().reshape(1, n_draws)
                continue

            ndim = len(first_idx)
            if any(idx is None for idx, _ in entries):
                raise ValueError(f"Mixed scalar/indexed columns for parameter '{base_name}'")
            if any(len(idx) != ndim for idx, _ in entries):
                raise ValueError(f"Inconsistent index dimensionality for parameter '{base_name}'")

            shape = tuple(max(idx[d] for idx, _ in entries) for d in range(ndim))
            arr = np.full((n_draws, *shape), np.nan, dtype=float)

            for idx, col in entries:
                arr[(slice(None),) + tuple(i - 1 for i in idx)] = draws_df[col].to_numpy()

            if np.isnan(arr).any():
                raise ValueError(f"Incomplete indexed draws for parameter '{base_name}'")

            posterior_dict[base_name] = arr[np.newaxis, ...]
            dims[base_name] = [f"{base_name}_dim_{i}" for i in range(ndim)]

        idata = az.from_dict(posterior=posterior_dict, dims=dims)
        
        print(f"Saving draws to: {draws_file}")
        idata.to_zarr(draws_file)
        
        return {
            'fit': fit,
            'draws': idata,
            'timing': timing_data
        }
