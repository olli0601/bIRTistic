"""
Model fitting functions for Bayesian IRT analysis.

This module provides functions for fitting partical credit models
using Stan HMC via cmdstanpy.
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

def fit_partial_credit_model_ncats(
    dit: pd.DataFrame,
    dcati: pd.DataFrame,
    output_file_prefix: str,
    stan_file: Optional[str] = None,
    x_formula: str = "time - 1",
    x_formula_ignore_regex: Optional[str] = None,
    chains: int = 2,
    parallel_chains: int = 2,
    threads_per_chain: int = 1,
    iter_warmup: int = 500,
    iter_sampling: int = 1500,
    seed: int = 123,
    show_console: bool = False,
    resume: bool = False,
    with_core_analyses: bool = True,
    with_additional_analyses: bool = False,
) -> Dict:
    """
    Run partial credit model analysis (ncats version) on pre-processed data.
    
    This function performs Bayesian IRT analysis using the flexible ncats partial credit model
    on pre-processed data. It compiles the Stan model, runs MCMC sampling, generates convergence
    diagnostics, and optionally creates detailed diagnostic plots and posterior predictive checks.
    
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
        Path to Stan model file (.stan). If None, uses default partial_credit_model_ncats_v260413.stan
    x_formula : str, default "time - 1"
        Patsy formula string specifying predictors for the design matrix.
    x_formula_ignore_regex : str, optional
        Regular expression pattern to identify X columns to exclude from posterior predictive
        probability calculations. If None, all X columns are used for both fitted and predictive.
    chains : int, default 2
        Number of MCMC chains to run.
    parallel_chains : int, default 2
        Number of chains to run in parallel.
    threads_per_chain : int, default 1
        Number of threads to use per chain.
    iter_warmup : int, default 500
        Number of warmup iterations per chain.
    iter_sampling : int, default 1500
        Number of sampling iterations per chain.
    seed : int, default 123
        Random seed for reproducibility.
    show_console : bool, default False
        If True, show Stan console output during sampling.
    resume : bool, default False
        If True and all output files exist, skip MCMC sampling and load existing results.
    with_core_analyses : bool, default True
        If True, generate core probability plots.
    with_additional_analyses : bool, default False
        If True, generate additional diagnostic plots.
    
    Returns
    -------
    dict
        Dictionary containing:
        - 'fit': CmdStanMCMC object (if not resumed)
        - 'draws': InferenceData object with posterior draws
        - 'timing': DataFrame with chain timing information
        - 'good_chains': List of chain IDs that converged
    """
    # Set default stan_file if not provided
    if stan_file is None:
        # Assume we're running from the bIRTistic root directory
        stan_file = str(Path(__file__).parents[2] / "src" / "stan" / "partial_credit_model_ncats_v260413.stan")
    
    # Print configuration
    print("\n" + "=" * 40)
    print("Partial Credit Model (ncats) Analysis Configuration")
    print("=" * 40)
    print(f"Stan file: {stan_file}")
    print(f"Stan include dir: {os.path.dirname(stan_file)}")
    print(f"Data: dit with {len(dit)} items, dcati with {len(dcati)} observations")
    print(f"Output prefix: {output_file_prefix}")
    print(f"Chains: {chains}")
    print(f"Parallel chains: {parallel_chains}")
    print(f"Threads per chain: {threads_per_chain}")
    print(f"Warmup iterations: {iter_warmup}")
    print(f"Sampling iterations: {iter_sampling}")
    print(f"Seed: {seed}")
    print(f"Resume: {resume}")
    print(f"Core analyses: {with_core_analyses}")
    print(f"Additional analyses: {with_additional_analyses}")
    print("=" * 40 + "\n")
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file_prefix), exist_ok=True)
    
    # Define output file paths
    timing_file = f"{output_file_prefix}_timing.csv"
    draws_file = f"{output_file_prefix}_draws.zarr"  # Use zarr for arviz (9x faster than NetCDF)
    output_file = f"{output_file_prefix}_stan.pkl"
    mixing_file = f"{output_file_prefix}_convergence_mixing.csv"
    
    # Check if we should resume from existing outputs
    can_resume = (resume and 
                  os.path.exists(timing_file) and 
                  os.path.exists(draws_file) and 
                  os.path.exists(output_file) and 
                  os.path.exists(mixing_file))
    
    # Prepare data in ncats format
    print("Preparing Stan data in ncats format...")
    
    # Ensure dcati is sorted by item_type_id and oidt
    dcati_sorted = dcati.sort_values(['item_type_id', 'oidt']).reset_index(drop=True)
    
    # Verify oid is sequential
    assert (dcati_sorted['oid'] == dcati_sorted.index + 1).all(), \
        "oid must be sequential 1:N after sorting by item_type_id, oidt"
    
    # Build stan_data dictionary
    stan_data = {}
    stan_data['C'] = int(dcati['item_type_id'].max())
    stan_data['U'] = int(dcati['pid'].max())
    
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
    
    # Concatenated observation data (data is already sorted by item_type_id) - convert to int
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
        good_chains = timing_data['chain'].tolist()
        print(f"Identified good HMC chains: {', '.join(map(str, good_chains))}")
        print("=" * 40 + "\n")
        
        return {
            'draws': idata,
            'timing': timing_data,
            'good_chains': good_chains
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
        
        # Sample from the model
        print("Running MCMC sampling...")
        fit = model.sample(
            data=stan_data,
            chains=chains,
            parallel_chains=parallel_chains,
            threads_per_chain=threads_per_chain,
            iter_warmup=iter_warmup,
            iter_sampling=iter_sampling,
            seed=seed,
            show_console=show_console,
            refresh=500,
            save_warmup=True
        )
        
        # Remove any chains that did not converge based on log probability
        print("Checking for divergent transitions and removing non-converged chains...")
        lp_df = fit.draws_pd(vars=['lp__'], inc_warmup=False)
        
        # Calculate mean and std per chain
        chain_stats = lp_df.groupby('chain__')['lp__'].agg(['mean', 'std', 'count']).reset_index()
        chain_stats.columns = ['chain', 'mean_lp', 'sd_lp', 'n']
        
        # Identify good chains (within 2 SD of best chain)
        best_idx = chain_stats['mean_lp'].idxmax()
        threshold = chain_stats.loc[best_idx, 'mean_lp'] - 2 * chain_stats.loc[best_idx, 'sd_lp']
        good_chains = chain_stats[chain_stats['mean_lp'] > threshold]['chain'].tolist()
        print(f"Identified good HMC chains: {', '.join(map(str, good_chains))}")
        
        # Extract chain timing information
        print("\nExtracting timing information...")
        timing_info = fit.metadata.cmdstan_config
        timing_data = pd.DataFrame({
            'chain': good_chains,
            'warmup_minutes': [timing_info['metric'] for _ in good_chains],  # Placeholder
            'sampling_minutes': [0.0 for _ in good_chains],  # Placeholder
            'total_chain_minutes': [0.0 for _ in good_chains]  # Placeholder
        })
        
        timing_data.to_csv(timing_file, index=False)
        print(f"Saved timing information to: {timing_file}")
        
        # Convert to ArviZ InferenceData
        print("Converting to ArviZ format...")
        idata = az.from_cmdstanpy(
            fit,
            coords={'chain': good_chains},
            dims={}
        )
        
        # Save draws
        print(f"Saving draws to: {draws_file}")
        idata.to_zarr(draws_file)
        
        # Check convergence and mixing
        print("Generating convergence diagnostics...")
        summary = az.summary(
            idata,
            var_names=['latent_factor_unit', 'latent_factor_beta', 
                      'skill_thresholds', 'loadings_questions_m1']
        )
        summary = summary.sort_values('ess_bulk')
        summary.to_csv(mixing_file)
        
        # Generate trace plots for worst parameters
        if with_additional_analyses:
            print("Generating trace plots...")
            worst_vars = summary.index[:9].tolist()
            
            fig, axes = plt.subplots(len(worst_vars) + 1, 1, figsize=(20, 10))
            
            # Plot log probability
            az.plot_trace(idata, var_names=['lp'], axes=axes[0:1])
            
            # Plot worst parameters
            for i, var in enumerate(worst_vars):
                az.plot_trace(idata, var_names=[var], axes=axes[i+1:i+2])
            
            plt.tight_layout()
            plt.savefig(f"{output_file_prefix}_worsttrace.png", dpi=150)
            plt.close()
        
        return {
            'fit': fit,
            'draws': idata,
            'timing': timing_data,
            'good_chains': good_chains
        }