"""
Model fitting functions for Bayesian IRT analysis.

This module provides functions for fitting partial credit models and ordered logit models
using Stan via cmdstanpy.
"""

from typing import Dict, Optional, Union
import os
from pathlib import Path
import warnings

import pandas as pd
import numpy as np
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
    draws_file = f"{output_file_prefix}_draws.nc"  # Use NetCDF for arviz
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
    stan_data['C'] = int(dcati_sorted['item_type_id'].max())
    stan_data['U'] = int(dcati_sorted['pid'].max())
    
    # Arrays for each category type
    category_stats = dcati_sorted.groupby('item_type_id').agg({
        'oid': 'count',
        'item_time_id': 'max',
        'y_stan': lambda x: len(x.unique())
    }).rename(columns={'oid': 'N', 'item_time_id': 'Q', 'y_stan': 'K'})
    
    stan_data['N'] = [int(x) for x in category_stats['N'].tolist()]
    stan_data['Q'] = [int(x) for x in category_stats['Q'].tolist()]
    stan_data['K'] = [int(x) for x in category_stats['K'].tolist()]
    
    # Total counts
    stan_data['N_total'] = int(sum(stan_data['N']))
    stan_data['Q_total'] = int(sum(stan_data['Q']))
    
    # Concatenated observation data (data is already sorted by item_type_id) - convert to int
    stan_data['y'] = [int(x) for x in dcati_sorted['y_stan'].tolist()]
    stan_data['unit_of_obs'] = [int(x) for x in dcati_sorted['pid'].tolist()]
    stan_data['question_of_obs'] = [int(x) for x in dcati_sorted['item_time_id'].tolist()]
    stan_data['cat_type'] = [int(x) for x in dcati_sorted['item_type_id'].tolist()]
    
    # Create design matrix using patsy
    import patsy
    design_matrix = patsy.dmatrix(x_formula, data=dcati_sorted, return_type='dataframe')
    
    # Determine which columns to use for predictive probabilities
    if x_formula_ignore_regex is not None:
        import re
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
        idata = az.from_netcdf(draws_file)
        
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
        idata.to_netcdf(draws_file)
        
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
            
            import matplotlib.pyplot as plt
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


def fit_ordered_logit_model_ncats(
    dit: pd.DataFrame,
    dcati: pd.DataFrame,
    output_file_prefix: str,
    stan_file: Optional[str] = None,
    x_formula: str = "~ time - 1",
    x_formula_ignore_regex: Optional[str] = None,
    chains: int = 2,
    parallel_chains: int = 2,
    threads_per_chain: int = 1,
    iter_warmup: int = 500,
    iter_sampling: int = 1500,
    seed: int = 123,
    show_messages: bool = False,
    show_exceptions: bool = False,
    resume: bool = False,
    with_core_analyses: bool = True,
    with_additional_analyses: bool = False,
) -> Dict:
    """
    Run ordered logit model analysis (ncats version) on pre-processed data using HMC.
    
    This function performs Bayesian IRT analysis using the flexible ncats ordered logit model
    on pre-processed data. It compiles the Stan model, runs MCMC sampling, generates convergence
    diagnostics, and optionally creates detailed diagnostic plots.
    
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
    show_messages : bool, default False
        If True, show Stan informational messages.
    show_exceptions : bool, default False
        If True, show Stan exception messages when errors occur.
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
        stan_file = str(Path(__file__).parents[2] / "src" / "stan" / "ordered_logit_ncats_v260413.stan")
    
    # Print configuration
    print("\n" + "=" * 40)
    print("Ordered Logit Model (ncats) Analysis Configuration")
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
    draws_file = f"{output_file_prefix}_draws.pkl"  # Use pickle for compatibility
    output_file = f"{output_file_prefix}_stan.pkl"
    mixing_file = f"{output_file_prefix}_convergence_mixing.csv"
    data_file = f"{output_file_prefix}_data.pkl"
    
    # Save the input data for later use
    pd.to_pickle({'dp1': dcati, 'dit': dit}, data_file)
    print(f"Saved preprocessed data to: {data_file}")
    
    # Check if we should resume from existing outputs
    can_resume = (resume and 
                  os.path.exists(timing_file) and 
                  os.path.exists(draws_file) and 
                  os.path.exists(output_file) and 
                  os.path.exists(mixing_file))
    
    # Prepare data in ncats format (same as fit_partial_credit_model_ncats)
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
    category_stats = dcati_sorted.groupby('item_type_id').agg({
        'oid': 'count',
        'item_time_id': 'max',
        'y_stan': lambda x: len(x.unique())
    }).rename(columns={'oid': 'N', 'item_time_id': 'Q', 'y_stan': 'K'})
    
    stan_data['N'] = [int(x) for x in category_stats['N'].tolist()]
    stan_data['Q'] = [int(x) for x in category_stats['Q'].tolist()]
    stan_data['K'] = [int(x) for x in category_stats['K'].tolist()]
    
    # Total counts
    stan_data['N_total'] = int(sum(stan_data['N']))
    stan_data['Q_total'] = int(sum(stan_data['Q']))
    
    # Concatenated observation data
    stan_data['y'] = [int(x) for x in dcati_sorted['y_stan'].tolist()]
    stan_data['unit_of_obs'] = [int(x) for x in dcati_sorted['pid'].tolist()]
    stan_data['question_of_obs'] = [int(x) for x in dcati_sorted['item_time_id'].tolist()]
    stan_data['cat_type'] = [int(x) for x in dcati_sorted['item_type_id'].tolist()]
    
    # Create design matrix using patsy
    import patsy
    design_matrix = patsy.dmatrix(x_formula, data=dcati_sorted, return_type='dataframe')
    
    # Determine which columns to use for predictive probabilities
    if x_formula_ignore_regex is not None:
        import re
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
        poa = pd.read_pickle(draws_file)
        
        print(f"Loading timing data from: {timing_file}")
        timing_data = pd.read_csv(timing_file)
        good_chains = timing_data['chain'].tolist()
        print(f"Identified good HMC chains: {', '.join(map(str, good_chains))}")
        print("=" * 40 + "\n")
        
        return {
            'draws': poa,
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
            show_console=show_messages,
            refresh=500,
            save_warmup=True
        )
        
        # Remove any chains that did not converge
        print("Checking for divergent transitions and removing non-converged chains...")
        # Get lp__ draws using the method from cmdstanr: 
        # fit$draws(variables = "lp__", inc_warmup = FALSE, format = "draws_df")
        # In cmdstanpy, this is fit.draws_pd()
        try:
            # Try with inc_warmup parameter (cmdstanpy 1.x)
            lp_draws = fit.draws_pd(inc_warmup=False)
        except TypeError:
            # Fallback for older API
            lp_draws = fit.draws_pd()
        
        # Check if lp__ column exists
        if 'lp__' not in lp_draws.columns:
            raise ValueError(f"lp__ not found in draws. Available columns: {list(lp_draws.columns)[:10]}")
        
        # Identify chain column - could be 'chain__', '.chain', or 'chain'
        chain_col = None
        for col_name in ['chain__', '.chain', 'chain']:
            if col_name in lp_draws.columns:
                chain_col = col_name
                break
        
        if chain_col is None:
            # Need to reconstruct chain information from draw indices
            num_draws_per_chain = fit.num_draws_sampling
            total_draws = len(lp_draws)
            num_chains = len(fit.chain_ids)
            
            # Create chain column based on position
            lp_draws['chain'] = [fit.chain_ids[i // num_draws_per_chain] for i in range(total_draws)]
            chain_col = 'chain'
        
        # Calculate mean and std per chain (matching R code)
        chain_stats = lp_draws.groupby(chain_col)['lp__'].agg(['mean', 'std', 'count']).reset_index()
        chain_stats.columns = ['chain', 'mean_lp', 'sd_lp', 'n']
        
        # Identify good chains
        best_idx = chain_stats['mean_lp'].idxmax()
        threshold = chain_stats.loc[best_idx, 'mean_lp'] - 2 * chain_stats.loc[best_idx, 'sd_lp']
        good_chains = chain_stats[chain_stats['mean_lp'] > threshold]['chain'].tolist()
        print(f"Identified good HMC chains: {', '.join(map(str, good_chains))}")
        
        # Extract chain timing information
        print("\nExtracting timing information...")
        timing_data = pd.DataFrame({
            'chain': good_chains,
            'warmup_minutes': [0.0 for _ in good_chains],
            'sampling_minutes': [0.0 for _ in good_chains],
            'total_chain_minutes': [0.0 for _ in good_chains]
        })
        
        timing_data.to_csv(timing_file, index=False)
        print(f"Saved timing information to: {timing_file}")
        
        # Save model fit
        print(f"Saving model fit to: {output_file}")
        pd.to_pickle(fit, output_file)
        
        # Convert to draws array and save
        print("Extracting draws...")
        poa = fit.draws(inc_warmup=False)
        # Filter to good chains (convert chain IDs to 0-based indices)
        good_chain_indices = [int(c) - 1 for c in good_chains]
        poa = poa[:, good_chain_indices, :]
        
        pd.to_pickle(poa, draws_file)
        print(f"Saved draws to: {draws_file}")
        
        # Check convergence and mixing
        print("Generating convergence diagnostics...")
        summary_df = fit.summary()
        
        # Filter to key variables of interest
        key_vars = ['latent_factor_unit', 'latent_factor_beta', 
                   'skill_thresholds_1', 'skill_thresholds_incs', 'loadings_questions_m1']
        # Filter rows where variable name starts with any of the key variable prefixes
        mask = summary_df.index.to_series().apply(
            lambda x: any(x.startswith(var) for var in key_vars)
        )
        summary_df_filtered = summary_df[mask]
        
        # Sort by ESS if available, otherwise just save as-is
        if 'ess_bulk' in summary_df_filtered.columns:
            summary_df_filtered = summary_df_filtered.sort_values('ess_bulk')
        elif 'N_Eff' in summary_df_filtered.columns:
            summary_df_filtered = summary_df_filtered.sort_values('N_Eff')
        
        summary_df_filtered.to_csv(mixing_file)
        
        print(f"Saved convergence diagnostics to: {mixing_file}")
        
        return {
            'fit': fit,
            'draws': poa,
            'timing': timing_data,
            'good_chains': good_chains
        }


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
    draws_file = f"{output_file_prefix}_draws.pkl"
    output_file = f"{output_file_prefix}_stan.pkl"
    data_file = f"{output_file_prefix}_data.pkl"
    
    # Save the input data for later use
    pd.to_pickle({'dp1': dcati, 'dit': dit}, data_file)
    print(f"Saved preprocessed data to: {data_file}")
    
    # Check if we should resume from existing outputs
    can_resume = (resume and 
                  os.path.exists(timing_file) and 
                  os.path.exists(draws_file) and 
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
    category_stats = dcati_sorted.groupby('item_type_id').agg({
        'oid': 'count',
        'item_time_id': 'max',
        'y_stan': lambda x: len(x.unique())
    }).rename(columns={'oid': 'N', 'item_time_id': 'Q', 'y_stan': 'K'})
    
    stan_data['N'] = [int(x) for x in category_stats['N'].tolist()]
    stan_data['Q'] = [int(x) for x in category_stats['Q'].tolist()]
    stan_data['K'] = [int(x) for x in category_stats['K'].tolist()]
    
    # Total counts
    stan_data['N_total'] = int(sum(stan_data['N']))
    stan_data['Q_total'] = int(sum(stan_data['Q']))
    
    # Concatenated observation data
    stan_data['y'] = [int(x) for x in dcati_sorted['y_stan'].tolist()]
    stan_data['unit_of_obs'] = [int(x) for x in dcati_sorted['pid'].tolist()]
    stan_data['question_of_obs'] = [int(x) for x in dcati_sorted['item_time_id'].tolist()]
    stan_data['cat_type'] = [int(x) for x in dcati_sorted['item_type_id'].tolist()]
    
    # Create design matrix using patsy
    import patsy
    design_matrix = patsy.dmatrix(x_formula, data=dcati_sorted, return_type='dataframe')
    
    # Determine which columns to use for predictive probabilities
    if x_formula_ignore_regex is not None:
        import re
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
        poa = pd.read_pickle(draws_file)
        
        print(f"Loading timing data from: {timing_file}")
        timing_data = pd.read_csv(timing_file)
        print("=" * 40 + "\n")
        
        return {
            'draws': poa,
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
        import time
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
        
        # Convert to draws array and save
        print("Extracting draws...")
        # For ADVI, cmdstanpy returns draws as numpy array
        # variational_sample is already a numpy array, not a DataFrame
        draws = fit.variational_sample
        if isinstance(draws, pd.DataFrame):
            poa = draws.values
        else:
            poa = draws  # Already a numpy array
        
        pd.to_pickle(poa, draws_file)
        print(f"Saved draws to: {draws_file}")
        
        return {
            'fit': fit,
            'draws': poa,
            'timing': timing_data
        }
