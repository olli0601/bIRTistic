"""
Model fitting functions for Bayesian IRT analysis.

This module provides functions for fitting ordered logit models
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
import ggsci
import matplotlib.colors as mcolors
import arviz as az
import matplotlib.pyplot as plt
from cmdstanpy import CmdStanModel
from plotnine import *
from plotnine import options as p9_options

from utils import _plot_ppcheck, _compute_ordinal_brier_scores, _plot_prob_barplots, _plot_worst_chain_traces, _fit_ordered_logit_make_stan_data

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
    #%%
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
    draws_file = f"{output_file_prefix}_draws.zarr"  # Use zarr for arviz (9x faster than NetCDF)
    output_file = f"{output_file_prefix}_stan.pkl"
    mixing_file = f"{output_file_prefix}_convergence_mixing.csv"
    data_file = f"{output_file_prefix}_data.csv"
    
    # Save the input data for later use
    # Save dp1 and dit metadata separately
    dcati.to_csv(data_file.replace('.csv', '_dp1.csv'), index=False)
    dit.to_csv(data_file.replace('.csv', '_dit.csv'), index=False)
    print(f"Saved preprocessed data to: {data_file.replace('.csv', '_dp1.csv')} and _dit.csv")
    
    # Check if we should resume from existing outputs
    can_resume = (resume and 
                  os.path.exists(timing_file) and 
                  os.path.exists(draws_file) and 
                  os.path.exists(output_file) and 
                  os.path.exists(mixing_file) and
                  os.path.exists(data_file.replace('.csv', '_dp1.csv')) and
                  os.path.exists(data_file.replace('.csv', '_dit.csv')))
    
    # Prepare data in ncats format (same as fit_partial_credit_model_ncats)
    print("Preparing Stan data in ncats format...")
    stan_data = _fit_ordered_logit_make_stan_data(
        dit=dit,
        dcati=dcati,
        x_formula=x_formula,
        x_formula_ignore_regex=x_formula_ignore_regex,
    )
    
    # Check if resuming from existing outputs
    if can_resume:
        print("\n" + "=" * 40)
        print("RESUMING from existing outputs")
        print("=" * 40)
        
        print(f"Loading fit from: {output_file}")        
        fit = pd.read_pickle(output_file)

        print(f"Loading draws from: {draws_file}")        
        idata = az.from_zarr(draws_file)
        
        print(f"Loading timing data from: {timing_file}")
        timing_data = pd.read_csv(timing_file)
        good_chains = timing_data['chain'].tolist()
        print(f"Identified good HMC chains: {', '.join(map(str, good_chains))}")
        print("=" * 40 + "\n")
        
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
        key_vars = ['latent_factor_unit', 'latent_factor_beta', 
               'skill_thresholds_1', 'skill_thresholds_incs', 'loadings_questions_m1']
        summary_df = az.summary(idata, var_names=key_vars)
        summary_df = summary_df.sort_values('ess_bulk')        
                
        print(f"Saved convergence diagnostics to: {mixing_file}")
        summary_df.to_csv(mixing_file)

        if with_additional_analyses and not os.path.exists(f"{output_file_prefix}_worsttrace.pdf"):
            print("Generating worst-chain trace plot...")

            worst_vars = summary_df.head(9).index.tolist()
            worst_vars = [sub(r"\[(\d+)\]", lambda match: f"[{int(match.group(1)) + 1}]", v)
                          for v in worst_vars]
            po = fit.draws_pd(inc_warmup=True)            
            assert all(v in po.columns for v in worst_vars), \
                f"Missing worst_vars in draws: {[v for v in worst_vars if v not in po.columns]}"        
            
            _plot_worst_chain_traces(
                po=po[['lp__', 'chain__', 'iter__', 'draw__'] + worst_vars],
                iter_warmup=iter_warmup,
                output_file_stem=f"{output_file_prefix}_worsttrace",
            )

    #%%    
    # Additional diagnostic analyses
    if with_additional_analyses and not os.path.exists(f"{output_file_prefix}_intervals.pdf"):
        print("\nRunning additional diagnostic analyses...")
        
        # 1. Generate parameter intervals plot
        print("Generating parameter plots...")
        key_vars_pattern = '|'.join(['latent_factor_unit', 'latent_factor_beta',
                                    'skill_thresholds_1', 'skill_thresholds_incs', 
                                    'loadings_questions_m1'])
        
        # Use ArviZ plot_forest for parameter intervals
        var_names = [v for v in idata.posterior.data_vars 
                    if any(pat in v for pat in key_vars_pattern.split('|'))]
        
        if var_names:
            p = az.plot_forest(
                idata,
                var_names=var_names,
                combined=True,
                hdi_prob=0.95,
                textsize=8,
                figsize=(18, min(160, max(30, len(var_names) * 1.1)))
            )

            # Increase label area and reduce y-label size for long parameter names.
            for ax in np.ravel(np.atleast_1d(p)):
                ax.tick_params(axis='y', labelsize=6)

            plt.subplots_adjust(left=0.42, right=0.98, top=0.98, bottom=0.02)
            plt.savefig(f"{output_file_prefix}_intervals.pdf", bbox_inches='tight')
            plt.close()
    #%%    
    # 2. Posterior predictive checks
    if with_additional_analyses and not os.path.exists(f"{output_file_prefix}_ppcheck.png"):
        print("Generating posterior predictive checks...")
        
        if 'ypred' in idata.posterior.data_vars:
            _plot_ppcheck(idata.posterior['ypred'].values, dcati, f"{output_file_prefix}_ppcheck")
    #%%    
    # 3. Compute ordinal Brier scores
    if with_additional_analyses and not os.path.exists(f"{output_file_prefix}_ordered_brierscore.csv"):
        print("Computing ordinal Brier scores by question...")
        
        if 'ordinal_brier_score' in idata.posterior.data_vars:
            _compute_ordinal_brier_scores(idata.posterior['ordinal_brier_score'].values, dcati, f"{output_file_prefix}_ordered_brierscore")
    #%%    
    # Core probability analyses
    if with_core_analyses and not os.path.exists(f"{output_file_prefix}_prob_by_question_fit.pdf"):
        print("\nGenerating fitted probability plots...")        
        if 'ordered_prob_by_cat_qu_fit' in idata.posterior.data_vars:
            _plot_prob_barplots(
                idata.posterior['ordered_prob_by_cat_qu_fit'].values,
                dcati,
                dit,
                f"{output_file_prefix}_prob_by_question_fit",
            )
                
    #%%
    if with_core_analyses and x_formula_ignore_regex is not None and not os.path.exists(f"{output_file_prefix}_prob_by_question_pr.pdf"):
        print("\nGenerating predictive probability plots with columns matching x_formula_ignore_regex removed from X...")        
        if 'ordered_prob_by_cat_qu_pr' in idata.posterior.data_vars:
            _plot_prob_barplots(
                idata.posterior['ordered_prob_by_cat_qu_pr'].values,
                dcati,
                dit,
                f"{output_file_prefix}_prob_by_question_pr",
            )
        
    #%% 
    return {
        'fit': fit,
        'draws': idata,
        'timing': timing_data,
        'good_chains': good_chains
    }
