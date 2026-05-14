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
import ggsci
import matplotlib.colors as mcolors
import arviz as az
import matplotlib.pyplot as plt
from cmdstanpy import CmdStanModel
from plotnine import *
from plotnine import options as p9_options

from utils import _plot_ppcheck, _compute_ordinal_brier_scores, _plot_prob_barplots, _fit_ordered_logit_make_stan_data, _make_idata_from_advi_fit

def fit_ordered_logit_model_ncats_stanadvi(
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
            'chain_id': [1],
            'n_warmup': [0],
            'n_sample': [output_samples],
            # CmdStan variational API does not expose per-step timings; store wall-clock as init.
            'mins_init': [elapsed_time / 60],
            'mins_sample': [0.0],
            'mins_generate_samples': [0.0],
            'mins_total': [elapsed_time / 60],
        })
        timing_cols = ['mins_init', 'mins_sample', 'mins_generate_samples', 'mins_total']
        timing_data[timing_cols] = timing_data[timing_cols].round(3)
        
        timing_data.to_csv(timing_file, index=False)
        print(f"Saved timing information to: {timing_file}")
        
        # Save model fit
        print(f"Saving model fit to: {output_file}")
        pd.to_pickle(fit, output_file)
        
        # Convert to ArviZ InferenceData and save
        print("Converting to ArviZ format...")
        idata = _make_idata_from_advi_fit(fit)
        
        print(f"Saving draws to: {draws_file}")
        idata.to_zarr(draws_file)

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
    if with_additional_analyses and not os.path.exists(f"{output_file_prefix}_ppcheck.pdf"):
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
        'timing': timing_data
    }
