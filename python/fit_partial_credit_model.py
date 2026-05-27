"""
Model fitting functions for Bayesian IRT analysis.

This module provides functions for fitting partial credit models
using Stochastic Variational Inference (SVI) with NumPyro.
"""

from typing import Dict, Optional, Callable
import os
from pathlib import Path
import time
import sys
import runpy
from functools import partial

import pandas as pd
import numpy as np
import arviz as az
import matplotlib.pyplot as plt

import jax
import jax.numpy as jnp
import numpyro
from numpyro.infer import SVI, Trace_ELBO
from numpyro.infer.autoguide import (
    AutoLaplaceApproximation,
    AutoMultivariateNormal,
    AutoLowRankMultivariateNormal,
    AutoDiagonalNormal,
    AutoIAFNormal,
)

from utils import (
    _plot_ppcheck,
    _compute_ordinal_brier_scores,
    _plot_prob_barplots,
    _fit_partial_credit_make_stan_data,
)


def _get_autoguide_factory(algorithm: str) -> Callable:
    """
    Get the autoguide factory for the specified algorithm.

    Parameters
    ----------
    algorithm : str
        One of 'AutoLaplaceApproximation', 'AutoMultivariateNormal',
        'AutoDiagonalNormal', 'AutoIAFNormal'.

    Returns
    -------
    Callable
        Autoguide factory class.

    Raises
    ------
    ValueError
        If algorithm is not recognized.
    """
    factories = {
        'AutoLaplaceApproximation': AutoLaplaceApproximation,
        'AutoMultivariateNormal': AutoMultivariateNormal,
        'AutoLowRankMultivariateNormal': AutoLowRankMultivariateNormal,
        'AutoDiagonalNormal': AutoDiagonalNormal,
        'AutoIAFNormal': AutoIAFNormal,
    }
    if algorithm not in factories:
        raise ValueError(
            f"Unknown algorithm '{algorithm}'. Must be one of: {', '.join(factories.keys())}"
        )
    return factories[algorithm]


def _make_idata_from_svi_posterior(stan_data: Dict, posterior_samples: Dict, predictions: Dict) -> object:
    """
    Build an ArviZ InferenceData object from SVI posterior samples.

    Parameters
    ----------
    stan_data : dict
        Stan data dict used for the model.
    posterior_samples : dict
        Dictionary of posterior samples from numpyro model, with shape (n_draws, *param_shape).
    predictions : dict
        Dictionary of generated quantity samples, with shape (n_draws, *var_shape).

    Returns
    -------
    arviz.InferenceData
        InferenceData containing posterior draws in standard Stan format.
    """
    # Keep a reference for future schema checks (expected by caller contract).
    _ = stan_data

    # Reshape all samples to (1, n_draws, *shape) for chain/draw dims
    posterior_dict = {}
    dims = {}

    for var_name, samples in posterior_samples.items():
        if var_name.startswith('_'):
            continue

        samples_arr = np.asarray(samples)
        if samples_arr.ndim == 1:
            # Scalar parameter
            posterior_dict[var_name] = samples_arr.reshape(1, -1)
        else:
            # Multi-dimensional parameter: (n_draws, *shape) -> (1, n_draws, *shape)
            posterior_dict[var_name] = samples_arr[np.newaxis, ...]
            ndim = samples_arr.ndim - 1
            dims[var_name] = [f"{var_name}_dim_{i}" for i in range(ndim)]

    # Include generated quantities in posterior to match existing ADVI/HMC zarr layout.
    for var_name, samples in predictions.items():
        samples_arr = np.asarray(samples)
        if samples_arr.ndim == 1:
            posterior_dict[var_name] = samples_arr.reshape(1, -1)
        else:
            posterior_dict[var_name] = samples_arr[np.newaxis, ...]
            ndim = samples_arr.ndim - 1
            dims[var_name] = [f"{var_name}_dim_{i}" for i in range(ndim)]

    return az.from_dict(posterior=posterior_dict, dims=dims)


def fit_partial_credit_model_ncats_pyrosvi(
    dit: pd.DataFrame,
    dcati: pd.DataFrame,
    output_file_prefix: str,
    x_formula: str = "~ time - 1",
    x_formula_ignore_regex: Optional[str] = None,
    algorithm: str = 'AutoDiagonalNormal',
    lr: float = 0.01,
    num_steps: int = 10000,
    output_samples: int = 4000,
    seed: int = 123,
    resume: bool = False,
    with_core_analyses: bool = True,
    with_additional_analyses: bool = False,
) -> Dict:
    """
    Run partial credit model analysis (ncats version) using SVI with NumPyro.

    This function performs Bayesian IRT analysis using the flexible ncats partial credit model
    with Stochastic Variational Inference (SVI) as the approximate inference backend.

    Parameters
    ----------
    dit : pd.DataFrame
        Item metadata table with item types and labels.
    dcati : pd.DataFrame
        Pre-processed data with observations for analysis. Must include
        item_type, y_stan, item_time_id, pid, and oidt columns.
    output_file_prefix : str
        Full path prefix for output files (without extension).
    x_formula : str, default "~ time - 1"
        Patsy formula string specifying predictors for the design matrix.
    x_formula_ignore_regex : str, optional
        Regular expression pattern to identify X columns to exclude from posterior predictive.
    algorithm : str, default 'AutoDiagonalNormal'
        Autoguide algorithm to use. One of:
        - 'AutoLaplaceApproximation': Laplace approximation
        - 'AutoMultivariateNormal': Full covariance Gaussian
        - 'AutoDiagonalNormal': Diagonal covariance Gaussian
        - 'AutoIAFNormal': Inverse autoregressive flow with Gaussian base
    lr : float, default 0.01
        Learning rate for the Adam optimizer.
    num_steps : int, default 10000
        Number of SVI optimization steps.
    output_samples : int, default 4000
        Number of samples to draw from the learned approximate posterior.
    seed : int, default 123
        Random seed for reproducibility.
    resume : bool, default False
        If True and all output files exist, skip SVI and load existing results.
    with_core_analyses : bool, default True
        If True, generate core probability plots.
    with_additional_analyses : bool, default False
        If True, generate additional diagnostic plots.

    Returns
    -------
    dict
        Dictionary containing:
        - 'posterior_samples': Dictionary of posterior samples from SVI
        - 'draws': InferenceData object with posterior draws
        - 'timing': DataFrame with timing information
        - 'algorithm': str, the algorithm name used
    """

    print("\n" + "=" * 40)
    print("Partial Credit Model (ncats) SVI Analysis Configuration")
    print("=" * 40)
    print(f"Data: dit with {len(dit)} items, dcati with {len(dcati)} observations")
    print(f"Output prefix: {output_file_prefix}")
    print(f"Algorithm: {algorithm}")
    print(f"Learning rate: {lr}")
    print(f"SVI steps: {num_steps}")
    print(f"Output samples: {output_samples}")
    print(f"Seed: {seed}")
    print(f"Resume: {resume}")
    print(f"Core analyses: {with_core_analyses}")
    print(f"Additional analyses: {with_additional_analyses}")
    print("=" * 40 + "\n")

    os.makedirs(os.path.dirname(output_file_prefix), exist_ok=True)

    timing_file = f"{output_file_prefix}_timing.csv"
    draws_file = f"{output_file_prefix}_draws.zarr"
    data_file = f"{output_file_prefix}_data.csv"

    dcati.to_csv(data_file.replace('.csv', '_dp1.csv'), index=False)
    dit.to_csv(data_file.replace('.csv', '_dit.csv'), index=False)
    print(f"Saved preprocessed data to: {data_file.replace('.csv', '_dp1.csv')} and _dit.csv")

    # Resume only needs the two artifacts that capture the fit itself: the
    # zarr-serialised draws (= posterior + generated quantities) and the timing
    # CSV. 
    can_resume = (
        resume
        and os.path.exists(timing_file)
        and os.path.exists(draws_file)
    )

    print("Preparing Stan data in ncats format...")
    stan_data = _fit_partial_credit_make_stan_data(
        dit=dit,
        dcati=dcati,
        x_formula=x_formula,
        x_formula_ignore_regex=x_formula_ignore_regex,
    )

    if can_resume:
        print("\n" + "=" * 40)
        print("RESUMING from existing outputs")
        print("=" * 40)

        print(f"Loading draws from: {draws_file}")
        idata = az.from_zarr(draws_file)

        posterior_samples = {}
        for var_name in idata.posterior.data_vars:
            arr = np.asarray(idata.posterior[var_name])  # (chain, draw, *shape)
            posterior_samples[var_name] = arr.reshape((-1,) + arr.shape[2:])

        print(f"Loading timing data from: {timing_file}")
        timing_data = pd.read_csv(timing_file)
        print("=" * 40 + "\n")
    else:
        if resume:
            print("\nNote: Resume requested but not all output files exist. Running full analysis.\n")

        # Import model from .pyro source file.
        print("Loading NumPyro model...")
        numpyro_model_file = Path(__file__).resolve().parents[1] / "src" / "numpyro" / "partial_credit_model_ncats_v260413.pyro"
        if not numpyro_model_file.exists():
            raise FileNotFoundError(f"NumPyro model file not found: {numpyro_model_file}")

        model_ns = runpy.run_path(str(numpyro_model_file))
        if "partial_credit_model_ncats" not in model_ns:
            raise AttributeError(
                f"Function 'partial_credit_model_ncats' not found in {numpyro_model_file}"
            )
        partial_credit_model_ncats = model_ns["partial_credit_model_ncats"]

        # Use a lean variant during SVI optimization:
        # - sample_ypred=False avoids discrete-latent funsor dependency
        # - compute_generated=False skips expensive generated-quantity tracing
        model_for_svi = partial(
            partial_credit_model_ncats,
            sample_ypred=False,
            compute_generated=False,
        )

        # Get autoguide factory
        guide_factory = _get_autoguide_factory(algorithm)

        # Initialize JAX RNG
        rng_key = jax.random.PRNGKey(seed)

        print("Running SVI with Adam optimizer...")
        print(f"  Algorithm: {algorithm}")

        # Create guide
        guide = guide_factory(model_for_svi)

        # Create SVI object
        svi = SVI(
            model_for_svi,
            guide,
            numpyro.optim.Adam(lr),
            Trace_ELBO(),
        )

        # Initialize SVI state
        print("Initializing SVI state...")
        t_init0 = time.time()
        rng_key, subkey = jax.random.split(rng_key)
        init_state = svi.init(subkey, stan_data)
        init_minutes = (time.time() - t_init0) / 60.0
        print(f"  SVI init completed in {init_minutes:.2f} minutes")

        # Run SVI optimization using the compiled run loop.
        print(f"Optimizing variational parameters for {num_steps} steps...")
        t_opt0 = time.time()
        rng_key, subkey = jax.random.split(rng_key)
        run_result = svi.run(
            subkey,
            num_steps,
            stan_data,
            init_state=init_state,
            progress_bar=False,
        )
        opt_minutes = (time.time() - t_opt0) / 60.0
        svi_state = run_result.state
        losses = np.asarray(run_result.losses)

        for i in range(999, len(losses), 1000):
            print(f"  Step {i + 1:5d} / {num_steps}: ELBO = {-losses[i]:,.1f}", flush=True)
        if len(losses) and (len(losses) % 1000) != 0:
            i = len(losses) - 1
            print(f"  Step {i + 1:5d} / {num_steps}: ELBO = {-losses[i]:,.1f}", flush=True)

        print(f"\nSVI optimization completed in {opt_minutes:.2f} minutes (after {init_minutes:.2f} min init)")

        # Extract posterior samples from learned guide
        print(f"Sampling {output_samples} draws from the learned approximate posterior...")
        t_post0 = time.time()
        rng_key, subkey = jax.random.split(rng_key)
        params = svi.get_params(svi_state)
        if hasattr(svi, "get_posterior"):
            posterior_samples = svi.get_posterior(svi_state).sample(
                subkey,
                sample_shape=(output_samples,),
            )
        else:
            posterior_samples = guide.sample_posterior(
                subkey,
                params,
                stan_data,
                sample_shape=(output_samples,),
            )
        posterior_samples_dict = {k: np.asarray(v) for k, v in posterior_samples.items()}
        post_minutes = (time.time() - t_post0) / 60.0
        print(f"  Posterior sampling completed in {post_minutes:.2f} minutes")

        # Re-run model with posterior samples to get generated quantities
        print("Generating posterior predictive samples...")

        # Use predictive to get generated quantities
        from numpyro.infer import Predictive

        predictive = Predictive(
            partial_credit_model_ncats,
            posterior_samples=posterior_samples_dict,
            return_sites=['log_lik', 'ypred', 'ordered_prob_by_cat_qu_fit',
                         'ordered_prob_by_cat_qu_pr', 'ordinal_brier_score'],
        )

        t_pred0 = time.time()
        rng_key, subkey = jax.random.split(rng_key)
        predictions = predictive(subkey, stan_data)
        pred_minutes = (time.time() - t_pred0) / 60.0
        print(f"  Posterior predictive generation completed in {pred_minutes:.2f} minutes")

        print("Converting to ArviZ format...")
        t_idata0 = time.time()
        idata = _make_idata_from_svi_posterior(stan_data, posterior_samples_dict, predictions)
        idata_minutes = (time.time() - t_idata0) / 60.0
        print(f"  ArviZ conversion completed in {idata_minutes:.2f} minutes")

        print("Extracting timing information...")
        mins_generate_samples = post_minutes + pred_minutes + idata_minutes
        mins_total = init_minutes + opt_minutes + mins_generate_samples
        timing_data = pd.DataFrame({
            'chain_id': [1],
            'n_warmup': [0],
            'n_sample': [output_samples],
            'mins_init': [init_minutes],
            'mins_sample': [opt_minutes],
            'mins_generate_samples': [mins_generate_samples],
            'mins_total': [mins_total],
        })
        timing_cols = ['mins_init', 'mins_sample', 'mins_generate_samples', 'mins_total']
        timing_data[timing_cols] = timing_data[timing_cols].round(3)
        timing_data.to_csv(timing_file, index=False)
        print(f"Saved timing information to: {timing_file}")

        print(f"Saving draws to: {draws_file}")
        idata.to_zarr(draws_file)

    if with_additional_analyses and not os.path.exists(f"{output_file_prefix}_intervals.pdf"):
        print("\nRunning additional diagnostic analyses...")
        print("Generating parameter plots...")

        key_vars_pattern = '|'.join([
            'latent_factor_unit',
            'latent_factor_beta',
            'skill_thresholds',
            'loadings_questions_m1',
        ])

        var_names = [
            v for v in idata.posterior.data_vars
            if any(pat in v for pat in key_vars_pattern.split('|'))
        ]

        if var_names:
            p = az.plot_forest(
                idata,
                var_names=var_names,
                combined=True,
                hdi_prob=0.95,
                textsize=8,
                figsize=(18, min(160, max(30, len(var_names) * 1.1))),
            )

            for ax in np.ravel(np.atleast_1d(p)):
                ax.tick_params(axis='y', labelsize=6)

            plt.subplots_adjust(left=0.42, right=0.98, top=0.98, bottom=0.02)
            plt.savefig(f"{output_file_prefix}_intervals.pdf", bbox_inches='tight')
            plt.close()

    if with_additional_analyses and not os.path.exists(f"{output_file_prefix}_ppcheck.pdf"):
        print("Generating posterior predictive checks...")
        if 'ypred' in idata.posterior.data_vars:
            _plot_ppcheck(idata.posterior['ypred'].values, dcati, f"{output_file_prefix}_ppcheck")

    if with_additional_analyses and not os.path.exists(f"{output_file_prefix}_ordered_brierscore.csv"):
        print("Computing ordinal Brier scores by question...")
        if 'ordinal_brier_score' in idata.posterior.data_vars:
            _compute_ordinal_brier_scores(
                idata.posterior['ordinal_brier_score'].values,
                dcati,
                f"{output_file_prefix}_ordered_brierscore",
            )

    if with_core_analyses and not os.path.exists(f"{output_file_prefix}_prob_by_question_fit.pdf"):
        print("\nGenerating fitted probability plots...")
        if 'ordered_prob_by_cat_qu_fit' in idata.posterior.data_vars:
            _plot_prob_barplots(
                idata.posterior['ordered_prob_by_cat_qu_fit'].values,
                dcati,
                dit,
                f"{output_file_prefix}_prob_by_question_fit",
            )

    if (
        with_core_analyses
        and x_formula_ignore_regex is not None
        and not os.path.exists(f"{output_file_prefix}_prob_by_question_pr.pdf")
    ):
        print("\nGenerating predictive probability plots with columns matching x_formula_ignore_regex removed from X...")
        if 'ordered_prob_by_cat_qu_pr' in idata.posterior.data_vars:
            _plot_prob_barplots(
                idata.posterior['ordered_prob_by_cat_qu_pr'].values,
                dcati,
                dit,
                f"{output_file_prefix}_prob_by_question_pr",
            )

    return {
        'posterior_samples': posterior_samples_dict if not can_resume else posterior_samples,
        'draws': idata,
        'timing': timing_data,
        'algorithm': algorithm,
    }
