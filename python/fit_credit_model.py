"""
Model fitting functions for the credit (ncats) Bayesian IRT model.

Three inference back-ends are provided as separate top-level functions:

- ``fit_credit_model_ncats_pyrosvi``: SVI with NumPyro.
- ``fit_credit_model_ncats_stanadvi``: ADVI with Stan (cmdstanpy).
- ``fit_credit_model_ncats_stanhmc``: HMC with Stan (cmdstanpy).

The private helper ``_fit_credit_make_stan_data`` (formerly in ``utils``) lives
here because all three back-ends use it and nothing else does.
"""

from typing import Dict, Optional
import os
from pathlib import Path
import time
import sys
import runpy
import warnings
import re
from functools import partial

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

import jax
import jax.numpy as jnp
import numpyro
from numpyro.infer import SVI, Trace_ELBO

from utils import (
    _plot_ppcheck,
    _compute_ordinal_brier_scores,
    _plot_prob_barplots,
    _plot_worst_chain_traces,
    _make_idata_from_advi_fit,
    _get_autoguide_factory,
    _make_idata_from_svi_posterior,
)


def _fit_credit_make_stan_data(
    dit: pd.DataFrame,
    dcati: pd.DataFrame,
    x_formula: str,
    x_formula_ignore_regex: Optional[str] = None,
) -> dict:
    """
    Build stan_data and design-matrix metadata for credit model ncats models.

    Parameters
    ----------
    dit : pd.DataFrame
        Item metadata table with item types and labels.
    dcati : pd.DataFrame
        Pre-processed data frame used for fitting.
    x_formula : str
        Patsy formula string specifying predictors for the design matrix.
    x_formula_ignore_regex : str | None, optional
        Regex pattern for columns to exclude from predictive probabilities.

    Returns
    -------
    dict
        Stan input dictionary for credit model ncats models.
    """
    # Ensure dcati is sorted by item_type_id and oidt.
    dcati_sorted = dcati.sort_values(['item_type_id', 'oidt']).reset_index(drop=True)

    # Verify oid is sequential.
    assert (dcati_sorted['oid'] == dcati_sorted.index + 1).all(), \
        "oid must be sequential 1:N after sorting by item_type_id, oidt"

    stan_data = {}
    stan_data['C'] = int(dcati_sorted['item_type_id'].max())
    stan_data['U'] = int(dcati_sorted['pid'].max())

    category_stats = dcati.groupby('item_type_id').agg({
        'oid': 'count',
        'item_time_id': 'max',
        'y_stan': lambda x: len(x.unique()),
    }).rename(columns={'oid': 'N', 'item_time_id': 'Q', 'y_stan': 'K'})

    stan_data['N'] = category_stats['N'].astype(int).tolist()
    stan_data['Q'] = category_stats['Q'].astype(int).tolist()
    stan_data['K'] = category_stats['K'].astype(int).tolist()

    stan_data['N_total'] = int(sum(stan_data['N']))
    stan_data['Q_total'] = int(sum(stan_data['Q']))

    stan_data['y'] = dcati['y_stan'].astype(int).tolist()
    stan_data['unit_of_obs'] = dcati['pid'].astype(int).tolist()
    stan_data['question_of_obs'] = dcati['item_time_id'].astype(int).tolist()
    stan_data['cat_type'] = dcati['item_type_id'].astype(int).tolist()

    design_matrix = patsy.dmatrix(x_formula, data=dcati, return_type='dataframe')

    if x_formula_ignore_regex is not None:
        xpr_id = [
            i + 1
            for i, col in enumerate(design_matrix.columns)
            if not re.search(x_formula_ignore_regex, col)
        ]
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

    return stan_data


def fit_credit_model_ncats_pyrosvi(
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
    Run credit model analysis (ncats version) using SVI with NumPyro.

    Mirrors ``fit_partial_credit_model_ncats_pyrosvi`` but loads the credit
    numpyro model and uses ``_fit_credit_make_stan_data``.
    """

    print("\n" + "=" * 40)
    print("Credit Model (ncats) SVI Analysis Configuration")
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
    # zarr-serialised draws and the timing CSV. posterior_samples is rebuilt
    # from the zarr below; auxiliary plots gated on `with_*_analyses` are
    # recomputed if their PDFs are missing, but never block resume.
    can_resume = (
        resume
        and os.path.exists(timing_file)
        and os.path.exists(draws_file)
    )

    print("Preparing Stan data in ncats format...")
    stan_data = _fit_credit_make_stan_data(
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
        # Reconstruct (n_draws, *shape) posterior_samples dict from idata so the
        # return contract matches a fresh fit.
        posterior_samples = {}
        for var_name in idata.posterior.data_vars:
            arr = np.asarray(idata.posterior[var_name])
            posterior_samples[var_name] = arr.reshape((-1,) + arr.shape[2:])
        print(f"Loading timing data from: {timing_file}")
        timing_data = pd.read_csv(timing_file)
        print("=" * 40 + "\n")
    else:
        if resume:
            print("\nNote: Resume requested but not all output files exist. Running full analysis.\n")

        print("Loading NumPyro model...")
        numpyro_model_file = Path(__file__).resolve().parents[1] / "src" / "numpyro" / "credit_model_ncats_v260413.pyro"
        if not numpyro_model_file.exists():
            raise FileNotFoundError(f"NumPyro model file not found: {numpyro_model_file}")

        model_ns = runpy.run_path(str(numpyro_model_file))
        if "credit_model_ncats" not in model_ns:
            raise AttributeError(
                f"Function 'credit_model_ncats' not found in {numpyro_model_file}"
            )
        credit_model_ncats = model_ns["credit_model_ncats"]

        model_for_svi = partial(
            credit_model_ncats,
            sample_ypred=False,
            compute_generated=False,
        )

        guide_factory = _get_autoguide_factory(algorithm)
        rng_key = jax.random.PRNGKey(seed)

        print("Running SVI with Adam optimizer...")
        print(f"  Algorithm: {algorithm}")

        guide = guide_factory(model_for_svi)
        svi = SVI(
            model_for_svi,
            guide,
            numpyro.optim.Adam(lr),
            Trace_ELBO(),
        )

        print("Initializing SVI state...")
        t_init0 = time.time()
        rng_key, subkey = jax.random.split(rng_key)
        init_state = svi.init(subkey, stan_data)
        init_minutes = (time.time() - t_init0) / 60.0
        print(f"  SVI init completed in {init_minutes:.2f} minutes")

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

        print("Generating posterior predictive samples...")
        from numpyro.infer import Predictive

        predictive = Predictive(
            credit_model_ncats,
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


def fit_credit_model_ncats_stanadvi(
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
    Run credit model analysis (ncats version) using ADVI.

    Mirrors ``fit_partial_credit_model_ncats_stanadvi`` with the credit
    Stan model and ``_fit_credit_make_stan_data``.
    """
    if stan_file is None:
        stan_file = str(Path(__file__).parents[2] / "src" / "stan" / "credit_model_ncats_v260413.stan")

    print("\n" + "=" * 40)
    print("Credit Model (ncats) ADVI Analysis Configuration")
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

    os.makedirs(os.path.dirname(output_file_prefix), exist_ok=True)

    timing_file = f"{output_file_prefix}_timing.csv"
    draws_file = f"{output_file_prefix}_draws.zarr"
    output_file = f"{output_file_prefix}_stan.pkl"
    data_file = f"{output_file_prefix}_data.csv"

    dcati.to_csv(data_file.replace('.csv', '_dp1.csv'), index=False)
    dit.to_csv(data_file.replace('.csv', '_dit.csv'), index=False)
    print(f"Saved preprocessed data to: {data_file.replace('.csv', '_dp1.csv')} and _dit.csv")

    can_resume = (
        resume
        and os.path.exists(timing_file)
        and os.path.exists(draws_file)
        and os.path.exists(data_file.replace('.csv', '_dp1.csv'))
        and os.path.exists(data_file.replace('.csv', '_dit.csv'))
        and os.path.exists(output_file)
    )

    print("Preparing Stan data in ncats format...")
    stan_data = _fit_credit_make_stan_data(
        dit=dit,
        dcati=dcati,
        x_formula=x_formula,
        x_formula_ignore_regex=x_formula_ignore_regex,
    )

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

        print("Compiling Stan model...")
        model = CmdStanModel(
            stan_file=stan_file,
            cpp_options={'STAN_THREADS': 'TRUE'},
        )

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
            show_console=show_messages,
        )

        elapsed_time = time.time() - start_time

        print("\nExtracting timing information...")
        timing_data = pd.DataFrame({
            'chain_id': [1],
            'n_warmup': [0],
            'n_sample': [output_samples],
            'mins_init': [elapsed_time / 60],
            'mins_sample': [0.0],
            'mins_generate_samples': [0.0],
            'mins_total': [elapsed_time / 60],
        })
        timing_cols = ['mins_init', 'mins_sample', 'mins_generate_samples', 'mins_total']
        timing_data[timing_cols] = timing_data[timing_cols].round(3)
        timing_data.to_csv(timing_file, index=False)
        print(f"Saved timing information to: {timing_file}")

        print(f"Saving model fit to: {output_file}")
        pd.to_pickle(fit, output_file)

        print("Converting to ArviZ format...")
        idata = _make_idata_from_advi_fit(fit)

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
        'fit': fit,
        'draws': idata,
        'timing': timing_data,
    }


def fit_credit_model_ncats_stanhmc(
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
    Run credit model analysis (ncats version) on pre-processed data using HMC.

    This function performs Bayesian IRT analysis using the flexible ncats credit model
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
        Path to Stan model file (.stan). If None, uses default credit_model_ncats_v260413.stan
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
    if stan_file is None:
        stan_file = str(Path(__file__).parents[2] / "src" / "stan" / "credit_model_ncats_v260413.stan")

    print("\n" + "=" * 40)
    print("Credit Model (ncats) Analysis Configuration")
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

    os.makedirs(os.path.dirname(output_file_prefix), exist_ok=True)

    timing_file = f"{output_file_prefix}_timing.csv"
    draws_file = f"{output_file_prefix}_draws.zarr"
    output_file = f"{output_file_prefix}_stan.pkl"
    mixing_file = f"{output_file_prefix}_convergence_mixing.csv"
    data_file = f"{output_file_prefix}_data.csv"

    dcati.to_csv(data_file.replace('.csv', '_dp1.csv'), index=False)
    dit.to_csv(data_file.replace('.csv', '_dit.csv'), index=False)
    print(f"Saved preprocessed data to: {data_file.replace('.csv', '_dp1.csv')} and _dit.csv")

    can_resume = (resume and
                  os.path.exists(timing_file) and
                  os.path.exists(draws_file) and
                  os.path.exists(output_file) and
                  os.path.exists(mixing_file) and
                  os.path.exists(data_file.replace('.csv', '_dp1.csv')) and
                  os.path.exists(data_file.replace('.csv', '_dit.csv')))

    print("Preparing Stan data in ncats format...")
    stan_data = _fit_credit_make_stan_data(
        dit=dit,
        dcati=dcati,
        x_formula=x_formula,
        x_formula_ignore_regex=x_formula_ignore_regex,
    )

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
        chain_col = 'chain_id' if 'chain_id' in timing_data.columns else 'chain'
        good_chains = timing_data[chain_col].astype(int).tolist()
        print(f"Identified good HMC chains: {', '.join(map(str, good_chains))}")
        print("=" * 40 + "\n")

    else:
        if resume:
            print("\nNote: Resume requested but not all output files exist. Running full analysis.\n")

        print("Compiling Stan model...")
        model = CmdStanModel(
            stan_file=stan_file,
            cpp_options={'STAN_THREADS': 'TRUE'}
        )

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

        print("Checking for divergent transitions and removing non-converged chains...")
        try:
            lp_draws = fit.draws_pd(inc_warmup=False)
        except TypeError:
            lp_draws = fit.draws_pd()

        if 'lp__' not in lp_draws.columns:
            raise ValueError(f"lp__ not found in draws. Available columns: {list(lp_draws.columns)[:10]}")

        chain_col = None
        for col_name in ['chain__', '.chain', 'chain']:
            if col_name in lp_draws.columns:
                chain_col = col_name
                break

        if chain_col is None:
            num_draws_per_chain = fit.num_draws_sampling
            total_draws = len(lp_draws)
            num_chains = len(fit.chain_ids)
            lp_draws['chain'] = [fit.chain_ids[i // num_draws_per_chain] for i in range(total_draws)]
            chain_col = 'chain'

        chain_stats = lp_draws.groupby(chain_col)['lp__'].agg(['mean', 'std', 'count']).reset_index()
        chain_stats.columns = ['chain', 'mean_lp', 'sd_lp', 'n']

        best_idx = chain_stats['mean_lp'].idxmax()
        threshold = chain_stats.loc[best_idx, 'mean_lp'] - 2 * chain_stats.loc[best_idx, 'sd_lp']
        good_chains = chain_stats[chain_stats['mean_lp'] > threshold]['chain'].astype(int).tolist()
        print(f"Identified good HMC chains: {', '.join(map(str, good_chains))}")

        print("\nExtracting timing information...")
        chain_times = fit.time
        timing_data = pd.DataFrame({
            'chain_id': good_chains,
            'n_warmup': [iter_warmup for _ in good_chains],
            'n_sample': [iter_sampling for _ in good_chains],
            'mins_init': [chain_times[int(chain_id) - 1]['warmup'] / 60.0 for chain_id in good_chains],
            'mins_sample': [chain_times[int(chain_id) - 1]['sampling'] / 60.0 for chain_id in good_chains],
            'mins_generate_samples': [0.0 for _ in good_chains],
            'mins_total': [chain_times[int(chain_id) - 1]['total'] / 60.0 for chain_id in good_chains],
        })
        timing_cols = ['mins_init', 'mins_sample', 'mins_generate_samples', 'mins_total']
        timing_data[timing_cols] = timing_data[timing_cols].round(3)

        timing_data.to_csv(timing_file, index=False)
        print(f"Saved timing information to: {timing_file}")

        print(f"Saving model fit to: {output_file}")
        pd.to_pickle(fit, output_file)

        print("Converting to ArviZ format...")
        idata = az.from_cmdstanpy(
            fit,
            coords={'chain': good_chains},
            dims={}
        )

        print(f"Saving draws to: {draws_file}")
        idata.to_zarr(draws_file)

        print("Generating convergence diagnostics...")
        key_vars = ['latent_factor_unit', 'latent_factor_beta',
                   'skill_thresholds', 'loadings_questions_m1']
        summary_df = az.summary(idata, var_names=key_vars)
        summary_df = summary_df.sort_values('ess_bulk')

        print(f"Saved convergence diagnostics to: {mixing_file}")
        summary_df.to_csv(mixing_file)

        if with_additional_analyses and not os.path.exists(f"{output_file_prefix}_worsttrace.pdf"):
            print("Generating worst-chain trace plot...")

            worst_vars = summary_df.head(9).index.tolist()
            worst_vars = [re.sub(r"\[(\d+)\]", lambda match: f"[{int(match.group(1)) + 1}]", v)
                          for v in worst_vars]
            po = fit.draws_pd(inc_warmup=True)
            assert all(v in po.columns for v in worst_vars), \
                f"Missing worst_vars in draws: {[v for v in worst_vars if v not in po.columns]}"

            _plot_worst_chain_traces(
                po=po[['lp__', 'chain__', 'iter__', 'draw__'] + worst_vars],
                iter_warmup=iter_warmup,
                output_file_stem=f"{output_file_prefix}_worsttrace",
            )

    if with_additional_analyses and not os.path.exists(f"{output_file_prefix}_intervals.pdf"):
        print("\nRunning additional diagnostic analyses...")

        print("Generating parameter plots...")
        key_vars_pattern = '|'.join(['latent_factor_unit', 'latent_factor_beta',
                                    'skill_thresholds', 'loadings_questions_m1'])

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

            for ax in np.ravel(np.atleast_1d(p)):
                ax.tick_params(axis='y', labelsize=6)

            plt.subplots_adjust(left=0.42, right=0.98, top=0.98, bottom=0.02)
            plt.savefig(f"{output_file_prefix}_intervals.pdf", bbox_inches='tight')
            plt.close()

    if with_additional_analyses and not os.path.exists(f"{output_file_prefix}_ppcheck.png"):
        print("Generating posterior predictive checks...")

        if 'ypred' in idata.posterior.data_vars:
            _plot_ppcheck(idata.posterior['ypred'].values, dcati, f"{output_file_prefix}_ppcheck")

    if with_additional_analyses and not os.path.exists(f"{output_file_prefix}_ordered_brierscore.csv"):
        print("Computing ordinal Brier scores by question...")

        if 'ordinal_brier_score' in idata.posterior.data_vars:
            _compute_ordinal_brier_scores(idata.posterior['ordinal_brier_score'].values, dcati, f"{output_file_prefix}_ordered_brierscore")

    if with_core_analyses and not os.path.exists(f"{output_file_prefix}_prob_by_question_fit.pdf"):
        print("\nGenerating fitted probability plots...")
        if 'ordered_prob_by_cat_qu_fit' in idata.posterior.data_vars:
            _plot_prob_barplots(
                idata.posterior['ordered_prob_by_cat_qu_fit'].values,
                dcati,
                dit,
                f"{output_file_prefix}_prob_by_question_fit",
            )

    if with_core_analyses and x_formula_ignore_regex is not None and not os.path.exists(f"{output_file_prefix}_prob_by_question_pr.pdf"):
        print("\nGenerating predictive probability plots with columns matching x_formula_ignore_regex removed from X...")
        if 'ordered_prob_by_cat_qu_pr' in idata.posterior.data_vars:
            _plot_prob_barplots(
                idata.posterior['ordered_prob_by_cat_qu_pr'].values,
                dcati,
                dit,
                f"{output_file_prefix}_prob_by_question_pr",
            )

    return {
        'fit': fit,
        'draws': idata,
        'timing': timing_data,
        'good_chains': good_chains
    }
