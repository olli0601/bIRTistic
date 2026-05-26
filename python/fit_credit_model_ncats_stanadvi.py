"""
Model fitting functions for Bayesian IRT analysis.

This module provides functions for fitting credit models (ncats)
using ADVI Stan via cmdstanpy.
"""

from typing import Dict, Optional
import os
from pathlib import Path
import time

import pandas as pd
import numpy as np
import arviz as az
import matplotlib.pyplot as plt
from cmdstanpy import CmdStanModel

from utils import (
    _plot_ppcheck,
    _compute_ordinal_brier_scores,
    _plot_prob_barplots,
    _fit_credit_make_stan_data,
    _make_idata_from_advi_fit,
)


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
