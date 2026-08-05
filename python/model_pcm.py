"""
:class:`PartialCreditModel` -- the IRT partial-credit model subclass of
:class:`model._IRTModel`. All PCM-specific math, the stan-data builder,
and the three Stan/Pyro fit drivers are now methods on this class
(previously free functions in ``fit_partial_credit_model.py``, deleted
after this consolidation).

Each method body matches the corresponding pre-refactor free function
verbatim; the only change at the method entry point is::

    dit = self.dit
    dcati = self.dcati
    if x_formula is None: x_formula = self.x_formula
    if seed is None: seed = self.seed

so the original bodies (data saving, logging, ELBO loops, ArviZ
conversion, plotting fan-out) compile unchanged.
"""

from typing import Dict, Optional
import os
from pathlib import Path
import time
import runpy
import re
from functools import partial

import pandas as pd
import numpy as np
import patsy
import arviz as az
import matplotlib.pyplot as plt
from cmdstanpy import CmdStanModel

import jax
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
    _pdf_or_parts_exist,
)
from model_irt import IRTModel


_MODEL_NS = None


def _model_namespace():
    """Lazily load + cache the NumPyro model namespace from the .pyro source."""
    global _MODEL_NS
    if _MODEL_NS is None:
        model_file = Path(__file__).resolve().parents[1] / "src" / "numpyro" / "partial_credit_model_ncats_v260413.pyro"
        if not model_file.exists():
            raise FileNotFoundError(f"NumPyro model file not found: {model_file}")
        _MODEL_NS = runpy.run_path(str(model_file))
    return _MODEL_NS


class PartialCreditModel(IRTModel):
    """Partial credit (ncats) IRT model. Implements every :class:`model.Model`
    blueprint method inline; ``_module`` is unset so the :class:`_IRTModel`
    mixin's dict dispatch is bypassed."""

    param_names = (
        'latent_factor_unit',
        'latent_factor_beta',
        'skill_thresholds',
        'skill_thresholds_1',
        'skill_thresholds_incs',
        'loadings_questions_m1',
    )
    positive_params = ('loadings_questions_m1',)

    # ------------------------------------------------------------------
    # Thin proxies to the .pyro namespace (verbatim from the old module).
    # ------------------------------------------------------------------

    @staticmethod
    def eval_loglik(data, params):
        """
        Pointwise log-likelihood ``array[N_total]`` for a single parameter set,
        given a stan_data dict ``data`` (as built by ``make_stan_data``)
        and ``params`` holding ``latent_factor_unit``, ``latent_factor_beta``,
        ``skill_thresholds``, ``loadings_questions_m1``. Delegates to the
        model's ``get_log_likelihood_of_partial_credit_model_ncats``.
        """
        return _model_namespace()['get_log_likelihood_of_partial_credit_model_ncats'](data, params)

    @staticmethod
    def eval_loglik_annealed(data, params):
        """
        Temperature-annealed pointwise log-likelihood
        ``params['temperature'] * eval_loglik(data, params)``. Carrying the
        temperature inside ``params`` (a traced value) lets an SMC move kernel
        compile once and be reused across the whole annealing schedule.
        """
        return _model_namespace()['get_log_likelihood_of_partial_credit_model_ncats_with_annealing'](data, params)

    @staticmethod
    def eval_log_prior(params):
        """Pointwise log p(theta) for the partial credit ncats priors."""
        return _model_namespace()['get_prior_of_partial_credit_model_ncats'](params)

    @staticmethod
    def eval_outcome_for_endpoint(data, params):
        """``ordered_prob_by_cat_qu_fit`` for one parameter set; feeds
        :func:`get_endpoints.get_endpoints_per_draw` so reweighted /
        moved particles can be scored for per-item improvement ratios."""
        return _model_namespace()['get_ordered_prob_of_partial_credit_model_ncats'](data, params)

    # ------------------------------------------------------------------
    # Stan-data builder. Body verbatim; only `dit = self.dit` is prepended.
    # ------------------------------------------------------------------

    def make_stan_data(self, dcati: pd.DataFrame, x_formula: str,
                       x_formula_ignore_regex: Optional[str] = None,
                       verbose: bool = False) -> dict:
        """Build stan_data and design-matrix metadata for partial credit
        ncats models. Reads ``self.dit`` plus the caller-supplied ``dcati``
        (so interim algorithms can pass a z-cohort frame without instantiating
        a fresh model)."""
        dit = self.dit
        vprint = print if verbose else (lambda *args, **kwargs: None)
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

        vprint(f"Design matrix number of predictors: P = {stan_data['P']}")
        vprint(f"Design matrix column names (for fitting): {', '.join(design_matrix.columns)}")
        vprint(f"Design matrix column names (for prediction): {', '.join([design_matrix.columns[i-1] for i in xpr_id])}")
        vprint(f"Number of category types: C = {stan_data['C']}")
        vprint(f"Category type labels: {', '.join(dit['item_type'].unique())}")
        vprint(f"N per category: {', '.join(map(str, stan_data['N']))}")
        vprint(f"Q per category: {', '.join(map(str, stan_data['Q']))}")
        vprint(f"K per category: {', '.join(map(str, stan_data['K']))}")

        return stan_data

    # ------------------------------------------------------------------
    # NumPyro SVI fit driver. Body verbatim; only the entry rebinds dit /
    # dcati / x_formula / seed from self.
    # ------------------------------------------------------------------

    def fit_pyro_svi(
        self,
        output_file_prefix: str,
        x_formula: Optional[str] = None,
        x_formula_ignore_regex: Optional[str] = None,
        algorithm: str = 'AutoDiagonalNormal',
        lr: float = 0.01,
        num_steps: int = 10000,
        output_samples: int = 4000,
        seed: Optional[int] = None,
        resume: bool = False,
        with_core_analyses: bool = True,
        with_additional_analyses: bool = False,
        save_to_file: bool = True,
        verbose: bool = True,
    ) -> Dict:
        """
        Run partial credit model analysis (ncats version) using SVI with NumPyro.

        When ``save_to_file`` is False the fit runs purely in memory: no data/timing
        /draws files are written, resume is disabled, and the core/additional
        analyses (which produce on-disk plots) are skipped. The returned dict still
        holds the in-memory ``draws`` idata, ``posterior_samples`` and ``timing``.
        """
        dit = self.dit
        dcati = self.dcati
        if x_formula is None:
            x_formula = self.x_formula
        if seed is None:
            seed = self.seed

        vprint = print if verbose else (lambda *args, **kwargs: None)

        vprint("\n" + "=" * 40)
        vprint("Partial Credit Model (ncats) SVI Analysis Configuration")
        vprint("=" * 40)
        vprint(f"Data: dit with {len(dit)} items, dcati with {len(dcati)} observations")
        vprint(f"Output prefix: {output_file_prefix}")
        vprint(f"Algorithm: {algorithm}")
        vprint(f"Learning rate: {lr}")
        vprint(f"SVI steps: {num_steps}")
        vprint(f"Output samples: {output_samples}")
        vprint(f"Seed: {seed}")
        vprint(f"Resume: {resume}")
        vprint(f"Core analyses: {with_core_analyses}")
        vprint(f"Additional analyses: {with_additional_analyses}")
        vprint(f"Save to file: {save_to_file}")
        vprint("=" * 40 + "\n")

        if save_to_file:
            os.makedirs(os.path.dirname(output_file_prefix), exist_ok=True)

        timing_file = f"{output_file_prefix}_timing.csv"
        draws_file = f"{output_file_prefix}_draws.zarr"

        # Resume only needs the two artifacts that capture the fit itself: the
        # zarr-serialised draws (= posterior + generated quantities) and the timing
        # CSV. Disabled when not persisting (nothing reliable to reuse).
        can_resume = (
            resume
            and save_to_file
            and os.path.exists(timing_file)
            and os.path.exists(draws_file)
        )

        if can_resume:
            vprint("\n" + "=" * 40)
            vprint("RESUMING from existing outputs")
            vprint("=" * 40)

            vprint(f"Loading draws from: {draws_file}")
            idata = az.from_zarr(draws_file)

            posterior_samples = {}
            for var_name in idata.posterior.data_vars:
                arr = np.asarray(idata.posterior[var_name])  # (chain, draw, *shape)
                posterior_samples[var_name] = arr.reshape((-1,) + arr.shape[2:])

            vprint(f"Loading timing data from: {timing_file}")
            timing_data = pd.read_csv(timing_file)
            vprint("=" * 40 + "\n")
        else:
            if resume:
                vprint("\nNote: Resume requested but not all output files exist. Running full analysis.\n")

            if save_to_file:
                data_file = f"{output_file_prefix}_data.csv"
                dcati.to_csv(data_file.replace('.csv', '_dp1.csv'), index=False)
                dit.to_csv(data_file.replace('.csv', '_dit.csv'), index=False)
                vprint(f"Saved preprocessed data to: {data_file.replace('.csv', '_dp1.csv')} and _dit.csv")

            vprint("Preparing Stan data in ncats format...")
            stan_data = self.make_stan_data(
                dcati=dcati,
                x_formula=x_formula,
                x_formula_ignore_regex=x_formula_ignore_regex,
                verbose=verbose,
            )

            vprint("Loading NumPyro model...")
            # Use the cached module-level _model_namespace() so repeated SVI
            # fits in the same Python session share the same partial_credit
            # _model_ncats function identity and therefore the JAX JIT cache.
            # Previously this was a fresh runpy.run_path per call, causing
            # 5-10x recompile cost on every algorithm in the SVI sweep.
            model_ns = _model_namespace()
            if "partial_credit_model_ncats" not in model_ns:
                raise AttributeError(
                    "Function 'partial_credit_model_ncats' not found in"
                    " partial_credit_model_ncats_v260413.pyro"
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

            guide_factory = _get_autoguide_factory(algorithm)
            rng_key = jax.random.PRNGKey(seed)

            vprint("Running SVI with Adam optimizer...")
            vprint(f"  Algorithm: {algorithm}")

            guide = guide_factory(model_for_svi)
            svi = SVI(
                model_for_svi,
                guide,
                numpyro.optim.Adam(lr),
                Trace_ELBO(),
            )

            vprint("Initializing SVI state...")
            t_init0 = time.time()
            rng_key, subkey = jax.random.split(rng_key)
            init_state = svi.init(subkey, stan_data)
            init_minutes = (time.time() - t_init0) / 60.0
            vprint(f"  SVI init completed in {init_minutes:.2f} minutes")

            vprint(f"Optimizing variational parameters for {num_steps} steps...")
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
                vprint(f"  Step {i + 1:5d} / {num_steps}: ELBO = {-losses[i]:,.1f}", flush=True)
            if len(losses) and (len(losses) % 1000) != 0:
                i = len(losses) - 1
                vprint(f"  Step {i + 1:5d} / {num_steps}: ELBO = {-losses[i]:,.1f}", flush=True)

            vprint(f"\nSVI optimization completed in {opt_minutes:.2f} minutes (after {init_minutes:.2f} min init)")

            vprint(f"Sampling {output_samples} draws from the learned approximate posterior...")
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
            vprint(f"  Posterior sampling completed in {post_minutes:.2f} minutes")

            vprint("Generating posterior predictive samples...")
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
            vprint(f"  Posterior predictive generation completed in {pred_minutes:.2f} minutes")

            vprint("Converting to ArviZ format...")
            t_idata0 = time.time()
            idata = _make_idata_from_svi_posterior(stan_data, posterior_samples_dict, predictions)
            idata_minutes = (time.time() - t_idata0) / 60.0
            vprint(f"  ArviZ conversion completed in {idata_minutes:.2f} minutes")

            vprint("Extracting timing information...")
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
            if save_to_file:
                timing_data.to_csv(timing_file, index=False)
                vprint(f"Saved timing information to: {timing_file}")

                vprint(f"Saving draws to: {draws_file}")
                idata.to_zarr(draws_file)

        if save_to_file and with_additional_analyses and not os.path.exists(f"{output_file_prefix}_intervals.pdf"):
            vprint("\nRunning additional diagnostic analyses...")
            vprint("Generating parameter plots...")

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

        if save_to_file and with_additional_analyses and not os.path.exists(f"{output_file_prefix}_ppcheck.pdf"):
            vprint("Generating posterior predictive checks...")
            if 'ypred' in idata.posterior.data_vars:
                _plot_ppcheck(idata.posterior['ypred'].values, dcati, f"{output_file_prefix}_ppcheck", verbose=verbose)

        if save_to_file and with_additional_analyses and not os.path.exists(f"{output_file_prefix}_ordered_brierscore.csv"):
            vprint("Computing ordinal Brier scores by question...")
            if 'ordinal_brier_score' in idata.posterior.data_vars:
                _compute_ordinal_brier_scores(
                    idata.posterior['ordinal_brier_score'].values,
                    dcati,
                    f"{output_file_prefix}_ordered_brierscore",
                    verbose=verbose,
                )

        if save_to_file and with_core_analyses and not _pdf_or_parts_exist(f"{output_file_prefix}_prob_by_question_fit"):
            vprint("\nGenerating fitted probability plots...")
            if 'ordered_prob_by_cat_qu_fit' in idata.posterior.data_vars:
                _plot_prob_barplots(
                    idata.posterior['ordered_prob_by_cat_qu_fit'].values,
                    dcati,
                    dit,
                    f"{output_file_prefix}_prob_by_question_fit",
                    verbose=verbose,
                )

        if (
            save_to_file
            and with_core_analyses
            and x_formula_ignore_regex is not None
            and not _pdf_or_parts_exist(f"{output_file_prefix}_prob_by_question_pr")
        ):
            vprint("\nGenerating predictive probability plots with columns matching x_formula_ignore_regex removed from X...")
            if 'ordered_prob_by_cat_qu_pr' in idata.posterior.data_vars:
                _plot_prob_barplots(
                    idata.posterior['ordered_prob_by_cat_qu_pr'].values,
                    dcati,
                    dit,
                    f"{output_file_prefix}_prob_by_question_pr",
                    verbose=verbose,
                )

        return {
            'posterior_samples': posterior_samples_dict if not can_resume else posterior_samples,
            'draws': idata,
            'timing': timing_data,
            'algorithm': algorithm,
        }

    # ------------------------------------------------------------------
    # Stan ADVI fit driver. Body verbatim; only the entry rebinds dit /
    # dcati / x_formula / seed.
    # ------------------------------------------------------------------

    def fit_stan_svi(
        self,
        output_file_prefix: str,
        stan_file: Optional[str] = None,
        x_formula: Optional[str] = None,
        x_formula_ignore_regex: Optional[str] = None,
        iter: int = 10000,
        grad_samples: int = 1,
        elbo_samples: int = 100,
        output_samples: int = 4000,
        seed: Optional[int] = None,
        show_messages: bool = False,
        show_exceptions: bool = False,
        resume: bool = False,
        with_core_analyses: bool = True,
        with_additional_analyses: bool = False,
        save_to_file: bool = True,
        verbose: bool = True,
    ) -> Dict:
        """
        Run partial credit model analysis (ncats version) using ADVI.

        When ``save_to_file`` is False the fit runs purely in memory: no data/timing
        /draws/pickle files are written, resume is disabled, and the core/additional
        analyses are skipped.
        """
        dit = self.dit
        dcati = self.dcati
        if x_formula is None:
            x_formula = self.x_formula
        if seed is None:
            seed = self.seed

        vprint = print if verbose else (lambda *args, **kwargs: None)
        if stan_file is None:
            stan_file = str(Path(__file__).parents[2] / "src" / "stan" / "partial_credit_model_ncats_v260413.stan")

        vprint("\n" + "=" * 40)
        vprint("Partial Credit Model (ncats) ADVI Analysis Configuration")
        vprint("=" * 40)
        vprint(f"Stan file: {stan_file}")
        vprint(f"Stan include dir: {os.path.dirname(stan_file)}")
        vprint(f"Data: dit with {len(dit)} items, dcati with {len(dcati)} observations")
        vprint(f"Output prefix: {output_file_prefix}")
        vprint(f"ADVI iterations: {iter}")
        vprint(f"Gradient samples: {grad_samples}")
        vprint(f"ELBO samples: {elbo_samples}")
        vprint(f"Output samples: {output_samples}")
        vprint(f"Seed: {seed}")
        vprint(f"Resume: {resume}")
        vprint(f"Core analyses: {with_core_analyses}")
        vprint(f"Additional analyses: {with_additional_analyses}")
        vprint(f"Save to file: {save_to_file}")
        vprint("=" * 40 + "\n")

        if save_to_file:
            os.makedirs(os.path.dirname(output_file_prefix), exist_ok=True)

        timing_file = f"{output_file_prefix}_timing.csv"
        draws_file = f"{output_file_prefix}_draws.zarr"
        output_file = f"{output_file_prefix}_stan.pkl"
        data_file = f"{output_file_prefix}_data.csv"

        can_resume = (
            resume
            and save_to_file
            and os.path.exists(timing_file)
            and os.path.exists(draws_file)
            and os.path.exists(data_file.replace('.csv', '_dp1.csv'))
            and os.path.exists(data_file.replace('.csv', '_dit.csv'))
            and os.path.exists(output_file)
        )

        if can_resume:
            vprint("\n" + "=" * 40)
            vprint("RESUMING from existing outputs")
            vprint("=" * 40)

            vprint(f"Loading fit from: {output_file}")
            fit = pd.read_pickle(output_file)

            vprint(f"Loading draws from: {draws_file}")
            idata = az.from_zarr(draws_file)

            vprint(f"Loading timing data from: {timing_file}")
            timing_data = pd.read_csv(timing_file)
            vprint("=" * 40 + "\n")
        else:
            if resume:
                vprint("\nNote: Resume requested but not all output files exist. Running full analysis.\n")

            if save_to_file:
                dcati.to_csv(data_file.replace('.csv', '_dp1.csv'), index=False)
                dit.to_csv(data_file.replace('.csv', '_dit.csv'), index=False)
                vprint(f"Saved preprocessed data to: {data_file.replace('.csv', '_dp1.csv')} and _dit.csv")

            vprint("Preparing Stan data in ncats format...")
            stan_data = self.make_stan_data(
                dcati=dcati,
                x_formula=x_formula,
                x_formula_ignore_regex=x_formula_ignore_regex,
                verbose=verbose,
            )

            vprint("Compiling Stan model...")
            model = CmdStanModel(
                stan_file=stan_file,
                cpp_options={'STAN_THREADS': 'TRUE'},
            )

            vprint("Running ADVI...")
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

            vprint("\nExtracting timing information...")
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

            vprint("Converting to ArviZ format...")
            idata = _make_idata_from_advi_fit(fit)

            if save_to_file:
                timing_data.to_csv(timing_file, index=False)
                vprint(f"Saved timing information to: {timing_file}")

                vprint(f"Saving model fit to: {output_file}")
                pd.to_pickle(fit, output_file)

                vprint(f"Saving draws to: {draws_file}")
                idata.to_zarr(draws_file)

        if save_to_file and with_additional_analyses and not os.path.exists(f"{output_file_prefix}_intervals.pdf"):
            vprint("\nRunning additional diagnostic analyses...")
            vprint("Generating parameter plots...")

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

        if save_to_file and with_additional_analyses and not os.path.exists(f"{output_file_prefix}_ppcheck.pdf"):
            vprint("Generating posterior predictive checks...")
            if 'ypred' in idata.posterior.data_vars:
                _plot_ppcheck(idata.posterior['ypred'].values, dcati, f"{output_file_prefix}_ppcheck", verbose=verbose)

        if save_to_file and with_additional_analyses and not os.path.exists(f"{output_file_prefix}_ordered_brierscore.csv"):
            vprint("Computing ordinal Brier scores by question...")
            if 'ordinal_brier_score' in idata.posterior.data_vars:
                _compute_ordinal_brier_scores(
                    idata.posterior['ordinal_brier_score'].values,
                    dcati,
                    f"{output_file_prefix}_ordered_brierscore",
                    verbose=verbose,
                )

        if save_to_file and with_core_analyses and not _pdf_or_parts_exist(f"{output_file_prefix}_prob_by_question_fit"):
            vprint("\nGenerating fitted probability plots...")
            if 'ordered_prob_by_cat_qu_fit' in idata.posterior.data_vars:
                _plot_prob_barplots(
                    idata.posterior['ordered_prob_by_cat_qu_fit'].values,
                    dcati,
                    dit,
                    f"{output_file_prefix}_prob_by_question_fit",
                    verbose=verbose,
                )

        if (
            save_to_file
            and with_core_analyses
            and x_formula_ignore_regex is not None
            and not _pdf_or_parts_exist(f"{output_file_prefix}_prob_by_question_pr")
        ):
            vprint("\nGenerating predictive probability plots with columns matching x_formula_ignore_regex removed from X...")
            if 'ordered_prob_by_cat_qu_pr' in idata.posterior.data_vars:
                _plot_prob_barplots(
                    idata.posterior['ordered_prob_by_cat_qu_pr'].values,
                    dcati,
                    dit,
                    f"{output_file_prefix}_prob_by_question_pr",
                    verbose=verbose,
                )

        return {
            'fit': fit,
            'draws': idata,
            'timing': timing_data,
        }

    # ------------------------------------------------------------------
    # Stan HMC fit driver. Body verbatim; only the entry rebinds dit /
    # dcati / x_formula / seed.
    # ------------------------------------------------------------------

    def fit_stan_hmc(
        self,
        output_file_prefix: str,
        stan_file: Optional[str] = None,
        x_formula: Optional[str] = None,
        x_formula_ignore_regex: Optional[str] = None,
        chains: int = 2,
        parallel_chains: int = 2,
        threads_per_chain: int = 1,
        iter_warmup: int = 500,
        iter_sampling: int = 1500,
        seed: Optional[int] = None,
        show_messages: bool = False,
        show_exceptions: bool = False,
        resume: bool = False,
        with_core_analyses: bool = True,
        with_additional_analyses: bool = False,
        save_to_file: bool = True,
        verbose: bool = True,
    ) -> Dict:
        """
        Run partial credit model analysis (ncats version) on pre-processed data using HMC.

        When ``save_to_file`` is False the fit runs purely in memory: no data/timing
        /draws/pickle/convergence files are written, resume is disabled, and the
        core/additional analyses are skipped.
        """
        dit = self.dit
        dcati = self.dcati
        if x_formula is None:
            x_formula = self.x_formula
        if seed is None:
            seed = self.seed

        vprint = print if verbose else (lambda *args, **kwargs: None)
        if stan_file is None:
            stan_file = str(Path(__file__).parents[2] / "src" / "stan" / "partial_credit_model_ncats_v260413.stan")

        vprint("\n" + "=" * 40)
        vprint("Partial Credit Model (ncats) Analysis Configuration")
        vprint("=" * 40)
        vprint(f"Stan file: {stan_file}")
        vprint(f"Stan include dir: {os.path.dirname(stan_file)}")
        vprint(f"Data: dit with {len(dit)} items, dcati with {len(dcati)} observations")
        vprint(f"Output prefix: {output_file_prefix}")
        vprint(f"Chains: {chains}")
        vprint(f"Parallel chains: {parallel_chains}")
        vprint(f"Threads per chain: {threads_per_chain}")
        vprint(f"Warmup iterations: {iter_warmup}")
        vprint(f"Sampling iterations: {iter_sampling}")
        vprint(f"Seed: {seed}")
        vprint(f"Resume: {resume}")
        vprint(f"Core analyses: {with_core_analyses}")
        vprint(f"Additional analyses: {with_additional_analyses}")
        vprint(f"Save to file: {save_to_file}")
        vprint("=" * 40 + "\n")

        if save_to_file:
            os.makedirs(os.path.dirname(output_file_prefix), exist_ok=True)

        timing_file = f"{output_file_prefix}_timing.csv"
        draws_file = f"{output_file_prefix}_draws.zarr"
        output_file = f"{output_file_prefix}_stan.pkl"
        mixing_file = f"{output_file_prefix}_convergence_mixing.csv"
        data_file = f"{output_file_prefix}_data.csv"

        can_resume = (resume and
                      save_to_file and
                      os.path.exists(timing_file) and
                      os.path.exists(draws_file) and
                      os.path.exists(output_file) and
                      os.path.exists(mixing_file) and
                      os.path.exists(data_file.replace('.csv', '_dp1.csv')) and
                      os.path.exists(data_file.replace('.csv', '_dit.csv')))

        if can_resume:
            vprint("\n" + "=" * 40)
            vprint("RESUMING from existing outputs")
            vprint("=" * 40)

            vprint(f"Loading fit from: {output_file}")
            fit = pd.read_pickle(output_file)

            vprint(f"Loading draws from: {draws_file}")
            idata = az.from_zarr(draws_file)

            vprint(f"Loading timing data from: {timing_file}")
            timing_data = pd.read_csv(timing_file)
            chain_col = 'chain_id' if 'chain_id' in timing_data.columns else 'chain'
            good_chains = timing_data[chain_col].astype(int).tolist()
            vprint(f"Identified good HMC chains: {', '.join(map(str, good_chains))}")
            vprint("=" * 40 + "\n")

        else:
            if resume:
                vprint("\nNote: Resume requested but not all output files exist. Running full analysis.\n")

            if save_to_file:
                dcati.to_csv(data_file.replace('.csv', '_dp1.csv'), index=False)
                dit.to_csv(data_file.replace('.csv', '_dit.csv'), index=False)
                vprint(f"Saved preprocessed data to: {data_file.replace('.csv', '_dp1.csv')} and _dit.csv")

            vprint("Preparing Stan data in ncats format...")
            stan_data = self.make_stan_data(
                dcati=dcati,
                x_formula=x_formula,
                x_formula_ignore_regex=x_formula_ignore_regex,
                verbose=verbose,
            )

            vprint("Compiling Stan model...")
            model = CmdStanModel(
                stan_file=stan_file,
                cpp_options={'STAN_THREADS': 'TRUE'}
            )

            vprint("Running MCMC sampling...")
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

            vprint("Checking for divergent transitions and removing non-converged chains...")
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
            vprint(f"Identified good HMC chains: {', '.join(map(str, good_chains))}")

            vprint("\nExtracting timing information...")
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

            vprint("Converting to ArviZ format...")
            # Convert with all chains first; arviz refuses a chain coord shorter
            # than the data. Then drop chains that failed the lp threshold above.
            idata = az.from_cmdstanpy(fit, dims={})
            good_chain_idx = [list(fit.chain_ids).index(int(c)) for c in good_chains]
            if len(good_chain_idx) < len(fit.chain_ids):
                idata = idata.isel(chain=good_chain_idx)

            if save_to_file:
                timing_data.to_csv(timing_file, index=False)
                vprint(f"Saved timing information to: {timing_file}")

                vprint(f"Saving model fit to: {output_file}")
                pd.to_pickle(fit, output_file)

                vprint(f"Saving draws to: {draws_file}")
                idata.to_zarr(draws_file)

                vprint("Generating convergence diagnostics...")
                key_vars = ['latent_factor_unit', 'latent_factor_beta',
                           'skill_thresholds', 'loadings_questions_m1']
                summary_df = az.summary(idata, var_names=key_vars)
                summary_df = summary_df.sort_values('ess_bulk')

                vprint(f"Saved convergence diagnostics to: {mixing_file}")
                summary_df.to_csv(mixing_file)

            if save_to_file and with_additional_analyses and not os.path.exists(f"{output_file_prefix}_worsttrace.pdf"):
                vprint("Generating worst-chain trace plot...")

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

        if save_to_file and with_additional_analyses and not os.path.exists(f"{output_file_prefix}_intervals.pdf"):
            vprint("\nRunning additional diagnostic analyses...")

            vprint("Generating parameter plots...")
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

        if save_to_file and with_additional_analyses and not os.path.exists(f"{output_file_prefix}_ppcheck.png"):
            vprint("Generating posterior predictive checks...")

            if 'ypred' in idata.posterior.data_vars:
                _plot_ppcheck(idata.posterior['ypred'].values, dcati, f"{output_file_prefix}_ppcheck", verbose=verbose)

        if save_to_file and with_additional_analyses and not os.path.exists(f"{output_file_prefix}_ordered_brierscore.csv"):
            vprint("Computing ordinal Brier scores by question...")

            if 'ordinal_brier_score' in idata.posterior.data_vars:
                _compute_ordinal_brier_scores(idata.posterior['ordinal_brier_score'].values, dcati, f"{output_file_prefix}_ordered_brierscore", verbose=verbose)

        if save_to_file and with_core_analyses and not _pdf_or_parts_exist(f"{output_file_prefix}_prob_by_question_fit"):
            vprint("\nGenerating fitted probability plots...")
            if 'ordered_prob_by_cat_qu_fit' in idata.posterior.data_vars:
                _plot_prob_barplots(
                    idata.posterior['ordered_prob_by_cat_qu_fit'].values,
                    dcati,
                    dit,
                    f"{output_file_prefix}_prob_by_question_fit",
                    verbose=verbose,
                )

        if save_to_file and with_core_analyses and x_formula_ignore_regex is not None and not _pdf_or_parts_exist(f"{output_file_prefix}_prob_by_question_pr"):
            vprint("\nGenerating predictive probability plots with columns matching x_formula_ignore_regex removed from X...")
            if 'ordered_prob_by_cat_qu_pr' in idata.posterior.data_vars:
                _plot_prob_barplots(
                    idata.posterior['ordered_prob_by_cat_qu_pr'].values,
                    dcati,
                    dit,
                    f"{output_file_prefix}_prob_by_question_pr",
                    verbose=verbose,
                )

        return {
            'fit': fit,
            'draws': idata,
            'timing': timing_data,
            'good_chains': good_chains
        }

    # ------------------------------------------------------------------
    # Prior-predictive sampler for the deepset amortisers (§8 + §14 of
    # dev/amortised_decision_making.md)
    # ------------------------------------------------------------------

    @staticmethod
    def make_training_data_with_participant_tokens_prior(
        rng: np.random.Generator,
        S: int,
        *,
        item_type_idx: np.ndarray,
        item_high_idx: np.ndarray,
        n_max: int,
        n_min: int = 2,
        K_out_of_7: int = 8,
        K_categorical: int = 4,
        categorical_threshold: int = 2,
        queries_per_sample: int = 4,
        target_clip: float = 20.0,
        threshold_scale: float = 3.5,
    ):
        """Prior-predictive sampler for the deepset{S,X}compAtt amortisers
        on the PCM (§8 algorithm, §14.1 dimensions).

        Per training sample ``s``:

        1. draw ``n ~ U{n_min, .., n_max - 1}``, ``m = n_max - n``;
        2. draw ``theta_i ~ N(0, 1)`` for all ``N = n_max`` participants
           (ZeroSumNormal marginal), ``beta_t ~ N(0, 1)`` for the two
           time-points (x_formula "~ time - 1");
        3. per (real item j, time t): thresholds
           ``tau_{j, t, .} ~ N(0, threshold_scale)`` (K_j - 1 of them)
           and loading ``lam_{j, t} ~ |StudentT(3)|`` with the first
           fake item per category type pinned to 1 (Stan convention);
        4. simulate ``y_{i, j, t} ~ CategoricalLogit(eta)`` with the PCM
           cumulative-increment logits;
        5. label ``rho_j(theta)``: population endpoint per real item,
           direction-aware ratio across the two time-points of
           ``w_{j, t} = mean_i E[y | theta_i]`` (out-of-7) or
           ``mean_i P(y >= categorical_threshold | theta_i)``
           (categorical), clipped at +-target_clip.

        Responses enter the batch rescaled to [0, 1] via
        ``(k - 1) / (K_j - 1)`` (1-indexed category k), matching the
        §14.1 deployment-side input rescale.

        Returns ``(batch, rho)`` in the deepset batch schema with
        ``R = 2`` (baseline, endline response) and
        ``item_metadata = (item_type_idx, item_high_idx)``, ``aux =
        (n / n_max, m / n_max)``; ``B = S * queries_per_sample``.
        """
        item_type_idx = np.asarray(item_type_idx, dtype=np.float32)
        item_high_idx = np.asarray(item_high_idx, dtype=np.float32)
        J = item_type_idx.shape[0]
        Q = int(queries_per_sample)
        N = int(n_max)

        type_groups = [
            (np.flatnonzero(item_type_idx == 0.0), K_out_of_7),
            (np.flatnonzero(item_type_idx == 1.0), K_categorical),
        ]

        n_arr = rng.integers(n_min, n_max, size=S).astype(np.int64)
        m_arr = (n_max - n_arr).astype(np.int64)

        Y_all = np.zeros((S, N, J, 2), dtype=np.float32)   # rescaled [0, 1]
        rho_all = np.zeros((S, J), dtype=np.float64)

        for s in range(S):
            theta = rng.standard_normal(N)
            beta = rng.standard_normal(2)
            for js, K in type_groups:
                if js.size == 0:
                    continue
                Jt = js.size
                Km1 = K - 1
                tau = threshold_scale * rng.standard_normal((Jt, 2, Km1))
                lam = np.abs(rng.standard_t(3.0, size=(Jt, 2)))
                lam[0, 0] = 1.0
                # eta increments: (N, Jt, 2, K-1)
                latent = theta[:, None, None] + beta[None, None, :]
                inc = lam[None, :, :, None] * (
                    latent[:, :, :, None] - tau[None, :, :, :]
                )
                eta = np.concatenate(
                    [np.zeros((N, Jt, 2, 1)), np.cumsum(inc, axis=-1)],
                    axis=-1,
                )                                           # (N, Jt, 2, K)
                eta -= eta.max(axis=-1, keepdims=True)
                p = np.exp(eta)
                p /= p.sum(axis=-1, keepdims=True)
                # sample category 1..K via inverse cdf
                cdf = np.cumsum(p, axis=-1)
                u = rng.random((N, Jt, 2))
                k_idx = (cdf < u[..., None]).sum(axis=-1)   # 0..K-1
                Y_all[s, :, js, :] = (
                    k_idx.astype(np.float32) / float(K - 1)
                ).transpose(1, 0, 2)
                # population endpoint per (j, t)
                p_bar = p.mean(axis=0)                       # (Jt, 2, K)
                if K == K_out_of_7 and js.size and item_type_idx[js[0]] == 0.0:
                    ks = np.arange(1, K + 1, dtype=np.float64)
                    w_jt = (p_bar * ks[None, None, :]).sum(axis=-1)  # (Jt, 2)
                else:
                    w_jt = p_bar[:, :, categorical_threshold - 1:].sum(axis=-1)
                w_base = np.maximum(w_jt[:, 0], 1e-3)
                w_end = w_jt[:, 1]
                high = item_high_idx[js] > 0.5
                ratio = np.where(high, w_end / w_base - 1.0,
                                 1.0 - w_end / w_base)
                rho_all[s, js] = np.clip(ratio, -target_clip, target_clip)

        col = np.arange(N)[None, :]
        mask_x = (col < n_arr[:, None]).astype(np.float32)   # (S, N)
        mask_z = (col >= n_arr[:, None]).astype(np.float32)  # (S, N)
        # x cohort = first n participants; z cohort = remaining m. Zero out
        # the complementary rows so padded slots carry no signal.
        Y_x = Y_all * mask_x[:, :, None, None]
        Y_z = Y_all * mask_z[:, :, None, None]

        j_star = rng.integers(0, J, size=(S, Q))
        B = S * Q
        x_rep = np.broadcast_to(
            Y_x[:, None], (S, Q, N, J, 2),
        ).reshape(B, N, J, 2).astype(np.float32)
        z_rep = np.broadcast_to(
            Y_z[:, None], (S, Q, N, J, 2),
        ).reshape(B, N, J, 2).astype(np.float32)
        mask_x_rep = np.broadcast_to(
            mask_x[:, None], (S, Q, N),
        ).reshape(B, N).astype(np.float32)
        mask_z_rep = np.broadcast_to(
            mask_z[:, None], (S, Q, N),
        ).reshape(B, N).astype(np.float32)
        meta = np.stack([item_type_idx, item_high_idx], axis=-1)  # (J, 2)
        meta_rep = np.broadcast_to(
            meta[None], (B, J, 2),
        ).astype(np.float32)
        sizes = np.stack([n_arr / n_max, m_arr / n_max], axis=-1)  # (S, 2)
        aux_rep = np.broadcast_to(
            sizes[:, None], (S, Q, 2),
        ).reshape(B, 2).astype(np.float32)
        s_idx = np.arange(S)[:, None]
        rho = rho_all[s_idx, j_star].reshape(B).astype(np.float32)
        return (
            {
                'x_responses':   x_rep,
                'mask_x':        mask_x_rep,
                'z_responses':   z_rep,
                'mask_z':        mask_z_rep,
                'item_metadata': meta_rep,
                'query_idx':     j_star.reshape(B).astype(np.int32),
                'aux':           aux_rep,
            },
            rho,
        )
