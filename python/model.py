"""
``Model`` blueprint for the interim PPS / label-estimation algorithms.

A ``Model`` instance bundles a generative model with the data it was
trained on (``dit``, ``dcati``, ``stan_data``) and exposes the
loglik / prior / endpoint callables that the IS / MM / SMC / regression
algorithms in ``fit_interim.py`` consume.

Subclasses fall into two camps:

* **IRT family** (PCM, credit, ordered-logit) — share a near-identical
  wrapper pattern, so the boilerplate lives in :class:`_IRTModel`. Each
  concrete IRT subclass only needs to declare ``param_names``,
  ``positive_params``, ``_module`` (the underlying numpyro namespace),
  and ``_make_stan_data`` (the model-specific stan-data builder).
* **Binomial** — different data shape and posterior; subclasses
  :class:`Model` directly.

See ``dev/interim_model_blueprint_plan.md`` for the full design.
"""

from abc import ABC, abstractmethod
from typing import Any, Dict

import arviz as az
import jax.numpy as jnp
import numpy as np
import pandas as pd


class Model(ABC):
    """Base blueprint. Holds (dit, dcati, stan_data) on the instance.
    Concrete subclasses must override the abstract methods; the three
    Stan/Pyro fitters all return the existing
    ``{'draws', 'timing', 'posterior_samples'}`` dict shape so downstream
    code is unchanged."""

    # ---- 2.1 Subclass-declared metadata --------------------------------
    param_names: tuple = ()
    """Posterior parameter names this model exposes. Used by
    :meth:`get_stacked_posterior` to flatten arviz draws into a
    ``(K, ...)`` jnp array dict."""

    positive_params: tuple = ()
    """Subset of ``param_names`` that live on the positive reals --
    log-transformed into the MM 'match space' (Paananen 2021)."""

    # ---- 2.2 Construction ----------------------------------------------
    def __init__(self, dit: pd.DataFrame, dcati: pd.DataFrame,
                 x_formula: str = "~ time - 1", *, seed: int = 123):
        self.dit = dit
        self.dcati = dcati
        self.x_formula = x_formula
        self.seed = seed
        self.stan_data = self.make_stan_data(dcati, x_formula)

    @abstractmethod
    def make_stan_data(self, dcati: pd.DataFrame, x_formula: str) -> dict:
        """Build the data dict consumed by :meth:`eval_loglik` /
        :meth:`eval_outcome_for_endpoint`. Reused by interim algorithms with
        a z-cohort ``dcati`` to produce a fresh stan_data dict without
        instantiating another :class:`Model`."""

    # ---- 2.3 Fitters ----------------------------------------------------
    @abstractmethod
    def fit_pyro_svi(self, output_file_prefix: str, *,
                     save_to_file: bool = True, resume: bool = False,
                     verbose: bool = True, **method_kwargs) -> Dict[str, Any]:
        """Pyro / NumPyro SVI fit. Method-specific tuning (lr, num_steps,
        algorithm) goes via ``**method_kwargs``."""

    @abstractmethod
    def fit_stan_svi(self, output_file_prefix: str, *,
                     save_to_file: bool = True, resume: bool = False,
                     verbose: bool = True, **method_kwargs) -> Dict[str, Any]:
        """Stan ADVI fit. ``**method_kwargs`` carries iter, grad_samples,
        elbo_samples, output_samples."""

    @abstractmethod
    def fit_stan_hmc(self, output_file_prefix: str, *,
                     save_to_file: bool = True, resume: bool = False,
                     verbose: bool = True, **method_kwargs) -> Dict[str, Any]:
        """Stan HMC / NUTS fit. ``**method_kwargs`` carries chains,
        iter_warmup, iter_sampling."""

    def fit_closed_form(self, output_file_prefix: str, *,
                        save_to_file: bool = True, verbose: bool = True,
                        **method_kwargs) -> Dict[str, Any]:
        """Closed-form posterior sampler. Defaults to ``NotImplementedError``
        for models without a closed-form posterior (PCM, credit,
        ordered-logit). The Binomial subclass overrides with a Beta-Binomial
        conjugate sampler."""
        raise NotImplementedError(
            f"{type(self).__name__} has no closed-form posterior."
        )

    # ---- 2.4 Posterior <-> params --------------------------------------
    def get_stacked_posterior(self, draws) -> dict:
        """Flatten the (chain, draw) dims of ``draws.posterior[name]`` for
        every ``name in self.param_names`` into a dict of ``(K, ...)`` jnp
        arrays. Only the parameters declared by the subclass are returned.
        """
        post = draws.posterior
        return {
            name: jnp.asarray(
                np.asarray(post[name].values).reshape(-1, *post[name].shape[2:])
            )
            for name in self.param_names if name in post
        }

    # ---- 2.5 Likelihood + prior ----------------------------------------
    @abstractmethod
    def eval_loglik(self, data: dict, params: dict) -> jnp.ndarray:
        """Pointwise log-likelihood vector; shape ``(N_total,)``."""

    @abstractmethod
    def eval_loglik_annealed(self, data: dict, params: dict) -> jnp.ndarray:
        """Annealed log-likelihood for SMC tempering;
        ``params['temperature']`` in ``[0, 1]`` scales the loglik."""

    @abstractmethod
    def eval_log_prior(self, params: dict) -> jnp.ndarray:
        """Scalar log-prior."""

    # ---- 2.6 Endpoint / outcome quantity per draw ----------------------
    @abstractmethod
    def eval_outcome_for_endpoint(self, data: dict, params: dict) -> jnp.ndarray:
        """Quantity the per-draw endpoint code consumes. PCM:
        ``ordered_prob_by_cat_qu_fit``. Binomial: ``p(theta)``."""

    @abstractmethod
    def get_endpoints_per_draw(self, draws, categorical_threshold: int,
                               endpoint_type: str = 'items') -> pd.DataFrame:
        """One row per (item, draw) with at least
        ``item_label, item_type, item_high_label, draw, ratio``. IRT
        subclasses inherit :class:`model_irt.IRTModel`'s implementation."""

    @abstractmethod
    def get_w(self, zi: pd.DataFrame) -> pd.DataFrame:
        """Strong-Oakley per-(item, draw) summary frame W(z^(s)) from
        the future-data block ``zi``. Must return at least columns
        ``item_label, item_type, item_high_label, draw, w_ratio`` so
        :func:`fit_interim.fit_interim_regress_endptx_on_wz` and the
        H1x sibling can consume it. IRT subclasses inherit
        :class:`model_irt.IRTModel`'s per-time pivot; Binomial overrides
        with a per-draw mean-of-ypred summary."""

    def get_p_h1(self, ratio_xz: pd.DataFrame,
                 pps_H1_def: float = 0.5) -> pd.DataFrame:
        """Per-item ``mean(endpt > pps_H1_def)`` over draws.

        Both IRT (where ``endpt`` is the directional improvement
        ``1 - p_endline/p_baseline`` for ``lower_is_better`` items) and
        Binomial (where ``endpt = 1 - p / p_0``) use the same
        right-tail rule. Returns one row per
        ``(item_label, item_type, item_high_label)`` with a ``p_h1``
        column. Callers rename to ``p_h1_x`` / ``p_h1_xz`` /
        ``pps_ProbH1_x`` etc per context."""
        return (
            ratio_xz
            .groupby(['item_label', 'item_type', 'item_high_label'])['ratio']
            .apply(lambda endpt: float((endpt > pps_H1_def).mean()))
            .reset_index(name='p_h1')
        )

    # ---- 2.6b IS / SMC scoring helpers ---------------------------------
    def make_stan_data_from_xi(self) -> dict:
        """Build stan_data for the x cohort (``self.dcati``). Default
        delegates to :meth:`make_stan_data`; IRT subclasses override to
        re-index oid/oidt/item_type sort first."""
        return self.make_stan_data(self.dcati, self.x_formula)

    def make_stan_data_from_zi(self, zi: pd.DataFrame, s_idx: int) -> dict:
        """Build stan_data for the future-data sample s drawn from ``zi``.
        Default raises ``NotImplementedError``; IRT subclasses do the
        ``ypred_s -> y_stan`` rename + index rebuild + sort + make_stan_data
        path. Binomial subclasses override per their data shape."""
        raise NotImplementedError(
            f"{type(self).__name__}.make_stan_data_from_zi"
        )

    def get_endpoints_per_draw_from_theta_batch(
        self, theta_batch, x_stan, endpoint_type: str = 'items',
    ) -> pd.DataFrame:
        """Score a batch of parameter dicts (each leaf shape ``(K, ...)``)
        on the x-cohort design ``x_stan`` and return the per-(item, draw)
        endpoint frame. Default raises ``NotImplementedError``; IRT
        subclasses vmap ``eval_outcome_for_endpoint`` -> arviz idata with
        ``ordered_prob_by_cat_qu_pr`` -> :meth:`get_endpoints_per_draw`.
        Binomial subclasses override per their posterior parameterisation."""
        raise NotImplementedError(
            f"{type(self).__name__}.get_endpoints_per_draw_from_theta_batch"
        )

    # ---- 2.7 Generic posterior-predictive utilities --------------------
    @staticmethod
    def _load_ypred(draws_file: str, pps_z_total: int, rng,
                    keep_order: bool = False) -> np.ndarray:
        """Load ``posterior['ypred']`` from a draws zarr, flatten the
        (chain, draw) dims, and return ``pps_z_total`` draws as
        ``(pps_z_total, N_total)``. ``keep_order=False`` samples without
        replacement via ``rng.choice``; ``keep_order=True`` returns the
        first ``pps_z_total`` draws in posterior order."""
        ypred = az.from_zarr(draws_file).posterior['ypred'].values
        ypred = ypred.reshape(-1, ypred.shape[-1])
        if keep_order:
            return ypred[:pps_z_total]
        return ypred[rng.choice(ypred.shape[0], size=pps_z_total, replace=False)]

    def get_interim_z_from_ypredi(self, draws_file: str, interim_m: int,
                                  pps_z_total: int = 10, seed: int = 123,
                                  keep_order: bool = False
                                  ) -> pd.DataFrame:
        """Build the future-data block z from this model's INTERIM fit
        predictions.

        Uses ``self.dcati`` as the interim cohort ``xi``. ``interim_m``
        new (shadow) participants are resampled WITH replacement from
        ``xi``, each assigned a fresh ``pid`` offset above ``xi``'s
        maximum so they don't collide on a later concat. Their predicted
        outcomes are taken from ``draws_file``'s ``posterior['ypred']``
        at the original ``xi`` observation positions; for sample ``s``
        all rows share the same posterior draw (y_s is linked to
        theta_s by the shared draw index).

        ``keep_order=False`` (default) randomly subsamples
        ``pps_z_total`` posterior draws. ``keep_order=True`` keeps
        posterior order so row ``s`` corresponds to posterior draw
        index ``s`` -- required when downstream code merges on the
        posterior-draw index (e.g. against ``get_endpoints_per_draw``).

        Requires ``self.dcati`` to carry ``pid`` and ``oid``."""
        xi = self.dcati
        rng = np.random.default_rng(seed)
        ypred = self._load_ypred(draws_file, pps_z_total, rng,
                                 keep_order=keep_order)

        src_pids = rng.choice(xi['pid'].unique(), size=interim_m, replace=True)
        tmp = pd.DataFrame({
            'src_pid': src_pids,
            'pid': np.arange(1, interim_m + 1) + int(xi['pid'].max()),
        })
        zi = tmp.merge(xi.rename(columns={'pid': 'src_pid'}), on='src_pid')
        zi.sort_values(['pid', 'oid'], inplace=True)
        zi.reset_index(drop=True, inplace=True)
        zi[[f'ypred_{s}' for s in range(pps_z_total)]] = (
            ypred[:, zi['oid'].to_numpy() - 1].T
        )
        return zi
