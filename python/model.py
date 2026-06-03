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
    :meth:`stack_posterior_theta` to flatten arviz draws into a
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
    def stack_posterior_theta(self, draws) -> dict:
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
