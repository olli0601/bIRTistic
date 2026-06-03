"""
:class:`BinomialModel` -- a Beta-Binomial generative model that
implements the :class:`model.Model` blueprint. The first non-IRT
subclass; proves the blueprint generalises beyond the partial-credit
family.

The model is the §3.1 example from ``dev/amortised_decision_making.md``:

    p ~ Beta(a, b)
    y_i | p ~ Bernoulli(p)        i = 1, ..., n

with a single "item" (no item-effect structure -- this is what makes the
endpoint code trivial). The closed-form posterior is Beta(a + k, b + n - k)
with k = sum_i y_i.

Implementation status:
- ``fit_closed_form``: direct Beta sampling (fully implemented).
- ``fit_pyro_svi`` / ``fit_stan_svi`` / ``fit_stan_hmc``: raise
  NotImplementedError for now. They will reuse the same Beta-Binomial
  numpyro / Stan program once written; the blueprint already accepts
  them as separate methods so plugging the implementations in later is a
  drop-in replacement.

The point of this first cut is to prove the architecture: the IS /
regression algorithms in ``fit_interim.py`` should run unchanged on a
BinomialModel instance because they only touch the methods declared in
the blueprint.
"""

import time
from typing import Any, Dict, Optional

import arviz as az
import jax.numpy as jnp
import jax.scipy.stats as jstats
import numpy as np
import pandas as pd
import xarray as xr

from model import Model


class BinomialModel(Model):
    """Beta-Binomial generative model. Single item; no item-effect
    structure. ``dcati`` is expected to be a long DataFrame with a ``y``
    column carrying 0/1 outcomes. ``dit`` is an item table with a single
    row whose ``item_label`` / ``item_type`` / ``item_high_label`` columns
    are echoed back by :meth:`endpoints_per_draw`."""

    param_names = ('p',)
    positive_params = ()  # p lives on (0, 1); MM would use a logit transform.

    def __init__(self, dit: pd.DataFrame, dcati: pd.DataFrame,
                 x_formula: str = "~ 1", *, seed: int = 123,
                 prior_a: float = 1.0, prior_b: float = 1.0):
        self.prior_a = float(prior_a)
        self.prior_b = float(prior_b)
        super().__init__(dit=dit, dcati=dcati, x_formula=x_formula, seed=seed)

    # ---- Stan data --------------------------------------------------------
    def make_stan_data(self, dcati: pd.DataFrame, x_formula: str) -> dict:
        y = dcati['y'].astype(int).to_numpy()
        return {
            'y': jnp.asarray(y, dtype=jnp.int32),
            'N_total': int(len(y)),
            'k': int(y.sum()),
            'a': self.prior_a,
            'b': self.prior_b,
        }

    # ---- Loglik / prior ---------------------------------------------------
    def eval_loglik(self, data: dict, params: dict) -> jnp.ndarray:
        """Pointwise log Bernoulli(y_i | p). Returns shape ``(N_total,)``."""
        p = params['p']
        y = data['y']
        return y * jnp.log(p) + (1 - y) * jnp.log1p(-p)

    def eval_loglik_annealed(self, data: dict, params: dict) -> jnp.ndarray:
        temperature = params.get('temperature', 1.0)
        return temperature * self.eval_loglik(data, params)

    def logprior(self, params: dict) -> jnp.ndarray:
        return jstats.beta.logpdf(params['p'], self.prior_a, self.prior_b)

    # ---- Endpoint / outcome -----------------------------------------------
    def eval_outcome_for_endpoint(self, data: dict, params: dict) -> jnp.ndarray:
        return params['p']

    def endpoints_per_draw(self, dcati: pd.DataFrame, draws,
                           categorical_threshold: int,
                           endpoint_type: str = 'items') -> pd.DataFrame:
        """One row per (item, draw) with ratio = p (the §3.1 example uses
        baseline p_0 = 0, so ratio = p / p_0 - 1 collapses to p with the
        convention used by the IS / regression algorithms)."""
        p_draws = np.asarray(draws.posterior['p'].values).reshape(-1)
        K = p_draws.size
        item_row = self.dit.iloc[0]
        return pd.DataFrame({
            'item_label':       [item_row['item_label']] * K,
            'item_type':        [item_row['item_type']] * K,
            'item_high_label':  [item_row['item_high_label']] * K,
            'draw':             np.arange(K),
            'ratio':            p_draws,
        })

    # ---- Fitters ----------------------------------------------------------
    def fit_pyro_svi(self, output_file_prefix: str, **kw) -> Dict[str, Any]:
        raise NotImplementedError(
            "BinomialModel.fit_pyro_svi is not implemented yet; use"
            " fit_closed_form instead. SVI implementation is a follow-up."
        )

    def fit_stan_svi(self, output_file_prefix: str, **kw) -> Dict[str, Any]:
        raise NotImplementedError(
            "BinomialModel.fit_stan_svi is not implemented yet; use"
            " fit_closed_form instead."
        )

    def fit_stan_hmc(self, output_file_prefix: str, **kw) -> Dict[str, Any]:
        raise NotImplementedError(
            "BinomialModel.fit_stan_hmc is not implemented yet; use"
            " fit_closed_form instead."
        )

    def fit_closed_form(self, output_file_prefix: Optional[str] = None, *,
                        save_to_file: bool = False, verbose: bool = True,
                        output_samples: int = 4000,
                        **method_kwargs) -> Dict[str, Any]:
        """Sample directly from the analytic Beta(a + k, b + n - k) posterior."""
        t0 = time.time()
        k = int(self.stan_data['k'])
        n = int(self.stan_data['N_total'])
        a_post = self.prior_a + k
        b_post = self.prior_b + (n - k)
        rng = np.random.default_rng(self.seed)
        p_samples = rng.beta(a_post, b_post, size=output_samples)
        # Wrap as a (chain=1, draw=output_samples) arviz idata for parity with
        # the other fitters.
        ds = xr.Dataset({'p': (('chain', 'draw'), p_samples[None, :])})
        draws = az.InferenceData(posterior=ds)
        timing = {'mins': (time.time() - t0) / 60.0}
        if verbose:
            print(f"[BinomialModel] closed-form posterior: Beta({a_post:.2f},"
                  f" {b_post:.2f}); drew {output_samples} samples in"
                  f" {timing['mins']:.4f} min")
        return {
            'draws': draws,
            'posterior_samples': {'p': p_samples},
            'timing': timing,
            'a_post': a_post,
            'b_post': b_post,
        }
