"""
:class:`BinomialModel` -- a Beta-Binomial generative model that
implements the :class:`model.Model` blueprint. The §3.1 example from
``dev/amortised_decision_making.md``:

    p ~ Beta(a, b)
    y_i | p ~ Bernoulli(p)        i = 1, ..., n

Closed-form posterior is Beta(a + k, b + n - k) with k = sum_i y_i.

Implementation:
- ``fit_closed_form_posterior``: direct Beta sampling.
- ``fit_closed_form_pps``: analytic Beta-Binomial PPS tail sum (§3.1).
- ``fit_pyro_svi``: NumPyro SVI with an Adam + Trace_ELBO loop on an
  inline Beta-Bernoulli program.
- ``fit_stan_svi`` / ``fit_stan_hmc``: cmdstanpy variational + NUTS on
  ``src/stan/binomial_v260603.stan``.
"""

from pathlib import Path
import time
from typing import Any, Dict, Optional

import arviz as az
import jax
import jax.numpy as jnp
import jax.scipy.stats as jstats
import numpy as np
import numpyro
import numpyro.distributions as dist
from numpyro.infer import SVI, Trace_ELBO
from numpyro.infer.autoguide import AutoDiagonalNormal
import pandas as pd
from scipy.stats import beta as _scipy_beta
from scipy.stats import betabinom as _scipy_betabinom
import xarray as xr

from model import Model


def _bernoulli_numpyro_model(N_total, y, a, b):
    """Inline Beta-Bernoulli program. Used by :meth:`BinomialModel.fit_pyro_svi`."""
    p = numpyro.sample('p', dist.Beta(a, b))
    with numpyro.plate('obs', N_total):
        numpyro.sample('y', dist.Bernoulli(probs=p), obs=y)


class BinomialModel(Model):
    """Beta-Binomial generative model. Single item; ``dcati`` is a long
    DataFrame with a ``y`` column carrying 0/1 outcomes. ``dit`` is an
    item table with a single row whose ``item_label`` / ``item_type`` /
    ``item_high_label`` columns are echoed back by
    :meth:`get_endpoints_per_draw`."""

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

    def eval_log_prior(self, params: dict) -> jnp.ndarray:
        return jstats.beta.logpdf(params['p'], self.prior_a, self.prior_b)

    # ---- Endpoint / outcome -----------------------------------------------
    def eval_outcome_for_endpoint(self, data: dict, params: dict) -> jnp.ndarray:
        return params['p']

    def get_endpoints_per_draw(self, draws,
                               categorical_threshold: int = 3,
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

    # ------------------------------------------------------------------
    # Closed-form posterior + closed-form PPS
    # ------------------------------------------------------------------

    def fit_closed_form_posterior(
        self,
        output_file_prefix: Optional[str] = None, *,
        save_to_file: bool = False,
        verbose: bool = True,
        output_samples: int = 4000,
        **method_kwargs,
    ) -> Dict[str, Any]:
        """Sample directly from the analytic Beta(a + k, b + n - k) posterior."""
        t0 = time.time()
        k = int(self.stan_data['k'])
        n = int(self.stan_data['N_total'])
        a_post = self.prior_a + k
        b_post = self.prior_b + (n - k)
        rng = np.random.default_rng(self.seed)
        p_samples = rng.beta(a_post, b_post, size=output_samples)
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

    def fit_closed_form_pps(
        self,
        m: int,
        *,
        pps_H1_def: float = 0.5,
        pps_ProbH1_thresh: float = 0.89,
    ) -> float:
        """Analytic Beta-Binomial PPS over ``m`` future trials (§3.1).

        Posterior after observing x: ``Beta(a_post, b_post)`` with
        ``a_post = prior_a + k``, ``b_post = prior_b + n - k``.

        For each possible future success count ``k_m in {0, ..., m}``:
            posterior after (x, z): ``Beta(a_post + k_m, b_post + m - k_m)``
            P(H_1 | x, z) = 1 - F_Beta(pps_H1_def; a_post + k_m, b_post + m - k_m)

        PPS is the Beta-Binomial(m, a_post, b_post) tail mass at the
        smallest ``k_m`` for which ``P(H_1 | x, z) > pps_ProbH1_thresh``
        (the decision is monotone increasing in ``k_m``).
        """
        if m < 0:
            raise ValueError("m must be non-negative.")
        k = int(self.stan_data['k'])
        n = int(self.stan_data['N_total'])
        a_post = self.prior_a + k
        b_post = self.prior_b + (n - k)

        k_m_grid = np.arange(m + 1)
        a_z = a_post + k_m_grid
        b_z = b_post + (m - k_m_grid)
        p_h1_given_z = 1.0 - _scipy_beta.cdf(pps_H1_def, a_z, b_z)
        crosses = p_h1_given_z > pps_ProbH1_thresh
        if not crosses.any():
            return 0.0
        k_star = int(k_m_grid[crosses][0])
        # Tail probability under Beta-Binomial(m, a_post, b_post).
        tail = _scipy_betabinom.sf(k_star - 1, m, a_post, b_post)
        return float(tail)

    # ------------------------------------------------------------------
    # NumPyro SVI
    # ------------------------------------------------------------------

    def fit_pyro_svi(
        self,
        output_file_prefix: Optional[str] = None,
        *,
        lr: float = 0.05,
        num_steps: int = 2000,
        output_samples: int = 4000,
        save_to_file: bool = False,
        resume: bool = False,
        verbose: bool = True,
        **method_kwargs,
    ) -> Dict[str, Any]:
        """NumPyro SVI fit on the inline Beta-Bernoulli program."""
        t0 = time.time()
        rng = jax.random.PRNGKey(self.seed)
        guide = AutoDiagonalNormal(_bernoulli_numpyro_model)
        svi = SVI(_bernoulli_numpyro_model, guide,
                  numpyro.optim.Adam(lr), Trace_ELBO())
        rng, sub = jax.random.split(rng)
        run = svi.run(
            sub, num_steps,
            int(self.stan_data['N_total']),
            jnp.asarray(self.stan_data['y']),
            float(self.stan_data['a']),
            float(self.stan_data['b']),
            progress_bar=False,
        )
        rng, sub = jax.random.split(rng)
        if hasattr(svi, 'get_posterior'):
            post = svi.get_posterior(run.state).sample(
                sub, sample_shape=(output_samples,),
            )
        else:
            params = svi.get_params(run.state)
            post = guide.sample_posterior(sub, params,
                                          sample_shape=(output_samples,))
        p_samples = np.asarray(post['p'])
        ds = xr.Dataset({'p': (('chain', 'draw'), p_samples[None, :])})
        draws = az.InferenceData(posterior=ds)
        timing = {'mins': (time.time() - t0) / 60.0,
                  'final_elbo': float(-run.losses[-1])}
        if verbose:
            print(f"[BinomialModel] pyro SVI: {num_steps} steps,"
                  f" final ELBO={timing['final_elbo']:.1f},"
                  f" {output_samples} posterior draws in"
                  f" {timing['mins']:.4f} min")
        return {
            'draws': draws,
            'posterior_samples': {'p': p_samples},
            'timing': timing,
        }

    # ------------------------------------------------------------------
    # Stan fitters
    # ------------------------------------------------------------------

    _STAN_FILE = (
        Path(__file__).resolve().parents[1]
        / 'src' / 'stan' / 'binomial_v260603.stan'
    )

    def _stan_data_native(self) -> dict:
        return {
            'N': int(self.stan_data['N_total']),
            'y': np.asarray(self.stan_data['y'], dtype=int).tolist(),
            'a': self.prior_a,
            'b': self.prior_b,
        }

    def fit_stan_hmc(
        self,
        output_file_prefix: Optional[str] = None,
        *,
        chains: int = 2,
        iter_warmup: int = 500,
        iter_sampling: int = 1500,
        stan_file: Optional[str] = None,
        save_to_file: bool = False,
        resume: bool = False,
        verbose: bool = True,
        **method_kwargs,
    ) -> Dict[str, Any]:
        """Stan NUTS via cmdstanpy on the Beta-Bernoulli model."""
        from cmdstanpy import CmdStanModel
        t0 = time.time()
        sf = stan_file or str(self._STAN_FILE)
        model = CmdStanModel(stan_file=sf)
        fit = model.sample(
            data=self._stan_data_native(),
            chains=chains,
            iter_warmup=iter_warmup,
            iter_sampling=iter_sampling,
            seed=self.seed,
            show_console=False,
            show_progress=False,
        )
        idata = az.from_cmdstanpy(fit)
        p_samples = np.asarray(idata.posterior['p'].values).reshape(-1)
        timing = {'mins': (time.time() - t0) / 60.0}
        if verbose:
            print(f"[BinomialModel] Stan HMC: {chains} chains x"
                  f" {iter_sampling} draws in {timing['mins']:.4f} min")
        return {
            'fit': fit,
            'draws': idata,
            'posterior_samples': {'p': p_samples},
            'timing': timing,
        }

    def fit_stan_svi(
        self,
        output_file_prefix: Optional[str] = None,
        *,
        iter: int = 5000,
        grad_samples: int = 1,
        elbo_samples: int = 100,
        output_samples: int = 4000,
        stan_file: Optional[str] = None,
        save_to_file: bool = False,
        resume: bool = False,
        verbose: bool = True,
        **method_kwargs,
    ) -> Dict[str, Any]:
        """Stan ADVI via cmdstanpy on the Beta-Bernoulli model."""
        from cmdstanpy import CmdStanModel
        t0 = time.time()
        sf = stan_file or str(self._STAN_FILE)
        model = CmdStanModel(stan_file=sf)
        fit = model.variational(
            data=self._stan_data_native(),
            algorithm='meanfield',
            iter=iter,
            grad_samples=grad_samples,
            elbo_samples=elbo_samples,
            output_samples=output_samples,
            seed=self.seed,
            show_console=False,
        )
        # cmdstanpy variational draws are accessed via variational_sample_pd.
        p_samples = np.asarray(fit.variational_sample_pd['p']).reshape(-1)
        ds = xr.Dataset({'p': (('chain', 'draw'), p_samples[None, :])})
        idata = az.InferenceData(posterior=ds)
        timing = {'mins': (time.time() - t0) / 60.0}
        if verbose:
            print(f"[BinomialModel] Stan ADVI: {iter} iter,"
                  f" {output_samples} draws in {timing['mins']:.4f} min")
        return {
            'fit': fit,
            'draws': idata,
            'posterior_samples': {'p': p_samples},
            'timing': timing,
        }
