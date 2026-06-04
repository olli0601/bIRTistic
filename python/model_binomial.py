"""
:class:`BinomialModel` -- Beta-Bernoulli generative model. §3.1 example.

    p ~ Beta(a, b)
    y_i | p ~ Bernoulli(p)       i = 1, ..., N

Closed-form posterior is Beta(a + k, b + n - k) with k = sum_i y_i.
All loglik / prior / endpoint kernels live in
``src/numpyro/binomial_v260603.pyro``; this class wraps them with the
``Model`` blueprint surface (``make_stan_data`` / fit drivers /
endpoint helpers).

Fitters: closed-form posterior, closed-form PPS (analytic Beta-Binomial
tail sum), NumPyro SVI + HMC, Stan ADVI + HMC.
"""

from functools import partial
import os
from pathlib import Path
import runpy
import time
from typing import Any, Dict, Optional

import arviz as az
import jax
import jax.numpy as jnp
import numpy as np
import numpyro
from numpyro.infer import MCMC, NUTS, SVI, Trace_ELBO
import pandas as pd
from scipy.stats import beta as _scipy_beta
from scipy.stats import betabinom as _scipy_betabinom
import xarray as xr

from model import Model
from utils import _get_autoguide_factory


_MODEL_NS = None


def _model_namespace():
    """Lazily load + cache the binomial numpyro model namespace."""
    global _MODEL_NS
    if _MODEL_NS is None:
        model_file = Path(__file__).resolve().parents[1] / 'src' / 'numpyro' / 'binomial_v260603.pyro'
        if not model_file.exists():
            raise FileNotFoundError(f"NumPyro model file not found: {model_file}")
        _MODEL_NS = runpy.run_path(str(model_file))
    return _MODEL_NS


class BinomialModel(Model):
    """Beta-Bernoulli generative model. ``dcati`` is a long DataFrame
    with a ``y`` column carrying 0/1 outcomes; ``dit`` is an item table
    with one row (single proportion). ``p_0`` is the null reference for
    the endpoint ratio ``ratio = p / p_0``."""

    param_names = ('p',)
    positive_params = ()  # p lives on (0, 1); MM would use a logit transform.

    def __init__(self, dit: pd.DataFrame, dcati: pd.DataFrame,
                 x_formula: str = "~ 1", *, seed: int = 123,
                 prior_a: float = 1.0, prior_b: float = 1.0,
                 p_0: float = 0.5):
        self.prior_a = float(prior_a)
        self.prior_b = float(prior_b)
        self.p_0 = float(p_0)
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

    # ---- Loglik / prior / outcome (delegate to the .pyro namespace) -------

    @staticmethod
    def eval_loglik(data, params):
        return _model_namespace()['get_log_likelihood_of_binomial'](data, params)

    @staticmethod
    def eval_loglik_annealed(data, params):
        return _model_namespace()['get_log_likelihood_of_binomial_with_annealing'](data, params)

    def eval_log_prior(self, params):
        """Log Beta(prior_a, prior_b) prior. The model wrapper passes the
        instance's ``prior_a`` / ``prior_b`` to the .pyro function (which
        defaults to Beta(1, 1) when called standalone)."""
        return _model_namespace()['get_log_prior_of_binomial'](
            params, a=self.prior_a, b=self.prior_b,
        )

    def eval_outcome_for_endpoint(self, data, params):
        """Endpoint scalar ``p / p_0`` per draw."""
        return _model_namespace()['get_endpoint_of_binomial'](params, p_0=self.p_0)

    # ---- Endpoint per draw ------------------------------------------------
    def get_endpoints_per_draw(self, draws,
                               categorical_threshold: int = 3,
                               endpoint_type: str = 'items') -> pd.DataFrame:
        """One row per (item, draw) with ``ratio = p / p_0``."""
        p_draws = np.asarray(draws.posterior['p'].values).reshape(-1)
        ratio = p_draws / self.p_0
        K = p_draws.size
        item_row = self.dit.iloc[0]
        return pd.DataFrame({
            'item_label':       [item_row['item_label']] * K,
            'item_type':        [item_row['item_type']] * K,
            'item_high_label':  [item_row['item_high_label']] * K,
            'draw':             np.arange(K),
            'p':                p_draws,
            'ratio':            ratio,
        })

    def get_endpoints(self, draws,
                      categorical_threshold: int = 3,
                      endpoint_type: str = 'items',
                      verbose: bool = False) -> pd.DataFrame:
        """Quantile-summarised endpoint at q=[2.5, 25, 50, 75, 97.5]% on
        both ``p`` and the ratio ``p / p_0``. Returns one row per
        (item, variable) with columns ``q_lower``, ``iqr_lower``,
        ``median``, ``iqr_upper``, ``q_upper``."""
        per_draw = self.get_endpoints_per_draw(
            draws=draws,
            categorical_threshold=categorical_threshold,
            endpoint_type=endpoint_type,
        )
        quantiles = [0.025, 0.25, 0.5, 0.75, 0.975]
        quantile_names = ['q_lower', 'iqr_lower', 'median', 'iqr_upper', 'q_upper']
        rows = []
        item_row = self.dit.iloc[0]
        for var in ('p', 'ratio'):
            q = per_draw[var].quantile(quantiles).to_list()
            rows.append({
                'item_label': item_row['item_label'],
                'item_type':  item_row['item_type'],
                'item_high_label': item_row['item_high_label'],
                'variable':   var,
                **dict(zip(quantile_names, q)),
            })
        return pd.DataFrame(rows)

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
        ``pps_H1_def`` is the threshold on ``p`` -- by default the same
        ``p_0`` used by the endpoint ratio. ``pps_ProbH1_thresh`` is the
        decision threshold on ``P(H_1 | x, z)``."""
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
        tail = _scipy_betabinom.sf(k_star - 1, m, a_post, b_post)
        return float(tail)

    # ------------------------------------------------------------------
    # NumPyro SVI + HMC
    # ------------------------------------------------------------------

    def _numpyro_model_fn(self):
        return _model_namespace()['bernoulli_model']

    def fit_pyro_svi(
        self,
        output_file_prefix: Optional[str] = None,
        *,
        algorithm: str = 'AutoDiagonalNormal',
        lr: float = 0.05,
        num_steps: int = 2000,
        output_samples: int = 4000,
        save_to_file: bool = False,
        resume: bool = False,
        verbose: bool = True,
        **method_kwargs,
    ) -> Dict[str, Any]:
        """NumPyro SVI fit. ``algorithm`` selects the autoguide
        (``AutoDiagonalNormal`` / ``AutoLowRankMultivariateNormal`` /
        etc). When ``save_to_file`` and ``resume`` are True and the
        ``{prefix}_draws.zarr`` + ``{prefix}_timing.csv`` artifacts
        exist, the cached posterior is returned instead of refitting."""
        t0 = time.time()
        draws_file = (f"{output_file_prefix}_draws.zarr"
                      if output_file_prefix is not None else None)
        timing_file = (f"{output_file_prefix}_timing.csv"
                       if output_file_prefix is not None else None)
        can_resume = (
            resume and save_to_file
            and draws_file is not None
            and os.path.exists(draws_file)
            and timing_file is not None and os.path.exists(timing_file)
        )
        if can_resume:
            if verbose:
                print(f"[BinomialModel] pyro SVI ({algorithm}) RESUME from"
                      f" {draws_file}")
            idata = az.from_zarr(draws_file)
            p_samples = np.asarray(idata.posterior['p'].values).reshape(-1)
            timing = pd.read_csv(timing_file).iloc[0].to_dict()
            return {
                'draws': idata,
                'posterior_samples': {'p': p_samples},
                'timing': timing,
            }

        model_fn = self._numpyro_model_fn()
        model_for_svi = partial(model_fn, sample_ypred=False,
                                compute_generated=False)

        rng = jax.random.PRNGKey(self.seed)
        guide_factory = _get_autoguide_factory(algorithm)
        guide = guide_factory(model_for_svi)
        svi = SVI(model_for_svi, guide,
                  numpyro.optim.Adam(lr), Trace_ELBO())
        rng, sub = jax.random.split(rng)
        run = svi.run(sub, num_steps, self.stan_data, progress_bar=False)
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
                  'final_elbo': float(-run.losses[-1]),
                  'algorithm': algorithm}
        if verbose:
            print(f"[BinomialModel] pyro SVI ({algorithm}): {num_steps}"
                  f" steps, final ELBO={timing['final_elbo']:.1f},"
                  f" {output_samples} posterior draws in"
                  f" {timing['mins']:.4f} min")
        if save_to_file and output_file_prefix is not None:
            os.makedirs(os.path.dirname(output_file_prefix) or '.',
                        exist_ok=True)
            draws.to_zarr(draws_file)
            pd.DataFrame([timing]).to_csv(timing_file, index=False)
        return {
            'draws': draws,
            'posterior_samples': {'p': p_samples},
            'timing': timing,
        }

    def fit_pyro_hmc(
        self,
        output_file_prefix: Optional[str] = None,
        *,
        chains: int = 2,
        iter_warmup: int = 500,
        iter_sampling: int = 1500,
        save_to_file: bool = False,
        resume: bool = False,
        verbose: bool = True,
        **method_kwargs,
    ) -> Dict[str, Any]:
        """NumPyro NUTS fit on the Beta-Bernoulli program. NUTS does not
        enumerate discrete sites so ``sample_ypred`` is forced to False
        (ypred is only useful for posterior-predictive simulation, not
        inference). ``resume=True`` returns the cached
        ``{prefix}_draws.zarr`` if it exists."""
        t0 = time.time()
        draws_file = (f"{output_file_prefix}_draws.zarr"
                      if output_file_prefix is not None else None)
        timing_file = (f"{output_file_prefix}_timing.csv"
                       if output_file_prefix is not None else None)
        can_resume = (
            resume and save_to_file
            and draws_file is not None and os.path.exists(draws_file)
            and timing_file is not None and os.path.exists(timing_file)
        )
        if can_resume:
            if verbose:
                print(f"[BinomialModel] pyro HMC RESUME from {draws_file}")
            idata = az.from_zarr(draws_file)
            p_samples = np.asarray(idata.posterior['p'].values).reshape(-1)
            timing = pd.read_csv(timing_file).iloc[0].to_dict()
            return {
                'draws': idata,
                'posterior_samples': {'p': p_samples},
                'timing': timing,
            }

        model_fn = self._numpyro_model_fn()
        model_for_hmc = partial(model_fn, sample_ypred=False,
                                compute_generated=True)
        rng = jax.random.PRNGKey(self.seed)
        kernel = NUTS(model_for_hmc)
        mcmc = MCMC(
            kernel,
            num_warmup=iter_warmup,
            num_samples=iter_sampling,
            num_chains=chains,
            progress_bar=False,
        )
        mcmc.run(rng, self.stan_data)
        samples = mcmc.get_samples(group_by_chain=True)
        p_samples_chained = np.asarray(samples['p'])  # (chains, draws)
        ds = xr.Dataset({'p': (('chain', 'draw'), p_samples_chained)})
        draws = az.InferenceData(posterior=ds)
        timing = {'mins': (time.time() - t0) / 60.0}
        if verbose:
            print(f"[BinomialModel] pyro HMC: {chains} chains x"
                  f" {iter_sampling} draws in {timing['mins']:.4f} min")
        if save_to_file and output_file_prefix is not None:
            os.makedirs(os.path.dirname(output_file_prefix) or '.',
                        exist_ok=True)
            draws.to_zarr(draws_file)
            pd.DataFrame([timing]).to_csv(timing_file, index=False)
        return {
            'mcmc': mcmc,
            'draws': draws,
            'posterior_samples': {'p': p_samples_chained.reshape(-1)},
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
        """Stan NUTS via cmdstanpy on the Beta-Bernoulli model.
        ``resume=True`` returns the cached ``{prefix}_draws.zarr`` if
        it exists."""
        from cmdstanpy import CmdStanModel
        t0 = time.time()
        draws_file = (f"{output_file_prefix}_draws.zarr"
                      if output_file_prefix is not None else None)
        timing_file = (f"{output_file_prefix}_timing.csv"
                       if output_file_prefix is not None else None)
        can_resume = (
            resume and save_to_file
            and draws_file is not None and os.path.exists(draws_file)
            and timing_file is not None and os.path.exists(timing_file)
        )
        if can_resume:
            if verbose:
                print(f"[BinomialModel] Stan HMC RESUME from {draws_file}")
            idata = az.from_zarr(draws_file)
            p_samples = np.asarray(idata.posterior['p'].values).reshape(-1)
            timing = pd.read_csv(timing_file).iloc[0].to_dict()
            return {
                'draws': idata,
                'posterior_samples': {'p': p_samples},
                'timing': timing,
            }
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
        if save_to_file and output_file_prefix is not None:
            os.makedirs(os.path.dirname(output_file_prefix) or '.',
                        exist_ok=True)
            idata.to_zarr(draws_file)
            pd.DataFrame([timing]).to_csv(timing_file, index=False)
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
