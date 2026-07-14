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
from numpyro.infer import MCMC, NUTS, Predictive, SVI, Trace_ELBO
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
        if 'oid' not in dcati.columns:
            raise RuntimeError(
                "BinomialModel.dcati must carry an 'oid' column "
                "(1-indexed observation id)."
            )
        super().__init__(dit=dit, dcati=dcati, x_formula=x_formula, seed=seed)

    # ---- data builder ---------------------------------------------------
    @staticmethod
    def get_interim_data_x(df: pd.DataFrame, interim_date=None) -> pd.DataFrame:
        """Trivial Binomial data builder. Returns ``df`` optionally
        sliced to ``submission_date <= interim_date`` with ``oid``
        reset to a contiguous 1..N index (matches the IRT cohort
        builder's invariant).

        If a ``y_stan`` column is present (rows arriving from the
        nested-MC ``zi`` block carry the sampled future outcomes here),
        it is coalesced into ``y`` so the augmented ``pd.concat([xi,
        zi])`` cohort (where xi has ``y`` only and zi has ``y_stan``
        only, with NaN in the missing column per row) feeds
        :meth:`make_stan_data` uniformly."""
        dcati = df
        if interim_date is not None and 'submission_date' in dcati.columns:
            dcati = dcati[dcati['submission_date'] <= interim_date]
        dcati = dcati.reset_index(drop=True)
        if 'y_stan' in dcati.columns:
            dcati = dcati.copy()
            if 'y' in dcati.columns:
                dcati['y'] = dcati['y'].where(
                    dcati['y'].notna(), dcati['y_stan'],
                )
            else:
                dcati['y'] = dcati['y_stan']
            dcati = dcati.drop(columns=['y_stan'])
        dcati['y'] = dcati['y'].astype(int)
        dcati['oid'] = np.arange(1, len(dcati) + 1, dtype=int)
        return dcati

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
        ratio = 1 - p_draws / self.p_0
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

    def get_interim_z_from_ypredi(
        self, draws_file: str, interim_m: int,
        pps_z_total: int = 10, seed: int = 123,
        keep_order: bool = False,
    ) -> pd.DataFrame:
        """Build ``zi`` for the Binomial model.

        Overrides the base class default so that ``m`` future outcomes
        are sampled as **fresh** Bernoulli(``p_s``) draws per posterior
        column, not by resampling from the interim cohort's ``n``
        observed ypred values. When ``m > n`` (small interim, big future
        cohort) the base-class resample-with-replacement scheme inflates
        the per-draw ``km`` variance by ``m^2 / n * E[p(1-p)]`` on top
        of the true Beta-Binomial variance; drawing fresh Bernoullis
        removes that inflation and matches the exact posterior-
        predictive marginal.

        Reads ``p`` samples from ``draws_file``'s
        ``posterior['p']`` and generates ``ypred`` at ``interim_m``
        new participants. ``keep_order`` and ``pps_z_total`` follow the
        base-class contract.
        """
        idata = az.from_zarr(draws_file)
        p_all = np.asarray(idata.posterior['p'].values).reshape(-1)  # (S_total,)
        S_total = p_all.shape[0]
        if pps_z_total > S_total:
            raise ValueError(
                f"pps_z_total={pps_z_total} > cached S={S_total}",
            )
        rng = np.random.default_rng(seed)
        if keep_order:
            sel = np.arange(pps_z_total)
        else:
            sel = rng.choice(S_total, size=pps_z_total, replace=False)
        p_sel = p_all[sel]                                          # (S,)
        # Fresh Bernoulli(p_s) at each of the interim_m future positions.
        ypred = rng.binomial(
            1, np.broadcast_to(p_sel[:, None], (pps_z_total, interim_m)),
        ).astype(np.int8).T                                         # (m, S)

        xi = self.dcati
        src_pids = rng.choice(xi['pid'].unique(), size=interim_m, replace=True)
        tmp = pd.DataFrame({
            'src_pid': src_pids,
            'pid': np.arange(1, interim_m + 1) + int(xi['pid'].max()),
        })
        zi = tmp.merge(xi.rename(columns={'pid': 'src_pid'}), on='src_pid')
        zi.sort_values(['pid', 'oid'], inplace=True)
        zi.reset_index(drop=True, inplace=True)
        ypred_df = pd.DataFrame(
            ypred, columns=[f'ypred_{s}' for s in range(pps_z_total)],
        )
        zi = pd.concat([zi.reset_index(drop=True), ypred_df], axis=1)
        return zi

    def get_w(self, zi: pd.DataFrame) -> pd.DataFrame:
        """Strong-Oakley per-draw summary W(z^(s)) for the single Binomial
        item.

        For each draw column ``ypred_s`` in ``zi``:
        - ``w`` = mean # successes per row of zi  (i.e. empirical Bernoulli
          rate of the m future trials for that draw),
        - ``w_ratio`` = ``1 - w / p_0`` (directional endpoint, same
          convention as :meth:`get_endpoints_per_draw`)."""
        ypred_cols = [c for c in zi.columns if c.startswith('ypred_')]
        item_row = self.dit.iloc[0]
        w_means = zi[ypred_cols].mean(axis=0).to_numpy()
        S = len(ypred_cols)
        df = pd.DataFrame({
            'item_label':       [item_row['item_label']] * S,
            'item_type':        [item_row['item_type']] * S,
            'item_high_label':  [item_row['item_high_label']] * S,
            'draw':             np.arange(S),
            'w':                w_means,
        })
        df['w_ratio'] = 1.0 - df['w'] / self.p_0
        return df

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
        resume: bool = False,
        verbose: bool = True,
        output_samples: int = 4000,
        **method_kwargs,
    ) -> Dict[str, Any]:
        """Sample directly from the analytic Beta(a + k, b + n - k) posterior.

        Also produces ``ypred`` (per-draw Bernoulli(p_s) at every
        observation) and analytic pointwise ``log_lik`` so the returned
        zarr matches the four numerical fit drivers.

        When ``save_to_file`` and ``resume`` are True and the
        ``{prefix}_draws.zarr`` + ``{prefix}_timing.csv`` artifacts
        exist, the cached posterior is returned instead of resampling.
        Matches the save/resume idiom of :meth:`fit_pyro_hmc` /
        :meth:`fit_pyro_svi` so downstream scripts can swap this in as
        a drop-in replacement.
        """
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
                print(f"[BinomialModel] closed-form posterior RESUME from"
                      f" {draws_file}")
            idata = az.from_zarr(draws_file)
            p_samples = np.asarray(idata.posterior['p'].values).reshape(-1)
            timing = pd.read_csv(timing_file).iloc[0].to_dict()
            k = int(self.stan_data['k'])
            n = int(self.stan_data['N_total'])
            return {
                'draws': idata,
                'posterior_samples': {'p': p_samples},
                'timing': timing,
                'a_post': self.prior_a + k,
                'b_post': self.prior_b + (n - k),
            }

        k = int(self.stan_data['k'])
        n = int(self.stan_data['N_total'])
        a_post = self.prior_a + k
        b_post = self.prior_b + (n - k)
        rng = np.random.default_rng(self.seed)
        p_samples = rng.beta(a_post, b_post, size=output_samples)

        y = np.asarray(self.stan_data['y'], dtype=np.int8)         # (N,)
        # ypred[s, i] ~ Bernoulli(p_s); log_lik[s, i] = y_i log p_s + (1-y_i) log(1-p_s).
        ypred = rng.binomial(
            1, np.broadcast_to(p_samples[:, None], (output_samples, n)),
        ).astype(np.int8)                                          # (S, N)
        log_p = np.log(p_samples)[:, None]                         # (S, 1)
        log_1mp = np.log1p(-p_samples)[:, None]
        log_lik = y[None, :] * log_p + (1 - y[None, :]) * log_1mp  # (S, N)

        ds = xr.Dataset({
            'p':       (('chain', 'draw'), p_samples[None, :]),
            'ypred':   (('chain', 'draw', 'obs'), ypred[None, :, :]),
            'log_lik': (('chain', 'draw', 'obs'), log_lik[None, :, :]),
        })
        draws = az.InferenceData(posterior=ds)
        timing = {'mins': (time.time() - t0) / 60.0}
        if verbose:
            print(f"[BinomialModel] closed-form posterior: Beta({a_post:.2f},"
                  f" {b_post:.2f}); drew {output_samples} samples in"
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
            'a_post': a_post,
            'b_post': b_post,
        }

    def fit_closed_form_pps(
        self,
        m: int,
        *,
        pps_H1_min_effect_size_thresh: float = 0.5,
        pps_ProbH1_target_lwr_quantile: float = 0.89,
    ) -> float:
        """Analytic Beta-Binomial PPS over ``m`` future trials (§3.1).

        H_1 is defined on the endpoint ``1 - p / p_0``:
        ``H_1: 1 - p / p_0 > pps_H1_min_effect_size_thresh``, equivalently
        ``p < (1 - pps_H1_min_effect_size_thresh) * p_0``.

        ``pps_ProbH1_target_lwr_quantile`` is the decision threshold on
        ``P(H_1 | x, z)``. Returns
        ``P( P(H_1 | x, z) > pps_ProbH1_target_lwr_quantile )`` integrating zi
        analytically over the Beta-Binomial(m, a_post, b_post)
        posterior-predictive.
        """
        if m < 0:
            raise ValueError("m must be non-negative.")
        k = int(self.stan_data['k'])
        n = int(self.stan_data['N_total'])
        a_post = self.prior_a + k
        b_post = self.prior_b + (n - k)

        # H_1: p < p_h1_threshold; P(H_1 | x, z) is decreasing in k_m
        # (more successes shift the posterior right, away from the threshold).
        p_h1_threshold = (1.0 - pps_H1_min_effect_size_thresh) * self.p_0
        k_m_grid = np.arange(m + 1)
        a_z = a_post + k_m_grid
        b_z = b_post + (m - k_m_grid)
        p_h1_given_z = _scipy_beta.cdf(p_h1_threshold, a_z, b_z)
        crosses = p_h1_given_z > pps_ProbH1_target_lwr_quantile
        if not crosses.any():
            return 0.0
        # Largest k_m where the decision still triggers; PPS = P(K_m <= k_star).
        k_star = int(k_m_grid[crosses][-1])
        return float(_scipy_betabinom.cdf(k_star, m, a_post, b_post))

    # ------------------------------------------------------------------
    # Amortiser training-data sampler
    # ------------------------------------------------------------------

    def _sample_binomial_prior_predictive(
        self, rng, S, n_max, n_min, fixed_total, m_min, m_max,
        n_dist='uniform',
    ):
        """Shared helper: draw ``(p, n, m, kn, km)`` per sample.

        ``n_dist='uniform'`` (default) samples ``n`` uniformly on
        ``{n_min, ..., n_max - 1}``. ``n_dist='log_uniform'`` samples
        ``n`` log-uniformly over the same range to oversample the small-
        ``n`` regime.
        """
        p = rng.beta(self.prior_a, self.prior_b, size=S)
        if n_dist == 'uniform':
            n = rng.integers(n_min, n_max, size=S)
        elif n_dist == 'log_uniform':
            log_lo = float(np.log(max(n_min, 1)))
            log_hi = float(np.log(n_max))
            n = np.exp(rng.uniform(log_lo, log_hi, size=S)).astype(np.int64)
            n = np.clip(n, n_min, n_max - 1)
        else:
            raise ValueError(f"unknown n_dist={n_dist!r}")
        if fixed_total:
            m = n_max - n
        else:
            if m_max is None:
                m_max = n_max
            m = rng.integers(m_min, m_max + 1, size=S)
        kn = rng.binomial(n, p)
        km = rng.binomial(m, p)
        return p, n, m, kn, km

    def make_training_data_with_features(
        self,
        rng: np.random.Generator,
        S: int,
        *,
        n_max: int,
        n_min: int = 1,
        fixed_total: bool = True,
        m_min: int = 0,
        m_max: int = None,
        n_dist: str = 'uniform',
    ):
        """Prior-predictive joint sampler for the features-fixed amortiser
        (§8 / §9.3).

        DeepSets pooling for hardcoded ``T(x_i) = x_i`` collapses to
        ``sum_T_x + sum_T_z = kn + km = k_total``. Features passed to
        the head are ``(k_total, n_total) / n_max`` (both normalised to
        ``[0, 1]``); this two-scalar input matches the true Binomial
        sufficient statistic (§3.1).

        Returns
        -------
        (batch, rho) : ``(dict, np.ndarray)``
            - ``batch['features']``: shape ``(S, 2)``, dtype ``float32``.
            - ``rho``: shape ``(S,)``, dtype ``float32``.
        """
        p, n, m, kn, km = self._sample_binomial_prior_predictive(
            rng, S, n_max, n_min, fixed_total, m_min, m_max, n_dist,
        )
        k_total = kn + km
        n_total = n + m
        features = np.stack(
            [k_total, n_total], axis=-1,
        ).astype(np.float32) / n_max
        rho = (1.0 - p / self.p_0).astype(np.float32)
        return {'features': features}, rho

    def make_training_data_with_raw_sequences(
        self,
        rng: np.random.Generator,
        S: int,
        *,
        n_max: int,
        n_min: int = 1,
        fixed_total: bool = True,
        m_min: int = 0,
        m_max: int = None,
        n_dist: str = 'uniform',
    ):
        """Prior-predictive joint sampler for the features-MLP amortiser.

        Emits raw padded item sequences instead of pooled sufficient
        statistics. Suitable for ``q_tau`` MLP encoders that must see
        individual ``x_i`` observations. Padding is right-aligned; the
        mask carries ``1`` for real items and ``0`` for padded slots.

        Returns
        -------
        (batch, rho) : ``(dict, np.ndarray)``
            - ``batch['x']``: shape ``(S, n_max, 1)``, dtype ``float32``,
              padded Bernoulli current-data items.
            - ``batch['mask_x']``: shape ``(S, n_max)``, dtype
              ``float32``, ``1`` for real, ``0`` for padded.
            - ``batch['z']``: shape ``(S, n_max, 1)``, dtype ``float32``,
              padded Bernoulli future-data items.
            - ``batch['mask_z']``: shape ``(S, n_max)``, dtype
              ``float32``.
            - ``batch['sizes']``: shape ``(S, 2)``, dtype ``float32``,
              cohort sizes ``(n / n_max, m / n_max)``.
            - ``rho``: shape ``(S,)``, dtype ``float32``.
        """
        p, n, m, kn, km = self._sample_binomial_prior_predictive(
            rng, S, n_max, n_min, fixed_total, m_min, m_max, n_dist,
        )
        # Build padded x. Slot 0..kn-1 = 1, slot kn..n-1 = 0, slot n..n_max-1 = padding.
        x = np.zeros((S, n_max, 1), dtype=np.float32)
        z = np.zeros((S, n_max, 1), dtype=np.float32)
        mask_x = np.zeros((S, n_max), dtype=np.float32)
        mask_z = np.zeros((S, n_max), dtype=np.float32)
        col_idx = np.arange(n_max)[None, :]                       # (1, n_max)
        # For x: 1 in slots [0, kn), 0 in [kn, n), padded in [n, n_max).
        x[..., 0] = (col_idx < kn[:, None]).astype(np.float32)
        mask_x[:] = (col_idx < n[:, None]).astype(np.float32)
        z[..., 0] = (col_idx < km[:, None]).astype(np.float32)
        mask_z[:] = (col_idx < m[:, None]).astype(np.float32)
        sizes = np.stack(
            [n / n_max, m / n_max], axis=-1,
        ).astype(np.float32)
        rho = (1.0 - p / self.p_0).astype(np.float32)
        return (
            {
                'x': x,
                'mask_x': mask_x,
                'z': z,
                'mask_z': mask_z,
                'sizes': sizes,
            },
            rho,
        )

    # ------------------------------------------------------------------
    # NumPyro SVI + HMC
    # ------------------------------------------------------------------

    def _numpyro_model_fn(self):
        return _model_namespace()['bernoulli_model']

    def _pyro_predictive(self, posterior_samples_dict, rng_key):
        """Run numpyro :class:`Predictive` to attach ``ypred`` + ``log_lik``
        per posterior draw. Returns ``{'ypred': (n, N_total),
        'log_lik': (n, N_total)}`` as plain numpy arrays."""
        pred_fn = partial(
            self._numpyro_model_fn(), compute_generated=True,
        )
        predictive = Predictive(
            pred_fn,
            posterior_samples=posterior_samples_dict,
            return_sites=['log_lik', 'ypred'],
        )
        predictions = predictive(rng_key, self.stan_data)
        return {
            'ypred': np.asarray(predictions['ypred']).astype(np.int8),
            'log_lik': np.asarray(predictions['log_lik']),
        }

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
        model_for_svi = partial(model_fn, compute_generated=False)

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
        rng, sub = jax.random.split(rng)
        predictions = self._pyro_predictive({'p': p_samples}, sub)
        N_total = int(self.stan_data['N_total'])
        n_draw = p_samples.shape[0]
        ds = xr.Dataset({
            'p':       (('chain', 'draw'), p_samples[None, :]),
            'ypred':   (('chain', 'draw', 'obs'),
                        predictions['ypred'].reshape(1, n_draw, N_total)),
            'log_lik': (('chain', 'draw', 'obs'),
                        predictions['log_lik'].reshape(1, n_draw, N_total)),
        })
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
        """NumPyro NUTS fit on the Beta-Bernoulli program. NUTS cannot
        enumerate the discrete ypred site, so ``compute_generated`` is
        forced to False during inference; ``ypred`` + ``log_lik`` are
        produced post-fit via :meth:`_pyro_predictive`. ``resume=True``
        returns the cached ``{prefix}_draws.zarr`` if it exists."""
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
        # NUTS cannot enumerate the discrete ypred site; disable generated
        # quantities during inference. ypred + log_lik come from Predictive below.
        model_for_hmc = partial(model_fn, compute_generated=False)
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
        rng, sub = jax.random.split(rng)
        predictions = self._pyro_predictive(
            {'p': p_samples_chained.reshape(-1)}, sub,
        )
        n_chain, n_draw = p_samples_chained.shape
        N_total = int(self.stan_data['N_total'])
        ds = xr.Dataset({
            'p':       (('chain', 'draw'), p_samples_chained),
            'ypred':   (('chain', 'draw', 'obs'),
                        predictions['ypred'].reshape(n_chain, n_draw, N_total)),
            'log_lik': (('chain', 'draw', 'obs'),
                        predictions['log_lik'].reshape(n_chain, n_draw, N_total)),
        })
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
        raw_idata = az.from_cmdstanpy(fit, posterior_predictive='y_rep')
        ds = xr.Dataset({
            'p':       (('chain', 'draw'), raw_idata.posterior['p'].values),
            'ypred':   (('chain', 'draw', 'obs'),
                        raw_idata.posterior_predictive['y_rep'].values.astype(np.int8)),
            'log_lik': (('chain', 'draw', 'obs'),
                        raw_idata.log_likelihood['log_lik'].values),
        })
        idata = az.InferenceData(posterior=ds)
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
        """Stan ADVI via cmdstanpy on the Beta-Bernoulli model. ``ypred`` +
        ``log_lik`` are attached post-fit by parsing the cmdstanpy
        variational sample (Stan's ``generated quantities``)."""
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
        df = fit.variational_sample_pd
        p_samples = np.asarray(df['p']).reshape(-1)
        N_total = int(self.stan_data['N_total'])
        yrep_cols = [c for c in df.columns if c.startswith('y_rep')]
        ll_cols = [c for c in df.columns if c.startswith('log_lik')]

        def _key(c):
            return int(c.split('.', 1)[1]) if '.' in c else int(
                c.split('[', 1)[1].rstrip(']')
            )

        yrep_cols.sort(key=_key)
        ll_cols.sort(key=_key)
        if len(yrep_cols) != N_total or len(ll_cols) != N_total:
            raise RuntimeError(
                f"[fit_stan_svi] expected {N_total} y_rep / log_lik columns, "
                f"got {len(yrep_cols)} / {len(ll_cols)}"
            )
        ypred = df[yrep_cols].to_numpy().astype(np.int8)   # (n_draw, N)
        log_lik = df[ll_cols].to_numpy()
        ds = xr.Dataset({
            'p':       (('chain', 'draw'), p_samples[None, :]),
            'ypred':   (('chain', 'draw', 'obs'), ypred[None, :, :]),
            'log_lik': (('chain', 'draw', 'obs'), log_lik[None, :, :]),
        })
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
