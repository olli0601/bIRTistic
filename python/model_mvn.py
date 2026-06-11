"""
:class:`MVNModel` -- multivariate Gaussian model with known correlation
structure R and known noise scale sigma. §3.3 of
``dev/amortised_decision_making.md``.

Generative model:
    K = R R^T            (known PSD, J x J)
    sigma > 0            (known)
    mu ~ MVN(prior_mu, prior_tau^2 * K)
    y_i | mu ~ MVN(mu, sigma^2 * K),     i = 1, ..., N

With the g-prior in the metric K and sigma known, the conjugate posterior
collapses to MVN and the per-component PPS reduces to a single Phi-tail per
(item, interim). This subclass supplies that closed-form benchmark plus a
NumPyro NUTS fitter so the nested-MC interim driver can be validated against
the analytic Phi-tail at J in {5, 20, 50, 100}.

Data layout. ``dit`` carries J rows, one per component j in {0, ..., J-1}
(``item_label = 'mu_<j>'`` etc). ``dcati`` is in long form with one row per
(pid, j): columns ``pid``, ``oid`` (1..N*J unique row id), ``j``,
``submission_date``, ``y``. This matches the Binomial-style ``oid`` indexing
so :class:`Model._load_ypred` and the nested-MC driver work unchanged.

Fitters: closed-form posterior, closed-form per-component PPS (Phi tail),
NumPyro NUTS. Stan + SVI raise NotImplementedError -- the analytic benchmark
makes them unnecessary for the MVN tests.
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
from numpyro.infer import MCMC, NUTS, Predictive
import pandas as pd
from scipy.stats import norm as _scipy_norm
import xarray as xr

from model import Model


_MODEL_NS = None


def _model_namespace():
    """Lazily load + cache the MVN numpyro model namespace."""
    global _MODEL_NS
    if _MODEL_NS is None:
        model_file = Path(__file__).resolve().parents[1] / 'src' / 'numpyro' / 'mvn_v260609.pyro'
        if not model_file.exists():
            raise FileNotFoundError(f"NumPyro model file not found: {model_file}")
        _MODEL_NS = runpy.run_path(str(model_file))
    return _MODEL_NS


# ----------------------------------------------------------------------
# R-matrix factories (cell A / B / C of §3.3.1)
# ----------------------------------------------------------------------


def make_ar1_chol(J: int, rho: float = 0.7) -> np.ndarray:
    """Cholesky factor R such that ``K = R R^T`` is the AR(1) correlation
    matrix with entries ``K_ij = rho^|i-j|`` (unit diagonal). Cell A of
    §3.3.1; well-conditioned for any rho in (-1, 1)."""
    idx = np.arange(J)
    K = rho ** np.abs(idx[:, None] - idx[None, :])
    return np.linalg.cholesky(K).astype(np.float64)


def make_block_equicorr_chol(J: int, block_size: int = 10,
                             rho_w: float = 0.8, rho_b: float = 0.1) -> np.ndarray:
    """Cell B: block-equicorrelation with unit diagonal. Cholesky factor
    of ``K`` where intra-block correlation is ``rho_w`` and inter-block is
    ``rho_b`` (PSD by construction when ``rho_w > rho_b >= 0``)."""
    if J % block_size != 0:
        raise ValueError(f"J ({J}) must be divisible by block_size ({block_size}).")
    B = J // block_size
    K = np.full((J, J), rho_b)
    for b in range(B):
        i0 = b * block_size
        i1 = i0 + block_size
        K[i0:i1, i0:i1] = rho_w
    np.fill_diagonal(K, 1.0)
    return np.linalg.cholesky(K).astype(np.float64)


def make_factor_chol(J: int, r: int = 5, psi: float = 0.1,
                     seed: int = 123) -> np.ndarray:
    """Cell C: low-rank-plus-diagonal factor structure. ``K = beta beta^T
    + psi * I``, rescaled to unit diagonal, then Cholesky."""
    rng = np.random.default_rng(seed)
    beta = rng.normal(scale=1.0 / np.sqrt(r), size=(J, r))
    K = beta @ beta.T + psi * np.eye(J)
    d = np.sqrt(np.diag(K))
    K = K / np.outer(d, d)
    return np.linalg.cholesky(K).astype(np.float64)


class MVNModel(Model):
    """Multivariate-normal model with known ``R`` (Cholesky of ``K``) and
    known noise scale ``sigma``. ``dit`` carries J item rows (one per
    component); ``dcati`` is in long form with one row per (pid, j) keyed
    by an integer 1-indexed ``oid``. ``prior_mu`` defaults to
    ``mu_0_baseline * 1_J`` and the g-prior covariance is
    ``prior_tau^2 * K``."""

    param_names = ('mu',)
    positive_params = ()

    def __init__(
        self,
        dit: pd.DataFrame,
        dcati: pd.DataFrame,
        x_formula: str = "~ 1",
        *,
        J: int,
        K_chol: np.ndarray,
        sigma: float = 1.0,
        prior_mu: Optional[np.ndarray] = None,
        prior_tau: float = 10.0,
        mu_0_baseline: float = 1.0,
        seed: int = 123,
    ):
        self.J = int(J)
        self.K_chol = np.asarray(K_chol, dtype=np.float64).copy()
        if self.K_chol.shape != (self.J, self.J):
            raise ValueError(
                f"K_chol must be ({J}, {J}); got {self.K_chol.shape}",
            )
        self.K = self.K_chol @ self.K_chol.T
        self.K_inv = np.linalg.inv(self.K)
        self.K_diag = np.diag(self.K).astype(np.float64).copy()
        self.sigma = float(sigma)
        self.sigma2 = float(sigma) ** 2
        self.mu_0_baseline = float(mu_0_baseline)
        if prior_mu is None:
            prior_mu = np.full(self.J, self.mu_0_baseline, dtype=np.float64)
        else:
            prior_mu = np.asarray(prior_mu, dtype=np.float64).reshape(-1)
            if prior_mu.shape != (self.J,):
                raise ValueError(
                    f"prior_mu must have shape ({J},); got {prior_mu.shape}",
                )
        self.prior_mu = prior_mu
        self.prior_tau = float(prior_tau)
        self.prior_tau2 = float(prior_tau) ** 2
        for col in ('pid', 'oid', 'j'):
            if col not in dcati.columns:
                raise RuntimeError(
                    f"MVNModel.dcati must carry a '{col}' column (long form).",
                )
        super().__init__(dit=dit, dcati=dcati, x_formula=x_formula, seed=seed)

    # ------------------------------------------------------------------
    # Data builder
    # ------------------------------------------------------------------

    @staticmethod
    def get_interim_data_x(df: pd.DataFrame,
                           interim_date=None) -> pd.DataFrame:
        """Slice ``df`` to ``submission_date <= interim_date`` and rebuild
        the long-form pid / oid / j indexing.

        Coalesces ``y_stan -> y`` (rows arriving from the nested-MC ``zi``
        block carry the sampled future outcomes in ``y_stan``), then
        re-indexes ``pid`` consecutively and resets ``oid`` to a
        contiguous 1..N*J row id."""
        dcati = df
        if interim_date is not None and 'submission_date' in dcati.columns:
            dcati = dcati[dcati['submission_date'] <= interim_date]
        if dcati.empty:
            return dcati.reset_index(drop=True)
        dcati = dcati.copy()
        if 'y_stan' in dcati.columns:
            if 'y' in dcati.columns:
                dcati['y'] = dcati['y'].where(
                    dcati['y'].notna(), dcati['y_stan'],
                )
            else:
                dcati['y'] = dcati['y_stan']
            dcati = dcati.drop(columns=['y_stan'])
        dcati['y'] = dcati['y'].astype(float)
        # Re-index pid consecutively (preserve relative order).
        pid_map = pd.DataFrame({'pid_orig': sorted(dcati['pid'].unique())})
        pid_map['pid_new'] = np.arange(1, len(pid_map) + 1, dtype=int)
        dcati = dcati.merge(pid_map, left_on='pid', right_on='pid_orig', how='left')
        dcati['pid'] = dcati['pid_new'].astype(int)
        dcati = dcati.drop(columns=['pid_orig', 'pid_new'])
        # Sort by (pid, j) lex so flat ypred index i*J + j matches.
        dcati = dcati.sort_values(['pid', 'j']).reset_index(drop=True)
        dcati['oid'] = np.arange(1, len(dcati) + 1, dtype=int)
        return dcati

    # ------------------------------------------------------------------
    # Stan / numpyro-data dict
    # ------------------------------------------------------------------

    def make_stan_data(self, dcati: pd.DataFrame, x_formula: str) -> dict:
        """Build the data dict consumed by the numpyro program. Pivots the
        long-form ``dcati`` to a (N, J) ``Y`` matrix sorted by (pid, j)."""
        Y = (
            dcati.pivot_table(index='pid', columns='j', values='y')
            .sort_index()
            .reindex(columns=range(self.J))
            .to_numpy()
            .astype(np.float64)
        )
        N = int(Y.shape[0])
        return {
            'N': N,
            'J': self.J,
            'Y': jnp.asarray(Y),
            'K_chol': jnp.asarray(self.K_chol),
            'sigma': self.sigma,
            'prior_mu': jnp.asarray(self.prior_mu),
            'prior_tau': self.prior_tau,
            'N_total': N * self.J,
        }

    # ------------------------------------------------------------------
    # Loglik / prior / endpoint
    # ------------------------------------------------------------------

    @staticmethod
    def eval_loglik(data, params):
        return _model_namespace()['get_log_likelihood_of_mvn'](data, params)

    @staticmethod
    def eval_loglik_annealed(data, params):
        return _model_namespace()['get_log_likelihood_of_mvn_with_annealing'](data, params)

    def eval_log_prior(self, params):
        return _model_namespace()['get_log_prior_of_mvn'](
            params,
            prior_mu=self.prior_mu,
            prior_tau=self.prior_tau,
            K_chol=self.K_chol,
        )

    def eval_outcome_for_endpoint(self, data, params):
        """Per-component endpoint mu_j - mu_0_baseline. H_1j: mu_j >
        mu_0_baseline  <=>  endpoint > 0."""
        return _model_namespace()['get_endpoint_of_mvn'](
            params, mu_0_baseline=self.mu_0_baseline,
        )

    # ------------------------------------------------------------------
    # Per-draw endpoint frame (one row per item j x draw)
    # ------------------------------------------------------------------

    def get_endpoints_per_draw(self, draws,
                               categorical_threshold: int = 3,
                               endpoint_type: str = 'items') -> pd.DataFrame:
        """One row per (item j, draw) with
        ``ratio = mu_j - mu_0_baseline`` and ``mu = mu_j``. H_1j: ratio > 0.
        """
        mu_draws = np.asarray(draws.posterior['mu'].values)
        mu_draws = mu_draws.reshape(-1, mu_draws.shape[-1])    # (K, J)
        K = mu_draws.shape[0]
        ratio = mu_draws - self.mu_0_baseline                  # (K, J)
        dit = self.dit.reset_index(drop=True)
        rows = []
        for j in range(self.J):
            item = dit.iloc[j]
            rows.append(pd.DataFrame({
                'j':                np.full(K, j, dtype=int),
                'item_label':       [item['item_label']] * K,
                'item_type':        [item['item_type']] * K,
                'item_high_label':  [item['item_high_label']] * K,
                'draw':             np.arange(K),
                'mu':               mu_draws[:, j],
                'ratio':            ratio[:, j],
            }))
        return pd.concat(rows, ignore_index=True)

    def get_w(self, zi: pd.DataFrame) -> pd.DataFrame:
        """Strong-Oakley per-(item j, draw) summary frame W(z^(s)).

        Per-component summary tied to the MVN exponential-family
        sufficient statistic ``T_1(z_{1:m}) = sum_i z_i``: the j-th
        component of ``T_1`` is ``sum_i z_{i,j}``, so taking the per-
        component mean ``w_j = (1/m) sum_i z_{i,j} = T_{1,j}/m``
        preserves all information about ``mu_j`` in the future cohort
        (see ``dev/amortised_decision_making.md`` §3.3 "Exponential family
        and sufficient statistics"). For each draw column ``ypred_s``
        (carrying the long-format flat outcomes for sample s), reshape
        to ``(m, J)`` by pivot on ``(pid, j)``. ``w_ratio = w_j -
        mu_0_baseline`` matches :meth:`get_endpoints_per_draw`'s
        ``ratio = mu_j - mu_0_baseline`` so the downstream Strong-Oakley
        regression of endpoint on w is on a common axis.
        """
        ypred_cols = [c for c in zi.columns if c.startswith('ypred_')]
        dit = self.dit.reset_index(drop=True)
        wide = (
            zi[['pid', 'j', *ypred_cols]]
            .melt(id_vars=['pid', 'j'], value_vars=ypred_cols,
                  var_name='s_col', value_name='val')
        )
        wide['draw'] = wide['s_col'].str.replace('ypred_', '').astype(int)
        means = (
            wide.groupby(['j', 'draw'])['val']
            .mean()
            .reset_index(name='w')
        )
        means = means.merge(
            dit[['item_label', 'item_type', 'item_high_label']]
               .reset_index().rename(columns={'index': 'j'}),
            on='j', how='left',
        )
        means['w_ratio'] = means['w'] - self.mu_0_baseline
        return means[['item_label', 'item_type', 'item_high_label',
                      'j', 'draw', 'w', 'w_ratio']]

    def get_endpoints(self, draws,
                      categorical_threshold: int = 3,
                      endpoint_type: str = 'items',
                      verbose: bool = False) -> pd.DataFrame:
        """Quantile summary q=[2.5, 25, 50, 75, 97.5]% on (mu, ratio) per
        item j."""
        per_draw = self.get_endpoints_per_draw(
            draws=draws,
            categorical_threshold=categorical_threshold,
            endpoint_type=endpoint_type,
        )
        quantiles = [0.025, 0.25, 0.5, 0.75, 0.975]
        names = ['q_lower', 'iqr_lower', 'median', 'iqr_upper', 'q_upper']
        rows = []
        for j, g in per_draw.groupby('j'):
            for var in ('mu', 'ratio'):
                q = g[var].quantile(quantiles).to_list()
                rows.append({
                    'item_label':      g.iloc[0]['item_label'],
                    'item_type':       g.iloc[0]['item_type'],
                    'item_high_label': g.iloc[0]['item_high_label'],
                    'j':               j,
                    'variable':        var,
                    **dict(zip(names, q)),
                })
        return pd.DataFrame(rows)

    # ------------------------------------------------------------------
    # Closed-form posterior + closed-form PPS (per-component Phi tail)
    # ------------------------------------------------------------------

    def _posterior_moments(self) -> Dict[str, np.ndarray]:
        """Closed-form NIG-with-sigma-known posterior moments. Returns
        ``mu_n`` (J,), ``cov_n`` = K / c_n (J, J), ``c_n``."""
        Y = np.asarray(self.stan_data['Y'])           # (N, J)
        N = int(Y.shape[0])
        c_n = 1.0 / self.prior_tau2 + N / self.sigma2
        mu_n = (
            self.prior_mu / self.prior_tau2
            + Y.sum(axis=0) / self.sigma2
        ) / c_n
        cov_n = self.K / c_n
        return {'mu_n': mu_n, 'cov_n': cov_n, 'c_n': float(c_n)}

    def fit_closed_form_posterior(
        self,
        output_file_prefix: Optional[str] = None, *,
        save_to_file: bool = False,
        verbose: bool = True,
        output_samples: int = 4000,
        **method_kwargs,
    ) -> Dict[str, Any]:
        """Sample directly from the analytic MVN posterior
        ``mu | x ~ MVN(mu_n, K / c_n)``. Also produces ``ypred`` (per-draw
        MVN(mu_s, sigma^2 K) at every obs) and analytic pointwise
        ``log_lik`` (per obs vector) so the returned zarr matches the
        numerical HMC driver's interface."""
        t0 = time.time()
        mom = self._posterior_moments()
        rng = np.random.default_rng(self.seed)
        # mu_s ~ MVN(mu_n, cov_n) drawn via Cholesky of cov_n = K/c_n.
        cov_chol = self.K_chol / np.sqrt(mom['c_n'])
        z = rng.standard_normal(size=(output_samples, self.J))
        mu_samples = mom['mu_n'][None, :] + z @ cov_chol.T   # (S, J)

        Y = np.asarray(self.stan_data['Y'])                  # (N, J)
        N = int(Y.shape[0])
        # ypred[s, i, :] = mu_s + sigma * R @ eps_i   (eps_i ~ N(0, I_J))
        eps = rng.standard_normal(size=(output_samples, N, self.J))
        ypred = (mu_samples[:, None, :]
                 + self.sigma * (eps @ self.K_chol.T))       # (S, N, J)
        ypred_flat = ypred.reshape(output_samples, N * self.J).astype(np.float64)

        # log_lik per obs vector: MVN(mu_s, sigma^2 K).log_prob(Y_i).
        log_det = 2.0 * float(np.sum(np.log(np.diag(self.K_chol))))
        L_obs = self.sigma * self.K_chol
        L_inv = np.linalg.solve(L_obs, np.eye(self.J))
        # quadratic form per (s, i): (Y_i - mu_s)^T (sigma^2 K)^{-1} (Y_i - mu_s)
        diff = Y[None, :, :] - mu_samples[:, None, :]        # (S, N, J)
        scaled = diff @ L_inv.T
        quad = (scaled ** 2).sum(axis=-1)                    # (S, N)
        norm_const = (
            self.J * np.log(2 * np.pi)
            + 2.0 * self.J * np.log(self.sigma)
            + log_det
        )
        log_lik = -0.5 * (quad + norm_const)

        ds = xr.Dataset({
            'mu':      (('chain', 'draw', 'j'), mu_samples[None, :, :]),
            'ypred':   (('chain', 'draw', 'obs'), ypred_flat[None, :, :]),
            'log_lik': (('chain', 'draw', 'obs_vec'), log_lik[None, :, :]),
        })
        draws = az.InferenceData(posterior=ds)
        timing = {'mins': (time.time() - t0) / 60.0}
        if verbose:
            print(f"[MVNModel] closed-form posterior: c_n={mom['c_n']:.3f};"
                  f" drew {output_samples} samples in"
                  f" {timing['mins']:.4f} min")
        if save_to_file and output_file_prefix is not None:
            os.makedirs(os.path.dirname(output_file_prefix) or '.',
                        exist_ok=True)
            draws.to_zarr(f"{output_file_prefix}_draws.zarr")
            pd.DataFrame([timing]).to_csv(
                f"{output_file_prefix}_timing.csv", index=False,
            )
        return {
            'draws': draws,
            'posterior_samples': {'mu': mu_samples},
            'timing': timing,
            'mu_n': mom['mu_n'],
            'cov_n': mom['cov_n'],
            'c_n': mom['c_n'],
        }

    def fit_closed_form_pH1(self, zi: pd.DataFrame) -> pd.DataFrame:
        """Analytic per-component ``P(H_{1j} | x, z_s)`` for each future-
        data sample ``s`` in ``zi`` (§3.3, lines 322-340 of
        ``dev/amortised_decision_making.md``).

        With sigma known and the g-prior, the posterior after the future
        cohort is Gaussian; for each sample ``s`` we read the per-component
        future cohort mean ``z_bar_s,j`` directly from ``ypred_s``,
        compute the updated posterior moments and return ``1 - Phi`` of
        the standardised distance to ``mu_0_baseline``.

        Updated posterior moments at ``s``:

            c_{n+m}      = 1/prior_tau^2 + (n + m)/sigma^2
            mu_{n+m,j}^s = (c_n * mu_n,j + (m/sigma^2) * z_bar_s,j) / c_{n+m}
            sd_{n+m,j}   = sqrt(K_jj / c_{n+m})

        and ``P(H_{1j} | x, z_s) = 1 - Phi((mu_0_baseline - mu_{n+m,j}^s)
        / sd_{n+m,j})``.

        Parameters
        ----------
        zi : pd.DataFrame
            Long-form future-data block from
            :meth:`get_interim_z_from_ypredi` carrying ``pid``, ``j`` and
            ``ypred_0 .. ypred_{S-1}`` columns.

        Returns
        -------
        pd.DataFrame
            Long form, one row per ``(item, sample)`` with columns
            ``item_label``, ``item_type``, ``item_high_label``, ``j``,
            ``s`` (1-indexed) and ``p_h1_xz``.
        """
        ypred_cols = [c for c in zi.columns if c.startswith('ypred_')]
        S = len(ypred_cols)
        if S == 0:
            return pd.DataFrame()

        zi_sorted = zi.sort_values(['pid', 'j']).reset_index(drop=True)
        m = int(zi_sorted['pid'].nunique())
        if m * self.J != len(zi_sorted):
            raise RuntimeError(
                f"zi has {len(zi_sorted)} rows but expected m*J = "
                f"{m}*{self.J} = {m * self.J}."
            )
        # (S, m, J) row-major over (s, pid, j).
        ypred_arr = (
            zi_sorted[ypred_cols].to_numpy().T.reshape(S, m, self.J)
        )
        z_bar = ypred_arr.mean(axis=1)                  # (S, J)

        mom = self._posterior_moments()
        N = int(self.stan_data['N'])
        c_n = mom['c_n']
        c_nm = 1.0 / self.prior_tau2 + (N + m) / self.sigma2
        mu_nm = (c_n * mom['mu_n'][None, :]
                 + (m / self.sigma2) * z_bar) / c_nm    # (S, J)
        sd = np.sqrt(self.K_diag / c_nm)                # (J,)
        p_h1_xz = 1.0 - _scipy_norm.cdf(
            (self.mu_0_baseline - mu_nm) / sd[None, :],
        )                                               # (S, J)

        dit = self.dit.reset_index(drop=True)
        s_idx = np.arange(1, S + 1)
        rows = pd.DataFrame({
            'item_label':       np.tile(dit['item_label'].to_numpy(), S),
            'item_type':        np.tile(dit['item_type'].to_numpy(), S),
            'item_high_label':  np.tile(dit['item_high_label'].to_numpy(), S),
            'j':                np.tile(np.arange(self.J), S),
            's':                np.repeat(s_idx, self.J),
            'p_h1_xz':          p_h1_xz.reshape(-1),
        })
        return rows

    def fit_closed_form_pps(
        self,
        m: int,
        *,
        pps_H1_def: float = 0.0,
        pps_ProbH1_thresh: float = 0.89,
    ) -> pd.DataFrame:
        """Analytic per-component Phi-tail PPS (§3.3.1).

        For each component j, with c_t = 1/prior_tau^2 + t/sigma^2 and
        K_jj on the diagonal (we report unit-diag K, so K_jj = 1):

            PPS_j(x) = Phi( sigma * sqrt(c_n * c_{n+m} / m) *
                            (mu_n,j - 1 - z_eta / sqrt(c_{n+m})) )

        where z_eta = Phi^{-1}(pps_ProbH1_thresh). Returns a DataFrame
        with one row per item j carrying ``pps``, ``mu_n``, ``c_n``,
        ``c_n_plus_m``, plus the dit metadata columns.

        ``pps_H1_def`` is accepted for interface parity but the H_1j
        threshold lives in ``mu_0_baseline`` (set on the model)."""
        if m < 0:
            raise ValueError("m must be non-negative.")
        if m == 0:
            return pd.DataFrame()
        mom = self._posterior_moments()
        N = int(self.stan_data['N'])
        c_n = mom['c_n']
        c_np = 1.0 / self.prior_tau2 + (N + m) / self.sigma2
        z_eta = float(_scipy_norm.ppf(pps_ProbH1_thresh))
        K_diag = self.K_diag
        # Decision rule: P(H_{1j} | x, z) > eta  <=>
        #     mu_{n+m,j} > mu_0_baseline + z_eta * sqrt(K_jj / c_{n+m}).
        # Marginal predictive of mu_{n+m,j} under z | x:
        #     E[mu_{n+m,j} | x] = mu_{n,j}
        #     Var[mu_{n+m,j} | x] = K_jj * m / (sigma^2 * c_n * c_{n+m}).
        # PPS_j = Phi( (mu_{n,j} - threshold) / sqrt(var_mu_{n+m,j} | x) ).
        threshold = self.mu_0_baseline + z_eta * np.sqrt(K_diag / c_np)
        var_pred = K_diag * m / (self.sigma2 * c_n * c_np)
        z = (mom['mu_n'] - threshold) / np.sqrt(var_pred)
        pps = _scipy_norm.cdf(z)
        dit = self.dit.reset_index(drop=True)
        return pd.DataFrame({
            'j':                np.arange(self.J),
            'item_label':       dit['item_label'].to_numpy(),
            'item_type':        dit['item_type'].to_numpy(),
            'item_high_label':  dit['item_high_label'].to_numpy(),
            'm':                int(m),
            'pps':              pps.astype(np.float64),
            'eta':              float(pps_ProbH1_thresh),
            'mu_n':             mom['mu_n'].astype(np.float64),
            'c_n':              float(c_n),
            'c_n_plus_m':       float(c_np),
        })

    # ------------------------------------------------------------------
    # NumPyro NUTS
    # ------------------------------------------------------------------

    def _numpyro_model_fn(self):
        return _model_namespace()['mvn_model']

    def _pyro_predictive(self, posterior_samples_dict, rng_key, N):
        """Predictive to attach ``ypred`` + ``log_lik`` per posterior draw.
        Returns ``{'ypred': (n_draw, N*J), 'log_lik': (n_draw, N)}``."""
        pred_fn = partial(self._numpyro_model_fn(), compute_generated=True)
        predictive = Predictive(
            pred_fn,
            posterior_samples=posterior_samples_dict,
            return_sites=['log_lik', 'ypred'],
        )
        predictions = predictive(rng_key, self.stan_data)
        return {
            'ypred':   np.asarray(predictions['ypred']).astype(np.float64),
            'log_lik': np.asarray(predictions['log_lik']).astype(np.float64),
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
        """NumPyro NUTS on the MVN program. ``compute_generated`` is False
        during inference; ``ypred`` + ``log_lik`` come from
        :meth:`_pyro_predictive` post-fit."""
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
                print(f"[MVNModel] pyro HMC RESUME from {draws_file}")
            idata = az.from_zarr(draws_file)
            mu_samples = np.asarray(idata.posterior['mu'].values).reshape(
                -1, self.J,
            )
            timing = pd.read_csv(timing_file).iloc[0].to_dict()
            return {
                'draws': idata,
                'posterior_samples': {'mu': mu_samples},
                'timing': timing,
            }

        model_fn = self._numpyro_model_fn()
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
        mu_chained = np.asarray(samples['mu'])      # (chains, draws, J)
        rng, sub = jax.random.split(rng)
        N = int(self.stan_data['N'])
        predictions = self._pyro_predictive(
            {'mu': mu_chained.reshape(-1, self.J)}, sub, N,
        )
        n_chain, n_draw, _ = mu_chained.shape
        ds = xr.Dataset({
            'mu':      (('chain', 'draw', 'j'), mu_chained),
            'ypred':   (('chain', 'draw', 'obs'),
                        predictions['ypred'].reshape(n_chain, n_draw, N * self.J)),
            'log_lik': (('chain', 'draw', 'obs_vec'),
                        predictions['log_lik'].reshape(n_chain, n_draw, N)),
        })
        draws = az.InferenceData(posterior=ds)
        timing = {'mins': (time.time() - t0) / 60.0}
        if verbose:
            print(f"[MVNModel] pyro HMC: {chains} chains x {iter_sampling}"
                  f" draws (J={self.J}, N={N}) in {timing['mins']:.4f} min")
        if save_to_file and output_file_prefix is not None:
            os.makedirs(os.path.dirname(output_file_prefix) or '.',
                        exist_ok=True)
            draws.to_zarr(draws_file)
            pd.DataFrame([timing]).to_csv(timing_file, index=False)
        return {
            'mcmc': mcmc,
            'draws': draws,
            'posterior_samples': {'mu': mu_chained.reshape(-1, self.J)},
            'timing': timing,
        }

    # ------------------------------------------------------------------
    # Abstract Model overrides for fitters we don't implement.
    # ------------------------------------------------------------------

    def fit_pyro_svi(self, output_file_prefix, *, save_to_file=True,
                     resume=False, verbose=True, **method_kwargs):
        raise NotImplementedError(
            "MVNModel: fit_pyro_svi not implemented (HMC + closed-form cover "
            "the validation surface used in the MVN benchmark)."
        )

    def fit_stan_svi(self, output_file_prefix, *, save_to_file=True,
                     resume=False, verbose=True, **method_kwargs):
        raise NotImplementedError(
            "MVNModel: fit_stan_svi not implemented."
        )

    def fit_stan_hmc(self, output_file_prefix, *, save_to_file=True,
                     resume=False, verbose=True, **method_kwargs):
        raise NotImplementedError(
            "MVNModel: fit_stan_hmc not implemented (use fit_pyro_hmc)."
        )
