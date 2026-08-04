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
from numpyro.infer import MCMC, NUTS, Predictive, SVI, Trace_ELBO
import pandas as pd
from scipy.stats import norm as _scipy_norm
import xarray as xr

from model import Model
from utils import _get_autoguide_factory


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


def sample_random_K_chol(rng, J, families=None):
    """Draw a random K_chol from a mixture of standard families.

    Used at TRAINING for the K-family-invariant amortiser (§12.11). Each
    call returns a fresh (J, J) Cholesky factor for a K matrix drawn
    from one of ``families``:

      - ``'identity'``           : K = I.
      - ``'ar1'``                : K_{ij} = rho^|i-j|, rho ~ U(-0.9, 0.9).
      - ``'block'``              : block-equicorrelation, block size
                                    drawn from divisors of J in {5, 10,
                                    20, 25, 50}, rho_w ~ U(0.3, 0.9),
                                    rho_b ~ U(0, min(0.5, rho_w - 0.1)).
      - ``'factor'``             : factor structure with r ~ {2, 5, 10},
                                    psi ~ U(0.05, 0.5).

    Default families: all four. All returned K have unit diagonal so
    per-component posterior variance is invariant across families.
    """
    families = families or ('identity', 'ar1', 'block', 'factor')
    fam = families[int(rng.integers(0, len(families)))]
    if fam == 'identity':
        K = np.eye(J, dtype=np.float64)
    elif fam == 'ar1':
        rho = float(rng.uniform(-0.9, 0.9))
        idx = np.arange(J)
        K = (rho ** np.abs(idx[:, None] - idx[None, :])).astype(np.float64)
    elif fam == 'block':
        candidates = [b for b in (5, 10, 20, 25, 50) if J % b == 0]
        if not candidates:
            candidates = [1]
        block_size = int(candidates[int(rng.integers(0, len(candidates)))])
        rho_w = float(rng.uniform(0.3, 0.9))
        rho_b = float(rng.uniform(0, max(0, min(0.5, rho_w - 0.1))))
        B = J // block_size
        K = np.full((J, J), rho_b, dtype=np.float64)
        for b in range(B):
            i0 = b * block_size
            i1 = i0 + block_size
            K[i0:i1, i0:i1] = rho_w
        np.fill_diagonal(K, 1.0)
    elif fam == 'factor':
        r = int(rng.choice([2, 5, 10]))
        psi = float(rng.uniform(0.05, 0.5))
        beta = rng.normal(scale=1.0 / np.sqrt(r), size=(J, r))
        K = beta @ beta.T + psi * np.eye(J)
        d = np.sqrt(np.diag(K))
        K = K / np.outer(d, d)
    else:
        raise ValueError(f"unknown K family {fam!r}")
    return np.linalg.cholesky(K).astype(np.float64), fam


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
        resume: bool = False,
        verbose: bool = True,
        output_samples: int = 4000,
        **method_kwargs,
    ) -> Dict[str, Any]:
        """Sample directly from the analytic MVN posterior
        ``mu | x ~ MVN(mu_n, K / c_n)``. Also produces ``ypred`` (per-draw
        MVN(mu_s, sigma^2 K) at every obs) and analytic pointwise
        ``log_lik`` (per obs vector) so the returned zarr matches the
        numerical HMC driver's interface.

        When ``save_to_file`` and ``resume`` are True and the
        ``{prefix}_draws.zarr`` + ``{prefix}_timing.csv`` artifacts
        exist, the cached posterior is returned instead of resampling.
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
                print(f"[MVNModel] closed-form posterior RESUME from"
                      f" {draws_file}")
            idata = az.from_zarr(draws_file)
            mu_samples = np.asarray(
                idata.posterior['mu'].values,
            ).reshape(-1, self.J)
            timing = pd.read_csv(timing_file).iloc[0].to_dict()
            mom_cached = self._posterior_moments()
            return {
                'draws': idata,
                'posterior_samples': {'mu': mu_samples},
                'timing': timing,
                'mu_n': mom_cached['mu_n'],
                'cov_n': mom_cached['cov_n'],
                'c_n': mom_cached['c_n'],
            }
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
            draws.to_zarr(draws_file)
            pd.DataFrame([timing]).to_csv(timing_file, index=False)
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
        pps_H1_min_effect_size_thresh: float = 0.0,
        pps_ProbH1_target_lwr_quantile: float = 0.89,
    ) -> pd.DataFrame:
        """Analytic per-component Phi-tail PPS (§3.3.1).

        For each component j, with c_t = 1/prior_tau^2 + t/sigma^2 and
        K_jj on the diagonal (we report unit-diag K, so K_jj = 1):

            PPS_j(x) = Phi( sigma * sqrt(c_n * c_{n+m} / m) *
                            (mu_n,j - 1 - z_eta / sqrt(c_{n+m})) )

        where z_eta = Phi^{-1}(pps_ProbH1_target_lwr_quantile). Returns a DataFrame
        with one row per item j carrying ``pps``, ``mu_n``, ``c_n``,
        ``c_n_plus_m``, plus the dit metadata columns.

        ``pps_H1_min_effect_size_thresh`` is accepted for interface parity but the H_1j
        threshold lives in ``mu_0_baseline`` (set on the model)."""
        if m < 0:
            raise ValueError("m must be non-negative.")
        if m == 0:
            return pd.DataFrame()
        mom = self._posterior_moments()
        N = int(self.stan_data['N'])
        c_n = mom['c_n']
        c_np = 1.0 / self.prior_tau2 + (N + m) / self.sigma2
        z_eta = float(_scipy_norm.ppf(pps_ProbH1_target_lwr_quantile))
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
            'eta':              float(pps_ProbH1_target_lwr_quantile),
            'mu_n':             mom['mu_n'].astype(np.float64),
            'c_n':              float(c_n),
            'c_n_plus_m':       float(c_np),
        })

    # ------------------------------------------------------------------
    # Amortiser training-data samplers
    # ------------------------------------------------------------------

    def _sample_mvn_prior_predictive_summary(
        self, rng, S, n_max, n_min, fixed_total, m_min, m_max,
        n_dist='uniform',
    ):
        """Fast per-component prior-predictive summary sampler.

        Draws ``S`` per-component training rows and returns the
        sufficient-statistic summary ``(mu, n, m, sum_y_x, sum_y_z)``
        WITHOUT materialising the raw ``(N_max, J)`` observation tensor:

          - Under the g-prior with unit-diagonal ``K``, the per-component
            distribution of ``y_{i, j} | mu`` is
            ``Normal(mu_j, sigma^2 * K_{jj})``, independent across i and
            j given ``mu``.
          - So ``sum_i y_{i, j} | mu`` is
            ``Normal(N * mu_j, N * sigma^2 * K_{jj})`` — a scalar draw per
            ``(sample, j)``, replacing the ``O(N_max * J^2)`` matmul in
            the raw sampler.

        Speedup on the MVN case study: 40000 training steps drop from
        ~1000 s to ~10 s.

        Returns
        -------
        mu        : (S, J) float64  -- prior draws MVN(prior_mu, prior_tau^2 K).
        n         : (S,)  int64     -- observed-cohort size per sample.
        m         : (S,)  int64     -- future-cohort size per sample.
        sum_y_x   : (S, J) float64  -- sum_i y_{i, j} over the observed cohort.
        sum_y_z   : (S, J) float64  -- sum_i y_{i, j} over the future cohort.
        """
        # Prior mu: MVN(prior_mu, prior_tau^2 K).
        eps_mu = rng.standard_normal(size=(S, self.J))
        mu = self.prior_mu[None, :] + self.prior_tau * (eps_mu @ self.K_chol.T)
        # Cohort sizes.
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
            m = (n_max - n).astype(np.int64)
        else:
            if m_max is None:
                m_max = n_max
            m = rng.integers(m_min, m_max + 1, size=S).astype(np.int64)
        n = n.astype(np.int64)
        # Per-component sum draws: sum_i y_{i, j} ~ N(N * mu_j,
        # N * sigma^2 * K_{jj}).
        n_col = n.astype(np.float64)[:, None]              # (S, 1)
        m_col = m.astype(np.float64)[:, None]
        std_x = self.sigma * np.sqrt(
            n_col * self.K_diag[None, :],
        )                                                   # (S, J)
        std_z = self.sigma * np.sqrt(
            m_col * self.K_diag[None, :],
        )
        eps_sx = rng.standard_normal(size=(S, self.J))
        eps_sz = rng.standard_normal(size=(S, self.J))
        sum_y_x = n_col * mu + std_x * eps_sx
        sum_y_z = m_col * mu + std_z * eps_sz
        return mu, n, m, sum_y_x, sum_y_z

    def _sample_mvn_prior_predictive_raw(
        self, rng, S, n_max, n_min, fixed_total, m_min, m_max,
        n_dist='uniform',
    ):
        """Prior-predictive sampler for the raw-sequences amortiser.

        Under the g-prior with unit-diagonal K, ``y_{i, j} | mu ~
        Normal(mu_j, sigma^2 * K_{jj})`` independent across i, j given
        mu. Emitting per-component scalar sequences avoids the
        ``(S, N_max, J) @ (J, J)`` matmul in the original raw sampler
        while preserving the full raw-sequence signal the MLP q_tau
        encoder consumes.

        Returns
        -------
        mu        : (S, J) float64          -- prior draws.
        n         : (S,)  int64             -- observed cohort size.
        m         : (S,)  int64             -- future cohort size.
        Y_x       : (S, J, n_max) float64   -- per-component observed sequences,
                                               zero-padded beyond n.
        Y_z       : (S, J, n_max) float64   -- per-component future sequences,
                                               zero-padded beyond m.
        mask_x    : (S, n_max) float32      -- shared mask.
        mask_z    : (S, n_max) float32
        """
        # Prior mu via K_chol.
        eps_mu = rng.standard_normal(size=(S, self.J))
        mu = self.prior_mu[None, :] + self.prior_tau * (eps_mu @ self.K_chol.T)
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
            m = (n_max - n).astype(np.int64)
        else:
            if m_max is None:
                m_max = n_max
            m = rng.integers(m_min, m_max + 1, size=S).astype(np.int64)
        n = n.astype(np.int64)
        std_per_j = self.sigma * np.sqrt(self.K_diag)      # (J,)
        eps_x = rng.standard_normal(size=(S, self.J, n_max))
        eps_z = rng.standard_normal(size=(S, self.J, n_max))
        Y_x_full = mu[..., None] + std_per_j[None, :, None] * eps_x
        Y_z_full = mu[..., None] + std_per_j[None, :, None] * eps_z
        col_idx = np.arange(n_max)[None, :]
        mask_x = (col_idx < n[:, None]).astype(np.float32)
        mask_z = (col_idx < m[:, None]).astype(np.float32)
        Y_x = Y_x_full * mask_x[:, None, :]
        Y_z = Y_z_full * mask_z[:, None, :]
        return mu, n, m, Y_x, Y_z, mask_x, mask_z

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
        """Prior-predictive joint sampler for the features-fixed MVN
        amortiser (§8 / §9.3, MVN counterpart).

        Under the conjugate MVN g-prior with known ``K``, ``sigma``,
        ``prior_tau``, the posterior of ``mu_j | (x, z)`` depends on the
        per-component sufficient statistic ``(sum_i y_{i,j}, n + m)``.
        DeepSets pooling for hardcoded ``T(y_i) = y_i`` collapses to
        ``sum_T_x + sum_T_z = sum_j y_j``. Features passed to the head
        are ``(sum_y_j, n_total) / n_max`` (both normalised to
        ``[0, 1]`` in expectation for ``sum_y_j`` when ``mu_j`` is near 1).
        The amortiser is J-invariant: one net trained on flattened
        per-component examples deploys across all J values.

        Returns
        -------
        (batch, rho) : ``(dict, np.ndarray)``
            - ``batch['features']``: shape ``(S * J, 2)`` float32.
            - ``rho``: shape ``(S * J,)`` float32, per-component
              ``mu_j - mu_0_baseline`` (matches
              :meth:`get_endpoints_per_draw`'s ``ratio``).
        """
        mu, n, m, sum_y_x, sum_y_z = self._sample_mvn_prior_predictive_summary(
            rng, S, n_max, n_min, fixed_total, m_min, m_max, n_dist,
        )
        sum_y = sum_y_x + sum_y_z                          # (S, J)
        n_total = (n + m).astype(np.float32)               # (S,)
        n_total_b = np.broadcast_to(n_total[:, None], (S, self.J))
        features = np.stack(
            [sum_y.astype(np.float32), n_total_b], axis=-1,
        ) / float(n_max)                                   # (S, J, 2)
        rho = (mu - self.mu_0_baseline).astype(np.float32) # (S, J)
        features_flat = features.reshape(S * self.J, 2)
        rho_flat = rho.reshape(S * self.J)
        return {'features': features_flat}, rho_flat

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
        """Prior-predictive joint sampler for the features-MLP MVN
        amortiser.

        Emits raw padded per-component sequences (``item_dim = 1`` scalar
        per obs). One prior draw contributes ``J`` per-component training
        examples, so the amortiser is J-invariant.

        Returns
        -------
        (batch, rho) : ``(dict, np.ndarray)``
            - ``batch['x']``: shape ``(S * J, n_max, 1)`` float32.
            - ``batch['mask_x']``: shape ``(S * J, n_max)`` float32.
            - ``batch['z']``: shape ``(S * J, n_max, 1)`` float32.
            - ``batch['mask_z']``: shape ``(S * J, n_max)`` float32.
            - ``batch['sizes']``: shape ``(S * J, 2)`` float32,
              ``(n / n_max, m / n_max)``.
            - ``rho``: shape ``(S * J,)`` float32.
        """
        mu, n, m, Y_x, Y_z, mask_x, mask_z = (
            self._sample_mvn_prior_predictive_raw(
                rng, S, n_max, n_min, fixed_total, m_min, m_max, n_dist,
            )
        )
        # Y_x, Y_z are (S, J, n_max) already; add the item_dim=1 axis.
        x = Y_x[..., None].astype(np.float32).reshape(S * self.J, n_max, 1)
        z = Y_z[..., None].astype(np.float32).reshape(S * self.J, n_max, 1)
        # Masks are shared across components -- tile per J.
        mask_x_tiled = np.broadcast_to(
            mask_x[:, None, :], (S, self.J, n_max),
        ).reshape(S * self.J, n_max).astype(np.float32)
        mask_z_tiled = np.broadcast_to(
            mask_z[:, None, :], (S, self.J, n_max),
        ).reshape(S * self.J, n_max).astype(np.float32)
        sizes_row = np.stack(
            [n / n_max, m / n_max], axis=-1,
        ).astype(np.float32)                               # (S, 2)
        sizes = np.broadcast_to(
            sizes_row[:, None, :], (S, self.J, 2),
        ).reshape(S * self.J, 2).astype(np.float32)
        rho = (mu - self.mu_0_baseline).astype(np.float32).reshape(-1)
        return (
            {
                'x': x,
                'mask_x': mask_x_tiled,
                'z': z,
                'mask_z': mask_z_tiled,
                'sizes': sizes,
            },
            rho,
        )

    def make_training_data_with_participant_sequences(
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
        queries_per_sample: int = 4,
    ):
        """Prior-predictive joint sampler for the features-MLP **xcomp**
        amortiser (§12.10).

        Emits FULL participant vectors ``y_i in R^J`` (rather than the
        per-component scalar sequences of
        :meth:`make_training_data_with_raw_sequences`), plus the
        queried component's K-row and K-diag entry so the amortiser can
        exploit the known cross-component correlation ``K``.

        For each of the ``S`` prior draws we draw ``queries_per_sample``
        random query components ``j*`` and emit ``S * queries_per_sample``
        training examples that share the same participant matrices but
        differ in ``(k_row, k_diag, rho)``. This amortises the
        expensive prior sampling across multiple queries per
        gradient step.

        Returns
        -------
        (batch, rho) : ``(dict, np.ndarray)``
            - ``batch['x']``:      ``(B, n_max, J)`` float32.
            - ``batch['mask_x']``: ``(B, n_max)``    float32.
            - ``batch['z']``:      ``(B, n_max, J)`` float32.
            - ``batch['mask_z']``: ``(B, n_max)``    float32.
            - ``batch['sizes']``:  ``(B, 2)``        float32.
            - ``batch['k_row']``:  ``(B, J)``        float32, ``K[j*, :]``.
            - ``batch['k_diag']``: ``(B, 1)``        float32, ``K[j*, j*]``.
            - ``rho``:             ``(B,)`` float32, ``mu_{j*} - mu_0``.

            ``B = S * queries_per_sample``.
        """
        # Full participant matrices via K_chol (cross-component correlated).
        eps_mu = rng.standard_normal(size=(S, self.J))
        mu = self.prior_mu[None, :] + self.prior_tau * (eps_mu @ self.K_chol.T)
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
            m = (n_max - n).astype(np.int64)
        else:
            if m_max is None:
                m_max = n_max
            m = rng.integers(m_min, m_max + 1, size=S).astype(np.int64)
        n = n.astype(np.int64)
        eps_x = rng.standard_normal(size=(S, n_max, self.J))
        eps_z = rng.standard_normal(size=(S, n_max, self.J))
        Y_x_full = mu[:, None, :] + self.sigma * (eps_x @ self.K_chol.T)
        Y_z_full = mu[:, None, :] + self.sigma * (eps_z @ self.K_chol.T)
        col_idx = np.arange(n_max)[None, :]
        mask_x = (col_idx < n[:, None]).astype(np.float32)   # (S, n_max)
        mask_z = (col_idx < m[:, None]).astype(np.float32)
        Y_x = (Y_x_full * mask_x[..., None]).astype(np.float32)
        Y_z = (Y_z_full * mask_z[..., None]).astype(np.float32)
        Q = int(queries_per_sample)
        j_star = rng.integers(0, self.J, size=(S, Q))        # (S, Q)
        B = S * Q
        # Broadcast participant matrices across the Q queries per sample.
        x_rep = np.broadcast_to(
            Y_x[:, None, :, :], (S, Q, n_max, self.J),
        ).reshape(B, n_max, self.J).astype(np.float32)
        z_rep = np.broadcast_to(
            Y_z[:, None, :, :], (S, Q, n_max, self.J),
        ).reshape(B, n_max, self.J).astype(np.float32)
        mask_x_rep = np.broadcast_to(
            mask_x[:, None, :], (S, Q, n_max),
        ).reshape(B, n_max).astype(np.float32)
        mask_z_rep = np.broadcast_to(
            mask_z[:, None, :], (S, Q, n_max),
        ).reshape(B, n_max).astype(np.float32)
        sizes_row = np.stack(
            [n / n_max, m / n_max], axis=-1,
        ).astype(np.float32)                                  # (S, 2)
        sizes = np.broadcast_to(
            sizes_row[:, None, :], (S, Q, 2),
        ).reshape(B, 2).astype(np.float32)
        # K-row and K-diag for the queried component.
        k_row = self.K[j_star].astype(np.float32)             # (S, Q, J)
        k_diag = self.K_diag[j_star].astype(np.float32)       # (S, Q)
        # Target rho_{j*} = mu[s, j_star[s, q]] - mu_0.
        s_idx = np.arange(S)[:, None]
        mu_query = mu[s_idx, j_star]                          # (S, Q)
        rho = (mu_query - self.mu_0_baseline).astype(np.float32).reshape(B)
        return (
            {
                'x': x_rep,
                'mask_x': mask_x_rep,
                'z': z_rep,
                'mask_z': mask_z_rep,
                'sizes': sizes,
                'k_row': k_row.reshape(B, self.J),
                'k_diag': k_diag.reshape(B, 1),
            },
            rho,
        )

    def make_training_data_with_component_summaries_random_K(
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
        queries_per_sample: int = 4,
        K_families=None,
    ):
        """Prior-predictive sampler for the K-family-invariant
        attention amortiser (§12.11 xcompAtt).

        Each of the ``S`` prior draws picks a random ``K_chol`` from
        ``sample_random_K_chol`` (identity / AR(1) / block / factor by
        default) so the amortiser learns to interpret ``K[j*, :]`` under
        multiple correlation families. Per-component summaries
        ``(sum_i y_{i, j}, sum_i z_{i, j})`` replace the raw participant
        matrices of §12.10 -- the MVN sufficient statistic for
        ``mu_j | (x, z)`` is preserved, and the encoder no longer needs
        to run per-participant tokens.

        Returns
        -------
        (batch, rho) : ``(dict, np.ndarray)``
            - ``batch['x_sum']``: ``(B, J)`` float32,
              ``sum_i y_{i, j} / n_max``.
            - ``batch['z_sum']``: ``(B, J)`` float32,
              ``sum_i z_{i, j} / n_max``.
            - ``batch['k_row']``:  ``(B, J)`` float32, ``K[j*, :]``.
            - ``batch['k_diag']``: ``(B, 1)`` float32.
            - ``batch['sizes']``:  ``(B, 2)`` float32,
              ``(n / n_max, m / n_max)``.
            - ``rho``: ``(B,)`` float32, ``mu_{j*} - mu_0``.

            ``B = S * queries_per_sample``.
        """
        Q = int(queries_per_sample)
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
            m = (n_max - n).astype(np.int64)
        else:
            if m_max is None:
                m_max = n_max
            m = rng.integers(m_min, m_max + 1, size=S).astype(np.int64)
        n = n.astype(np.int64)

        # Loop per sample: draw K, then mu, sum_x, sum_z per component.
        # sum_i y_{i, j} ~ Normal(n * mu_j, n * sigma^2 K_{jj}); since
        # all K in the mixture have unit diagonal, K_{jj} = 1 always.
        sum_x = np.empty((S, self.J), dtype=np.float64)
        sum_z = np.empty((S, self.J), dtype=np.float64)
        mu = np.empty((S, self.J), dtype=np.float64)
        k_rows = np.empty((S, Q, self.J), dtype=np.float64)
        k_diags = np.empty((S, Q), dtype=np.float64)
        rho_out = np.empty((S, Q), dtype=np.float64)
        j_star = rng.integers(0, self.J, size=(S, Q))
        for s in range(S):
            K_chol_s, _fam = sample_random_K_chol(rng, self.J, K_families)
            K_s = K_chol_s @ K_chol_s.T
            # mu ~ MVN(prior_mu, tau^2 K)
            eps_mu = rng.standard_normal(self.J)
            mu[s] = self.prior_mu + self.prior_tau * (K_chol_s @ eps_mu)
            # Direct per-component sums (unit-diagonal K).
            eps_x = rng.standard_normal(self.J)
            eps_z = rng.standard_normal(self.J)
            sum_x[s] = n[s] * mu[s] + self.sigma * np.sqrt(n[s]) * eps_x
            sum_z[s] = m[s] * mu[s] + self.sigma * np.sqrt(m[s]) * eps_z
            k_rows[s] = K_s[j_star[s]]
            k_diags[s] = np.diag(K_s)[j_star[s]]
            rho_out[s] = mu[s, j_star[s]] - self.mu_0_baseline

        B = S * Q
        # Broadcast x_sum, z_sum, sizes across the Q queries per sample.
        x_sum_rep = np.broadcast_to(
            sum_x[:, None, :], (S, Q, self.J),
        ).reshape(B, self.J).astype(np.float32) / float(n_max)
        z_sum_rep = np.broadcast_to(
            sum_z[:, None, :], (S, Q, self.J),
        ).reshape(B, self.J).astype(np.float32) / float(n_max)
        sizes_row = np.stack(
            [n / n_max, m / n_max], axis=-1,
        ).astype(np.float32)
        sizes = np.broadcast_to(
            sizes_row[:, None, :], (S, Q, 2),
        ).reshape(B, 2).astype(np.float32)
        k_row_flat = k_rows.reshape(B, self.J).astype(np.float32)
        k_diag_flat = k_diags.reshape(B, 1).astype(np.float32)
        # Per (b, j) token: (x_sum[j], z_sum[j], K[j*, j]). The
        # itemXcompAtt class handles query identification via a learned
        # attention-score bias (see ``query_bias_alpha`` in the class)
        # driven by ``query_idx``, so no explicit is_query token
        # feature is needed here.
        tokens = np.stack(
            [x_sum_rep, z_sum_rep, k_row_flat], axis=-1,
        ).astype(np.float32)                                    # (B, J, 3)
        mask = np.ones((B, self.J), dtype=np.float32)
        query_idx = j_star.reshape(B).astype(np.int32)
        aux = np.concatenate(
            [sizes, k_diag_flat], axis=-1,
        ).astype(np.float32)                                    # (B, 3)
        return (
            {
                'tokens':    tokens,
                'mask':      mask,
                'query_idx': query_idx,
                'aux':       aux,
            },
            rho_out.reshape(B).astype(np.float32),
        )

    def make_training_data_with_participant_tokens_random_K(
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
        queries_per_sample: int = 4,
        K_families=None,
    ):
        """Prior-predictive sampler for the deepset{S,X}compAtt
        amortisers (§12.14).

        Emits RAW per-(participant, item) response scalars plus K-row
        per-(item, query) as item metadata. K-family-invariant training
        (random K per prior draw).

        Returns
        -------
        (batch, rho) : ``(dict, np.ndarray)`` with
            ``batch['x_responses']``:   ``(B, n_max, J, 1)`` float32.
            ``batch['mask_x']``:        ``(B, n_max)``       float32.
            ``batch['z_responses']``:   ``(B, n_max, J, 1)`` float32.
            ``batch['mask_z']``:        ``(B, n_max)``       float32.
            ``batch['item_metadata']``: ``(B, J, 1)`` float32, ``K[j*, :]``.
            ``batch['query_idx']``:     ``(B,)`` int32.
            ``batch['aux']``:           ``(B, 3)`` float32,
                                        ``(n/N, m/N, K[j*, j*])``.
            ``rho``:                    ``(B,)`` float32.
            ``B = S * queries_per_sample``.
        """
        Q = int(queries_per_sample)
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
            m = (n_max - n).astype(np.int64)
        else:
            if m_max is None:
                m_max = n_max
            m = rng.integers(m_min, m_max + 1, size=S).astype(np.int64)
        n = n.astype(np.int64)

        Y_x_full = np.zeros((S, n_max, self.J), dtype=np.float64)
        Y_z_full = np.zeros((S, n_max, self.J), dtype=np.float64)
        mu_all = np.empty((S, self.J), dtype=np.float64)
        K_all = np.empty((S, self.J, self.J), dtype=np.float64)
        for s in range(S):
            K_chol_s, _fam = sample_random_K_chol(rng, self.J, K_families)
            K_all[s] = K_chol_s @ K_chol_s.T
            eps_mu = rng.standard_normal(self.J)
            mu_all[s] = self.prior_mu + self.prior_tau * (K_chol_s @ eps_mu)
            eps_x = rng.standard_normal((n_max, self.J))
            eps_z = rng.standard_normal((n_max, self.J))
            Y_x_full[s] = mu_all[s][None, :] + self.sigma * (eps_x @ K_chol_s.T)
            Y_z_full[s] = mu_all[s][None, :] + self.sigma * (eps_z @ K_chol_s.T)

        col_idx = np.arange(n_max)[None, :]
        mask_x = (col_idx < n[:, None]).astype(np.float32)
        mask_z = (col_idx < m[:, None]).astype(np.float32)
        Y_x = (Y_x_full * mask_x[..., None]).astype(np.float32)
        Y_z = (Y_z_full * mask_z[..., None]).astype(np.float32)
        j_star = rng.integers(0, self.J, size=(S, Q))
        B = S * Q
        # Broadcast (S, n_max, J) to (S, Q, n_max, J) then reshape.
        x_rep = np.broadcast_to(
            Y_x[:, None, :, :], (S, Q, n_max, self.J),
        ).reshape(B, n_max, self.J).astype(np.float32)
        z_rep = np.broadcast_to(
            Y_z[:, None, :, :], (S, Q, n_max, self.J),
        ).reshape(B, n_max, self.J).astype(np.float32)
        mask_x_rep = np.broadcast_to(
            mask_x[:, None, :], (S, Q, n_max),
        ).reshape(B, n_max).astype(np.float32)
        mask_z_rep = np.broadcast_to(
            mask_z[:, None, :], (S, Q, n_max),
        ).reshape(B, n_max).astype(np.float32)
        sizes_row = np.stack(
            [n / n_max, m / n_max], axis=-1,
        ).astype(np.float32)
        sizes = np.broadcast_to(
            sizes_row[:, None, :], (S, Q, 2),
        ).reshape(B, 2).astype(np.float32)
        # K-row per (S, Q) queried component: K_all[s, j_star[s, q], :].
        s_idx = np.arange(S)[:, None]
        k_row = K_all[s_idx, j_star]                    # (S, Q, J)
        k_diag = np.take_along_axis(
            np.diagonal(K_all, axis1=-2, axis2=-1)[:, None, :],
            j_star[..., None], axis=-1,
        ).squeeze(-1)                                    # (S, Q)
        mu_query = mu_all[s_idx, j_star]                # (S, Q)
        rho = (mu_query - self.mu_0_baseline).astype(np.float32).reshape(B)
        return (
            {
                'x_responses':   x_rep.reshape(B, n_max, self.J, 1),
                'mask_x':        mask_x_rep,
                'z_responses':   z_rep.reshape(B, n_max, self.J, 1),
                'mask_z':        mask_z_rep,
                'item_metadata': k_row.reshape(B, self.J, 1).astype(np.float32),
                'query_idx':     j_star.reshape(B).astype(np.int32),
                'aux':           np.concatenate(
                    [sizes, k_diag.reshape(B, 1).astype(np.float32)],
                    axis=-1,
                ).astype(np.float32),
            },
            rho,
        )

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

    def fit_pyro_svi(
        self,
        output_file_prefix: Optional[str] = None,
        *,
        algorithm: str = 'AutoMultivariateNormal',
        lr: float = 0.05,
        num_steps: int = 2000,
        output_samples: int = 4000,
        save_to_file: bool = False,
        resume: bool = False,
        verbose: bool = True,
        **method_kwargs,
    ) -> Dict[str, Any]:
        """NumPyro SVI on the MVN program. ``algorithm`` selects the
        autoguide (default ``AutoMultivariateNormal``; pass
        ``AutoDiagonalNormal`` / ``AutoLowRankMultivariateNormal`` etc).
        Mirrors :meth:`BinomialModel.fit_pyro_svi`. When ``save_to_file``
        and ``resume`` are True and the ``{prefix}_draws.zarr`` +
        ``{prefix}_timing.csv`` artifacts exist, the cached posterior is
        returned instead of refitting."""
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
                print(f"[MVNModel] pyro SVI ({algorithm}) RESUME from"
                      f" {draws_file}")
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
        mu_samples = np.asarray(post['mu'])
        rng, sub = jax.random.split(rng)
        N = int(self.stan_data['N'])
        predictions = self._pyro_predictive(
            {'mu': mu_samples}, sub, N,
        )
        n_draw = mu_samples.shape[0]
        ds = xr.Dataset({
            'mu':      (('chain', 'draw', 'j'), mu_samples[None, :, :]),
            'ypred':   (('chain', 'draw', 'obs'),
                        predictions['ypred'].reshape(1, n_draw, N * self.J)),
            'log_lik': (('chain', 'draw', 'obs_vec'),
                        predictions['log_lik'].reshape(1, n_draw, N)),
        })
        draws = az.InferenceData(posterior=ds)
        timing = {'mins': (time.time() - t0) / 60.0,
                  'final_elbo': float(-run.losses[-1]),
                  'algorithm': algorithm}
        if verbose:
            print(f"[MVNModel] pyro SVI ({algorithm}): {num_steps}"
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
            'posterior_samples': {'mu': mu_samples},
            'timing': timing,
        }

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
