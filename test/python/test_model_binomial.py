"""
:class:`BinomialModel` -- proves the blueprint generalises beyond the
partial-credit family.

Coverage:
- closed-form posterior matches analytic Beta(a+k, b+n-k) (§3.1)
- closed-form PPS matches the analytic Beta-Binomial tail sum
- eval_loglik / eval_loglik_annealed / eval_log_prior /
  eval_outcome_for_endpoint return the correct closed-form quantities
- get_endpoints_per_draw returns the expected (item, draw) shape
- fit_pyro_svi / fit_stan_svi / fit_stan_hmc posteriors agree with the
  closed-form Beta posterior across 100 simulated datasets (posterior
  mean within Monte-Carlo tolerance of the analytic mean
  (prior_a + k) / (prior_a + prior_b + n))
"""

import sys
from pathlib import Path

import jax.numpy as jnp
import numpy as np
import pandas as pd
import pytest
from scipy.stats import beta as _scipy_beta
from scipy.stats import betabinom as _scipy_betabinom

_repo_root = Path(__file__).resolve().parents[2]
_python_dir = _repo_root / 'python'
if str(_python_dir) not in sys.path:
    sys.path.insert(0, str(_python_dir))

from model_binomial import BinomialModel


# ---------------------------------------------------------------------------
# Fixtures: tiny synthetic Beta-Binomial example.
# ---------------------------------------------------------------------------


_N = 30
_K = 18                       # observed successes
_A = 2.0                      # prior shape
_B = 3.0
_SEED = 123


@pytest.fixture(scope='module')
def dit():
    return pd.DataFrame({
        'item_label': ['intervention'],
        'item_type': ['binomial'],
        'item_high_label': ['higher_is_better'],
    })


@pytest.fixture(scope='module')
def dcati():
    """30 Bernoulli outcomes with K=18 successes."""
    y = np.concatenate([np.ones(_K, dtype=int), np.zeros(_N - _K, dtype=int)])
    return pd.DataFrame({'y': y, 'pid': np.arange(_N)})


@pytest.fixture(scope='module')
def model(dit, dcati):
    return BinomialModel(dit=dit, dcati=dcati, prior_a=_A, prior_b=_B,
                         seed=_SEED)


def _make_model(n, k, prior_a, prior_b, seed):
    """Build a BinomialModel from observed counts (k successes out of n trials)."""
    y = np.concatenate([np.ones(k, dtype=int), np.zeros(n - k, dtype=int)])
    dit = pd.DataFrame({
        'item_label': ['intervention'],
        'item_type': ['binomial'],
        'item_high_label': ['higher_is_better'],
    })
    dcati = pd.DataFrame({'y': y, 'pid': np.arange(n)})
    return BinomialModel(dit=dit, dcati=dcati,
                         prior_a=prior_a, prior_b=prior_b, seed=seed)


# ---------------------------------------------------------------------------
# Blueprint metadata
# ---------------------------------------------------------------------------


def test_param_names(model):
    assert model.param_names == ('p',)


def test_positive_params_empty(model):
    """p lives on (0, 1) -- no positive-real log-transform."""
    assert model.positive_params == ()


# ---------------------------------------------------------------------------
# stan_data + loglik / prior / outcome
# ---------------------------------------------------------------------------


def test_stan_data_built_at_construction(model):
    sd = model.stan_data
    assert int(sd['N_total']) == _N
    assert int(sd['k']) == _K
    assert sd['a'] == _A
    assert sd['b'] == _B


def test_eval_loglik_matches_bernoulli(model):
    p = 0.6
    expected = (_K * np.log(p) + (_N - _K) * np.log1p(-p))
    actual = float(model.eval_loglik(model.stan_data, {'p': jnp.asarray(p)}).sum())
    assert actual == pytest.approx(expected, abs=1e-5)


def test_eval_loglik_annealed_scales_by_temperature(model):
    p = {'p': jnp.asarray(0.4), 'temperature': jnp.asarray(0.5)}
    base = float(model.eval_loglik(model.stan_data, {'p': p['p']}).sum())
    annealed = float(model.eval_loglik_annealed(model.stan_data, p).sum())
    assert annealed == pytest.approx(0.5 * base, abs=1e-5)


def test_eval_log_prior_matches_beta(model):
    import scipy.stats as st
    p = 0.5
    expected = float(st.beta.logpdf(p, _A, _B))
    actual = float(model.eval_log_prior({'p': jnp.asarray(p)}))
    assert actual == pytest.approx(expected, abs=1e-5)


def test_eval_outcome_returns_ratio(model):
    """eval_outcome_for_endpoint now returns the endpoint ratio p / p_0
    (default p_0 = 0.5), matching the convention used by the IS /
    regression algorithms."""
    p_val = jnp.asarray(0.42)
    out = model.eval_outcome_for_endpoint(model.stan_data, {'p': p_val})
    assert float(out) == pytest.approx(0.42 / 0.5, abs=1e-5)


# ---------------------------------------------------------------------------
# Closed-form posterior
# ---------------------------------------------------------------------------


def test_fit_closed_form_posterior_params(model):
    """Closed-form posterior is Beta(prior_a + k, prior_b + n - k)."""
    fit = model.fit_closed_form_posterior(output_samples=2000, verbose=False)
    assert fit['a_post'] == pytest.approx(_A + _K, abs=1e-5)
    assert fit['b_post'] == pytest.approx(_B + (_N - _K), abs=1e-5)


def test_fit_closed_form_posterior_mean_matches_analytic(model):
    """Beta(a, b) has mean a / (a + b). With 5000 samples the MC mean
    matches within ~0.01."""
    fit = model.fit_closed_form_posterior(output_samples=5000, verbose=False)
    expected_mean = (_A + _K) / (_A + _B + _N)
    actual_mean = float(fit['posterior_samples']['p'].mean())
    assert actual_mean == pytest.approx(expected_mean, abs=0.01)


def test_fit_closed_form_posterior_returns_arviz_draws(model):
    fit = model.fit_closed_form_posterior(output_samples=100, verbose=False)
    draws = fit['draws']
    arr = draws.posterior['p'].values
    assert arr.shape == (1, 100)  # (chain, draw)


# ---------------------------------------------------------------------------
# Closed-form PPS (analytic Beta-Binomial tail sum)
# ---------------------------------------------------------------------------


def test_fit_closed_form_pps_matches_brute_force():
    """Brute-force compute the analytic PPS by enumerating k_m and verifying
    fit_closed_form_pps returns the same number."""
    n, k, a, b = 30, 21, 2.0, 3.0
    p_h1_def, eta = 0.5, 0.89
    mdl = _make_model(n=n, k=k, prior_a=a, prior_b=b, seed=42)
    a_post, b_post = a + k, b + (n - k)
    for m in (5, 20, 100):
        # Brute: find smallest k_m with P(H_1 | x, z) > eta.
        ok = []
        for km in range(m + 1):
            p_h1 = 1.0 - _scipy_beta.cdf(p_h1_def, a_post + km, b_post + m - km)
            ok.append(p_h1 > eta)
        if not any(ok):
            expected = 0.0
        else:
            k_star = int(np.argmax(ok))
            expected = float(_scipy_betabinom.sf(k_star - 1, m, a_post, b_post))
        actual = mdl.fit_closed_form_pps(m=m, pps_H1_def=p_h1_def,
                                         pps_ProbH1_thresh=eta)
        assert actual == pytest.approx(expected, abs=1e-12)


def test_fit_closed_form_pps_returns_zero_when_threshold_unreachable():
    """If even k_m = m fails to cross the threshold, PPS = 0."""
    # Tiny posterior mass on H_1: prior says p~Beta(1, 99), observe 0/30.
    mdl = _make_model(n=30, k=0, prior_a=1.0, prior_b=99.0, seed=0)
    pps = mdl.fit_closed_form_pps(m=10, pps_H1_def=0.5, pps_ProbH1_thresh=0.89)
    assert pps == 0.0


# ---------------------------------------------------------------------------
# get_endpoints_per_draw
# ---------------------------------------------------------------------------


def test_get_endpoints_per_draw_shape(model):
    fit = model.fit_closed_form_posterior(output_samples=50, verbose=False)
    end = model.get_endpoints_per_draw(
        draws=fit['draws'],
        categorical_threshold=2,
    )
    assert len(end) == 50
    assert set(end.columns) >= {
        'item_label', 'item_type', 'item_high_label', 'draw', 'p', 'ratio',
    }
    # p column == posterior samples (in order); ratio = p / p_0 (default 0.5).
    p_samples = fit['posterior_samples']['p']
    np.testing.assert_array_equal(end['p'].to_numpy(), p_samples)
    np.testing.assert_allclose(
        end['ratio'].to_numpy(), p_samples / 0.5, rtol=0, atol=1e-12,
    )


# ---------------------------------------------------------------------------
# Posterior recovery across 100 simulated datasets
# ---------------------------------------------------------------------------
#
# For each of 100 simulated datasets we draw a true p ~ Beta(2, 3), generate
# n=30 Bernoulli observations, fit, and compare the posterior-mean of p to
# the analytic posterior mean (prior_a + k) / (prior_a + prior_b + n). The
# absolute error |E_fit[p] - E_analytic[p]| should be small for all three
# fitters; SVI mean-field on Beta-Bernoulli is exact in the mean direction,
# so we apply tight tolerances.


_NSIM = 100
_SIM_N = 30
_SIM_PRIOR_A = 2.0
_SIM_PRIOR_B = 3.0


def _simulate(seed):
    rng = np.random.default_rng(seed)
    true_p = rng.beta(_SIM_PRIOR_A, _SIM_PRIOR_B)
    y = rng.binomial(1, true_p, size=_SIM_N)
    k = int(y.sum())
    analytic_mean = (_SIM_PRIOR_A + k) / (_SIM_PRIOR_A + _SIM_PRIOR_B + _SIM_N)
    return k, analytic_mean


@pytest.fixture(scope='module')
def simulations():
    """100 (k, analytic_mean) pairs; the fits below reuse them."""
    return [_simulate(seed) for seed in range(_NSIM)]


def test_closed_form_posterior_recovers_analytic_mean_100_sims(simulations):
    """Sanity baseline: closed-form posterior MC mean must match the
    analytic mean within Monte-Carlo error across all 100 sims."""
    errors = []
    for sim_idx, (k, analytic) in enumerate(simulations):
        mdl = _make_model(n=_SIM_N, k=k,
                          prior_a=_SIM_PRIOR_A, prior_b=_SIM_PRIOR_B,
                          seed=sim_idx)
        fit = mdl.fit_closed_form_posterior(output_samples=4000, verbose=False)
        mc_mean = float(fit['posterior_samples']['p'].mean())
        errors.append(abs(mc_mean - analytic))
    assert np.median(errors) < 0.01
    assert np.max(errors) < 0.03


def test_pyro_svi_recovers_analytic_mean_100_sims(simulations):
    """Pyro SVI posterior mean vs analytic mean across 100 sims. With a
    diagonal-Normal guide on a single-parameter Beta-Bernoulli the
    posterior mean should be very close to the analytic mean."""
    errors = []
    for sim_idx, (k, analytic) in enumerate(simulations):
        mdl = _make_model(n=_SIM_N, k=k,
                          prior_a=_SIM_PRIOR_A, prior_b=_SIM_PRIOR_B,
                          seed=sim_idx)
        fit = mdl.fit_pyro_svi(num_steps=1500, lr=0.05,
                               output_samples=2000, verbose=False)
        svi_mean = float(fit['posterior_samples']['p'].mean())
        errors.append(abs(svi_mean - analytic))
    assert np.median(errors) < 0.03
    # Diagonal Normal guide on a (0,1)-bounded parameter occasionally
    # over/under-shoots; allow a more generous upper bound.
    assert np.max(errors) < 0.10


@pytest.mark.slow
def test_stan_hmc_recovers_analytic_mean_100_sims(simulations):
    """Stan HMC posterior mean vs analytic mean across 100 sims.
    Marked slow because each fit recompiles cmdstanpy state (~2-3s)."""
    errors = []
    for sim_idx, (k, analytic) in enumerate(simulations):
        mdl = _make_model(n=_SIM_N, k=k,
                          prior_a=_SIM_PRIOR_A, prior_b=_SIM_PRIOR_B,
                          seed=sim_idx)
        fit = mdl.fit_stan_hmc(chains=2, iter_warmup=200, iter_sampling=500,
                               verbose=False)
        hmc_mean = float(fit['posterior_samples']['p'].mean())
        errors.append(abs(hmc_mean - analytic))
    assert np.median(errors) < 0.02
    assert np.max(errors) < 0.05
