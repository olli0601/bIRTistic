"""
Step 7 of the OO-port refactor: :class:`BinomialModel` -- the first non-
IRT subclass of :class:`model.Model`. Proves the blueprint generalises
beyond the partial-credit family.

Tests pin down:
- closed-form posterior matches the analytic Beta(a+k, b+n-k) (§3.1 of
  ``dev/amortised_decision_making.md``)
- eval_loglik, eval_loglik_annealed, logprior, eval_outcome_for_endpoint
  compute the correct closed-form quantities
- endpoints_per_draw returns the expected (item, draw) shape with
  ratio = p
- fit_pyro_svi / fit_stan_svi / fit_stan_hmc raise NotImplementedError
  (future work; the blueprint accepts them as separate methods so
  plugging them in later is a drop-in replacement)
"""

import sys
from pathlib import Path

import jax.numpy as jnp
import numpy as np
import pandas as pd
import pytest

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


def test_logprior_matches_beta(model):
    import scipy.stats as st
    p = 0.5
    expected = float(st.beta.logpdf(p, _A, _B))
    actual = float(model.logprior({'p': jnp.asarray(p)}))
    assert actual == pytest.approx(expected, abs=1e-5)


def test_eval_outcome_returns_p(model):
    p_val = jnp.asarray(0.42)
    out = model.eval_outcome_for_endpoint(model.stan_data, {'p': p_val})
    assert float(out) == pytest.approx(0.42, abs=1e-5)


# ---------------------------------------------------------------------------
# Closed-form fit
# ---------------------------------------------------------------------------


def test_fit_closed_form_posterior_params(model):
    """Closed-form posterior is Beta(prior_a + k, prior_b + n - k)."""
    fit = model.fit_closed_form(output_samples=2000, verbose=False)
    assert fit['a_post'] == pytest.approx(_A + _K, abs=1e-5)
    assert fit['b_post'] == pytest.approx(_B + (_N - _K), abs=1e-5)


def test_fit_closed_form_posterior_mean_matches_analytic(model):
    """Beta(a, b) has mean a / (a + b). With 2000 samples the MC mean
    matches within ~0.01."""
    fit = model.fit_closed_form(output_samples=5000, verbose=False)
    expected_mean = (_A + _K) / (_A + _B + _N)
    actual_mean = float(fit['posterior_samples']['p'].mean())
    assert actual_mean == pytest.approx(expected_mean, abs=0.01)


def test_fit_closed_form_returns_arviz_draws(model):
    fit = model.fit_closed_form(output_samples=100, verbose=False)
    draws = fit['draws']
    arr = draws.posterior['p'].values
    assert arr.shape == (1, 100)  # (chain, draw)


# ---------------------------------------------------------------------------
# endpoints_per_draw
# ---------------------------------------------------------------------------


def test_endpoints_per_draw_shape(model):
    fit = model.fit_closed_form(output_samples=50, verbose=False)
    end = model.endpoints_per_draw(
        dcati=model.dcati, draws=fit['draws'],
        categorical_threshold=2,
    )
    assert len(end) == 50
    assert set(end.columns) >= {
        'item_label', 'item_type', 'item_high_label', 'draw', 'ratio',
    }
    # ratio column == p samples (in order)
    np.testing.assert_array_equal(
        end['ratio'].to_numpy(),
        fit['posterior_samples']['p'],
    )


# ---------------------------------------------------------------------------
# Unimplemented fitters raise loudly
# ---------------------------------------------------------------------------


@pytest.mark.parametrize('method', ['fit_pyro_svi', 'fit_stan_svi', 'fit_stan_hmc'])
def test_unimplemented_fitters_raise(model, method):
    with pytest.raises(NotImplementedError):
        getattr(model, method)('prefix')
