"""
Step 3 of the OO-port refactor: :class:`PartialCreditModelNCats` must
return bit-identical outputs to the existing free-function pipeline in
:mod:`fit_partial_credit_model` for every blueprint method.

These tests load the fixture cohort in ``test/test_data/`` and pair each
``model.X(...)`` call with the corresponding free-function call,
asserting the two agree.
"""

import sys
from pathlib import Path

import arviz as az
import jax.numpy as jnp
import numpy as np
import pandas as pd
import pytest

_repo_root = Path(__file__).resolve().parents[2]
_python_dir = _repo_root / 'python'
if str(_python_dir) not in sys.path:
    sys.path.insert(0, str(_python_dir))

from get_endpoints import get_endpoints_per_draw
from model_pcm import PartialCreditModelNCats, _model_namespace
from fit_interim import _stack_posterior_theta  # legacy version, hard-coded

_TEST_DATA = _repo_root / 'test' / 'test_data'
_INTERIM_STEM = _TEST_DATA / 'pcm_1_interim_1'
_DRAWS_FILE = f"{_INTERIM_STEM}_draws.zarr"


@pytest.fixture(scope='module')
def dit():
    return pd.read_csv(f"{_INTERIM_STEM}_data_dit.csv")


@pytest.fixture(scope='module')
def xi():
    return pd.read_csv(f"{_INTERIM_STEM}_data_dp1.csv")


@pytest.fixture(scope='module')
def draws():
    return az.from_zarr(_DRAWS_FILE)


@pytest.fixture(scope='module')
def model(dit, xi):
    return PartialCreditModelNCats(dit=dit, dcati=xi,
                                   x_formula="~ time - 1", seed=123)


@pytest.fixture(scope='module')
def stan_data(model):
    return model.stan_data


@pytest.fixture(scope='module')
def params_one(model, draws):
    """A single posterior-draw parameter dict (the first draw) for use
    with eval_loglik / logprior / eval_outcome."""
    theta = model.stack_posterior_theta(draws)
    return {k: v[0] for k, v in theta.items()}


# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------


def test_param_names_includes_known_pcm_params(model):
    expected = {
        'latent_factor_unit', 'latent_factor_beta',
        'loadings_questions_m1',
    }
    assert expected.issubset(set(model.param_names))


def test_positive_params_is_subset_of_param_names(model):
    assert set(model.positive_params).issubset(set(model.param_names))


# ---------------------------------------------------------------------------
# Stan data
# ---------------------------------------------------------------------------


def test_stan_data_built_at_construction(model, stan_data):
    assert isinstance(stan_data, dict)
    assert 'N_total' in stan_data
    assert int(stan_data['N_total']) == len(model.dcati)


def test_make_stan_data_shape_invariants(model, stan_data, dit, xi):
    """make_stan_data must produce the standard partial-credit stan_data
    keys + shapes against the fixture cohort."""
    expected_keys = {
        'C', 'U', 'N', 'Q', 'K', 'N_total', 'Q_total', 'y', 'unit_of_obs',
        'question_of_obs', 'cat_type', 'X', 'Xpr_id', 'P', 'P_pr',
    }
    assert expected_keys.issubset(set(stan_data.keys()))
    assert int(stan_data['N_total']) == len(xi)
    assert int(stan_data['C']) == int(dit['item_type_id'].max())
    # Call directly: model.make_stan_data(dcati, x_formula) returns the same.
    again = model.make_stan_data(xi, "~ time - 1")
    assert set(again.keys()) == set(stan_data.keys())


# ---------------------------------------------------------------------------
# Loglik / prior / outcome -- the @staticmethod proxies must agree with the
# underlying .pyro namespace they delegate to.
# ---------------------------------------------------------------------------


def test_eval_loglik_matches_pyro_namespace(model, stan_data, params_one):
    ns = _model_namespace()
    direct = ns['get_log_likelihood_of_partial_credit_model_ncats'](stan_data, params_one)
    by_method = model.eval_loglik(stan_data, params_one)
    np.testing.assert_array_equal(np.asarray(direct), np.asarray(by_method))


def test_eval_loglik_annealed_matches_pyro_namespace(model, stan_data, params_one):
    p = dict(params_one)
    p['temperature'] = jnp.asarray(0.5)
    ns = _model_namespace()
    direct = ns['get_log_likelihood_of_partial_credit_model_ncats_with_annealing'](stan_data, p)
    by_method = model.eval_loglik_annealed(stan_data, p)
    np.testing.assert_array_equal(np.asarray(direct), np.asarray(by_method))


def test_logprior_matches_pyro_namespace(model, params_one):
    direct = _model_namespace()['get_prior_of_partial_credit_model_ncats'](params_one)
    by_method = model.logprior(params_one)
    assert float(direct) == float(by_method)


def test_eval_outcome_matches_pyro_namespace(model, stan_data, params_one):
    direct = _model_namespace()['get_ordered_prob_of_partial_credit_model_ncats'](stan_data, params_one)
    by_method = model.eval_outcome_for_endpoint(stan_data, params_one)
    np.testing.assert_array_equal(np.asarray(direct), np.asarray(by_method))


# ---------------------------------------------------------------------------
# Endpoints per draw
# ---------------------------------------------------------------------------


def test_endpoints_per_draw_matches_free_function(model, dit, xi, draws):
    free = get_endpoints_per_draw(
        dcati=xi, dit=dit, draws=draws,
        categorical_threshold=2,
        endpoint_type='items', param_name='ordered_prob_by_cat_qu_fit',
        verbose=False,
    )
    by_method = model.endpoints_per_draw(
        dcati=xi, draws=draws, categorical_threshold=2,
        endpoint_type='items',
    )
    pd.testing.assert_frame_equal(
        free.reset_index(drop=True),
        by_method.reset_index(drop=True),
        check_dtype=True, check_exact=True,
    )


# ---------------------------------------------------------------------------
# Posterior stacking
# ---------------------------------------------------------------------------


def test_stack_posterior_theta_matches_legacy(model, draws):
    """The new param-name-driven stacker must produce identical arrays
    to the hard-coded :func:`fit_interim._stack_posterior_theta`."""
    legacy = _stack_posterior_theta(draws)
    by_method = model.stack_posterior_theta(draws)
    assert set(legacy.keys()) == set(by_method.keys())
    for k in legacy:
        np.testing.assert_array_equal(np.asarray(legacy[k]),
                                      np.asarray(by_method[k]))
