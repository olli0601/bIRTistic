"""
:class:`CreditModel` and :class:`OrderedLogit` -- the credit and
ordered-logit subclasses of :class:`model_irt.IRTModel`. Both inline
their fit drivers + loglik / prior proxies as class methods
(model_credit.py / model_ordered_logit.py); the legacy fit_credit_model.py
and fit_ordered_logit_model.py modules are deleted.

Coverage:
- structural invariants: param_names match the .pyro sample sites;
  proxies dispatch to functions present in the .pyro namespace;
  eval_outcome_for_endpoint raises NotImplementedError.
- make_stan_data: returns the expected stan-data dict against the PCM
  fixture cohort (the dp1 / dit columns are identical across the three
  IRT subclasses).
- @staticmethod proxy equivalence: ``model.eval_loglik(stan_data, params)``
  matches ``_model_namespace()[...](stan_data, params)`` for synthetic
  parameter dicts whose shapes are derived from stan_data (we cannot use
  cached PCM draws because credit / OL do not share the same param-name
  set).
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

from model_credit import CreditModel, _model_namespace as _credit_ns
from model_ordered_logit import (
    OrderedLogit, _model_namespace as _ol_ns,
)


_TEST_DATA = _repo_root / 'test' / 'test_data'
_INTERIM_STEM = _TEST_DATA / 'pcm_1_interim_1'


# ---------------------------------------------------------------------------
# Fixtures (shared with the PCM tests; same dp1/dit schema works for credit
# and OL even though the underlying numpyro models differ).
# ---------------------------------------------------------------------------


@pytest.fixture(scope='module')
def dit():
    return pd.read_csv(f"{_INTERIM_STEM}_data_dit.csv")


@pytest.fixture(scope='module')
def xi():
    return pd.read_csv(f"{_INTERIM_STEM}_data_dp1.csv")


@pytest.fixture(scope='module')
def credit_model(dit, xi):
    return CreditModel(dit=dit, dcati=xi, x_formula="~ time - 1", seed=123)


@pytest.fixture(scope='module')
def ol_model(dit, xi):
    return OrderedLogit(dit=dit, dcati=xi, x_formula="~ time - 1", seed=123)


def _synthetic_params(stan_data, param_shapes, seed):
    """Build a parameter dict with random arrays matching the requested
    shapes. ``loadings_questions_m1`` and ``skill_thresholds_incs`` are
    forced positive (the model expects positive loadings + monotone-
    increasing thresholds)."""
    rng = np.random.default_rng(seed)
    out = {}
    for name, shape in param_shapes.items():
        arr = rng.standard_normal(shape).astype(np.float32)
        if name in ('loadings_questions_m1', 'skill_thresholds_incs'):
            arr = np.abs(arr) + 0.1
        out[name] = jnp.asarray(arr)
    return out


def _credit_param_shapes(stan_data) -> dict:
    """Shapes for credit subclass parameters derived from stan_data."""
    U = int(stan_data['U'])
    P = int(stan_data['P'])
    Q_arr = stan_data['Q']
    K_arr = stan_data['K']
    skill_thresholds = sum(int(Q_arr[c]) * (int(K_arr[c]) - 1)
                           for c in range(len(Q_arr)))
    loadings_questions_m1 = sum(max(int(Q_arr[c]) - 1, 0)
                                for c in range(len(Q_arr)))
    return {
        'latent_factor_unit': (U,),
        'latent_factor_beta': (P,),
        'skill_thresholds': (skill_thresholds,),
        'loadings_questions_m1': (loadings_questions_m1,),
    }


def _ol_param_shapes(stan_data) -> dict:
    """Shapes for ordered-logit subclass parameters."""
    U = int(stan_data['U'])
    P = int(stan_data['P'])
    Q_arr = stan_data['Q']
    K_arr = stan_data['K']
    skill_thresholds_1 = sum(int(Q_arr[c]) for c in range(len(Q_arr)))
    skill_thresholds_incs = sum(
        int(Q_arr[c]) * max(int(K_arr[c]) - 2, 0)
        for c in range(len(Q_arr))
    )
    loadings_questions_m1 = sum(max(int(Q_arr[c]) - 1, 0)
                                for c in range(len(Q_arr)))
    return {
        'latent_factor_unit': (U,),
        'latent_factor_beta': (P,),
        'skill_thresholds_1': (skill_thresholds_1,),
        'skill_thresholds_incs': (skill_thresholds_incs,),
        'loadings_questions_m1': (loadings_questions_m1,),
    }


# ---------------------------------------------------------------------------
# Credit
# ---------------------------------------------------------------------------


def test_credit_param_names_match_pyro_model():
    """Sample sites in src/numpyro/credit_model_ncats_v260413.pyro are
    latent_factor_unit / latent_factor_beta / skill_thresholds /
    loadings_questions_m1 -- declared in the subclass."""
    assert set(CreditModel.param_names) == {
        'latent_factor_unit', 'latent_factor_beta',
        'skill_thresholds', 'loadings_questions_m1',
    }
    assert CreditModel.positive_params == ('loadings_questions_m1',)


def test_credit_proxies_resolve_in_pyro_namespace():
    ns = _credit_ns()
    assert 'get_log_likelihood_of_credit_model_ncats' in ns
    assert 'get_log_likelihood_of_credit_model_ncats_with_annealing' in ns
    assert 'get_prior_of_credit_model_ncats' in ns


def test_credit_eval_outcome_raises_not_implemented():
    """eval_outcome_for_endpoint must raise so MM / SMC fail loudly
    rather than silently scoring with wrong outputs."""
    with pytest.raises(NotImplementedError):
        CreditModel.eval_outcome_for_endpoint(data={}, params={})


def test_credit_make_stan_data_shape_invariants(credit_model, dit, xi):
    sd = credit_model.stan_data
    expected_keys = {
        'C', 'U', 'N', 'Q', 'K', 'N_total', 'Q_total', 'y', 'unit_of_obs',
        'question_of_obs', 'cat_type', 'X', 'Xpr_id', 'P', 'P_pr',
    }
    assert expected_keys.issubset(set(sd.keys()))
    assert int(sd['N_total']) == len(xi)
    assert int(sd['C']) == int(dit['item_type_id'].max())


def test_credit_eval_loglik_matches_pyro_namespace(credit_model):
    """The @staticmethod proxy must agree with the underlying .pyro call
    on synthetic parameters of the right shape."""
    sd = credit_model.stan_data
    params = _synthetic_params(sd, _credit_param_shapes(sd), seed=0)
    ns = _credit_ns()
    direct = ns['get_log_likelihood_of_credit_model_ncats'](sd, params)
    by_method = credit_model.eval_loglik(sd, params)
    np.testing.assert_array_equal(np.asarray(direct), np.asarray(by_method))


def test_credit_eval_loglik_annealed_matches_pyro_namespace(credit_model):
    sd = credit_model.stan_data
    params = _synthetic_params(sd, _credit_param_shapes(sd), seed=1)
    params['temperature'] = jnp.asarray(0.5)
    ns = _credit_ns()
    direct = ns['get_log_likelihood_of_credit_model_ncats_with_annealing'](sd, params)
    by_method = credit_model.eval_loglik_annealed(sd, params)
    np.testing.assert_array_equal(np.asarray(direct), np.asarray(by_method))


def test_credit_eval_log_prior_matches_pyro_namespace(credit_model):
    sd = credit_model.stan_data
    params = _synthetic_params(sd, _credit_param_shapes(sd), seed=2)
    ns = _credit_ns()
    direct = ns['get_prior_of_credit_model_ncats'](params)
    by_method = credit_model.eval_log_prior(params)
    assert float(direct) == float(by_method)


# ---------------------------------------------------------------------------
# OrderedLogit
# ---------------------------------------------------------------------------


def test_ordered_logit_param_names_match_pyro_model():
    """Sample sites in src/numpyro/ordered_logit_ncats_v260413.pyro are
    latent_factor_unit / latent_factor_beta / skill_thresholds_1 /
    skill_thresholds_incs / loadings_questions_m1."""
    assert set(OrderedLogit.param_names) == {
        'latent_factor_unit', 'latent_factor_beta',
        'skill_thresholds_1', 'skill_thresholds_incs',
        'loadings_questions_m1',
    }
    assert OrderedLogit.positive_params == ('loadings_questions_m1',)


def test_ordered_logit_proxies_resolve_in_pyro_namespace():
    ns = _ol_ns()
    assert 'get_log_likelihood_of_ordered_logit_ncats' in ns
    assert 'get_log_likelihood_of_ordered_logit_ncats_with_annealing' in ns
    assert 'get_prior_of_ordered_logit_ncats' in ns


def test_ordered_logit_eval_outcome_raises_not_implemented():
    with pytest.raises(NotImplementedError):
        OrderedLogit.eval_outcome_for_endpoint(data={}, params={})


def test_ordered_logit_make_stan_data_shape_invariants(ol_model, dit, xi):
    sd = ol_model.stan_data
    expected_keys = {
        'C', 'U', 'N', 'Q', 'K', 'N_total', 'Q_total', 'y', 'unit_of_obs',
        'question_of_obs', 'cat_type', 'X', 'Xpr_id', 'P', 'P_pr',
    }
    assert expected_keys.issubset(set(sd.keys()))
    assert int(sd['N_total']) == len(xi)
    assert int(sd['C']) == int(dit['item_type_id'].max())


def test_ordered_logit_eval_loglik_matches_pyro_namespace(ol_model):
    sd = ol_model.stan_data
    params = _synthetic_params(sd, _ol_param_shapes(sd), seed=10)
    ns = _ol_ns()
    direct = ns['get_log_likelihood_of_ordered_logit_ncats'](sd, params)
    by_method = ol_model.eval_loglik(sd, params)
    np.testing.assert_array_equal(np.asarray(direct), np.asarray(by_method))


def test_ordered_logit_eval_loglik_annealed_matches_pyro_namespace(ol_model):
    sd = ol_model.stan_data
    params = _synthetic_params(sd, _ol_param_shapes(sd), seed=11)
    params['temperature'] = jnp.asarray(0.5)
    ns = _ol_ns()
    direct = ns['get_log_likelihood_of_ordered_logit_ncats_with_annealing'](sd, params)
    by_method = ol_model.eval_loglik_annealed(sd, params)
    np.testing.assert_array_equal(np.asarray(direct), np.asarray(by_method))


def test_ordered_logit_eval_log_prior_matches_pyro_namespace(ol_model):
    sd = ol_model.stan_data
    params = _synthetic_params(sd, _ol_param_shapes(sd), seed=12)
    ns = _ol_ns()
    direct = ns['get_prior_of_ordered_logit_ncats'](params)
    by_method = ol_model.eval_log_prior(params)
    assert float(direct) == float(by_method)


# ---------------------------------------------------------------------------
# Cross-subclass sanity: both inherit the abstract contract via IRTModel
# ---------------------------------------------------------------------------


@pytest.mark.parametrize('cls', [CreditModel, OrderedLogit])
def test_subclass_passes_abc_contract(cls):
    """Constructor signature: (dit, dcati, x_formula='~ time - 1', *, seed=123)."""
    from model import Model
    from model_irt import IRTModel
    assert issubclass(cls, IRTModel)
    assert issubclass(cls, Model)
    import inspect
    sig = inspect.signature(cls.__init__)
    assert 'dit' in sig.parameters
    assert 'dcati' in sig.parameters
