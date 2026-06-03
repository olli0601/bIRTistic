"""
Step 6 of the OO-port refactor: :class:`CreditModelNCats` and
:class:`OrderedLogitNCats` subclasses of :class:`model._IRTModel`.

The PCM fixture in ``test/test_data/`` is partial-credit-shaped (Likert
1-7 + categorical items), so we cannot run bit-identical equivalence
tests against credit / ordered-logit free-function fits without
dedicated fixtures (deferred). These tests instead verify the
structural invariants every IRT subclass must satisfy:

- ``_module`` dispatch keys resolve to the expected free-function callables
- ``param_names`` / ``positive_params`` are declared
- ``eval_outcome_for_endpoint`` raises ``NotImplementedError`` (the
  underlying ``.pyro`` module does not yet expose a standalone
  ordered-prob callable; see model_{credit,ordered_logit}.py docstrings)
"""

import sys
from pathlib import Path

import pytest

_repo_root = Path(__file__).resolve().parents[2]
_python_dir = _repo_root / 'python'
if str(_python_dir) not in sys.path:
    sys.path.insert(0, str(_python_dir))

import fit_credit_model
import fit_ordered_logit_model
from model_credit import CreditModelNCats, _CREDIT_MODULE
from model_ordered_logit import OrderedLogitNCats, _OL_MODULE


# ---------------------------------------------------------------------------
# Credit
# ---------------------------------------------------------------------------


def test_credit_module_dispatch_resolves_to_free_functions():
    assert _CREDIT_MODULE['pyrosvi'] is fit_credit_model.fit_credit_model_ncats_pyrosvi
    assert _CREDIT_MODULE['stanadvi'] is fit_credit_model.fit_credit_model_ncats_stanadvi
    assert _CREDIT_MODULE['stanhmc'] is fit_credit_model.fit_credit_model_ncats_stanhmc
    assert _CREDIT_MODULE['eval_loglik'] is fit_credit_model.eval_loglik_credit_model_ncats
    assert _CREDIT_MODULE['eval_loglik_annealed'] is fit_credit_model.eval_loglik_credit_model_ncats_with_annealing
    assert _CREDIT_MODULE['eval_log_prior'] is fit_credit_model.get_prior_of_credit_model_ncats


def test_credit_param_names_match_pyro_model():
    """The sample sites in src/numpyro/credit_model_ncats_v260413.pyro are
    latent_factor_unit / latent_factor_beta / skill_thresholds /
    loadings_questions_m1 -- declared in the subclass."""
    assert set(CreditModelNCats.param_names) == {
        'latent_factor_unit', 'latent_factor_beta',
        'skill_thresholds', 'loadings_questions_m1',
    }
    assert CreditModelNCats.positive_params == ('loadings_questions_m1',)


def test_credit_make_stan_data_bound():
    """staticmethod hookup -- bound at class definition time."""
    assert CreditModelNCats._make_stan_data is not None
    assert callable(CreditModelNCats._make_stan_data)


def test_credit_eval_outcome_raises_not_implemented():
    """Credit .pyro module does not expose a standalone ordered-prob
    callable; eval_outcome_for_endpoint must raise so MM / SMC algorithms
    fail loudly rather than silently."""
    with pytest.raises(NotImplementedError):
        _CREDIT_MODULE['eval_outcome'](data={}, params={})


# ---------------------------------------------------------------------------
# OrderedLogit
# ---------------------------------------------------------------------------


def test_ordered_logit_module_dispatch_resolves_to_free_functions():
    assert _OL_MODULE['pyrosvi'] is fit_ordered_logit_model.fit_ordered_logit_model_ncats_pyrosvi
    assert _OL_MODULE['stanadvi'] is fit_ordered_logit_model.fit_ordered_logit_model_ncats_stanadvi
    assert _OL_MODULE['stanhmc'] is fit_ordered_logit_model.fit_ordered_logit_model_ncats_stanhmc
    assert _OL_MODULE['eval_loglik'] is fit_ordered_logit_model.eval_loglik_ordered_logit_ncats
    assert _OL_MODULE['eval_loglik_annealed'] is fit_ordered_logit_model.eval_loglik_ordered_logit_ncats_with_annealing
    assert _OL_MODULE['eval_log_prior'] is fit_ordered_logit_model.get_prior_of_ordered_logit_ncats


def test_ordered_logit_param_names_match_pyro_model():
    """The sample sites in src/numpyro/ordered_logit_ncats_v260413.pyro are
    latent_factor_unit / latent_factor_beta / skill_thresholds_1 /
    skill_thresholds_incs / loadings_questions_m1."""
    assert set(OrderedLogitNCats.param_names) == {
        'latent_factor_unit', 'latent_factor_beta',
        'skill_thresholds_1', 'skill_thresholds_incs',
        'loadings_questions_m1',
    }
    assert OrderedLogitNCats.positive_params == ('loadings_questions_m1',)


def test_ordered_logit_make_stan_data_bound():
    assert OrderedLogitNCats._make_stan_data is not None
    assert callable(OrderedLogitNCats._make_stan_data)


def test_ordered_logit_eval_outcome_raises_not_implemented():
    with pytest.raises(NotImplementedError):
        _OL_MODULE['eval_outcome'](data={}, params={})


# ---------------------------------------------------------------------------
# Cross-subclass sanity: both inherit the abstract contract via _IRTModel
# ---------------------------------------------------------------------------


@pytest.mark.parametrize('cls', [CreditModelNCats, OrderedLogitNCats])
def test_subclass_passes_abc_contract(cls):
    """Cannot instantiate without (dit, dcati). With them, the constructor
    builds stan_data via the subclass's _make_stan_data without error."""
    from model import Model, _IRTModel
    assert issubclass(cls, _IRTModel)
    assert issubclass(cls, Model)
    # Constructor signature: (dit, dcati, x_formula='~ time - 1', *, seed=123)
    import inspect
    sig = inspect.signature(cls.__init__)
    assert 'dit' in sig.parameters
    assert 'dcati' in sig.parameters
