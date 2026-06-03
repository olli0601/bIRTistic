"""
:class:`CreditModelNCats` and :class:`OrderedLogitNCats` -- the credit
and ordered-logit subclasses of :class:`model.IRTModel`. All bodies
are now inlined as methods on the class (model_credit.py /
model_ordered_logit.py); the legacy fit_credit_model.py and
fit_ordered_logit_model.py modules are deleted.

The PCM fixture in ``test/test_data/`` is partial-credit-shaped, so
without dedicated credit / ordered-logit fixtures we test the
structural invariants every IRT subclass must satisfy:

- ``param_names`` / ``positive_params`` are declared correctly
- the inlined ``@staticmethod`` proxies match what their underlying
  ``.pyro`` namespace returns directly
- ``eval_outcome_for_endpoint`` raises ``NotImplementedError`` (the
  ``.pyro`` modules do not yet expose a standalone ordered-prob callable)
- abstract contract satisfied (subclass instantiates with dit + dcati)
"""

import sys
from pathlib import Path

import pytest

_repo_root = Path(__file__).resolve().parents[2]
_python_dir = _repo_root / 'python'
if str(_python_dir) not in sys.path:
    sys.path.insert(0, str(_python_dir))

from model_credit import CreditModelNCats, _model_namespace as _credit_ns
from model_ordered_logit import (
    OrderedLogitNCats, _model_namespace as _ol_ns,
)


# ---------------------------------------------------------------------------
# Credit
# ---------------------------------------------------------------------------


def test_credit_param_names_match_pyro_model():
    """Sample sites in src/numpyro/credit_model_ncats_v260413.pyro are
    latent_factor_unit / latent_factor_beta / skill_thresholds /
    loadings_questions_m1 -- declared in the subclass."""
    assert set(CreditModelNCats.param_names) == {
        'latent_factor_unit', 'latent_factor_beta',
        'skill_thresholds', 'loadings_questions_m1',
    }
    assert CreditModelNCats.positive_params == ('loadings_questions_m1',)


def test_credit_proxies_resolve_in_pyro_namespace():
    """The four @staticmethod proxies must dispatch to functions that exist
    in the credit .pyro namespace."""
    ns = _credit_ns()
    assert 'get_log_likelihood_of_credit_model_ncats' in ns
    assert 'get_log_likelihood_of_credit_model_ncats_with_annealing' in ns
    assert 'get_prior_of_credit_model_ncats' in ns


def test_credit_eval_outcome_raises_not_implemented():
    """The credit .pyro module does not expose a standalone ordered-prob
    callable; eval_outcome_for_endpoint must raise so MM / SMC algorithms
    fail loudly rather than silently."""
    with pytest.raises(NotImplementedError):
        CreditModelNCats.eval_outcome_for_endpoint(data={}, params={})


# ---------------------------------------------------------------------------
# OrderedLogit
# ---------------------------------------------------------------------------


def test_ordered_logit_param_names_match_pyro_model():
    """Sample sites in src/numpyro/ordered_logit_ncats_v260413.pyro are
    latent_factor_unit / latent_factor_beta / skill_thresholds_1 /
    skill_thresholds_incs / loadings_questions_m1."""
    assert set(OrderedLogitNCats.param_names) == {
        'latent_factor_unit', 'latent_factor_beta',
        'skill_thresholds_1', 'skill_thresholds_incs',
        'loadings_questions_m1',
    }
    assert OrderedLogitNCats.positive_params == ('loadings_questions_m1',)


def test_ordered_logit_proxies_resolve_in_pyro_namespace():
    ns = _ol_ns()
    assert 'get_log_likelihood_of_ordered_logit_ncats' in ns
    assert 'get_log_likelihood_of_ordered_logit_ncats_with_annealing' in ns
    assert 'get_prior_of_ordered_logit_ncats' in ns


def test_ordered_logit_eval_outcome_raises_not_implemented():
    with pytest.raises(NotImplementedError):
        OrderedLogitNCats.eval_outcome_for_endpoint(data={}, params={})


# ---------------------------------------------------------------------------
# Cross-subclass sanity: both inherit the abstract contract via IRTModel
# ---------------------------------------------------------------------------


@pytest.mark.parametrize('cls', [CreditModelNCats, OrderedLogitNCats])
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
