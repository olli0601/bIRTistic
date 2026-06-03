"""
:class:`CreditModelNCats` -- the IRT credit-model subclass of
:class:`model._IRTModel`. Delegates every blueprint method to the
existing free-function pipeline in :mod:`fit_credit_model` and the
underlying ``.pyro`` namespace, so behaviour matches the pre-refactor
code bit-for-bit.

Note: the credit-model numpyro module does NOT expose a standalone
``get_ordered_prob_of_credit_model_ncats`` callable (the ordered
probabilities are recorded as ``ordered_prob_by_cat_qu_fit`` during the
fit and consumed by :func:`get_endpoints_per_draw`). Calling
``model.eval_outcome_for_endpoint(...)`` therefore raises
``NotImplementedError`` -- the moment-match and SMC algorithms that need
it currently target PCM only. If we extend MM / SMC to credit, expose
the function in the .pyro module and remove the override here.
"""

from fit_credit_model import (
    _fit_credit_make_stan_data,
    eval_loglik_credit_model_ncats,
    eval_loglik_credit_model_ncats_with_annealing,
    get_prior_of_credit_model_ncats,
    fit_credit_model_ncats_pyrosvi,
    fit_credit_model_ncats_stanadvi,
    fit_credit_model_ncats_stanhmc,
)
from model import _IRTModel


def _credit_eval_outcome_not_implemented(data, params):
    raise NotImplementedError(
        "CreditModelNCats does not yet expose get_ordered_prob_of_credit_model_ncats."
        " Add it to credit_model_ncats_v260413.pyro and update _CREDIT_MODULE if"
        " moment-matching or SMC algorithms need to score shifted draws."
    )


_CREDIT_MODULE = {
    'pyrosvi':              fit_credit_model_ncats_pyrosvi,
    'stanadvi':             fit_credit_model_ncats_stanadvi,
    'stanhmc':              fit_credit_model_ncats_stanhmc,
    'eval_loglik':          eval_loglik_credit_model_ncats,
    'eval_loglik_annealed': eval_loglik_credit_model_ncats_with_annealing,
    'logprior':             get_prior_of_credit_model_ncats,
    'eval_outcome':         _credit_eval_outcome_not_implemented,
}


class CreditModelNCats(_IRTModel):
    """Credit model (n-categories) wrapper around :mod:`fit_credit_model`."""

    param_names = (
        'latent_factor_unit',
        'latent_factor_beta',
        'skill_thresholds',
        'loadings_questions_m1',
    )
    positive_params = ('loadings_questions_m1',)

    _module = _CREDIT_MODULE
    _make_stan_data = staticmethod(_fit_credit_make_stan_data)
