"""
:class:`PartialCreditModelNCats` -- the IRT partial-credit model subclass of
:class:`model._IRTModel`. Currently delegates every blueprint method to
the existing module-level functions in :mod:`fit_partial_credit_model`
and the underlying ``.pyro`` namespace, so the runtime behaviour matches
the pre-refactor code bit-for-bit.

The actual numpyro program / Stan code is unchanged; only the call site
changes from free-function imports to method dispatch via the Model
instance.
"""

from fit_partial_credit_model import (
    _fit_partial_credit_make_stan_data,
    eval_loglik_partial_credit_model_ncats,
    eval_loglik_partial_credit_model_ncats_with_annealing,
    get_prior_of_partial_credit_model_ncats,
    get_ordered_prob_of_partial_credit_model_ncats,
    fit_partial_credit_model_ncats_pyrosvi,
    fit_partial_credit_model_ncats_stanadvi,
    fit_partial_credit_model_ncats_stanhmc,
)
from model import _IRTModel

# ``_module`` is a flat dict keyed by the generic names the _IRTModel
# mixin expects (pyrosvi / stanadvi / stanhmc / eval_loglik / ...). It is
# evaluated once at import time -- the underlying functions are stable.
_PCM_MODULE = {
    'pyrosvi':            fit_partial_credit_model_ncats_pyrosvi,
    'stanadvi':           fit_partial_credit_model_ncats_stanadvi,
    'stanhmc':            fit_partial_credit_model_ncats_stanhmc,
    'eval_loglik':        eval_loglik_partial_credit_model_ncats,
    'eval_loglik_annealed': eval_loglik_partial_credit_model_ncats_with_annealing,
    'logprior':           get_prior_of_partial_credit_model_ncats,
    'eval_outcome':       get_ordered_prob_of_partial_credit_model_ncats,
}


class PartialCreditModelNCats(_IRTModel):
    """Partial credit model (n-categories) wrapper around the existing
    free-function pipeline in :mod:`fit_partial_credit_model`."""

    param_names = (
        'latent_factor_unit',
        'latent_factor_beta',
        'skill_thresholds',
        'skill_thresholds_1',
        'skill_thresholds_incs',
        'loadings_questions_m1',
    )
    positive_params = ('loadings_questions_m1',)

    _module = _PCM_MODULE
    _make_stan_data = staticmethod(_fit_partial_credit_make_stan_data)
