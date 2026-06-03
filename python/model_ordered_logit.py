"""
:class:`OrderedLogitNCats` -- the IRT ordered-logit subclass of
:class:`model._IRTModel`. Delegates every blueprint method to the
existing free-function pipeline in :mod:`fit_ordered_logit_model` and
the underlying ``.pyro`` namespace.

Same caveat as :mod:`model_credit`: the ordered-logit numpyro module
does not expose a standalone ``get_ordered_prob_of_ordered_logit_ncats``
callable, so :meth:`OrderedLogitNCats.eval_outcome_for_endpoint` raises
``NotImplementedError``. Endpoint summaries via :meth:`endpoints_per_draw`
work because they read the precomputed ``ordered_prob_by_cat_qu_fit``
posterior variable rather than re-evaluating it.
"""

from fit_ordered_logit_model import (
    _fit_ordered_logit_make_stan_data,
    eval_loglik_ordered_logit_ncats,
    eval_loglik_ordered_logit_ncats_with_annealing,
    get_prior_of_ordered_logit_ncats,
    fit_ordered_logit_model_ncats_pyrosvi,
    fit_ordered_logit_model_ncats_stanadvi,
    fit_ordered_logit_model_ncats_stanhmc,
)
from model import _IRTModel


def _ordered_logit_eval_outcome_not_implemented(data, params):
    raise NotImplementedError(
        "OrderedLogitNCats does not yet expose"
        " get_ordered_prob_of_ordered_logit_ncats. Add it to"
        " ordered_logit_ncats_v260413.pyro and update _OL_MODULE if"
        " moment-matching or SMC algorithms need to score shifted draws."
    )


_OL_MODULE = {
    'pyrosvi':              fit_ordered_logit_model_ncats_pyrosvi,
    'stanadvi':             fit_ordered_logit_model_ncats_stanadvi,
    'stanhmc':              fit_ordered_logit_model_ncats_stanhmc,
    'eval_loglik':          eval_loglik_ordered_logit_ncats,
    'eval_loglik_annealed': eval_loglik_ordered_logit_ncats_with_annealing,
    'eval_log_prior':             get_prior_of_ordered_logit_ncats,
    'eval_outcome':         _ordered_logit_eval_outcome_not_implemented,
}


class OrderedLogitNCats(_IRTModel):
    """Ordered-logit (n-categories) wrapper around
    :mod:`fit_ordered_logit_model`."""

    param_names = (
        'latent_factor_unit',
        'latent_factor_beta',
        'skill_thresholds_1',
        'skill_thresholds_incs',
        'loadings_questions_m1',
    )
    positive_params = ('loadings_questions_m1',)

    _module = _OL_MODULE
    _make_stan_data = staticmethod(_fit_ordered_logit_make_stan_data)
