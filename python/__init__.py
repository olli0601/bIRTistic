"""
bIRTistic: Bayesian Item Response Theory with NumPyro and JAX

Python implementation of bIRTistic with support for:
- HMC (via NumPyro and cmdstanpy)
- SVI (via NumPyro)
- Amortised inference (via JAX/Flax)

Ported from the original R/Stan implementation.
"""

__version__ = "0.1.0"
__author__ = "bIRTistic Team"

from .data_loading import read_data_colombia, read_data_ukraine
from .fit_partial_credit_model import (
    fit_partial_credit_model_ncats_stanhmc,
    fit_partial_credit_model_ncats_stanadvi,
)
from .fit_ordered_logit_model import (
    fit_ordered_logit_model_ncats_stanhmc,
    fit_ordered_logit_model_ncats_stanadvi,
)
from .get_endpoints import get_endpoints

__all__ = [
    'read_data_colombia',
    'read_data_ukraine',
    'fit_partial_credit_model_ncats_stanhmc',
    'fit_ordered_logit_model_ncats_stanhmc',
    'fit_ordered_logit_model_ncats_stanadvi',
    'fit_partial_credit_model_ncats_stanadvi',
    'get_endpoints'
]

__all__ = ['read_data_colombia', 'read_data_ukraine']

# Package will be populated as modules are ported from R
