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
from .model_pcm import PartialCreditModel
from .model_credit import CreditModel
from .model_ordered_logit import OrderedLogit
from .model_binomial import BinomialModel

__all__ = [
    'read_data_colombia',
    'read_data_ukraine',
    'PartialCreditModel',
    'CreditModel',
    'OrderedLogit',
    'BinomialModel',
]

# Package will be populated as modules are ported from R
