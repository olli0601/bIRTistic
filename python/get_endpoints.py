"""
Back-compat shim. All endpoint logic now lives on
:class:`model_irt.IRTModel` as methods. This module exposes the same
free-function signatures the original ``get_endpoints.py`` had so
existing scripts and tests keep working.

Each free function builds a temporary :class:`PartialCreditModelNCats`
from the supplied ``dp1`` / ``dit`` (the IRT subclass chosen here is
arbitrary -- all three IRT subclasses share the same endpoint logic via
:class:`IRTModel`) and dispatches to the method.
"""

from typing import Literal, Optional

import pandas as pd


def get_endpoints(
    dp1: pd.DataFrame,
    dit: pd.DataFrame,
    draws_file: Optional[str] = None,
    draws=None,
    categorical_threshold: int = 3,
    endpoint_type: Literal["items", "item_groups"] = "items",
    param_name: str = "ordered_prob_by_cat_qu_pr",
    verbose: bool = True,
) -> pd.DataFrame:
    """Back-compat: dispatch to :meth:`IRTModel.get_endpoints`."""
    from model_pcm import PartialCreditModelNCats
    model = PartialCreditModelNCats(dit=dit, dcati=dp1)
    return model.get_endpoints(
        draws=draws, draws_file=draws_file,
        categorical_threshold=categorical_threshold,
        endpoint_type=endpoint_type,
        param_name=param_name, verbose=verbose,
    )


def get_endpoints_per_draw(
    dcati: pd.DataFrame,
    dit: pd.DataFrame,
    draws_file: Optional[str] = None,
    draws=None,
    categorical_threshold: int = 3,
    endpoint_type: Literal["items", "item_groups"] = "items",
    param_name: str = "ordered_prob_by_cat_qu_pr",
    verbose: bool = False,
) -> pd.DataFrame:
    """Back-compat: dispatch to :meth:`IRTModel.get_endpoints_per_draw`."""
    from model_pcm import PartialCreditModelNCats
    model = PartialCreditModelNCats(dit=dit, dcati=dcati)
    return model.get_endpoints_per_draw(
        draws=draws, draws_file=draws_file,
        categorical_threshold=categorical_threshold,
        endpoint_type=endpoint_type,
        param_name=param_name, verbose=verbose,
    )
