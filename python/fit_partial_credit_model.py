"""
Back-compat shim for :func:`fit_interim.fit_interim_MC_of_posterior_xz`'s
``fitting_method=`` multiprocessing default. All real PCM logic now lives
in :mod:`model_pcm` as methods on :class:`PartialCreditModelNCats`.

The MC path predates the OO refactor and ships its ``fitting_method``
callable into spawn workers (which re-import this module). To keep that
path working without refactoring its picklable callable interface, the
free function below adapts the legacy ``(dit, dcati, ...)`` signature
into a fresh :class:`PartialCreditModelNCats` and dispatches
``fit_pyro_svi``.
"""

from typing import Dict

import pandas as pd

from model_pcm import PartialCreditModelNCats


def fit_partial_credit_model_ncats_pyrosvi(
    dit: pd.DataFrame,
    dcati: pd.DataFrame,
    output_file_prefix: str,
    x_formula: str = "~ time - 1",
    seed: int = 123,
    **kwargs,
) -> Dict:
    """Shim: instantiate :class:`PartialCreditModelNCats` and forward to
    ``fit_pyro_svi``. Used as the multiprocessing ``fitting_method=``
    default in :func:`fit_interim.fit_interim_MC_of_posterior_xz`."""
    model = PartialCreditModelNCats(
        dit=dit, dcati=dcati, x_formula=x_formula, seed=seed,
    )
    return model.fit_pyro_svi(output_file_prefix=output_file_prefix, **kwargs)
