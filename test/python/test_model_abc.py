"""
The ``Model`` ABC + the ``IRTModel`` mixin. No concrete subclass yet --
these tests only verify the abstract contract.
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

_repo_root = Path(__file__).resolve().parents[2]
_python_dir = _repo_root / 'python'
if str(_python_dir) not in sys.path:
    sys.path.insert(0, str(_python_dir))

from model import Model  # noqa: E402
from model_irt import IRTModel  # noqa: E402


def test_model_is_abstract():
    """Cannot instantiate the bare base class -- every abstract method
    must be overridden."""
    with pytest.raises(TypeError):
        Model(dit=pd.DataFrame(), dcati=pd.DataFrame())


def test_irt_mixin_is_abstract():
    """IRTModel is still abstract until a concrete IRT subclass supplies
    the remaining abstract methods (make_stan_data, eval_loglik, etc.)."""
    with pytest.raises(TypeError):
        IRTModel(dit=pd.DataFrame(), dcati=pd.DataFrame())


def test_model_param_names_default_empty():
    assert Model.param_names == ()
    assert Model.positive_params == ()


def test_model_fit_closed_form_raises_by_default():
    """Default behaviour: NotImplementedError. Binomial subclass overrides."""

    class _Stub(Model):
        # Provide every abstract method with a no-op so we can instantiate.
        def make_stan_data(self, dcati, x_formula): return {}
        def fit_pyro_svi(self, output_file_prefix, **kw): return {}
        def fit_stan_svi(self, output_file_prefix, **kw): return {}
        def fit_stan_hmc(self, output_file_prefix, **kw): return {}
        def eval_loglik(self, data, params): return None
        def eval_loglik_annealed(self, data, params): return None
        def eval_log_prior(self, params): return None
        def eval_outcome_for_endpoint(self, data, params): return None
        def get_endpoints_per_draw(self, draws, categorical_threshold,
                                   endpoint_type='items'): return None
        def get_w(self, zi): return None

    stub = _Stub(dit=pd.DataFrame(), dcati=pd.DataFrame())
    with pytest.raises(NotImplementedError):
        stub.fit_closed_form('out_prefix')


def test_get_stacked_posterior_uses_param_names():
    """Base-class :meth:`get_stacked_posterior` reads from
    ``self.param_names`` -- subclasses only need to declare the tuple."""
    import numpy as np
    import xarray as xr

    class _Stub(Model):
        param_names = ('alpha',)
        def make_stan_data(self, dcati, x_formula): return {}
        def fit_pyro_svi(self, output_file_prefix, **kw): return {}
        def fit_stan_svi(self, output_file_prefix, **kw): return {}
        def fit_stan_hmc(self, output_file_prefix, **kw): return {}
        def eval_loglik(self, data, params): return None
        def eval_loglik_annealed(self, data, params): return None
        def eval_log_prior(self, params): return None
        def eval_outcome_for_endpoint(self, data, params): return None
        def get_endpoints_per_draw(self, draws, categorical_threshold,
                                   endpoint_type='items'): return None
        def get_w(self, zi): return None

    # Build a tiny fake arviz idata-like object.
    arr = np.arange(12, dtype=float).reshape(2, 3, 2)   # (chain, draw, d)
    ds = xr.Dataset({'alpha': (('chain', 'draw', 'd'), arr)})
    class _FakeDraws:
        posterior = ds

    stub = _Stub(dit=pd.DataFrame(), dcati=pd.DataFrame())
    theta = stub.get_stacked_posterior(_FakeDraws())
    assert set(theta) == {'alpha'}
    assert theta['alpha'].shape == (6, 2)
