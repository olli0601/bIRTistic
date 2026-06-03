"""
Step 2 of the OO-port refactor: the ``Model`` ABC and the ``_IRTModel``
mixin. No concrete subclass yet -- these tests only verify the abstract
contract.
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

_repo_root = Path(__file__).resolve().parents[2]
_python_dir = _repo_root / 'python'
if str(_python_dir) not in sys.path:
    sys.path.insert(0, str(_python_dir))

from model import Model, _IRTModel  # noqa: E402


def test_model_is_abstract():
    """Cannot instantiate the bare base class -- every abstract method
    must be overridden."""
    with pytest.raises(TypeError):
        Model(dit=pd.DataFrame(), dcati=pd.DataFrame())


def test_irt_mixin_is_abstract():
    """The mixin is still abstract until ``_module`` / ``_make_stan_data``
    are supplied by a concrete IRT subclass."""
    with pytest.raises(TypeError):
        _IRTModel(dit=pd.DataFrame(), dcati=pd.DataFrame())


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
        def logprior(self, params): return None
        def eval_outcome_for_endpoint(self, data, params): return None
        def endpoints_per_draw(self, dcati, draws, categorical_threshold,
                               endpoint_type='items'): return None

    stub = _Stub(dit=pd.DataFrame(), dcati=pd.DataFrame())
    with pytest.raises(NotImplementedError):
        stub.fit_closed_form('out_prefix')


def test_stack_posterior_theta_uses_param_names():
    """Base-class :meth:`stack_posterior_theta` reads from
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
        def logprior(self, params): return None
        def eval_outcome_for_endpoint(self, data, params): return None
        def endpoints_per_draw(self, dcati, draws, categorical_threshold,
                               endpoint_type='items'): return None

    # Build a tiny fake arviz idata-like object.
    arr = np.arange(12, dtype=float).reshape(2, 3, 2)   # (chain, draw, d)
    ds = xr.Dataset({'alpha': (('chain', 'draw', 'd'), arr)})
    class _FakeDraws:
        posterior = ds

    stub = _Stub(dit=pd.DataFrame(), dcati=pd.DataFrame())
    theta = stub.stack_posterior_theta(_FakeDraws())
    assert set(theta) == {'alpha'}
    assert theta['alpha'].shape == (6, 2)
