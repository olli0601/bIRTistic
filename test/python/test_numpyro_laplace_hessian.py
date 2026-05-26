"""
Hessian-rank tests for the three NumPyro ncats models under
``AutoLaplaceApproximation``.

For each model (partial credit / credit / ordered logit), this test:
1. Builds a small synthetic dataset using the same shape conventions as the
   Stan parity tests in this directory.
2. Runs a short SVI optimization with ``AutoLaplaceApproximation`` to land
   near the MAP.
3. Asks the guide for its posterior, which under Laplace is
   ``MultivariateNormal(loc=MAP, scale_tril=chol(inv(Hessian)))``.
4. Reconstructs the Hessian as ``inv(scale_tril @ scale_tril.T)`` (i.e. the
   precision matrix of the Laplace approximation, which equals the Hessian
   of the negative log joint at the MAP in the unconstrained space).
5. Asserts that every eigenvalue of that Hessian is strictly positive
   (above a small numerical tolerance) — i.e. no flat directions / no
   zero eigenvalues / no rank deficiency that would signal a
   non-identifiable model.
"""

import importlib.util
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path
from functools import partial

import jax
import jax.numpy as jnp
import numpy as np
import pytest
import numpyro
from numpyro.infer import SVI, Trace_ELBO
from numpyro.infer.autoguide import AutoLaplaceApproximation


project_root = Path(__file__).parent.parent.parent
numpyro_path = str(project_root / 'src' / 'numpyro')
if numpyro_path not in sys.path:
    sys.path.insert(0, numpyro_path)


def _load_pyro_module(pyro_path: Path, module_name: str):
    loader = SourceFileLoader(module_name, str(pyro_path))
    spec = importlib.util.spec_from_loader(module_name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def _build_pcm_credit_test_case():
    """Shared dataset for partial credit + credit models (categorical likelihood)."""
    C = 3
    P = 3
    U = 3

    N = np.array([4, 3, 5], dtype=np.int32)
    Q = np.array([2, 2, 3], dtype=np.int32)
    K = np.array([3, 4, 3], dtype=np.int32)

    N_total = int(N.sum())
    Q_total = int(Q.sum())

    cat_type = np.concatenate([
        np.full(N[0], 1, dtype=np.int32),
        np.full(N[1], 2, dtype=np.int32),
        np.full(N[2], 3, dtype=np.int32),
    ])
    question_of_obs = np.array([1, 2, 1, 2, 1, 2, 1, 1, 2, 3, 1, 2], dtype=np.int32)
    unit_of_obs = np.array([1, 2, 3, 1, 2, 3, 1, 1, 2, 3, 1, 2], dtype=np.int32)
    y = np.array([1, 2, 3, 1, 2, 4, 1, 3, 2, 1, 2, 3], dtype=np.int32)

    X = np.array([
        [0.5, -0.2, 1.0],
        [-0.1, 0.3, 0.2],
        [0.2, 0.1, -0.5],
        [0.0, -0.4, 0.7],
        [0.3, 0.2, -0.1],
        [-0.2, 0.1, 0.4],
        [0.1, -0.3, 0.6],
        [0.7, -0.1, 0.3],
        [0.2, 0.4, -0.2],
        [-0.3, 0.5, 0.1],
        [0.4, -0.2, -0.3],
        [0.6, 0.0, 0.2],
    ], dtype=float)
    Xpr_id = np.array([1, 3], dtype=np.int32)

    return {
        'C': C, 'P': P, 'U': U,
        'P_pr': int(Xpr_id.size),
        'N': N, 'Q': Q, 'K': K,
        'N_total': N_total, 'Q_total': Q_total,
        'y': y,
        'unit_of_obs': unit_of_obs,
        'question_of_obs': question_of_obs,
        'cat_type': cat_type,
        'X': X,
        'Xpr_id': Xpr_id,
    }


def _build_ordered_logit_test_case():
    """Dataset for the ordered logit model. Same shape as the PCM/credit case."""
    return _build_pcm_credit_test_case()


def _hessian_min_eigval(model_callable, stan_data, num_steps=400, lr=0.05, seed=0):
    """Train AutoLaplaceApproximation, extract Hessian = precision, return min eigenvalue."""
    # Run on CPU for stable Hessian extraction.
    jax.config.update('jax_platform_name', 'cpu')

    # Lean SVI variant: skip discrete ypred and expensive generated quantities.
    model_for_svi = partial(model_callable, sample_ypred=False, compute_generated=False)

    guide = AutoLaplaceApproximation(model_for_svi)
    svi = SVI(model_for_svi, guide, numpyro.optim.Adam(lr), Trace_ELBO())

    rng_key = jax.random.PRNGKey(seed)
    rng_key, subkey = jax.random.split(rng_key)
    init_state = svi.init(subkey, stan_data)
    rng_key, subkey = jax.random.split(rng_key)
    run_result = svi.run(subkey, num_steps, stan_data, init_state=init_state, progress_bar=False)
    params = svi.get_params(run_result.state)

    # AutoLaplaceApproximation.get_posterior returns MultivariateNormal whose
    # covariance equals inv(Hessian) at the MAP. Precision = Hessian.
    # In numpyro 0.21 the model args used during init are baked into the guide,
    # so get_posterior takes only the learned params.
    posterior = guide.get_posterior(params)
    scale_tril = np.asarray(posterior.scale_tril)
    cov = scale_tril @ scale_tril.T

    # Hessian eigenvalues are 1 / cov eigenvalues. We just need the minimum.
    cov_eigvals = np.linalg.eigvalsh(cov)
    cov_eigvals = np.sort(cov_eigvals)
    # Guard against numerical noise producing tiny negative eigenvalues.
    cov_max = float(cov_eigvals.max())
    if cov_max <= 0.0:
        raise AssertionError("Covariance has no positive eigenvalues — degenerate guide.")
    hessian_min_eigval = 1.0 / cov_max  # smallest Hessian eigenvalue = 1 / largest cov eigenvalue
    return hessian_min_eigval, cov_eigvals


@pytest.mark.integration
def test_partial_credit_model_hessian_has_nonzero_eigenvalues():
    """Hessian of pcm_ncats under AutoLaplaceApproximation has no zero eigenvalues."""
    pyro_path = project_root / 'src' / 'numpyro' / 'partial_credit_model_ncats_v260413.pyro'
    mod = _load_pyro_module(pyro_path, 'pcm_ncats_laplace')
    data = _build_pcm_credit_test_case()
    min_eig, cov_eigvals = _hessian_min_eigval(mod.partial_credit_model_ncats, data)
    assert min_eig > 1e-6, (
        f"Partial credit Hessian smallest eigenvalue {min_eig:.3e} below tolerance "
        f"(largest covariance eigenvalue {cov_eigvals.max():.3e}). Model likely "
        f"non-identifiable in some direction."
    )


@pytest.mark.integration
def test_credit_model_hessian_has_nonzero_eigenvalues():
    """Hessian of credit_model_ncats under AutoLaplaceApproximation has no zero eigenvalues."""
    pyro_path = project_root / 'src' / 'numpyro' / 'credit_model_ncats_v260413.pyro'
    mod = _load_pyro_module(pyro_path, 'cm_ncats_laplace')
    data = _build_pcm_credit_test_case()
    min_eig, cov_eigvals = _hessian_min_eigval(mod.credit_model_ncats, data)
    assert min_eig > 1e-6, (
        f"Credit model Hessian smallest eigenvalue {min_eig:.3e} below tolerance "
        f"(largest covariance eigenvalue {cov_eigvals.max():.3e}). Model likely "
        f"non-identifiable in some direction."
    )


@pytest.mark.integration
def test_ordered_logit_model_hessian_has_nonzero_eigenvalues():
    """Hessian of ordered_logit_ncats under AutoLaplaceApproximation has no zero eigenvalues."""
    pyro_path = project_root / 'src' / 'numpyro' / 'ordered_logit_ncats_v260413.pyro'
    mod = _load_pyro_module(pyro_path, 'ol_ncats_laplace')
    data = _build_ordered_logit_test_case()
    min_eig, cov_eigvals = _hessian_min_eigval(mod.ordered_logit_ncats, data)
    assert min_eig > 1e-6, (
        f"Ordered logit Hessian smallest eigenvalue {min_eig:.3e} below tolerance "
        f"(largest covariance eigenvalue {cov_eigvals.max():.3e}). Model likely "
        f"non-identifiable in some direction."
    )
