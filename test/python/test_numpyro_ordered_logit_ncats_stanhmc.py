"""
Parity tests for NumPyro ncats ordered logit model against Stan.

Compares deterministic generated quantities between:
- Stan: src/stan/ordered_logit_ncats_v260413.stan
- NumPyro: src/numpyro/ordered_logit_ncats_v260413.pyro

The test uses fixed parameters and validates both values and flattened formats.
"""

import importlib.util
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np
import pytest
from cmdstanpy import CmdStanModel
import numpyro
from numpyro.handlers import seed, substitute, trace


project_root = Path(__file__).parent.parent.parent
numpyro_path = str(project_root / 'src' / 'numpyro')
if numpyro_path not in sys.path:
    sys.path.insert(0, numpyro_path)

pyro_file = project_root / 'src' / 'numpyro' / 'ordered_logit_ncats_v260413.pyro'
loader = SourceFileLoader("ol_ncats_stanhmc", str(pyro_file))
spec = importlib.util.spec_from_loader("ol_ncats_stanhmc", loader)
ol_ncats_stanhmc = importlib.util.module_from_spec(spec)
loader.exec_module(ol_ncats_stanhmc)


def _squeeze_draw(x):
    arr = np.asarray(x)
    if arr.shape[0] == 1:
        return arr[0]
    return arr


def _build_test_case():
    # Three category types with heterogeneous K (3 and 4) and Q.
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

    question_of_obs = np.array([
        1, 2, 1, 2,
        1, 2, 1,
        1, 2, 3, 1, 2,
    ], dtype=np.int32)

    unit_of_obs = np.array([
        1, 2, 3, 1,
        2, 3, 1,
        1, 2, 3, 1, 2,
    ], dtype=np.int32)

    y = np.array([
        1, 2, 3, 1,
        2, 4, 1,
        3, 2, 1, 2, 3,
    ], dtype=np.int32)

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

    stan_data = {
        'C': C,
        'P': P,
        'U': U,
        'P_pr': int(Xpr_id.size),
        'N': N,
        'Q': Q,
        'K': K,
        'N_total': N_total,
        'Q_total': Q_total,
        'y': y,
        'unit_of_obs': unit_of_obs,
        'question_of_obs': question_of_obs,
        'cat_type': cat_type,
        'X': X,
        'Xpr_id': Xpr_id,
    }

    # skill_thresholds_1: one per question, Q_total = 7.
    skill_thresholds_1 = np.array([
        -0.3, 0.4,   # c1: Q=2
        -0.5, 0.1,   # c2: Q=2
        -0.6, 0.2, 0.5,  # c3: Q=3
    ], dtype=float)

    # skill_thresholds_incs: positive increments, len = sum(Q*(K-2)).
    # c1: 2*(3-2)=2, c2: 2*(4-2)=4, c3: 3*(3-2)=3 -> total 9.
    skill_thresholds_incs = np.array([
        0.6, 0.7,           # c1
        0.5, 0.6, 0.3, 0.7,  # c2
        0.4, 0.5, 0.8,       # c3
    ], dtype=float)

    loadings_questions_m1 = np.array([
        1.2,        # c1
        0.8,        # c2
        1.1, 0.9,   # c3
    ], dtype=float)

    latent_factor_unit = np.array([-0.6, 0.1, 0.5], dtype=float)
    latent_factor_beta = np.array([0.7, -0.4, 0.2], dtype=float)

    stan_params = {
        'latent_factor_unit': latent_factor_unit,
        'latent_factor_beta': latent_factor_beta,
        'skill_thresholds_1': skill_thresholds_1,
        'skill_thresholds_incs': skill_thresholds_incs,
        'loadings_questions_m1': loadings_questions_m1,
    }

    numpyro_substitute = {
        'latent_factor_unit': latent_factor_unit,
        'latent_factor_beta': latent_factor_beta,
        'skill_thresholds_1': skill_thresholds_1,
        'skill_thresholds_incs': skill_thresholds_incs,
        'loadings_questions_m1': loadings_questions_m1,
        'ypred_raw': np.zeros((N_total,), dtype=np.int32),
    }

    return stan_data, stan_params, numpyro_substitute


@pytest.mark.integration
def test_numpyro_ordered_logit_ncats_stanhmc_parity_with_stan_generated_quantities():
    """Parity for deterministic generated quantities and schema-compatible outputs."""
    np.random.seed(123)
    jax.config.update('jax_platform_name', 'cpu')

    stan_data, stan_params, numpyro_sub = _build_test_case()

    stan_model = CmdStanModel(
        stan_file=str(project_root / 'src' / 'stan' / 'ordered_logit_ncats_v260413.stan'),
        stanc_options={'include-paths': str(project_root / 'src' / 'stan')},
    )

    stan_fit = stan_model.sample(
        data=stan_data,
        inits=stan_params,
        chains=1,
        iter_warmup=0,
        iter_sampling=1,
        fixed_param=True,
        adapt_engaged=False,
        show_console=False,
    )

    stan_log_lik = _squeeze_draw(stan_fit.stan_variable('log_lik'))
    stan_fit_prob = _squeeze_draw(stan_fit.stan_variable('ordered_prob_by_cat_qu_fit'))
    stan_pr_prob = _squeeze_draw(stan_fit.stan_variable('ordered_prob_by_cat_qu_pr'))
    stan_brier = _squeeze_draw(stan_fit.stan_variable('ordinal_brier_score'))
    stan_ypred = _squeeze_draw(stan_fit.stan_variable('ypred'))
    stan_cutpoints = _squeeze_draw(stan_fit.stan_variable('cutpoints'))

    tr = trace(substitute(ol_ncats_stanhmc.ordered_logit_ncats, data=numpyro_sub)).get_trace(stan_data)

    npy_log_lik = np.asarray(tr['log_lik']['value'])
    npy_fit_prob = np.asarray(tr['ordered_prob_by_cat_qu_fit']['value'])
    npy_pr_prob = np.asarray(tr['ordered_prob_by_cat_qu_pr']['value'])
    npy_brier = np.asarray(tr['ordinal_brier_score']['value'])
    npy_ypred = np.asarray(tr['ypred']['value'])
    npy_cutpoints = np.asarray(tr['cutpoints']['value'])

    assert stan_log_lik.shape == npy_log_lik.shape
    assert stan_fit_prob.shape == npy_fit_prob.shape
    assert stan_pr_prob.shape == npy_pr_prob.shape
    assert stan_brier.shape == npy_brier.shape
    assert stan_ypred.shape == npy_ypred.shape
    assert stan_cutpoints.shape == npy_cutpoints.shape

    np.testing.assert_allclose(stan_cutpoints, npy_cutpoints, rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(stan_log_lik, npy_log_lik, rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(stan_fit_prob, npy_fit_prob, rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(stan_pr_prob, npy_pr_prob, rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(stan_brier, npy_brier, rtol=1e-6, atol=1e-6)

    assert np.issubdtype(npy_ypred.dtype, np.integer)
    assert np.all(npy_ypred >= 1)

    N = stan_data['N']
    K = stan_data['K']
    start = 0
    for c, Nc in enumerate(N):
        end = start + int(Nc)
        assert np.all(npy_ypred[start:end] <= int(K[c]))
        start = end
