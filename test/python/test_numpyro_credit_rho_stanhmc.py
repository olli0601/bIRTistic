"""
Parity test for the NumPyro ncats CREDIT-model reductions (rho_by_cat_qu) against Stan.

- Stan: src/stan/credit_model_ncats_v261012.stan  (compute_rho generated quantities)
- NumPyro: src/numpyro/credit_model_ncats_v261012.pyro  (get_rho_by_cat_qu)

One identical parameter draw is fed into both the Stan generated-quantities block and the pyro
get_rho_by_cat_qu; the R x Qcmp rho matrices (_fit full eta, _pr Xpr-only eta) are asserted equal
per-draw. Same observation layout + reduction manifest as the PCM rho parity test; the credit model
uses softmax(eta) category probs (same parameter set as the PCM model).
"""

import importlib.util
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

import jax
import numpy as np
import pytest
from cmdstanpy import CmdStanModel


project_root = Path(__file__).parent.parent.parent
for p in (str(project_root / 'src' / 'numpyro'), str(project_root / 'python')):
    if p not in sys.path:
        sys.path.insert(0, p)

pyro_file = project_root / 'src' / 'numpyro' / 'credit_model_ncats_v261012.pyro'
loader = SourceFileLoader("credit_ncats_rho", str(pyro_file))
spec = importlib.util.spec_from_loader("credit_ncats_rho", loader)
credit_ncats_rho = importlib.util.module_from_spec(spec)
loader.exec_module(credit_ncats_rho)


def _squeeze_draw(x):
    arr = np.asarray(x)
    if arr.ndim >= 1 and arr.shape[0] == 1:
        return arr[0]
    return arr


def _build_test_case():
    C, P, U = 3, 3, 3
    N = np.array([4, 3, 5], dtype=np.int32)
    Q = np.array([2, 2, 3], dtype=np.int32)
    K = np.array([3, 4, 3], dtype=np.int32)
    N_total, Q_total = int(N.sum()), int(Q.sum())

    cat_type = np.concatenate([np.full(N[0], 1), np.full(N[1], 2), np.full(N[2], 3)]).astype(np.int32)
    question_of_obs = np.array([1, 2, 1, 2, 1, 2, 1, 1, 2, 3, 1, 2], dtype=np.int32)
    unit_of_obs = np.array([1, 2, 3, 1, 2, 3, 1, 1, 2, 3, 1, 2], dtype=np.int32)
    y = np.array([1, 2, 3, 1, 2, 4, 1, 3, 2, 1, 2, 3], dtype=np.int32)
    X = np.array([
        [0.5, -0.2, 1.0], [-0.1, 0.3, 0.2], [0.2, 0.1, -0.5], [0.0, -0.4, 0.7],
        [0.3, 0.2, -0.1], [-0.2, 0.1, 0.4], [0.1, -0.3, 0.6], [0.7, -0.1, 0.3],
        [0.2, 0.4, -0.2], [-0.3, 0.5, 0.1], [0.4, -0.2, -0.3], [0.6, 0.0, 0.2],
    ], dtype=float)
    Xpr_id = np.array([1, 3], dtype=np.int32)

    Scmp, Qcmp = 2, 3
    cmp_cat_type = np.array([1, 2, 3], dtype=np.int32)
    question_cmp = np.array([
        [[1, 0], [2, 0]],
        [[1, 0], [2, 0]],
        [[1, 2], [3, 0]],
    ], dtype=np.int32)
    red_fn = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=np.int32)
    R = int(red_fn.size)
    cmp_reduce_thr = np.zeros((R, Qcmp), dtype=np.int32)
    cmp_reduce_thr[2, :] = 1
    cmp_reduce_thr[5, :] = 1

    stan_data = {
        'C': C, 'P': P, 'U': U, 'P_pr': int(Xpr_id.size),
        'N': N, 'Q': Q, 'K': K, 'N_total': N_total, 'Q_total': Q_total,
        'y': y, 'unit_of_obs': unit_of_obs, 'question_of_obs': question_of_obs,
        'cat_type': cat_type, 'X': X, 'Xpr_id': Xpr_id,
        'R': R, 'Qcmp': Qcmp, 'Scmp': Scmp, 'cmp_cat_type': cmp_cat_type,
        'question_cmp': question_cmp, 'red_fn': red_fn, 'cmp_reduce_thr': cmp_reduce_thr,
    }

    skill_thresholds = np.array([
        -0.3, 0.4, 0.2, 0.9,
        -0.7, -0.1, 0.3, 0.0, 0.5, 1.1,
        -0.5, 0.2, -0.1, 0.6, 0.3, 1.0,
    ], dtype=float)
    loadings_questions_m1 = np.array([1.2, 0.8, 1.1, 0.9], dtype=float)
    latent_factor_unit = np.array([-0.6, 0.1, 0.5], dtype=float)
    latent_factor_beta = np.array([0.7, -0.4, 0.2], dtype=float)

    stan_params = {
        'latent_factor_unit': latent_factor_unit,
        'latent_factor_beta': latent_factor_beta,
        'skill_thresholds': skill_thresholds,
        'loadings_questions_m1': loadings_questions_m1,
    }
    return stan_data, stan_params


@pytest.mark.integration
def test_numpyro_credit_rho_parity_with_stan():
    np.random.seed(123)
    jax.config.update('jax_platform_name', 'cpu')
    stan_data, stan_params = _build_test_case()

    stan_model = CmdStanModel(
        stan_file=str(project_root / 'src' / 'stan' / 'credit_model_ncats_v261012.stan'),
        stanc_options={'include-paths': str(project_root / 'src' / 'stan')},
    )
    stan_fit = stan_model.sample(
        data=stan_data, inits=stan_params, chains=1, iter_warmup=0, iter_sampling=1,
        fixed_param=True, adapt_engaged=False, show_console=False,
    )
    stan_rho_fit = _squeeze_draw(stan_fit.stan_variable('rho_by_cat_qu_fit'))
    stan_rho_pr = _squeeze_draw(stan_fit.stan_variable('rho_by_cat_qu_pr'))

    rho_fit, rho_pr = credit_ncats_rho.get_rho_by_cat_qu(stan_data, stan_params)

    assert stan_rho_fit.shape == rho_fit.shape == (stan_data['R'], stan_data['Qcmp'])
    assert stan_rho_pr.shape == rho_pr.shape == (stan_data['R'], stan_data['Qcmp'])
    # Pyro computes the per-obs probs in JAX float32 (CPU default) whereas Stan uses float64; the
    # ratio / 2^x reductions (rel%, GMFR) amplify that rounding, so parity holds to ~1e-5, not 1e-6.
    np.testing.assert_allclose(stan_rho_fit, rho_fit, rtol=1e-5, atol=1e-5)
    np.testing.assert_allclose(stan_rho_pr, rho_pr, rtol=1e-5, atol=1e-5)
