"""
Parity test for the NumPyro ncats PCM REDUCTIONS (rho_by_cat_qu) against Stan.

Compares the reduced effect sizes between:
- Stan: src/stan/partial_credit_model_ncats_v261012.stan  (compute_rho generated quantities)
- NumPyro: src/numpyro/partial_credit_model_ncats_v261012.pyro  (get_rho_by_cat_qu)

Reductions are deterministic given the parameters, so -- exactly like the R
test_get_endpoints_equivalence.R -- we feed ONE identical parameter draw into both the
Stan generated-quantities block and the pyro get_rho_by_cat_qu, and assert the resulting
R x Qcmp rho matrices match per-draw (both _fit, full eta, and _pr, Xpr-only eta).

The manifest exercises the whole red_fn catalog (per-participant 1-4 + group-mean 5-9),
single items (Scmp=1 padded) and an item-group (two baseline questions pooled), and both
the E[y] (thr=0) and P(y>=thr) (thr=1) summary families.
"""

import importlib.util
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

import jax
import numpy as np
import pytest
from cmdstanpy import CmdStanModel
from numpyro.handlers import substitute, trace  # noqa: F401  (kept for parity-test parallels)


project_root = Path(__file__).parent.parent.parent
numpyro_path = str(project_root / 'src' / 'numpyro')
python_path = str(project_root / 'python')
for p in (numpyro_path, python_path):
    if p not in sys.path:
        sys.path.insert(0, p)

pyro_file = project_root / 'src' / 'numpyro' / 'partial_credit_model_ncats_v261012.pyro'
loader = SourceFileLoader("pcm_ncats_rho", str(pyro_file))
spec = importlib.util.spec_from_loader("pcm_ncats_rho", loader)
pcm_ncats_rho = importlib.util.module_from_spec(spec)
loader.exec_module(pcm_ncats_rho)


def _squeeze_draw(x):
    arr = np.asarray(x)
    if arr.ndim >= 1 and arr.shape[0] == 1:
        return arr[0]
    return arr


def _build_test_case():
    """Same observation layout as the base parity test, plus a reduction manifest."""
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

    # Force P_pr < P to exercise the _pr (Xpr-only) reduction recomputation.
    Xpr_id = np.array([1, 3], dtype=np.int32)

    # --- reduction manifest ---------------------------------------------------------------
    # 3 comparisons: single item 1-vs-2 in cat1, cat2, and an ITEM-GROUP {1,2}-vs-{3} in cat3.
    Scmp = 2
    Qcmp = 3
    cmp_cat_type = np.array([1, 2, 3], dtype=np.int32)
    question_cmp = np.array([
        [[1, 0], [2, 0]],   # cat1: baseline q1, endline q2
        [[1, 0], [2, 0]],   # cat2: baseline q1, endline q2
        [[1, 2], [3, 0]],   # cat3: baseline pooled {q1,q2}, endline q3
    ], dtype=np.int32)

    # Whole catalog: per-participant 1-4, group-mean 5-9.
    red_fn = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=np.int32)
    R = int(red_fn.size)
    # thr per (r, jc): mostly E[y] (0); exercise P(y>=1) on rel%(r=3) and endline level(r=6).
    cmp_reduce_thr = np.zeros((R, Qcmp), dtype=np.int32)
    cmp_reduce_thr[2, :] = 1   # red_fn=3 (per-participant rel%) -> P(y>=1)
    cmp_reduce_thr[5, :] = 1   # red_fn=6 (group endline level)  -> P(y>=1)

    stan_data = {
        'C': C, 'P': P, 'U': U, 'P_pr': int(Xpr_id.size),
        'N': N, 'Q': Q, 'K': K,
        'N_total': N_total, 'Q_total': Q_total,
        'y': y,
        'unit_of_obs': unit_of_obs,
        'question_of_obs': question_of_obs,
        'cat_type': cat_type,
        'X': X,
        'Xpr_id': Xpr_id,
        'R': R, 'Qcmp': Qcmp, 'Scmp': Scmp,
        'cmp_cat_type': cmp_cat_type,
        'question_cmp': question_cmp,
        'red_fn': red_fn,
        'cmp_reduce_thr': cmp_reduce_thr,
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
def test_numpyro_pcm_rho_parity_with_stan():
    """Pyro get_rho_by_cat_qu == Stan compute_rho rho_by_cat_qu_fit/_pr for an identical draw."""
    np.random.seed(123)
    jax.config.update('jax_platform_name', 'cpu')

    stan_data, stan_params = _build_test_case()

    stan_model = CmdStanModel(
        stan_file=str(project_root / 'src' / 'stan' / 'partial_credit_model_ncats_v261012.stan'),
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

    stan_rho_fit = _squeeze_draw(stan_fit.stan_variable('rho_by_cat_qu_fit'))  # (R, Qcmp)
    stan_rho_pr = _squeeze_draw(stan_fit.stan_variable('rho_by_cat_qu_pr'))

    rho_fit, rho_pr = pcm_ncats_rho.get_rho_by_cat_qu(stan_data, stan_params)

    assert stan_rho_fit.shape == rho_fit.shape == (stan_data['R'], stan_data['Qcmp'])
    assert stan_rho_pr.shape == rho_pr.shape == (stan_data['R'], stan_data['Qcmp'])

    # Pyro computes the per-obs probs in JAX float32 (CPU default) whereas Stan uses float64; the
    # ratio / 2^x reductions (rel%, GMFR) amplify that rounding, so parity holds to ~1e-5, not 1e-6.
    np.testing.assert_allclose(stan_rho_fit, rho_fit, rtol=1e-5, atol=1e-5)
    np.testing.assert_allclose(stan_rho_pr, rho_pr, rtol=1e-5, atol=1e-5)
