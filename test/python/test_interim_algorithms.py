"""
Step 4 of the OO-port refactor: each interim algorithm in
``fit_interim.py`` gains a model-aware form and a back-compat shim. The
two forms must agree bit-for-bit on the fixture cohort.

Algorithms covered (incrementally as step 4 progresses):
- ``fit_interim_IS_reweight``  (vs ``fit_interim_importance_sampling_of_posterior_xz_from_x``)
- TODO: ``fit_interim_IS_moment_matching``,
        ``fit_interim_SMC_resample``,
        ``fit_interim_SMC_PPS``,
        ``get_interim_endpt_and_w_from_poi``
"""

import sys
from pathlib import Path

import arviz as az
import pandas as pd
import pytest

_repo_root = Path(__file__).resolve().parents[2]
_python_dir = _repo_root / 'python'
if str(_python_dir) not in sys.path:
    sys.path.insert(0, str(_python_dir))

from fit_interim import (
    fit_interim_IS_reweight,
    fit_interim_importance_sampling_of_posterior_xz_from_x,
    fit_interim_IS_moment_matching,
    fit_interim_IS_moment_matching_of_posterior_xz_from_x,
    fit_interim_SMC_resample,
    fit_interim_SMC_resample_of_posterior_xz_from_x,
    fit_interim_SMC_PPS,
    fit_interim_SMC_PPS_of_posterior_xz_from_x,
    get_interim_endpt_and_w_from_poi,
    get_interim_z_from_ypredi,
)
from model_pcm import PartialCreditModelNCats

_TEST_DATA = _repo_root / 'test' / 'test_data'
_INTERIM_STEM = _TEST_DATA / 'pcm_1_interim_1'
_DRAWS_FILE = f"{_INTERIM_STEM}_draws.zarr"

_SEED = 123
_PPS_Z_TOTAL = 8
_INTERIM_M = 12


@pytest.fixture(scope='module')
def dit():
    return pd.read_csv(f"{_INTERIM_STEM}_data_dit.csv")


@pytest.fixture(scope='module')
def xi():
    return pd.read_csv(f"{_INTERIM_STEM}_data_dp1.csv")


@pytest.fixture(scope='module')
def draws():
    return az.from_zarr(_DRAWS_FILE)


@pytest.fixture(scope='module')
def zi(xi):
    return get_interim_z_from_ypredi(
        xi, _DRAWS_FILE, _INTERIM_M,
        pps_z_total=_PPS_Z_TOTAL, seed=_SEED,
    )


@pytest.fixture(scope='module')
def pcm_model(dit, xi):
    return PartialCreditModelNCats(dit=dit, dcati=xi,
                                   x_formula="~ time - 1", seed=_SEED)


# ---------------------------------------------------------------------------
# fit_interim_IS_reweight: model-aware vs shim
# ---------------------------------------------------------------------------


def test_IS_reweight_model_matches_shim(pcm_model, dit, xi, zi, draws):
    p_model, perf_model = fit_interim_IS_reweight(
        model=pcm_model, zi=zi, draws=draws,
        pps_z_total=_PPS_Z_TOTAL, categorical_threshold=2,
        save_to_file=False, verbose=False,
    )
    p_shim, perf_shim = fit_interim_importance_sampling_of_posterior_xz_from_x(
        xi=xi, zi=zi, dit=dit, draws=draws,
        pps_z_total=_PPS_Z_TOTAL, categorical_threshold=2,
        save_to_file=False, verbose=False,
    )
    pd.testing.assert_frame_equal(p_model, p_shim, check_exact=True)
    pd.testing.assert_frame_equal(perf_model, perf_shim, check_exact=True)


# ---------------------------------------------------------------------------
# fit_interim_IS_moment_matching: model-aware vs shim
# ---------------------------------------------------------------------------


def test_IS_moment_matching_model_matches_shim(pcm_model, dit, xi, zi, draws):
    p_model, mm_model = fit_interim_IS_moment_matching(
        model=pcm_model, zi=zi, draws=draws,
        fitting_method_args={},
        pps_z_total=_PPS_Z_TOTAL, categorical_threshold=2,
        save_to_file=False, verbose=False,
    )
    p_shim, mm_shim = fit_interim_IS_moment_matching_of_posterior_xz_from_x(
        xi=xi, zi=zi, dit=dit, draws=draws,
        fitting_method_args={},
        pps_z_total=_PPS_Z_TOTAL, categorical_threshold=2,
        save_to_file=False, verbose=False,
    )
    # mm has a non-deterministic 'mins' wall-time field; drop before comparing.
    pd.testing.assert_frame_equal(p_model, p_shim, check_exact=True)
    pd.testing.assert_frame_equal(
        mm_model.drop(columns=['mins']), mm_shim.drop(columns=['mins']),
        check_exact=True,
    )


# ---------------------------------------------------------------------------
# fit_interim_SMC_resample: model-aware vs shim
# ---------------------------------------------------------------------------


def test_SMC_resample_model_matches_shim(pcm_model, dit, xi, zi, draws):
    fma = {
        's_idx': 0, 'n_particles': 8, 'max_temps': 3,
        'n_move_steps': 2, 'seed': 123,
    }
    sched_model, parts_model = fit_interim_SMC_resample(
        model=pcm_model, zi=zi, draws=draws,
        fitting_method_args=fma,
        save_to_file=False, verbose=False,
    )
    sched_shim, parts_shim = fit_interim_SMC_resample_of_posterior_xz_from_x(
        xi=xi, zi=zi, dit=dit, draws=draws,
        fitting_method_args=fma,
        save_to_file=False, verbose=False,
    )
    # 'mins_total' and 'move_secs' columns are wall-clock; drop before comparing.
    drop_cols = ['mins_total', 'move_secs']
    pd.testing.assert_frame_equal(
        sched_model.drop(columns=drop_cols),
        sched_shim.drop(columns=drop_cols),
        check_exact=True,
    )
    assert set(parts_model) == set(parts_shim)
    import numpy as np
    for k in parts_model:
        np.testing.assert_array_equal(
            np.asarray(parts_model[k]), np.asarray(parts_shim[k]),
        )


# ---------------------------------------------------------------------------
# fit_interim_SMC_PPS: model-aware wrapper vs free function
# ---------------------------------------------------------------------------


def test_get_interim_endpt_and_w_model_matches_shim(pcm_model, dit, xi, draws):
    """The model-aware call path (model=...) must produce the same wa frame
    as the legacy positional (xi=..., dit=...) path."""
    kw = dict(
        draws=draws, draws_file=_DRAWS_FILE,
        interim_m=_INTERIM_M, pps_z_total=_PPS_Z_TOTAL,
        pps_H1_def=0.5, pps_ProbH1_thresh=0.89,
        categorical_threshold=2, seed=_SEED,
    )
    wa_model = get_interim_endpt_and_w_from_poi(model=pcm_model, **kw)
    wa_shim = get_interim_endpt_and_w_from_poi(xi=xi, dit=dit, **kw)
    pd.testing.assert_frame_equal(wa_model, wa_shim, check_exact=True)


def test_SMC_PPS_model_matches_shim(pcm_model, dit, xi, zi):
    """The model-aware wrapper unpacks model into the legacy free function;
    output must be bit-identical."""
    fma = {
        'n_particles': 8, 'max_temps': 3, 'n_move_steps': 2, 'seed': 123,
    }
    p_model, summ_model = fit_interim_SMC_PPS(
        model=pcm_model, zi=zi,
        fitting_method_args=fma, draws_file=_DRAWS_FILE,
        pps_z_total=2, categorical_threshold=2, cpu_n=1,
        save_to_file=False, verbose=False,
    )
    p_shim, summ_shim = fit_interim_SMC_PPS_of_posterior_xz_from_x(
        xi=xi, zi=zi, dit=dit,
        fitting_method_args=fma, draws_file=_DRAWS_FILE,
        pps_z_total=2, categorical_threshold=2, cpu_n=1,
        save_to_file=False, verbose=False,
    )
    pd.testing.assert_frame_equal(p_model, p_shim, check_exact=True)
    pd.testing.assert_frame_equal(
        summ_model.drop(columns=['mins']), summ_shim.drop(columns=['mins']),
        check_exact=True,
    )
