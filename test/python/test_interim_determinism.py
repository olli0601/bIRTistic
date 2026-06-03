"""
Determinism tests for the interim helpers in python/fit_interim.py.

These pin down that any function carrying a ``seed`` argument (or a closed-over
RNG) produces bit-identical outputs across two invocations. Regressions to
catch include:
- non-deterministic shuffles inside ``_load_ypred`` / ``get_interim_z_from_ypredi``
- forgotten ``rng.choice`` calls in the wa construction
- accidental Python ``set`` / dict iteration order leaking into outputs
"""

import sys
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd
import pytest

_repo_root = Path(__file__).resolve().parents[2]
_python_dir = _repo_root / 'python'
if str(_python_dir) not in sys.path:
    sys.path.insert(0, str(_python_dir))

from interim_helpers import _load_ypred, get_interim_z_from_ypredi
from fit_interim import (
    get_interim_endpt_and_w_from_poi,
    fit_interim_regress_H1x_on_wz,
    fit_interim_regress_endptx_on_wz,
    fit_interim_IS_reweight,
)
from model_pcm import PartialCreditModelNCats

_TEST_DATA = _repo_root / 'test' / 'test_data'
_INTERIM_STEM = _TEST_DATA / 'pcm_1_interim_1'
_DRAWS_FILE = f"{_INTERIM_STEM}_draws.zarr"

_SEED = 123
_PPS_Z_TOTAL = 8
_INTERIM_M = 12


@pytest.fixture(scope='module')
def xi():
    return pd.read_csv(f"{_INTERIM_STEM}_data_dp1.csv")


@pytest.fixture(scope='module')
def dit():
    return pd.read_csv(f"{_INTERIM_STEM}_data_dit.csv")


@pytest.fixture(scope='module')
def draws():
    return az.from_zarr(_DRAWS_FILE)


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------


def test_load_ypred_random_subset_deterministic():
    a = _load_ypred(_DRAWS_FILE, _PPS_Z_TOTAL,
                    np.random.default_rng(_SEED), keep_order=False)
    b = _load_ypred(_DRAWS_FILE, _PPS_Z_TOTAL,
                    np.random.default_rng(_SEED), keep_order=False)
    np.testing.assert_array_equal(a, b)


def test_get_interim_z_from_ypredi_deterministic(xi):
    zi_a = get_interim_z_from_ypredi(
        xi, _DRAWS_FILE, _INTERIM_M, pps_z_total=_PPS_Z_TOTAL, seed=_SEED,
    )
    zi_b = get_interim_z_from_ypredi(
        xi, _DRAWS_FILE, _INTERIM_M, pps_z_total=_PPS_Z_TOTAL, seed=_SEED,
    )
    pd.testing.assert_frame_equal(zi_a, zi_b, check_dtype=True, check_exact=True)


def test_get_interim_z_from_ypredi_changes_with_seed(xi):
    zi_a = get_interim_z_from_ypredi(
        xi, _DRAWS_FILE, _INTERIM_M, pps_z_total=_PPS_Z_TOTAL, seed=_SEED,
    )
    zi_b = get_interim_z_from_ypredi(
        xi, _DRAWS_FILE, _INTERIM_M, pps_z_total=_PPS_Z_TOTAL, seed=_SEED + 1,
    )
    # Resampled src_pids should differ with a different seed (probabilistic
    # but with interim_m=12 and a different seed essentially certain).
    assert not zi_a['src_pid'].reset_index(drop=True).equals(
        zi_b['src_pid'].reset_index(drop=True)
    )


# ---------------------------------------------------------------------------
# wa construction
# ---------------------------------------------------------------------------


def _build_wa(xi, dit, draws, seed):
    model = PartialCreditModelNCats(dit=dit, dcati=xi,
                                    x_formula="~ time - 1", seed=seed)
    return get_interim_endpt_and_w_from_poi(
        model=model, draws=draws, draws_file=_DRAWS_FILE,
        interim_m=_INTERIM_M, pps_z_total=_PPS_Z_TOTAL,
        pps_H1_def=0.5, pps_ProbH1_thresh=0.89,
        categorical_threshold=2, seed=seed,
    )


def test_get_interim_endpt_and_w_deterministic(xi, dit, draws):
    a = _build_wa(xi, dit, draws, _SEED)
    b = _build_wa(xi, dit, draws, _SEED)
    pd.testing.assert_frame_equal(a, b, check_dtype=True, check_exact=True)


def test_get_interim_endpt_and_w_changes_with_seed(xi, dit, draws):
    a = _build_wa(xi, dit, draws, _SEED)
    b = _build_wa(xi, dit, draws, _SEED + 1)
    # w_baseline / w_endline depend on the resampled z block -> must differ.
    assert not a['w_baseline'].equals(b['w_baseline'])
    # But pps_ratio_x is computed from the x-posterior alone -> identical
    # across z resamples for the same draw index.
    merged = a[['item_label', 'draw', 'pps_ratio_x']].merge(
        b[['item_label', 'draw', 'pps_ratio_x']],
        on=['item_label', 'draw'], suffixes=('_a', '_b'),
    )
    np.testing.assert_array_equal(
        merged['pps_ratio_x_a'].to_numpy(), merged['pps_ratio_x_b'].to_numpy(),
    )


# ---------------------------------------------------------------------------
# Regression label estimators (pure functions of wa)
# ---------------------------------------------------------------------------


def test_fit_interim_regress_H1x_on_wz_deterministic(xi, dit, draws):
    wa = _build_wa(xi, dit, draws, _SEED)
    p_a, perf_a = fit_interim_regress_H1x_on_wz(wa)
    p_b, perf_b = fit_interim_regress_H1x_on_wz(wa)
    pd.testing.assert_frame_equal(p_a, p_b, check_dtype=True, check_exact=True)
    pd.testing.assert_frame_equal(perf_a, perf_b, check_dtype=True, check_exact=True)


def test_fit_interim_regress_endptx_on_wz_deterministic(xi, dit, draws):
    wa = _build_wa(xi, dit, draws, _SEED)
    p_a, perf_a = fit_interim_regress_endptx_on_wz(wa, pps_H1_def=0.5)
    p_b, perf_b = fit_interim_regress_endptx_on_wz(wa, pps_H1_def=0.5)
    pd.testing.assert_frame_equal(p_a, p_b, check_dtype=True, check_exact=True)
    pd.testing.assert_frame_equal(perf_a, perf_b, check_dtype=True, check_exact=True)


# ---------------------------------------------------------------------------
# IS reweighting
# ---------------------------------------------------------------------------


def test_fit_interim_IS_reweight_deterministic(xi, dit, draws):
    zi = get_interim_z_from_ypredi(
        xi, _DRAWS_FILE, _INTERIM_M, pps_z_total=_PPS_Z_TOTAL, seed=_SEED,
    )
    model = PartialCreditModelNCats(dit=dit, dcati=xi,
                                    x_formula="~ time - 1", seed=_SEED)
    p_a, _ = fit_interim_IS_reweight(
        model=model, zi=zi, draws=draws,
        pps_z_total=_PPS_Z_TOTAL, categorical_threshold=2,
        save_to_file=False, verbose=False,
    )
    p_b, _ = fit_interim_IS_reweight(
        model=model, zi=zi, draws=draws,
        pps_z_total=_PPS_Z_TOTAL, categorical_threshold=2,
        save_to_file=False, verbose=False,
    )
    pd.testing.assert_frame_equal(p_a, p_b, check_dtype=True, check_exact=False,
                                  atol=1e-12, rtol=0)
