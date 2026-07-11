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

from model import Model
from fit_interim import (
    fit_interim_regress_H1x_on_wz,
    fit_interim_regress_endptx_on_wz_with_Gaussian_approx,
    fit_interim_posterior_xz_from_z_with_IS_reweight,
)
from model_pcm import PartialCreditModel

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


def _make_model(xi, dit, seed):
    return PartialCreditModel(dit=dit, dcati=xi,
                              x_formula="~ time - 1", seed=seed)


def test_load_ypred_random_subset_deterministic():
    a = Model._load_ypred(_DRAWS_FILE, _PPS_Z_TOTAL,
                          np.random.default_rng(_SEED), keep_order=False)
    b = Model._load_ypred(_DRAWS_FILE, _PPS_Z_TOTAL,
                          np.random.default_rng(_SEED), keep_order=False)
    np.testing.assert_array_equal(a, b)


def test_get_interim_z_from_ypredi_deterministic(xi, dit):
    m = _make_model(xi, dit, _SEED)
    zi_a = m.get_interim_z_from_ypredi(
        _DRAWS_FILE, _INTERIM_M, pps_z_total=_PPS_Z_TOTAL, seed=_SEED,
    )
    zi_b = m.get_interim_z_from_ypredi(
        _DRAWS_FILE, _INTERIM_M, pps_z_total=_PPS_Z_TOTAL, seed=_SEED,
    )
    pd.testing.assert_frame_equal(zi_a, zi_b, check_dtype=True, check_exact=True)


def test_get_interim_z_from_ypredi_changes_with_seed(xi, dit):
    m = _make_model(xi, dit, _SEED)
    zi_a = m.get_interim_z_from_ypredi(
        _DRAWS_FILE, _INTERIM_M, pps_z_total=_PPS_Z_TOTAL, seed=_SEED,
    )
    zi_b = m.get_interim_z_from_ypredi(
        _DRAWS_FILE, _INTERIM_M, pps_z_total=_PPS_Z_TOTAL, seed=_SEED + 1,
    )
    assert not zi_a['src_pid'].reset_index(drop=True).equals(
        zi_b['src_pid'].reset_index(drop=True)
    )


# ---------------------------------------------------------------------------
# wa construction
# ---------------------------------------------------------------------------


def _build_wa(xi, dit, draws, seed):
    model = PartialCreditModel(dit=dit, dcati=xi,
                                    x_formula="~ time - 1", seed=seed)
    pps_H1_min_effect_size_thresh = 0.5
    pps_ProbH1_target_lwr_quantile = 0.89
    zi = model.get_interim_z_from_ypredi(
        _DRAWS_FILE, _INTERIM_M,
        pps_z_total=_PPS_Z_TOTAL, seed=seed, keep_order=True,
    )
    wa = model.get_w(zi, categorical_threshold=2)
    x_ratio = model.get_endpoints_per_draw(
        draws=draws, categorical_threshold=2, endpoint_type='items',
    )
    tmp = (
        model.get_p_h1(x_ratio, pps_H1_min_effect_size_thresh=pps_H1_min_effect_size_thresh)
        .rename(columns={'p_h1': 'pps_ProbH1_x'})
    )
    tmp['pps_H1_yes'] = (tmp['pps_ProbH1_x'] > pps_ProbH1_target_lwr_quantile).astype(int)
    x_ratio = x_ratio.rename(columns={'ratio': 'pps_ratio_x'})
    x_ratio['pps_H1_x'] = (x_ratio['pps_ratio_x'] > pps_H1_min_effect_size_thresh).astype(int)
    x_ratio = x_ratio.merge(
        tmp[['item_label', 'pps_ProbH1_x', 'pps_H1_yes']],
        on='item_label', how='left',
    )
    return wa.merge(
        x_ratio[['draw', 'item_label', 'item_type', 'pps_ratio_x', 'pps_H1_x']],
        on=['draw', 'item_label', 'item_type'], how='inner',
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
    p_a, perf_a = fit_interim_regress_endptx_on_wz_with_Gaussian_approx(wa, pps_H1_min_effect_size_thresh=0.5)
    p_b, perf_b = fit_interim_regress_endptx_on_wz_with_Gaussian_approx(wa, pps_H1_min_effect_size_thresh=0.5)
    pd.testing.assert_frame_equal(p_a, p_b, check_dtype=True, check_exact=True)
    pd.testing.assert_frame_equal(perf_a, perf_b, check_dtype=True, check_exact=True)


# ---------------------------------------------------------------------------
# IS reweighting
# ---------------------------------------------------------------------------


def test_fit_interim_IS_reweight_deterministic(xi, dit, draws):
    model = PartialCreditModel(dit=dit, dcati=xi,
                                    x_formula="~ time - 1", seed=_SEED,
                                    categorical_threshold=2)
    zi = model.get_interim_z_from_ypredi(
        _DRAWS_FILE, _INTERIM_M, pps_z_total=_PPS_Z_TOTAL, seed=_SEED,
    )
    is_kwargs = dict(
        interim_method_args={
            'pps_z_total': _PPS_Z_TOTAL,
            'pps_H1_min_effect_size_thresh': 0.5,
            'pps_ProbH1_target_lwr_quantile': 0.89,
            'output_file_prefix': None,
            'save_to_file': False,
            'verbose': False,
        },
    )
    p_a, _ = fit_interim_posterior_xz_from_z_with_IS_reweight(
        model, zi, draws=draws, **is_kwargs,
    )
    p_b, _ = fit_interim_posterior_xz_from_z_with_IS_reweight(
        model, zi, draws=draws, **is_kwargs,
    )
    pd.testing.assert_frame_equal(p_a, p_b, check_dtype=True, check_exact=False,
                                  atol=1e-12, rtol=0)
