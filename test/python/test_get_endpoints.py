"""Tests for python/get_endpoints.py."""

import sys
from pathlib import Path

import pandas as pd
import pytest

_repo_root = Path(__file__).resolve().parents[2]
_python_dir = _repo_root / 'python'
if str(_python_dir) not in sys.path:
    sys.path.insert(0, str(_python_dir))

from get_endpoints import get_endpoints

_TEST_DATA = _repo_root / 'test' / 'test_data'
_INTERIM_STEM = _TEST_DATA / 'pcm_1_interim_1'


@pytest.fixture(scope='module')
def interim_inputs():
    dcati = pd.read_csv(f"{_INTERIM_STEM}_data_dp1.csv")
    dit = pd.read_csv(f"{_INTERIM_STEM}_data_dit.csv")
    zarr_path = f"{_INTERIM_STEM}_draws.zarr"
    return dcati, dit, zarr_path


@pytest.fixture(scope='module')
def reference_endpoints():
    return pd.read_csv(f"{_INTERIM_STEM}_endpoints_items.csv")


def test_get_endpoints_matches_reference(interim_inputs, reference_endpoints):
    dcati, dit, zarr_path = interim_inputs
    pos = get_endpoints(
        dp1=dcati,
        dit=dit,
        draws_file=zarr_path,
        categorical_threshold=2,
        endpoint_type='items',
        param_name='ordered_prob_by_cat_qu_fit',
    )

    sort_cols = ['item_type', 'item_label', 'variable']
    actual = pos.sort_values(sort_cols).reset_index(drop=True)
    expected = reference_endpoints.sort_values(sort_cols).reset_index(drop=True)

    assert list(actual.columns) == list(expected.columns)
    assert len(actual) == len(expected)

    pd.testing.assert_frame_equal(
        actual, expected,
        check_dtype=False,
        check_exact=False,
        atol=1e-8,
        rtol=1e-6,
    )
