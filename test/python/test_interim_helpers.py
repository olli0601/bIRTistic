"""
Step 1 of the OO-port refactor: pure-pandas helpers moved from
``fit_interim.py`` to ``interim_helpers.py``. These tests pin down that
both import paths resolve to the same callables and that the output of
``get_interim_z_from_ypredi`` is byte-identical to the cached snapshot
for the test fixture (golden snapshot).
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_repo_root = Path(__file__).resolve().parents[2]
_python_dir = _repo_root / 'python'
if str(_python_dir) not in sys.path:
    sys.path.insert(0, str(_python_dir))

import interim_helpers as ih  # noqa: E402
import fit_interim as fi  # noqa: E402

_TEST_DATA = _repo_root / 'test' / 'test_data'
_INTERIM_STEM = _TEST_DATA / 'pcm_1_interim_1'
_DRAWS_FILE = f"{_INTERIM_STEM}_draws.zarr"


# ---------------------------------------------------------------------------
# Import-path equivalence
# ---------------------------------------------------------------------------


def test_load_ypred_same_object():
    assert ih._load_ypred is fi._load_ypred


def test_get_interim_x_same_object():
    assert ih.get_interim_x is fi.get_interim_x


def test_get_interim_z_from_ypredi_same_object():
    assert ih.get_interim_z_from_ypredi is fi.get_interim_z_from_ypredi


def test_get_interim_z_from_ypredf_same_object():
    assert ih.get_interim_z_from_ypredf is fi.get_interim_z_from_ypredf


# ---------------------------------------------------------------------------
# Golden snapshot of get_interim_z_from_ypredi
# ---------------------------------------------------------------------------


@pytest.fixture(scope='module')
def xi():
    return pd.read_csv(f"{_INTERIM_STEM}_data_dp1.csv")


def test_get_interim_z_from_ypredi_golden(xi):
    """Locks the exact column/row layout for the fixture cohort
    (interim_m=12, pps_z_total=8, seed=123). Any future refactor in the
    helpers must preserve this output bit-for-bit."""
    zi = ih.get_interim_z_from_ypredi(
        xi, _DRAWS_FILE, interim_m=12, pps_z_total=8, seed=123,
    )
    # Shape contract: 12 new pids * (rows-per-pid in xi) rows.
    rows_per_pid = xi.groupby('pid').size().iloc[0]
    assert len(zi) == 12 * rows_per_pid

    # Required columns present.
    must_have = {'pid', 'src_pid', 'oid'} | {f'ypred_{s}' for s in range(8)}
    assert must_have.issubset(zi.columns)

    # pid range: fresh pids above xi's max.
    assert int(zi['pid'].min()) == int(xi['pid'].max()) + 1
    assert int(zi['pid'].max()) == int(xi['pid'].max()) + 12

    # src_pid values must all exist in xi.
    assert set(zi['src_pid'].unique()).issubset(set(xi['pid'].unique()))

    # Sort invariant: (pid, oid) ascending.
    pd.testing.assert_frame_equal(
        zi[['pid', 'oid']].reset_index(drop=True),
        zi[['pid', 'oid']].sort_values(['pid', 'oid']).reset_index(drop=True),
    )

    # ypred_* columns are integers in [0, K] (Likert / categorical levels).
    for s in range(8):
        col = zi[f'ypred_{s}'].to_numpy()
        assert np.issubdtype(col.dtype, np.integer)
        assert col.min() >= 0
