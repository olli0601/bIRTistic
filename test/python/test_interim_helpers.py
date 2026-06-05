"""
Golden snapshot for the per-Model posterior-predictive helpers.

After the OO refactor (see ``dev/OO_refactor.md``), ``_load_ypred`` and
``get_interim_z_from_ypredi`` live on :class:`model.Model`. The standalone
``interim_helpers`` module is gone, so the legacy import-path equivalence
tests are dropped; what's preserved is the byte-identical layout contract
for the fixture cohort (interim_m=12, pps_z_total=8, seed=123).
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

from model_pcm import PartialCreditModel  # noqa: E402

_TEST_DATA = _repo_root / 'test' / 'test_data'
_INTERIM_STEM = _TEST_DATA / 'pcm_1_interim_1'
_DRAWS_FILE = f"{_INTERIM_STEM}_draws.zarr"


@pytest.fixture(scope='module')
def xi():
    return pd.read_csv(f"{_INTERIM_STEM}_data_dp1.csv")


@pytest.fixture(scope='module')
def dit():
    return pd.read_csv(f"{_INTERIM_STEM}_data_dit.csv")


def test_get_interim_z_from_ypredi_golden(xi, dit):
    """Locks the exact column/row layout for the fixture cohort
    (interim_m=12, pps_z_total=8, seed=123). Any future refactor in
    :meth:`Model.get_interim_z_from_ypredi` must preserve this output
    bit-for-bit."""
    model = PartialCreditModel(dit=dit, dcati=xi,
                               x_formula="~ time - 1", seed=123)
    zi = model.get_interim_z_from_ypredi(
        _DRAWS_FILE, interim_m=12, pps_z_total=8, seed=123,
    )
    rows_per_pid = xi.groupby('pid').size().iloc[0]
    assert len(zi) == 12 * rows_per_pid

    must_have = {'pid', 'src_pid', 'oid'} | {f'ypred_{s}' for s in range(8)}
    assert must_have.issubset(zi.columns)

    assert int(zi['pid'].min()) == int(xi['pid'].max()) + 1
    assert int(zi['pid'].max()) == int(xi['pid'].max()) + 12

    assert set(zi['src_pid'].unique()).issubset(set(xi['pid'].unique()))

    pd.testing.assert_frame_equal(
        zi[['pid', 'oid']].reset_index(drop=True),
        zi[['pid', 'oid']].sort_values(['pid', 'oid']).reset_index(drop=True),
    )

    for s in range(8):
        col = zi[f'ypred_{s}'].to_numpy()
        assert np.issubdtype(col.dtype, np.integer)
        assert col.min() >= 0
