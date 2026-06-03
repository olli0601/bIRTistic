"""
Draw-index alignment tests for the interim helpers in python/fit_interim.py.

These pin down two invariants the production scripts now rely on:

1. ``_load_ypred(..., keep_order=True)`` returns row ``s`` == posterior draw
   ``s`` (no shuffling). Regression test for the IS-from-x ρ-collapse bug.
2. ``get_interim_endpt_and_w_from_poi`` joins z-summaries and x-endpoint
   ratios on the same posterior-draw index: per (item, draw) the wa frame's
   ``pps_ratio_x`` equals the ``get_endpoints_per_draw`` ratio for the
   same draw index.
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

from interim_helpers import _load_ypred
from fit_interim import get_interim_endpt_and_w_from_poi
from model_pcm import PartialCreditModel

_TEST_DATA = _repo_root / 'test' / 'test_data'
_INTERIM_STEM = _TEST_DATA / 'pcm_1_interim_1'
_DRAWS_FILE = f"{_INTERIM_STEM}_draws.zarr"


@pytest.fixture(scope='module')
def xi():
    return pd.read_csv(f"{_INTERIM_STEM}_data_dp1.csv")


@pytest.fixture(scope='module')
def dit():
    return pd.read_csv(f"{_INTERIM_STEM}_data_dit.csv")


@pytest.fixture(scope='module')
def draws():
    return az.from_zarr(_DRAWS_FILE)


@pytest.fixture(scope='module')
def raw_ypred():
    """The (n_draw, N_total) ypred array in posterior order, flattened over
    (chain, draw). Source of truth for the alignment tests."""
    yp = az.from_zarr(_DRAWS_FILE).posterior['ypred'].values
    return yp.reshape(-1, yp.shape[-1])


# ---------------------------------------------------------------------------
# _load_ypred keep_order behaviour
# ---------------------------------------------------------------------------


def test_load_ypred_keep_order_returns_posterior_prefix(raw_ypred):
    S = 8
    rng = np.random.default_rng(0)
    out = _load_ypred(_DRAWS_FILE, S, rng, keep_order=True)
    assert out.shape == (S, raw_ypred.shape[1])
    np.testing.assert_array_equal(out, raw_ypred[:S])


def test_load_ypred_keep_order_independent_of_rng(raw_ypred):
    """keep_order=True must ignore the rng — different seeds give the same
    rows (the first S posterior draws)."""
    S = 5
    a = _load_ypred(_DRAWS_FILE, S, np.random.default_rng(1), keep_order=True)
    b = _load_ypred(_DRAWS_FILE, S, np.random.default_rng(999), keep_order=True)
    np.testing.assert_array_equal(a, b)
    np.testing.assert_array_equal(a, raw_ypred[:S])


def test_load_ypred_random_subset_uses_rng(raw_ypred):
    """keep_order=False permutes via rng.choice without replacement: rows are
    a SUBSET of raw_ypred, and different seeds give different selections."""
    S = 10
    a = _load_ypred(_DRAWS_FILE, S, np.random.default_rng(1), keep_order=False)
    b = _load_ypred(_DRAWS_FILE, S, np.random.default_rng(2), keep_order=False)
    assert a.shape == (S, raw_ypred.shape[1])
    # Each row of a must appear somewhere in raw_ypred (no fabricated rows).
    for row in a:
        assert any(np.array_equal(row, raw_ypred[k]) for k in range(raw_ypred.shape[0]))
    # Different seeds yield different selections (probabilistic, but with
    # S=10 out of n_draw collisions are astronomically unlikely).
    assert not np.array_equal(a, b)


# ---------------------------------------------------------------------------
# get_interim_endpt_and_w_from_poi draw-index alignment
# ---------------------------------------------------------------------------


@pytest.fixture(scope='module')
def wa(xi, dit, draws):
    model = PartialCreditModel(dit=dit, dcati=xi,
                                    x_formula="~ time - 1", seed=123)
    return get_interim_endpt_and_w_from_poi(
        model=model, draws=draws, draws_file=_DRAWS_FILE,
        interim_m=10, pps_z_total=20,
        pps_H1_def=0.5, pps_ProbH1_thresh=0.89,
        categorical_threshold=2, seed=123,
    )


def test_wa_pps_ratio_x_matches_get_endpoints_per_draw(wa, xi, dit, draws):
    """For every (item, draw) the wa frame's pps_ratio_x must equal the
    independently-computed get_endpoints_per_draw ratio at the same draw
    index. Catches any permutation of the posterior-draw axis between the
    two pipelines."""
    _ref_model = PartialCreditModel(dit=dit, dcati=xi)
    x_ratio = _ref_model.get_endpoints_per_draw(
        draws=draws, categorical_threshold=2,
        endpoint_type='items', param_name='ordered_prob_by_cat_qu_fit',
        verbose=False,
    ).rename(columns={'ratio': 'pps_ratio_x_truth'})

    merged = wa[['item_label', 'draw', 'pps_ratio_x']].merge(
        x_ratio[['item_label', 'draw', 'pps_ratio_x_truth']],
        on=['item_label', 'draw'], how='inner',
    )
    assert len(merged) == len(wa), "wa has rows missing from get_endpoints_per_draw"
    np.testing.assert_allclose(
        merged['pps_ratio_x'].to_numpy(),
        merged['pps_ratio_x_truth'].to_numpy(),
        rtol=0, atol=0,
    )


def test_wa_draw_index_is_dense_posterior_prefix(wa):
    """Every draw in [0, pps_z_total) must appear once per item; no shuffled
    or skipped indices."""
    S = 20
    counts_per_item = wa.groupby('item_label')['draw'].apply(
        lambda d: sorted(d.tolist())
    )
    expected = list(range(S))
    for item, draws_seen in counts_per_item.items():
        assert draws_seen == expected, f"item {item} has draws {draws_seen}"
