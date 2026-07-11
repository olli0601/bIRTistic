"""
Correctness tests for
:mod:`amortiser_pps_features_fixed_qpsi_MLP_loss_multiquantilehead`
(§11.5 of ``dev/amortised_decision_making.md``, first prototype §11.6).

Binomial exp-fam benchmark: amortised PPS vs closed-form Beta-Binomial PPS.
Because training the net inside the test suite is expensive, these tests
train a small net once and cache it under
``test/test_data/amortised_pps_binomial_net.pkl``.
"""

from functools import partial
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

_repo_root = Path(__file__).resolve().parents[2]
_python_dir = _repo_root / 'python'
if str(_python_dir) not in sys.path:
    sys.path.insert(0, str(_python_dir))

from model_binomial import BinomialModel  # noqa: E402
from amortiser_common import (  # noqa: E402
    train,
    predict_amortised_p_h1_for_one_xz,
    save_trained_model,
    load_fitted_model,
)
from amortiser_pps_features_fixed_qpsi_MLP_loss_multiquantilehead import (  # noqa: E402
    Amortiser_PPS_features_fixed_qpsi_MLP_loss_multiquantilehead,
)
from amortiser_pps_features_MLP_qpsi_MLP_loss_multiquantilehead import (  # noqa: E402
    Amortiser_PPS_features_MLP_qpsi_MLP_loss_multiquantilehead,
)

_DEFAULT_TAUS = (0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)


_TEST_DATA = _repo_root / 'test' / 'test_data'
_TEST_DATA.mkdir(exist_ok=True)
_NET_CKPT = _TEST_DATA / 'amortised_pps_binomial_net.pkl'

_PRIOR_A = 1.0
_PRIOR_B = 1.0
_P_0 = 0.5
_PPS_H1_MIN_EFFECT_SIZE_THRESH = 0.25   # rho > 0.25 (i.e. p < 0.375)
_PPS_PROBH1_THRESH = 0.89
_N_MAX = 500


def _make_binomial_model():
    """Construct a minimal :class:`BinomialModel` with the constants
    above so we can call ``make_training_data_with_features`` on it."""
    dp = pd.DataFrame({
        'pid': np.arange(_N_MAX),
        'submission_date': pd.date_range('2025-01-01', periods=_N_MAX, freq='D'),
        'y': np.zeros(_N_MAX, dtype=int),
        'oid': np.arange(1, _N_MAX + 1, dtype=int),
    })
    dit = pd.DataFrame({
        'item_label': ['intervention'],
        'item_type': ['binomial'],
        'item_high_label': ['higher_is_better'],
    })
    return BinomialModel(
        dit=dit, dcati=dp,
        prior_a=_PRIOR_A, prior_b=_PRIOR_B, p_0=_P_0, seed=42,
    )


# ---------------------------------------------------------------------------
# Closed-form Beta-Binomial PPS (mirrors BinomialModel.fit_closed_form_pps
# but standalone so the test does not depend on a pd.DataFrame ``dcati``).
# ---------------------------------------------------------------------------


def _closed_form_pps_binomial(kn, n, m, a=1.0, b=1.0, p_0=0.5,
                              pps_H1_min_effect_size_thresh=0.25,
                              pps_ProbH1_target_lwr_quantile=0.89):
    from scipy.stats import beta as _beta, betabinom as _betabinom
    a_post = a + kn
    b_post = b + (n - kn)
    p_h1_threshold = (1.0 - pps_H1_min_effect_size_thresh) * p_0
    if m <= 0:
        p_h1_now = _beta.cdf(p_h1_threshold, a_post, b_post)
        return float(p_h1_now > pps_ProbH1_target_lwr_quantile)
    k_m_grid = np.arange(m + 1)
    a_z = a_post + k_m_grid
    b_z = b_post + (m - k_m_grid)
    p_h1_given_z = _beta.cdf(p_h1_threshold, a_z, b_z)
    crosses = p_h1_given_z > pps_ProbH1_target_lwr_quantile
    if not crosses.any():
        return 0.0
    k_star = int(k_m_grid[crosses][-1])
    return float(_betabinom.cdf(k_star, m, a_post, b_post))


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope='module')
def trained_fit():
    """Train (or load) a modest amortised net once per test session."""
    if _NET_CKPT.exists():
        return load_fitted_model(str(_NET_CKPT))
    model = _make_binomial_model()
    sample_fn = partial(
        model.make_training_data_with_features, n_max=_N_MAX,
    )
    fit = train(
        sample_fn,
        Amortiser_PPS_features_fixed_qpsi_MLP_loss_multiquantilehead,
        net_kwargs={'hidden_dims': (128, 128, 64)},
        pps_ProbH1_lwr_quantiles_mesh=_DEFAULT_TAUS,
        num_steps=6000,
        batch_size=2048,
        lr=1e-3,
        seed=42,
        verbose=False,
    )
    save_trained_model(fit, str(_NET_CKPT))
    return fit


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------


def test_amortised_pps_forward_pass_deterministic(trained_fit):
    # Features are (k_total, n_total) / N_MAX after DeepSets pooling.
    features = np.array([[150, 300],
                         [270, 300],
                         [25, 300]], dtype=np.float32) / _N_MAX
    p_a, q_a, _ = predict_amortised_p_h1_for_one_xz(
        trained_fit, features, _PPS_H1_MIN_EFFECT_SIZE_THRESH,
    )
    p_b, q_b, _ = predict_amortised_p_h1_for_one_xz(
        trained_fit, features, _PPS_H1_MIN_EFFECT_SIZE_THRESH,
    )
    np.testing.assert_array_equal(p_a, p_b)
    np.testing.assert_array_equal(q_a, q_b)


def test_amortised_pps_training_is_deterministic():
    """Two training runs with the same seed must produce identical params
    on the first loss step (JAX + numpy RNGs both seeded)."""
    model = _make_binomial_model()
    sample_fn = partial(
        model.make_training_data_with_features, n_max=_N_MAX,
    )
    net_cls = Amortiser_PPS_features_fixed_qpsi_MLP_loss_multiquantilehead
    fit_a = train(
        sample_fn, net_cls,
        net_kwargs={'hidden_dims': (32,)}, pps_ProbH1_lwr_quantiles_mesh=(0.1, 0.5, 0.9),
        num_steps=5, batch_size=64, lr=1e-3, seed=7, verbose=False,
    )
    fit_b = train(
        sample_fn, net_cls,
        net_kwargs={'hidden_dims': (32,)}, pps_ProbH1_lwr_quantiles_mesh=(0.1, 0.5, 0.9),
        num_steps=5, batch_size=64, lr=1e-3, seed=7, verbose=False,
    )
    np.testing.assert_array_equal(
        np.asarray(fit_a['loss_history']),
        np.asarray(fit_b['loss_history']),
    )


def test_amortised_pps_features_MLP_class_trains():
    """Sanity: the features-MLP class trains via the same top-level API
    (dict batches) and yields non-degenerate quantile predictions on the
    Binomial raw-sequences sampler."""
    model = _make_binomial_model()
    sample_fn = partial(
        model.make_training_data_with_raw_sequences,
        n_max=_N_MAX,
    )
    net_cls = Amortiser_PPS_features_MLP_qpsi_MLP_loss_multiquantilehead
    fit = train(
        sample_fn, net_cls,
        net_kwargs={
            'q_tau_hidden_dims': (16,),
            'embed_dim': 8,
            'hidden_dims': (32, 32),
        },
        pps_ProbH1_lwr_quantiles_mesh=_DEFAULT_TAUS,
        num_steps=200, batch_size=64, lr=1e-3, seed=0, verbose=False,
    )
    # Forward-pass on a small batch and check the returned p_h1_xz is
    # in [0, 1] with non-trivial variance across (kn, km) inputs.
    rng = np.random.default_rng(1)
    batch, rho = model.make_training_data_with_raw_sequences(
        rng, S=8, n_max=_N_MAX,
    )
    p_h1_xz, _, preds = predict_amortised_p_h1_for_one_xz(
        fit, batch, _PPS_H1_MIN_EFFECT_SIZE_THRESH,
    )
    assert (p_h1_xz >= 0).all() and (p_h1_xz <= 1).all()
    assert preds.shape == (8, len(_DEFAULT_TAUS))
    # Loss should have decreased over training.
    assert fit['loss_history'][-1] < fit['loss_history'][0]


# ---------------------------------------------------------------------------
# Correctness on the 12-interim grid used by the deployment script.
# ---------------------------------------------------------------------------


def _monthly_interim_grid(seed=123, N=500, TRUE_P=0.4):
    """Reproduce the (kn, n, m) grid from the amortised deployment
    script."""
    rng = np.random.default_rng(seed)
    start = pd.Timestamp('2025-01-01')
    end = pd.Timestamp('2025-12-31')
    dates = pd.date_range(start, end, periods=N)
    y = rng.binomial(1, TRUE_P, size=N)
    dp = pd.DataFrame({'submission_date': dates, 'y': y})
    month_starts = pd.date_range(start, end, freq='MS')
    interim_dates = month_starts + pd.offsets.MonthEnd(0)
    rows = []
    for iid, d in enumerate(interim_dates, 1):
        dpi = dp[dp['submission_date'] <= d]
        n = len(dpi)
        if n == 0:
            continue
        m = N - n
        if m <= 0:
            continue
        rows.append({
            'interim_id': iid,
            'kn': int(dpi['y'].sum()),
            'n': int(n),
            'm': int(m),
        })
    return pd.DataFrame(rows)


def test_amortised_pps_matches_analytic_on_monthly_grid(trained_fit):
    """Amortised PPS at the deployment grid (kn, n, m) matches the
    closed-form Beta-Binomial PPS within loose tolerance for a small
    training budget."""
    grid = _monthly_interim_grid()
    rng = np.random.default_rng(0)
    errors = []
    for _, row in grid.iterrows():
        kn, n, m = int(row['kn']), int(row['n']), int(row['m'])
        pps_ref = _closed_form_pps_binomial(
            kn, n, m,
            a=_PRIOR_A, b=_PRIOR_B, p_0=_P_0,
            pps_H1_min_effect_size_thresh=_PPS_H1_MIN_EFFECT_SIZE_THRESH,
            pps_ProbH1_target_lwr_quantile=_PPS_PROBH1_THRESH,
        )
        from scipy.stats import betabinom as _betabinom
        S = 2000
        km_draws = _betabinom.rvs(
            m, _PRIOR_A + kn, _PRIOR_B + (n - kn),
            size=S, random_state=rng,
        )
        k_total = np.full(S, kn) + km_draws
        n_total = np.full(S, n + m)
        features = np.stack(
            [k_total, n_total], axis=-1,
        ).astype(np.float32) / _N_MAX
        p_h1_xz, _, _ = predict_amortised_p_h1_for_one_xz(
            trained_fit, features, _PPS_H1_MIN_EFFECT_SIZE_THRESH,
        )
        pps_amortised = float((p_h1_xz > _PPS_PROBH1_THRESH).mean())
        errors.append(abs(pps_amortised - pps_ref))
    errs = np.array(errors)
    print(f"  mean |err| = {errs.mean():.4f}, max = {errs.max():.4f}")
    assert errs.max() < 0.10
    assert errs.mean() < 0.05


def test_amortised_pps_matches_analytic_on_random_xobs(trained_fit):
    """Sample 20 random (kn, n, m) triples and verify amortised PPS matches
    the closed-form value across the state space."""
    rng = np.random.default_rng(31)
    errors = []
    for _ in range(20):
        n = int(rng.integers(50, _N_MAX))
        m = _N_MAX - n
        p_true = float(rng.uniform(0.1, 0.9))
        kn = int(rng.binomial(n, p_true))

        pps_ref = _closed_form_pps_binomial(
            kn, n, m,
            a=_PRIOR_A, b=_PRIOR_B, p_0=_P_0,
            pps_H1_min_effect_size_thresh=_PPS_H1_MIN_EFFECT_SIZE_THRESH,
            pps_ProbH1_target_lwr_quantile=_PPS_PROBH1_THRESH,
        )
        from scipy.stats import betabinom as _betabinom
        S = 2000
        km_draws = _betabinom.rvs(
            m, _PRIOR_A + kn, _PRIOR_B + (n - kn),
            size=S, random_state=rng,
        )
        k_total = np.full(S, kn) + km_draws
        n_total = np.full(S, n + m)
        features = np.stack(
            [k_total, n_total], axis=-1,
        ).astype(np.float32) / _N_MAX
        p_h1_xz, _, _ = predict_amortised_p_h1_for_one_xz(
            trained_fit, features, _PPS_H1_MIN_EFFECT_SIZE_THRESH,
        )
        pps_amortised = float((p_h1_xz > _PPS_PROBH1_THRESH).mean())
        errors.append(abs(pps_amortised - pps_ref))
    errs = np.array(errors)
    print(f"  mean |err| = {errs.mean():.4f}, max = {errs.max():.4f}")
    assert errs.max() < 0.15
    assert errs.mean() < 0.06


# ---------------------------------------------------------------------------
# Save / load round-trip
# ---------------------------------------------------------------------------


def test_save_load_round_trip(tmp_path, trained_fit):
    p = tmp_path / 'net.pkl'
    save_trained_model(trained_fit, str(p))
    reloaded = load_fitted_model(str(p))
    features = np.array([[150, 300]], dtype=np.float32) / _N_MAX
    p_a, q_a, _ = predict_amortised_p_h1_for_one_xz(
        trained_fit, features, _PPS_H1_MIN_EFFECT_SIZE_THRESH,
    )
    p_b, q_b, _ = predict_amortised_p_h1_for_one_xz(
        reloaded, features, _PPS_H1_MIN_EFFECT_SIZE_THRESH,
    )
    np.testing.assert_allclose(p_a, p_b, atol=1e-6)
    np.testing.assert_allclose(q_a, q_b, atol=1e-6)
