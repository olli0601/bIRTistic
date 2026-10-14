"""
Step 3 of the OO-port refactor: :class:`PartialCreditModel` must
return bit-identical outputs to the existing free-function pipeline in
:mod:`fit_partial_credit_model` for every blueprint method.

These tests load the fixture cohort in ``test/test_data/`` and pair each
``model.X(...)`` call with the corresponding free-function call,
asserting the two agree.
"""

import sys
from pathlib import Path

import arviz as az
import jax.numpy as jnp
import numpy as np
import pandas as pd
import pytest

_repo_root = Path(__file__).resolve().parents[2]
_python_dir = _repo_root / 'python'
if str(_python_dir) not in sys.path:
    sys.path.insert(0, str(_python_dir))

from model_pcm import PartialCreditModel, _model_namespace

_TEST_DATA = _repo_root / 'test' / 'test_data'
_INTERIM_STEM = _TEST_DATA / 'pcm_1_interim_1'
_DRAWS_FILE = f"{_INTERIM_STEM}_draws.zarr"


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
def model(dit, xi):
    return PartialCreditModel(dit=dit, dcati=xi,
                                   x_formula="~ time - 1", seed=123)


@pytest.fixture(scope='module')
def stan_data(model):
    return model.stan_data


@pytest.fixture(scope='module')
def params_one(model, draws):
    """A single posterior-draw parameter dict (the first draw) for use
    with eval_loglik / eval_log_prior / eval_outcome."""
    theta = model.get_stacked_posterior(draws)
    return {k: v[0] for k, v in theta.items()}


# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------


def test_param_names_includes_known_pcm_params(model):
    expected = {
        'latent_factor_unit', 'latent_factor_beta',
        'loadings_questions_m1',
    }
    assert expected.issubset(set(model.param_names))


def test_positive_params_is_subset_of_param_names(model):
    assert set(model.positive_params).issubset(set(model.param_names))


# ---------------------------------------------------------------------------
# Stan data
# ---------------------------------------------------------------------------


def test_stan_data_built_at_construction(model, stan_data):
    assert isinstance(stan_data, dict)
    assert 'N_total' in stan_data
    assert int(stan_data['N_total']) == len(model.dcati)


def test_make_stan_data_shape_invariants(model, stan_data, dit, xi):
    """make_stan_data must produce the standard partial-credit stan_data
    keys + shapes against the fixture cohort."""
    expected_keys = {
        'C', 'U', 'N', 'Q', 'K', 'N_total', 'Q_total', 'y', 'unit_of_obs',
        'question_of_obs', 'cat_type', 'X', 'Xpr_id', 'P', 'P_pr',
    }
    assert expected_keys.issubset(set(stan_data.keys()))
    assert int(stan_data['N_total']) == len(xi)
    assert int(stan_data['C']) == int(dit['item_type_id'].max())
    # Call directly: model.make_stan_data(dcati, x_formula) returns the same.
    again = model.make_stan_data(xi, "~ time - 1")
    assert set(again.keys()) == set(stan_data.keys())


# ---------------------------------------------------------------------------
# Loglik / prior / outcome -- the @staticmethod proxies must agree with the
# underlying .pyro namespace they delegate to.
# ---------------------------------------------------------------------------


def test_eval_loglik_matches_pyro_namespace(model, stan_data, params_one):
    ns = _model_namespace()
    direct = ns['get_log_likelihood_of_partial_credit_model_ncats'](stan_data, params_one)
    by_method = model.eval_loglik(stan_data, params_one)
    np.testing.assert_array_equal(np.asarray(direct), np.asarray(by_method))


def test_eval_loglik_annealed_matches_pyro_namespace(model, stan_data, params_one):
    p = dict(params_one)
    p['temperature'] = jnp.asarray(0.5)
    ns = _model_namespace()
    direct = ns['get_log_likelihood_of_partial_credit_model_ncats_with_annealing'](stan_data, p)
    by_method = model.eval_loglik_annealed(stan_data, p)
    np.testing.assert_array_equal(np.asarray(direct), np.asarray(by_method))


def test_eval_log_prior_matches_pyro_namespace(model, params_one):
    direct = _model_namespace()['get_prior_of_partial_credit_model_ncats'](params_one)
    by_method = model.eval_log_prior(params_one)
    assert float(direct) == float(by_method)


def test_eval_outcome_matches_pyro_namespace(model, stan_data, params_one):
    direct = _model_namespace()['get_ordered_prob_of_partial_credit_model_ncats'](stan_data, params_one)
    by_method = model.eval_outcome_for_endpoint(stan_data, params_one)
    np.testing.assert_array_equal(np.asarray(direct), np.asarray(by_method))


# ---------------------------------------------------------------------------
# Endpoints per draw
# ---------------------------------------------------------------------------


def test_get_endpoints_per_draw_shape(model, draws):
    """Output schema invariants for the model-aware endpoint method."""
    out = model.get_endpoints_per_draw(
        draws=draws, categorical_threshold=2,
        endpoint_type='items', param_name='ordered_prob_by_cat_qu_fit',
    )
    # ordered-prob path emits the per-condition levels as group1/group2 (the rho_by_cat_qu
    # path labels them identically, so get_endpoints consumes either unchanged).
    expected_cols = {'item_label', 'item_type', 'item_high_label',
                     'draw', 'ratio', 'group1', 'group2', 'diff'}
    assert expected_cols.issubset(set(out.columns))
    assert out['ratio'].notna().any()


def test_get_endpoints_per_draw_rho_by_cat_qu_decodes_labels(model):
    """The rho_by_cat_qu reader (param_name='rho_by_cat_qu_pr') label-decodes the in-model
    reductions via self.drho and applies the R direction rule: group1/group2 are the raw
    Baseline/Endline levels; diff/ratio keep the Baseline-Endline orientation for categorical OR
    lower-is-better items and flip for out-of-7 higher-is-better; diff2/ratio2 are the negatives;
    ratio is a percent in the matrix and a fraction out."""
    drho = model.drho
    R = int(drho['r'].nunique())
    Qcmp = int(drho['jc'].nunique())
    r_of = {lbl: int(drho.loc[drho['rho'] == lbl, 'r'].iloc[0]) for lbl in ('Baseline', 'Endline', 'ratio', 'diff')}

    # Deterministic, distinct per-comparison Baseline/Endline levels (B > 0), same across 3 draws.
    jc = np.arange(1, Qcmp + 1, dtype=float)
    B = 1.0 + 0.1 * jc
    E = 0.5 + 0.05 * jc
    ndraws = 3
    arr = np.zeros((1, ndraws, R, Qcmp))
    arr[:, :, r_of['Baseline'] - 1, :] = B
    arr[:, :, r_of['Endline'] - 1, :] = E
    arr[:, :, r_of['diff'] - 1, :] = B - E
    arr[:, :, r_of['ratio'] - 1, :] = (B - E) / B * 100.0        # model emits percent
    idata = az.from_dict(posterior={'rho_by_cat_qu_pr': arr})

    out = model.get_endpoints_per_draw(
        draws=idata, endpoint_type='items', param_name='rho_by_cat_qu_pr',
    )
    assert {'draw', 'item_type_id', 'item_type', 'item_label', 'construct', 'item_high_label',
            'group1', 'group2', 'diff', 'ratio', 'diff2', 'ratio2'}.issubset(out.columns)

    # Expected per item-scope comparison, keyed by (item_type_id, item_label).
    items = drho[drho['scope'] == 'item'].drop_duplicates('jc').copy()
    items['B'] = 1.0 + 0.1 * items['jc'].astype(float)
    items['E'] = 0.5 + 0.05 * items['jc'].astype(float)
    items['dir_pos'] = (items['item_type'] == 'categorical') | (items['item_high_label'] == 'lower_is_better')
    exp = items[['item_type_id', 'item_label', 'B', 'E', 'dir_pos']]

    chk = out.merge(exp, on=['item_type_id', 'item_label'])
    assert len(chk) == len(out)                                  # every output row mapped to a comparison
    raw_diff = chk['B'] - chk['E']
    raw_ratio = raw_diff / chk['B']
    np.testing.assert_allclose(chk['group1'], chk['B'], rtol=1e-9)
    np.testing.assert_allclose(chk['group2'], chk['E'], rtol=1e-9)
    np.testing.assert_allclose(chk['diff'], np.where(chk['dir_pos'], raw_diff, -raw_diff), rtol=1e-6)
    np.testing.assert_allclose(chk['ratio'], np.where(chk['dir_pos'], raw_ratio, -raw_ratio), rtol=1e-6)
    np.testing.assert_allclose(chk['diff2'], -raw_diff, rtol=1e-6)
    np.testing.assert_allclose(chk['ratio2'], -raw_ratio, rtol=1e-6)

    # item_groups scope: schema + non-empty (participant-weighted group reductions).
    grp = model.get_endpoints_per_draw(
        draws=idata, endpoint_type='item_groups', param_name='rho_by_cat_qu_pr',
    )
    assert {'draw', 'item_type', 'construct', 'group1', 'group2', 'diff', 'ratio'}.issubset(grp.columns)
    assert len(grp) == ndraws * drho[drho['scope'] == 'item_group']['jc'].nunique()


