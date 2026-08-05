#!/usr/bin/env python3
"""
MVN interim analysis: **amortised** PPS via a DeepSets multi-quantile
network with hardcoded per-component identity encoder (§13.1 of
``dev/amortised_decision_making.md``).

**idcomp variant** -- independent-components approximation. The
amortiser sees ONE component's summary at a time, so the joint
posterior over ``mu`` is factorised across components. The known
cross-component correlation ``K`` is used at TRAINING (to draw
``mu ~ MVN(prior_mu, tau^2 K)`` correctly) but discarded at DEPLOYMENT
(each ``mu_j`` is predicted marginally, so downstream utilities that
need the joint -- multivariate stopping rules, family-wise error --
would have to reconstruct it from something else). The full-K variant
that exploits the known correlation at deployment is future work.

Under the MVN g-prior with known ``K``, ``sigma``, ``prior_tau``:

  - Posterior of ``mu_j | y_{1:N}`` depends only on the per-component
    sufficient statistic ``(sum_i y_{i,j}, N)``. Block-equicorr K has
    ``K_{jj} = 1`` for all j, so per-component posterior distributions
    are J-invariant.
  - We train ONE amortiser on prior draws from an MVN with any J in
    the grid (per-component decomposition flattens (S, J) into S*J
    training examples) and deploy across all J values in
    ``simu_params['J_grid']``.

Prior-predictive training sampler:
    mu ~ MVN(prior_mu, prior_tau^2 K);
    n ~ U(1, N_full - 1);  m = N_full - n;
    y ~ MVN(mu, sigma^2 K), size N_full;
    features_j = (sum_i y_{i,j}, N_full) / N_full;
    rho_j = mu_j - mu_0_baseline.

Multi-quantile head at 11 tau levels; monotone-corrected CDF interp at
``pps_H1_min_effect_size_thresh = 0.0`` for continuous
``P(H_1j | x, z)``.

Outputs per J (``_RGEA_`` suffix; matches the MVN compare-methods
loader):

  - ``mvn_J{J}_pps_RGEA_p_h1_xz.pkl``
  - ``mvn_J{J}_pps_RGEA.csv``
  - ``mvn_J{J}_pps_RGEA_timing.csv``
  - ``mvn_J{J}_pps_RGEA_perf.csv``

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/MVN_interim_analysis_amortise_endptx_on_wz_with_features_fixed_idcomp_qpsi_MLP_loss_multiquantilehead.py
"""

# %%

import os
import sys
import time
import warnings
from functools import partial
from pathlib import Path

try:
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
except NameError:
    project_root = Path.cwd()
    if project_root.name == 'scripts-py':
        project_root = project_root.parent

python_path = str(project_root / 'python')
if python_path not in sys.path:
    sys.path.insert(0, python_path)

import numpy as np
import pandas as pd
import arviz as az
import xarray as xr

warnings.filterwarnings('ignore')

from model_mvn import MVNModel
from amortiser_common import (
    train,
    save_trained_model,
    load_fitted_model,
)
from amortiser_pps_features_fixed_qpsi_MLP_loss_multiquantilehead import (
    Amortiser_PPS_features_fixed_qpsi_MLP_loss_multiquantilehead,
    predict_amortised_p_h1_for_many_xz,
)

print("Imports successful")

# %%

# =============================================================================
# Configuration
# =============================================================================

DIR_SIM = os.path.join(
    "/Users/or105/sandbox/bIRTistic", "py-mvn-interim-simulations-260609",
)
_ROOT_DIR = (
    "/Users/or105/sandbox/bIRTistic/"
    "py-mvn-interim-amortise-endptx-on-wz-with-features-fixed-"
    "idcomp-qpsi-MLP-loss-multiquantilehead"
)
dir_out = f"{_ROOT_DIR}-260714"
os.makedirs(dir_out, exist_ok=True)
print(f"Output dir: {dir_out}")

PPS_Z_TOTAL = 4000

# Amortised net training config. Match the Binomial §12.1 default
# (64x64 head, 11 tau grid, adam warmup-cosine LR). The
# ``TRAIN_J`` picks which cell of ``sim_data['cells']`` to use for the
# prior K_chol -- the trained net is J-invariant under the per-component
# decomposition, so this is only a batching choice.
TRAIN_J = 20
NET_HIDDEN = (64, 64)
NET_TAUS = (0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
NET_STEPS = 40000
NET_BATCH = 512               # S samples -> S * TRAIN_J examples per step.
NET_LR = 1e-3
N_DIST = 'uniform'

file_prefix = "mvn_interim"
NET_CKPT = os.path.join(dir_out, f"{file_prefix}_amortised_pps_net.pkl")

# %%

# =============================================================================
# Load upstream simulation data (cohort + interim cache + closed-form pps).
# =============================================================================

sim_data = pd.read_pickle(os.path.join(DIR_SIM, 'mvn_sim_data.pkl'))
simu_params = sim_data['simu_params']
di = sim_data['interim_grid']
J_GRID = list(simu_params['J_grid'])
seed = simu_params['seed']
NET_SEED = seed
NET_N_MAX = simu_params['N_full']
pps_H1_min_effect_size_thresh = float(
    simu_params.get(
        'pps_H1_min_effect_size_thresh',
        simu_params.get('pps_H1_def', 0.0),
    ),
)
pps_ProbH1_target_lwr_quantile = float(
    simu_params.get(
        'pps_ProbH1_target_lwr_quantile',
        simu_params.get('pps_ProbH1_thresh', 0.89),
    ),
)

pps_cf = pd.read_pickle(os.path.join(DIR_SIM, 'mvn_pps_closed_form.pkl'))['pps_cf']

interim_data_by_J = {}
for J in J_GRID:
    pkl_path = os.path.join(DIR_SIM, f'mvn_J{J}_interim_data.pkl')
    if not os.path.exists(pkl_path):
        raise FileNotFoundError(
            f"{pkl_path} missing; "
            f"run MVN_interim_analyses_make_interim_data.py first."
        )
    interim_data_by_J[J] = pd.read_pickle(pkl_path)

print(f"Loaded sim_data + {len(J_GRID)} interim_data J cells + pps_cf "
      f"({len(pps_cf)} rows). N_full = {NET_N_MAX}.")

# %%

# =============================================================================
# Train amortised net once (or load).
# =============================================================================

if os.path.exists(NET_CKPT):
    print(f"\nLoading cached amortised net from {NET_CKPT}")
    fit = load_fitted_model(NET_CKPT)
else:
    print(f"\nTraining amortised net "
          f"({NET_STEPS} steps, batch {NET_BATCH}, TRAIN_J={TRAIN_J}) ...")

    cell_train = sim_data['cells'][TRAIN_J]
    _amortiser_model = MVNModel(
        dit=cell_train['dit'],
        dcati=MVNModel.get_interim_data_x(cell_train['dp']),
        seed=seed,
        J=TRAIN_J, K_chol=cell_train['R_chol'],
        sigma=simu_params['sigma'],
        prior_tau=simu_params['prior_tau'],
        mu_0_baseline=simu_params['mu_0_baseline'],
    )
    sample_fn = partial(
        _amortiser_model.make_training_data_with_features,
        n_max=NET_N_MAX, n_dist=N_DIST,
    )
    fit = train(
        sample_fn,
        Amortiser_PPS_features_fixed_qpsi_MLP_loss_multiquantilehead,
        net_kwargs={'hidden_dims': NET_HIDDEN},
        pps_ProbH1_lwr_quantiles_mesh=NET_TAUS,
        num_steps=NET_STEPS,
        batch_size=NET_BATCH,
        lr=NET_LR,
        seed=NET_SEED,
        verbose=True,
    )
    save_trained_model(fit, NET_CKPT)
    print(f"Saved amortised net to {NET_CKPT}")
mins_train = float(fit.get('training_mins', float('nan')))
print(f"Amortised net trained in {mins_train:.2f} min "
      f"(excludes JIT compilation, data-loading, and load-cache time).")

# %%

# =============================================================================
# Per-J deployment loop. Reuses ``interim_data_by_J`` (mu draws + zi
# already cached upstream) so the amortised path only computes features
# and a forward pass.
# =============================================================================

for J in J_GRID:
    cell = sim_data['cells'][J]
    dit, dp, R_chol = cell['dit'], cell['dp'], cell['R_chol']
    print(f"\n{'#' * 70}\n# J = {J} (n_items = {J})\n{'#' * 70}")

    p_h1_xz_rows = []
    pps_timing_rows = []
    perf_rows = []
    t_J0 = time.time()
    for interim_id, blk in interim_data_by_J[J].items():
        t0 = time.time()
        dpi = blk['dpi']
        zi_full = blk['zi']
        mu_draws_all = blk['mu_draws']                # (S_cache, J)
        interim_date = blk['interim_date']
        interim_month_year = blk['interim_month_year']
        interim_m = blk['interim_m']
        n_obs = blk['n_obs']

        ypred_cols = sorted(
            [c for c in zi_full.columns if c.startswith('ypred_')],
            key=lambda c: int(c.split('_')[1]),
        )
        if PPS_Z_TOTAL > len(ypred_cols):
            raise ValueError(
                f"PPS_Z_TOTAL={PPS_Z_TOTAL} > cached S={len(ypred_cols)}."
            )
        keep_cols = ypred_cols[:PPS_Z_TOTAL]
        meta_cols = [c for c in zi_full.columns if not c.startswith('ypred_')]
        zi = zi_full[meta_cols + keep_cols].copy()

        print(f"\n{'=' * 70}\nJ={J} RGEA interim {interim_id}"
              f" ({interim_date.date()}) n={n_obs} m={interim_m}"
              f" S={PPS_Z_TOTAL}\n{'=' * 70}")

        model_x = MVNModel(
            dit=dit, dcati=dpi, seed=seed,
            J=J, K_chol=R_chol, sigma=simu_params['sigma'],
            prior_tau=simu_params['prior_tau'],
            mu_0_baseline=simu_params['mu_0_baseline'],
        )

        # ---- Per-(j, draw) w summary from cached zi (mean over m rows).
        wa = model_x.get_w(zi)                             # (J*S rows)

        # ---- Per-(j) sum of observed x column.
        sum_x_j = (
            dpi.groupby('j')['y'].sum()
            .reindex(range(J), fill_value=0.0)
            .to_numpy(dtype=np.float64)
        )                                                  # (J,)

        # ---- DeepSets pooling: sum_T_x + sum_T_z = sum_j y_j over ALL rows.
        # wa['w'] = (1/m) sum_i z_{i, j}^(s) -> sum_z = w * interim_m.
        wa['sum_x'] = sum_x_j[wa['j'].to_numpy(dtype=int)]
        wa['sum_z'] = wa['w'].to_numpy() * float(interim_m)
        wa['n_total'] = float(n_obs + interim_m)
        scale = float(NET_N_MAX)
        wa['k_total'] = (wa['sum_x'] + wa['sum_z']) / scale
        wa['n_total'] = wa['n_total'] / scale

        p_h1_xz_interim, perf_interim = predict_amortised_p_h1_for_many_xz(
            wa, fit,
            pps_H1_min_effect_size_thresh=pps_H1_min_effect_size_thresh,
            pps_ProbH1_target_lwr_quantile=pps_ProbH1_target_lwr_quantile,
            feature_cols=('k_total', 'n_total'),
        )
        p_h1_xz_interim['J'] = J
        p_h1_xz_interim['interim_id'] = interim_id
        p_h1_xz_interim['interim_date'] = interim_date
        p_h1_xz_interim['interim_month_year'] = interim_month_year
        p_h1_xz_rows.append(p_h1_xz_interim)

        perf_interim['J'] = J
        perf_interim['interim_id'] = interim_id
        perf_interim['interim_date'] = interim_date
        perf_interim['interim_month_year'] = interim_month_year
        perf_rows.append(perf_interim)

        mins_interim = (time.time() - t0) / 60.0
        pps_timing_rows.append({
            'J':                  J,
            'interim_id':         interim_id,
            'interim_date':       interim_date,
            'interim_month_year': interim_month_year,
            'n_obs':              n_obs,
            'm_future':           interim_m,
            'mins_interim_id':    round(mins_interim, 3),
        })
        print(f"  J={J} interim {interim_id} amortised PPS done in"
              f" {mins_interim:.3f} min")

    if not p_h1_xz_rows:
        print(f"J={J}: no interims with m > 0; skipping save.")
        continue

    dp_h1_xz = pd.concat(p_h1_xz_rows, ignore_index=True)
    dp_h1_xz['pps_H1_min_effect_size_thresh'] = pps_H1_min_effect_size_thresh
    dp_h1_xz['pps_ProbH1_target_lwr_quantile'] = pps_ProbH1_target_lwr_quantile
    dp_h1_xz['S'] = PPS_Z_TOTAL
    dp_h1_xz['j'] = dp_h1_xz['item_label'].str.replace('mu_', '').astype(int)
    perf_all = pd.concat(perf_rows, ignore_index=True)
    perf_all['j'] = perf_all['item_label'].str.replace('mu_', '').astype(int)
    timing_all = pd.DataFrame(pps_timing_rows)
    timing_all['mins_J_total'] = round((time.time() - t_J0) / 60.0, 3)

    pps_df = (
        dp_h1_xz.groupby(['J', 'interim_id', 'interim_date',
                          'interim_month_year', 'item_label', 'item_type',
                          'item_high_label', 'j'], observed=True)['p_h1_xz']
        .apply(lambda p: float(
            (p > pps_ProbH1_target_lwr_quantile).mean()
        ))
        .reset_index(name='pps')
    )
    pps_df['eta'] = pps_ProbH1_target_lwr_quantile
    pps_df['S'] = PPS_Z_TOTAL

    dp_h1_xz.to_pickle(
        os.path.join(dir_out, f'mvn_J{J}_pps_RGEA_p_h1_xz.pkl'),
    )
    perf_all.to_csv(
        os.path.join(dir_out, f'mvn_J{J}_pps_RGEA_perf.csv'), index=False,
    )
    timing_all.to_csv(
        os.path.join(dir_out, f'mvn_J{J}_pps_RGEA_timing.csv'), index=False,
    )
    pps_df.to_csv(
        os.path.join(dir_out, f'mvn_J{J}_pps_RGEA.csv'), index=False,
    )

    # ---- Compare vs analytic pps_cf (this J) ----
    tab = pps_df.merge(
        pps_cf[pps_cf['J'] == J][['interim_id', 'j', 'pps']].rename(
            columns={'pps': 'pps_analytic'},
        ),
        on=['interim_id', 'j'], how='left',
    )
    tab['abs_err'] = (tab['pps'] - tab['pps_analytic']).abs()
    print(f"J={J}: RGEA vs analytic PPS: max abs {tab['abs_err'].max():.4f},"
          f" mean {tab['abs_err'].mean():.4f} over {len(tab)} (j, interim)"
          f" cells; total J time {(time.time() - t_J0) / 60.0:.2f} min.")

print("\nAll J done -- MVN features-fixed amortised PPS complete.")
