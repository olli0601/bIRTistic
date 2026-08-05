#!/usr/bin/env python3
"""
MVN interim analysis: **amortised** PPS via a nested-DeepSets encoder
over (participant, component) tokens + K-row query
(§13.4 of ``dev/amortised_decision_making.md``).

**xcomp variant** -- exploits the known cross-component covariance
``K``. Each per-token feature is the two scalars ``(K[j*, j],
y_{i, j})``, so the amortiser is told, per token, how correlated the
token's component is with the queried component. Independent-
components idcomp (§13.1 / §13.2) is the special case ``K = I``.

Outputs per J (``_RGEC_`` suffix; picked up by the MVN compare-methods
loader alongside the idcomp variants):

  - ``mvn_J{J}_pps_RGEC_p_h1_xz.pkl``
  - ``mvn_J{J}_pps_RGEC.csv``
  - ``mvn_J{J}_pps_RGEC_timing.csv``
  - ``mvn_J{J}_pps_RGEC_perf.csv``

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/MVN_interim_analysis_amortise_endptx_on_wz_with_features_MLP_xcomp_qpsi_MLP_loss_multiquantilehead.py
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

warnings.filterwarnings('ignore')

from model_mvn import MVNModel
from amortiser_common import (
    train,
    save_trained_model,
    load_fitted_model,
    predict_amortised_p_h1_for_one_xz,
)
from amortiser_pps_features_MLP_xcomp_qpsi_MLP_loss_multiquantilehead import (
    Amortiser_PPS_features_MLP_xcomp_qpsi_MLP_loss_multiquantilehead,
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
    "py-mvn-interim-amortise-endptx-on-wz-with-features-MLP-"
    "xcomp-qpsi-MLP-loss-multiquantilehead"
)
# 15k-step overnight run (see §13.4 follow-ups). Preserve the 4k
# checkpoint dir at _260715 for reference; write the longer run to a
# separate dir so both remain reproducible.
_DIR_VARIANT = os.environ.get('XCOMP_DIR_VARIANT', '_15k_260716')
dir_out = f"{_ROOT_DIR}{_DIR_VARIANT}"
os.makedirs(dir_out, exist_ok=True)
print(f"Output dir: {dir_out}")

PPS_Z_TOTAL = 500                # matches RGEB; each posterior draw
                                 # triggers a J-way batched forward pass
                                 # on (J, N_max, J) tokens.

TRAIN_J = 20                     # train on the smallest cell in the grid;
                                 # net is J-invariant via the K-row query.
NET_Q_TAU_INNER_HIDDEN = (32, 32)
NET_EMBED_INNER = 32
NET_HIDDEN = (64, 64)
NET_TAUS = (0.05, 0.25, 0.5, 0.75, 0.95)
NET_STEPS = int(os.environ.get('XCOMP_STEPS', 15000))
NET_BATCH = 16                   # S=16 samples * Q=4 queries -> B=64 rows.
NET_QUERIES = 4
NET_LR = 1e-3
N_DIST = 'uniform'

file_prefix = "mvn_interim"
NET_CKPT = os.path.join(dir_out, f"{file_prefix}_amortised_pps_net.pkl")

# %%

# =============================================================================
# Load upstream simulation data.
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
    print(f"\nTraining xcomp amortised net "
          f"({NET_STEPS} steps, S={NET_BATCH}, Q={NET_QUERIES}, "
          f"TRAIN_J={TRAIN_J}) ...")

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
        _amortiser_model.make_training_data_with_participant_sequences,
        n_max=NET_N_MAX, n_dist=N_DIST,
        queries_per_sample=NET_QUERIES,
    )
    fit = train(
        sample_fn,
        Amortiser_PPS_features_MLP_xcomp_qpsi_MLP_loss_multiquantilehead,
        net_kwargs={
            'q_tau_inner_hidden': NET_Q_TAU_INNER_HIDDEN,
            'embed_inner': NET_EMBED_INNER,
            'hidden_dims': NET_HIDDEN,
        },
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
      f"(excludes JIT + data-loading + load-cache).")

# %%

# =============================================================================
# Per-J deployment loop. Each interim posterior draw s becomes one
# batch of J elements (one per query component). K_row / K_diag per
# element carry the cross-component correlation profile.
# =============================================================================

for J in J_GRID:
    cell = sim_data['cells'][J]
    dit, dp, R_chol = cell['dit'], cell['dp'], cell['R_chol']
    K = R_chol @ R_chol.T
    K_diag = np.diag(K)
    print(f"\n{'#' * 70}\n# J = {J} (n_items = {J})\n{'#' * 70}")

    p_h1_xz_rows = []
    pps_timing_rows = []
    perf_rows = []
    t_J0 = time.time()
    for interim_id, blk in interim_data_by_J[J].items():
        t0 = time.time()
        dpi = blk['dpi']
        zi_full = blk['zi']
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

        print(f"\n{'=' * 70}\nJ={J} RGEC interim {interim_id}"
              f" ({interim_date.date()}) n={n_obs} m={interim_m}"
              f" S={PPS_Z_TOTAL}\n{'=' * 70}")

        # ---- Observed cohort as (n_obs, J) matrix.
        x_wide = (
            dpi.pivot_table(index='pid', columns='j', values='y')
            .sort_index()
            .reindex(columns=range(J))
            .to_numpy(dtype=np.float64)
        )
        # ---- Future cohort as (S, m, J) tensor.
        zi_sorted = zi_full.sort_values(['pid', 'j']).reset_index(drop=True)
        ypred_arr = (
            zi_sorted[keep_cols].to_numpy().T.reshape(
                PPS_Z_TOTAL, interim_m, J,
            )
        )
        # ---- Pad both to N_max on the participant axis.
        x_pad = np.zeros((NET_N_MAX, J), dtype=np.float32)
        x_pad[:n_obs] = x_wide.astype(np.float32)
        mask_x = np.zeros((NET_N_MAX,), dtype=np.float32)
        mask_x[:n_obs] = 1.0
        # z padded shape per draw (S, N_max, J).
        z_pad = np.zeros((PPS_Z_TOTAL, NET_N_MAX, J), dtype=np.float32)
        z_pad[:, :interim_m] = ypred_arr.astype(np.float32)
        mask_z = np.zeros((NET_N_MAX,), dtype=np.float32)
        mask_z[:interim_m] = 1.0
        sizes_row = np.array(
            [n_obs / NET_N_MAX, interim_m / NET_N_MAX], dtype=np.float32,
        )

        # ---- For each posterior draw s, batch over J query components.
        # Batch element (s, j*) : x = shared x_pad, z = z_pad[s],
        # k_row = K[j*], k_diag = K[j*, j*].
        S = PPS_Z_TOTAL
        # (J, N_max, J) shared x (broadcast per query), (J, N_max) shared mask.
        x_batch = np.broadcast_to(
            x_pad[None, :, :], (J, NET_N_MAX, J),
        ).astype(np.float32)
        mask_x_batch = np.broadcast_to(
            mask_x[None, :], (J, NET_N_MAX),
        ).astype(np.float32)
        mask_z_batch = np.broadcast_to(
            mask_z[None, :], (J, NET_N_MAX),
        ).astype(np.float32)
        sizes_batch = np.broadcast_to(
            sizes_row[None, :], (J, 2),
        ).astype(np.float32)
        k_row_batch = K.astype(np.float32)                    # (J, J)
        k_diag_batch = K_diag[:, None].astype(np.float32)     # (J, 1)

        preds_per_s = []  # each element (J, num_quantiles)
        for s in range(S):
            batch = {
                'x':       x_batch,
                'mask_x':  mask_x_batch,
                'z':       z_pad[s][None].repeat(J, axis=0),
                'mask_z':  mask_z_batch,
                'sizes':   sizes_batch,
                'k_row':   k_row_batch,
                'k_diag':  k_diag_batch,
            }
            p_h1_xz_s, _q_hat_s, preds_s = predict_amortised_p_h1_for_one_xz(
                fit, batch, pps_H1_min_effect_size_thresh,
            )
            preds_per_s.append((p_h1_xz_s, preds_s))
        # Stack: p_h1_xz shape (S, J), preds shape (S, J, num_quantiles).
        p_h1_xz = np.stack([p for p, _ in preds_per_s], axis=0)
        q_hat_all = np.stack(
            [pr[:, len(NET_TAUS) // 2] for _, pr in preds_per_s], axis=0,
        )

        # ---- Rebuild long-form (item_label = mu_j, s, p_h1_xz).
        item_labels = dit['item_label'].to_numpy()
        item_types  = dit['item_type'].to_numpy()
        item_highs  = dit['item_high_label'].to_numpy()
        rows = []
        for j in range(J):
            rows.append(pd.DataFrame({
                'item_label':      [item_labels[j]] * S,
                'item_type':       [item_types[j]] * S,
                'item_high_label': [item_highs[j]] * S,
                'j':               np.full(S, j, dtype=int),
                's':               np.arange(1, S + 1, dtype=int),
                'p_h1_xz':         p_h1_xz[:, j].astype(np.float64),
                'q_hat':           q_hat_all[:, j].astype(np.float64),
            }))
        df = pd.concat(rows, ignore_index=True)
        df['J'] = J
        df['interim_id'] = interim_id
        df['interim_date'] = interim_date
        df['interim_month_year'] = interim_month_year
        p_h1_xz_rows.append(df)

        for j in range(J):
            perf_rows.append({
                'J':                  J,
                'interim_id':         interim_id,
                'interim_date':       interim_date,
                'interim_month_year': interim_month_year,
                'j':                  j,
                'item_label':         item_labels[j],
                'item_type':          item_types[j],
                'rho':                float('nan'),
                'r2':                 float('nan'),
            })

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
        print(f"  J={J} interim {interim_id} xcomp amortised done in"
              f" {mins_interim:.3f} min")

    if not p_h1_xz_rows:
        print(f"J={J}: no interims with m > 0; skipping save.")
        continue

    dp_h1_xz = pd.concat(p_h1_xz_rows, ignore_index=True)
    dp_h1_xz['pps_H1_min_effect_size_thresh'] = pps_H1_min_effect_size_thresh
    dp_h1_xz['pps_ProbH1_target_lwr_quantile'] = pps_ProbH1_target_lwr_quantile
    dp_h1_xz['S'] = PPS_Z_TOTAL
    perf_all = pd.DataFrame(perf_rows)
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
        os.path.join(dir_out, f'mvn_J{J}_pps_RGEC_p_h1_xz.pkl'),
    )
    perf_all.to_csv(
        os.path.join(dir_out, f'mvn_J{J}_pps_RGEC_perf.csv'), index=False,
    )
    timing_all.to_csv(
        os.path.join(dir_out, f'mvn_J{J}_pps_RGEC_timing.csv'), index=False,
    )
    pps_df.to_csv(
        os.path.join(dir_out, f'mvn_J{J}_pps_RGEC.csv'), index=False,
    )

    tab = pps_df.merge(
        pps_cf[pps_cf['J'] == J][['interim_id', 'j', 'pps']].rename(
            columns={'pps': 'pps_analytic'},
        ),
        on=['interim_id', 'j'], how='left',
    )
    tab['abs_err'] = (tab['pps'] - tab['pps_analytic']).abs()
    print(f"J={J}: RGEC vs analytic PPS: max abs {tab['abs_err'].max():.4f},"
          f" mean {tab['abs_err'].mean():.4f} over {len(tab)} (j, interim)"
          f" cells; total J time {(time.time() - t_J0) / 60.0:.2f} min.")

print("\nAll J done -- MVN xcomp amortised PPS complete.")
