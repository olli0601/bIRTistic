#!/usr/bin/env python3
"""
MVN interim analysis: **amortised** PPS via attention over the
component axis + K-row query + K-family-invariant training
(§12.11 of ``dev/amortised_decision_making.md``).

**xcompAtt** replaces the inner participant-level sum-pool of the
`xcomp` variant (§12.10) with a single-head cross-attention block over
the J components. Component summaries `(sum_i y_{i, j}, sum_i z_{i, j},
K[j*, j])` are pre-pooled per component -- MVN sufficient statistic --
so the encoder no longer runs on `(participant, component)` tokens.

Training draws `K` from a mixture of standard families per prior draw
(identity / AR(1) / block-equicorrelation / factor) so the trained net
does NOT overfit to the block-equicorr K family used at deployment.

Outputs per J (``_RGEF_`` suffix; picked up by the MVN compare-methods
loader alongside the idcomp + xcomp variants):

  - ``mvn_J{J}_pps_RGEF_p_h1_xz.pkl``
  - ``mvn_J{J}_pps_RGEF.csv``
  - ``mvn_J{J}_pps_RGEF_timing.csv``
  - ``mvn_J{J}_pps_RGEF_perf.csv``

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/MVN_interim_analysis_amortise_endptx_on_wz_with_features_MLP_xcompAtt_qpsi_MLP_loss_multiquantilehead.py
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
from amortiser_pps_features_itemXcompAtt_qpsi_MLP_loss_multiquantilehead import (
    Amortiser_PPS_features_itemXcompAtt_qpsi_MLP_loss_multiquantilehead,
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
    "py-mvn-interim-amortise-endptx-on-wz-with-features-"
    "itemXcompAtt-qpsi-MLP-loss-multiquantilehead"
)
_DIR_VARIANT = os.environ.get('XCOMPATT_DIR_VARIANT', '_260716')
dir_out = f"{_ROOT_DIR}{_DIR_VARIANT}"
os.makedirs(dir_out, exist_ok=True)
print(f"Output dir: {dir_out}")

PPS_Z_TOTAL = 500

TRAIN_J = 20
NET_Q_TOK_HIDDEN = (32, 32)
NET_EMBED_DIM = 32
NET_HIDDEN = (64, 64)
NET_TAUS = (0.05, 0.25, 0.5, 0.75, 0.95)
NET_STEPS = int(os.environ.get('XCOMPATT_STEPS', 3000))
NET_BATCH = 32                   # S=32 samples * Q=4 queries -> B=128.
NET_QUERIES = 4
NET_LR = 1e-3
N_DIST = 'uniform'
# Train on the full mixture by default so the amortiser is
# K-family-invariant.
K_FAMILIES = ('identity', 'ar1', 'block', 'factor')

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
    print(f"\nTraining xcompAtt amortised net "
          f"({NET_STEPS} steps, S={NET_BATCH}, Q={NET_QUERIES}, "
          f"TRAIN_J={TRAIN_J}, K_families={K_FAMILIES}) ...")

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
        _amortiser_model.make_training_data_with_component_summaries_random_K,
        n_max=NET_N_MAX, n_dist=N_DIST,
        queries_per_sample=NET_QUERIES,
        K_families=K_FAMILIES,
    )
    fit = train(
        sample_fn,
        Amortiser_PPS_features_itemXcompAtt_qpsi_MLP_loss_multiquantilehead,
        net_kwargs={
            'q_tok_hidden': NET_Q_TOK_HIDDEN,
            'embed_dim': NET_EMBED_DIM,
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
# Per-J deployment loop. Encoder input is per-component summary
# ``(x_sum, z_sum, K[j*, :])`` -- no participant-level tokens -- so one
# posterior draw becomes one J-way batched forward pass.
# =============================================================================

for J in J_GRID:
    cell = sim_data['cells'][J]
    dit, dp, R_chol = cell['dit'], cell['dp'], cell['R_chol']
    K = (R_chol @ R_chol.T).astype(np.float64)
    K_diag = np.diag(K).astype(np.float64)
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

        print(f"\n{'=' * 70}\nJ={J} RGEF interim {interim_id}"
              f" ({interim_date.date()}) n={n_obs} m={interim_m}"
              f" S={PPS_Z_TOTAL}\n{'=' * 70}")

        # ---- Per-component observed sum sum_i y_{i, j}.
        x_sum_obs = (
            dpi.pivot_table(index='pid', columns='j', values='y')
            .sort_index()
            .reindex(columns=range(J))
            .to_numpy(dtype=np.float64)
            .sum(axis=0)                                # (J,)
        )
        # ---- Per-component future sum per posterior draw.
        zi_sorted = zi_full.sort_values(['pid', 'j']).reset_index(drop=True)
        ypred_arr = (
            zi_sorted[keep_cols].to_numpy().T.reshape(
                PPS_Z_TOTAL, interim_m, J,
            )
        )
        z_sum_per_s = ypred_arr.sum(axis=1)             # (S, J)

        sizes_row = np.array(
            [n_obs / NET_N_MAX, interim_m / NET_N_MAX], dtype=np.float32,
        )
        scale = float(NET_N_MAX)

        # ---- Batch of J elements per posterior draw, converted to the
        # unified itemXcompAtt schema (tokens, mask, query_idx, aux).
        S = PPS_Z_TOTAL
        item_labels = dit['item_label'].to_numpy()
        item_types  = dit['item_type'].to_numpy()
        item_highs  = dit['item_high_label'].to_numpy()
        x_sum_batch = np.broadcast_to(
            (x_sum_obs / scale).astype(np.float32)[None, :], (J, J),
        ).astype(np.float32)                                    # (J_query, J_tok)
        k_row_batch = K.astype(np.float32)                      # (J_query, J_tok)
        k_diag_batch = K_diag[:, None].astype(np.float32)       # (J_query, 1)
        sizes_batch = np.broadcast_to(
            sizes_row[None, :], (J, 2),
        ).astype(np.float32)
        mask_batch = np.ones((J, J), dtype=np.float32)
        query_idx_batch = np.arange(J, dtype=np.int32)          # each batch elem
                                                                 # queries its own j*
        aux_batch = np.concatenate(
            [sizes_batch, k_diag_batch], axis=-1,
        ).astype(np.float32)                                    # (J_query, 3)
        # Query identity handled inside the class via query_bias_alpha
        # (learned scalar bias on attention scores at query_idx). No
        # explicit is_query token feature needed at deployment.

        preds_p = np.empty((S, J), dtype=np.float64)
        preds_q = np.empty((S, J), dtype=np.float64)
        med_idx = len(NET_TAUS) // 2
        for s in range(S):
            z_sum_batch = np.broadcast_to(
                (z_sum_per_s[s] / scale).astype(np.float32)[None, :],
                (J, J),
            ).astype(np.float32)
            tokens = np.stack(
                [x_sum_batch, z_sum_batch, k_row_batch], axis=-1,
            ).astype(np.float32)                                # (J_query, J_tok, 3)
            batch = {
                'tokens':    tokens,
                'mask':      mask_batch,
                'query_idx': query_idx_batch,
                'aux':       aux_batch,
            }
            p_h1_xz_s, _q_s, preds_s = predict_amortised_p_h1_for_one_xz(
                fit, batch, pps_H1_min_effect_size_thresh,
            )
            preds_p[s] = p_h1_xz_s
            preds_q[s] = preds_s[:, med_idx]

        # ---- Rebuild long form.
        rows = []
        for j in range(J):
            rows.append(pd.DataFrame({
                'item_label':      [item_labels[j]] * S,
                'item_type':       [item_types[j]] * S,
                'item_high_label': [item_highs[j]] * S,
                'j':               np.full(S, j, dtype=int),
                's':               np.arange(1, S + 1, dtype=int),
                'p_h1_xz':         preds_p[:, j].astype(np.float64),
                'q_hat':           preds_q[:, j].astype(np.float64),
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
        print(f"  J={J} interim {interim_id} xcompAtt amortised done in"
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
        os.path.join(dir_out, f'mvn_J{J}_pps_RGEF_p_h1_xz.pkl'),
    )
    perf_all.to_csv(
        os.path.join(dir_out, f'mvn_J{J}_pps_RGEF_perf.csv'), index=False,
    )
    timing_all.to_csv(
        os.path.join(dir_out, f'mvn_J{J}_pps_RGEF_timing.csv'), index=False,
    )
    pps_df.to_csv(
        os.path.join(dir_out, f'mvn_J{J}_pps_RGEF.csv'), index=False,
    )

    tab = pps_df.merge(
        pps_cf[pps_cf['J'] == J][['interim_id', 'j', 'pps']].rename(
            columns={'pps': 'pps_analytic'},
        ),
        on=['interim_id', 'j'], how='left',
    )
    tab['abs_err'] = (tab['pps'] - tab['pps_analytic']).abs()
    mse = float(((tab['pps'] - tab['pps_analytic']) ** 2).mean())
    print(f"J={J}: RGEF vs analytic PPS: max abs {tab['abs_err'].max():.4f},"
          f" mean {tab['abs_err'].mean():.4f}, MSE {mse:.5f}"
          f" over {len(tab)} (j, interim) cells;"
          f" total J time {(time.time() - t_J0) / 60.0:.2f} min.")

print("\nAll J done -- MVN xcompAtt amortised PPS complete.")
