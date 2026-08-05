#!/usr/bin/env python3
"""
MVN interim analysis: **amortised** PPS via a DeepSets multi-quantile
network with a **learnable** per-item encoder ``q_tau`` (§13.2 of
``dev/amortised_decision_making.md``).

**idcomp variant** -- independent-components approximation. Same
per-component decomposition as the features-fixed sibling (§13.1):
the amortiser sees ONE component's raw sequence at a time and
predicts ``mu_j`` marginally; cross-component correlation ``K`` is
used at TRAINING (mu prior) but discarded at DEPLOYMENT. A full-K
variant that exploits the known cross-component correlation at
deployment is future work.

Parallel to the features-fixed sibling: the target distribution of
``mu_j | (x, z)`` is J-invariant under the MVN g-prior with unit-
diagonal ``K``, so a single amortiser trained on prior draws from
any J in ``simu_params['J_grid']`` deploys across all J values.

Difference from the features-fixed variant: raw padded per-component
scalar sequences (item_dim = 1) are consumed by a shared
per-item MLP encoder, and the pooled embeddings feed the MLP head.

Outputs per J (``_RGEB_`` suffix; matches the MVN compare-methods
loader):

  - ``mvn_J{J}_pps_RGEB_p_h1_xz.pkl``
  - ``mvn_J{J}_pps_RGEB.csv``
  - ``mvn_J{J}_pps_RGEB_timing.csv``
  - ``mvn_J{J}_pps_RGEB_perf.csv``

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/MVN_interim_analysis_amortise_endptx_on_wz_with_features_MLP_idcomp_qpsi_MLP_loss_multiquantilehead.py
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
from amortiser_pps_features_MLP_qpsi_MLP_loss_multiquantilehead import (
    Amortiser_PPS_features_MLP_qpsi_MLP_loss_multiquantilehead,
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
    "idcomp-qpsi-MLP-loss-multiquantilehead"
)
dir_out = f"{_ROOT_DIR}-260714"
os.makedirs(dir_out, exist_ok=True)
print(f"Output dir: {dir_out}")

PPS_Z_TOTAL = 500                 # smaller than RGEA (S = 4000): each MLP
                                  # forward is on a (S, N_max, 1) tensor per
                                  # (j, interim); the JAX apply_fn is the
                                  # deployment bottleneck, so we cap S here
                                  # at the S = 500 MC-noise floor of ~ 0.02.

TRAIN_J = 20
NET_Q_TAU_HIDDEN = (32, 32)
NET_EMBED_DIM = 16
NET_HIDDEN = (64, 64)             # matches Binomial §12.2 combo default
NET_TAUS = (0.05, 0.25, 0.5, 0.75, 0.95)
NET_STEPS = 3000                  # heavy per-step cost (~350 ms/step on 6-CPU)
                                  # over a (S*J=2560, N_max=1050, 1) MLP
                                  # forward. Loss plateaus after ~2k steps
                                  # for this net width; keep the wall-clock
                                  # in the 15-20 min band.
NET_BATCH = 128                   # per-sample tensor (N_max, 1) * TRAIN_J.
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
    print(f"\nTraining features-MLP amortised net "
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
        _amortiser_model.make_training_data_with_raw_sequences,
        n_max=NET_N_MAX, n_dist=N_DIST,
    )
    fit = train(
        sample_fn,
        Amortiser_PPS_features_MLP_qpsi_MLP_loss_multiquantilehead,
        net_kwargs={
            'q_tau_hidden_dims': NET_Q_TAU_HIDDEN,
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
# Per-J deployment loop.
# =============================================================================


def _padded_scalar_seq(vals: np.ndarray, n_max: int):
    """Return ``x`` shape ``(n_max, 1)`` and ``mask`` shape ``(n_max,)``
    with the ``len(vals)`` real scalars filled in slots ``[0, len(vals))``,
    zeros padded in slots ``[len(vals), n_max)``, and the mask carrying
    1s over real slots."""
    total = int(vals.shape[0])
    x = np.zeros((n_max, 1), dtype=np.float32)
    mask = np.zeros((n_max,), dtype=np.float32)
    x[:total, 0] = vals.astype(np.float32)
    mask[:total] = 1.0
    return x, mask


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

        print(f"\n{'=' * 70}\nJ={J} RGEB interim {interim_id}"
              f" ({interim_date.date()}) n={n_obs} m={interim_m}"
              f" S={PPS_Z_TOTAL}\n{'=' * 70}")

        # ---- Per-(j) observed sequence: sorted by pid, then rows for j.
        # ``dpi`` is long-form (pid, j) with contiguous 1..N*J oid. Pivot
        # to (n_obs, J) so column j has the observed scalars for that
        # component in pid order.
        x_wide = (
            dpi.pivot_table(index='pid', columns='j', values='y')
            .sort_index()
            .reindex(columns=range(J))
            .to_numpy(dtype=np.float64)
        )                                                  # (n_obs, J)

        # ---- Per-(j, s) future sequence from ``zi``: pivot ypred_s
        # so row = pid, col = j gives z_{i, j}^(s) in the same layout.
        zi_sorted = zi_full.sort_values(['pid', 'j']).reset_index(drop=True)
        # (S, m, J) row-major over (s, pid, j).
        ypred_arr = (
            zi_sorted[keep_cols].to_numpy().T.reshape(
                PPS_Z_TOTAL, interim_m, J,
            )
        )

        # ---- Loop per j; each j has its own batch of S entries. Keeps
        # the peak tensor size to S * N_max * 1 (fits in RAM for all J
        # in the grid).
        S = PPS_Z_TOTAL
        sizes_row = np.array(
            [n_obs / NET_N_MAX, interim_m / NET_N_MAX], dtype=np.float32,
        )
        item_labels = dit['item_label'].to_numpy()
        item_types  = dit['item_type'].to_numpy()
        item_highs  = dit['item_high_label'].to_numpy()
        per_j_rows = []
        for j in range(J):
            x_j, mask_x_j = _padded_scalar_seq(x_wide[:, j], NET_N_MAX)
            x_arr = np.broadcast_to(
                x_j[None, :, :], (S, NET_N_MAX, 1),
            ).astype(np.float32)
            mask_x_arr = np.broadcast_to(
                mask_x_j[None, :], (S, NET_N_MAX),
            ).astype(np.float32)
            z_arr = np.zeros((S, NET_N_MAX, 1), dtype=np.float32)
            mask_z_arr = np.zeros((S, NET_N_MAX), dtype=np.float32)
            for s in range(S):
                z_seq, mask_z_seq = _padded_scalar_seq(
                    ypred_arr[s, :, j], NET_N_MAX,
                )
                z_arr[s] = z_seq
                mask_z_arr[s] = mask_z_seq
            sizes_arr = np.broadcast_to(
                sizes_row[None, :], (S, 2),
            ).astype(np.float32)
            batch = {
                'x':       x_arr,
                'mask_x':  mask_x_arr,
                'z':       z_arr,
                'mask_z':  mask_z_arr,
                'sizes':   sizes_arr,
            }
            p_h1_xz_j, q_hat_j, _ = predict_amortised_p_h1_for_one_xz(
                fit, batch, pps_H1_min_effect_size_thresh,
            )
            per_j_rows.append(pd.DataFrame({
                'item_label':      [item_labels[j]] * S,
                'item_type':       [item_types[j]] * S,
                'item_high_label': [item_highs[j]] * S,
                'j':               np.full(S, j, dtype=int),
                's':               np.arange(1, S + 1, dtype=int),
                'p_h1_xz':         p_h1_xz_j.astype(np.float64),
                'q_hat':           q_hat_j.astype(np.float64),
            }))
        df = pd.concat(per_j_rows, ignore_index=True)
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
        print(f"  J={J} interim {interim_id} MLP amortised done in"
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
        os.path.join(dir_out, f'mvn_J{J}_pps_RGEB_p_h1_xz.pkl'),
    )
    perf_all.to_csv(
        os.path.join(dir_out, f'mvn_J{J}_pps_RGEB_perf.csv'), index=False,
    )
    timing_all.to_csv(
        os.path.join(dir_out, f'mvn_J{J}_pps_RGEB_timing.csv'), index=False,
    )
    pps_df.to_csv(
        os.path.join(dir_out, f'mvn_J{J}_pps_RGEB.csv'), index=False,
    )

    tab = pps_df.merge(
        pps_cf[pps_cf['J'] == J][['interim_id', 'j', 'pps']].rename(
            columns={'pps': 'pps_analytic'},
        ),
        on=['interim_id', 'j'], how='left',
    )
    tab['abs_err'] = (tab['pps'] - tab['pps_analytic']).abs()
    print(f"J={J}: RGEB vs analytic PPS: max abs {tab['abs_err'].max():.4f},"
          f" mean {tab['abs_err'].mean():.4f} over {len(tab)} (j, interim)"
          f" cells; total J time {(time.time() - t_J0) / 60.0:.2f} min.")

print("\nAll J done -- MVN features-MLP amortised PPS complete.")
