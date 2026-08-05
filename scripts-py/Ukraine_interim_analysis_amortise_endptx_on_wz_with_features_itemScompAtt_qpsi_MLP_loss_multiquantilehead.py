#!/usr/bin/env python3
"""
Ukraine interim analysis: **amortised** PPS via self-attention over the
20 PCM items + multi-quantile head (§14.1 of
``dev/amortised_decision_making.md``).

Adapts the MVN `xcompAtt` (§13.5) to the Ukraine partial-credit model.
Training data are the cached Strong-Oakley `wa` frames from
``Ukraine_interim_analysis_regression_endptx_on_wz_with_gauss_approx.py``:
each row already carries per-(item, draw) ``(w_baseline, w_endline,
w_ratio)`` and target ``pps_ratio_x``. Amortiser learns
``pps_ratio_x ~ (all items' summaries)`` in one shot, sharing weights
across items and letting attention over the 20 items propagate
cross-item information -- the analog of the K-row query in MVN,
replaced here by fully-learned attention because PCM has no closed-form
cross-item covariance.

Training uses a **leave-one-interim-out** protocol: for each held-out
interim, the amortiser is trained on the other seven interims' wa
frames and then applied to the held-out interim's rows for PPS
evaluation. This mirrors the "amortise across cohorts" story: one net
trained once, deployed at any new interim.

Outputs (``_RGEE_`` suffix; matches the Ukraine compare-methods
loader):

  - ``pcm_1_interim_pps_RGEE_p_h1_xz.pkl``
  - ``pcm_1_interim_pps_RGEE.csv``
  - ``pcm_1_interim_pps_RGEE_timing.csv``
  - ``pcm_1_interim_pps_RGEE_perf_long.pkl``

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/Ukraine_interim_analysis_amortise_endptx_on_wz_with_xcompAtt_qpsi_MLP_loss_multiquantilehead.py
"""

# %%

import os
import sys
import time
import warnings
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
from plotnine import (
    ggplot, aes, geom_boxplot, geom_density, geom_hline, geom_point,
    geom_smooth, geom_text, geom_vline,
    facet_wrap, position_dodge,
    scale_colour_manual, scale_fill_manual, scale_y_continuous,
    theme_bw, theme, element_blank, element_text, element_rect, labs,
)

warnings.filterwarnings('ignore')

from utils import _futurama_palette
from amortiser_common import (
    train as amortiser_train,
    save_trained_model,
    load_fitted_model,
    predict_amortised_p_h1_for_one_xz,
)
from amortiser_pps_features_itemScompAtt_qpsi_MLP_loss_multiquantilehead import (
    Amortiser_PPS_features_itemScompAtt_qpsi_MLP_loss_multiquantilehead,
)
from model_pcm import PartialCreditModel

print("Imports successful")

# %%

# =============================================================================
# Configuration
# =============================================================================

DIR_RGE = os.path.join(
    "/Users/or105/sandbox/bIRTistic",
    "py-ukraine-interim-with-regression-on-endptx-wz-260601",
)
dir_out = os.path.join(
    "/Users/or105/sandbox/bIRTistic",
    "py-ukraine-interim-amortise-endptx-on-wz-with-features-itemScompAtt-qpsi-MLP-loss-multiquantilehead-260716",
)
os.makedirs(dir_out, exist_ok=True)
print(f"Output dir: {dir_out}")

file_prefix = "pcm_1_interim"

pps_H1_min_effect_size_thresh = 0.5
pps_ProbH1_target_lwr_quantile = 0.89

NET_STEPS = int(os.environ.get('UKR_XCOMPATT_STEPS', 6000))
NET_BATCH = 256
NET_LR = 1e-3
NET_TAUS = (0.05, 0.25, 0.5, 0.75, 0.95)
NET_TOK_HIDDEN = (32, 32)
NET_EMBED_DIM = 32
NET_HIDDEN = (64, 64)
NET_SEED = 123
# Clip the pps_ratio_x target so the single numerical-blowup outlier at
# CG-VIO_ph-punish interim 1 (value -2.3e12) doesn't dominate the
# pinball loss. Legitimate targets reach [-22, +17] across items (99th
# percentile ~7); a symmetric clip at 20 removes only the blowup while
# preserving the true target spread so the quantile head can learn to
# predict values in the legitimate tail.
TARGET_CLIP = 20.0

# %%

# =============================================================================
# Load all interims' wa. Reshape to a (draws, items, features) tensor
# per interim keyed by interim_id.
# =============================================================================

interim_wa = {}
for i in range(1, 12):
    p = os.path.join(DIR_RGE, f"{file_prefix}_i{i}_regression_training.pkl")
    if not os.path.exists(p):
        continue
    interim_wa[i] = pd.read_pickle(p)
INTERIM_IDS = sorted(interim_wa.keys())
print(f"Loaded interims: {INTERIM_IDS}")

# Stable item order across interims.
items_df = (
    interim_wa[INTERIM_IDS[0]][['item_label', 'item_type', 'item_high_label']]
    .drop_duplicates()
    .sort_values(['item_type', 'item_label'])
    .reset_index(drop=True)
)
items_df['item_pos'] = np.arange(len(items_df), dtype=int)
J = len(items_df)
item_pos = dict(zip(items_df['item_label'], items_df['item_pos']))
item_type_idx = items_df['item_type'].map(
    {'out-of-7': 0.0, 'categorical': 1.0}
).to_numpy(dtype=np.float32)
item_high_idx = items_df['item_high_label'].map(
    {'higher_is_better': 1.0, 'lower_is_better': 0.0}
).to_numpy(dtype=np.float32)
print(f"J = {J} items:\n{items_df}")

# Full dit frame (needed to rebuild the PCM model for the z^(s) block,
# and the per-item response-level counts).
dit_df = pd.read_csv(os.path.join(DIR_RGE, f"{file_prefix}_1_data_dit.csv"))
items_df = items_df.merge(
    dit_df[['item_label', 'cat_length']].drop_duplicates(),
    on='item_label', how='left',
)
# Per-item response-level count from dit['cat_length'] (Ukraine: 8 for
# out-of-7 with raw y in 0..7, 4 for categorical with raw y in 0..3).
# Study-specific -- Colombia has different level counts (and hence a
# different caseness threshold c, which enters the TARGET only).
item_klevels = items_df['cat_length'].to_numpy(dtype=np.float32)
assert not np.isnan(item_klevels).any(), "cat_length missing for some items"


N_FULL = 503                      # total trial participants (n + m)
KHAT_SHRINK_N0 = 50.0             # shrinkage prior weight for empirical K
F_TOK = 8                         # per-item token feature dim (see below)


def _rescale_w(w):
    """wa-frame out-of-7 w's are means of 1-INDEXED categories in
    [1, K_j]; map to [0, 1] via (w - 1) / (K_j - 1) with K_j =
    dit['cat_length'] per item. Categorical columns are overwritten
    with mean-response summaries downstream."""
    o7 = (item_type_idx == 0.0)[None, :]
    return np.where(
        o7, (w - 1.0) / (item_klevels[None, :] - 1.0), w,
    ).astype(np.float32)


def _interim_x_features(interim_id):
    """Per-interim x-side (observed-cohort) features from the cached xi
    frame: (w_x [J, 2] in [0, 1], khat [J, J], n_pid).

    w_x: per-item baseline/endline group-mean RESPONSE over the n
    observed participants, rescaled by the item type's level count:
    mean(y) / (K_j - 1) with K_j = dit['cat_length'] per item and raw
    y in 0..K_j-1. No input-side thresholding (the caseness threshold
    lives in the target only).
    khat: empirical cross-item correlation of per-participant CHANGE
    scores d_ij = y_end - y_base (Spearman; robust to the ordinal
    scale), shrunk toward identity with weight n / (n + n0) so early
    interims rely less on the noisy estimate.
    """
    xi = pd.read_csv(
        os.path.join(DIR_RGE, f"{file_prefix}_{interim_id}_data_dp1.csv"),
    )
    piv = xi.pivot_table(index='pid', columns=['item_label', 'time'],
                         values='y')
    n_pid = piv.index.size
    w_x = np.zeros((J, 2), dtype=np.float32)
    chg = np.zeros((n_pid, J), dtype=np.float64)
    for j, lbl in enumerate(items_df['item_label']):
        y0 = piv[(lbl, 0)].to_numpy(dtype=np.float64)
        y1 = piv[(lbl, 1)].to_numpy(dtype=np.float64)
        kmax = item_klevels[j] - 1.0
        w_x[j, 0] = np.nanmean(y0) / kmax
        w_x[j, 1] = np.nanmean(y1) / kmax
        chg[:, j] = y1 - y0
    rho = pd.DataFrame(chg).corr(method='spearman').to_numpy()
    rho = np.nan_to_num(rho, nan=0.0)
    np.fill_diagonal(rho, 1.0)
    lam = n_pid / (n_pid + KHAT_SHRINK_N0)
    khat = (lam * rho + (1.0 - lam) * np.eye(J)).astype(np.float32)
    return w_x, khat, n_pid


def _interim_z_cat_stats(interim_id, D):
    """Per-draw mean-response z-side summaries for the CATEGORICAL
    items, which the cached wa frames only carry as thresholded
    proportions. Rebuild the posterior-predictive block z^(s) from the
    cached SVI draws (same seed + keep_order as the RGE script, so the
    shadow cohort and draw alignment match the wa targets) and compute
    mean(ypred - 1) / (K_j - 1) in [0, 1] per (draw, item, time), with
    K_j = dit['cat_length'] per item.

    Returns (D, n_cat, 2) for the items where item_type_idx == 1, in
    items_df order."""
    cat_j = [j for j in range(J) if item_type_idx[j] == 1.0]
    cat_labels = [items_df['item_label'].iloc[j] for j in cat_j]
    if not cat_labels:
        return np.zeros((D, 0, 2), dtype=np.float32)
    xi = pd.read_csv(
        os.path.join(DIR_RGE, f"{file_prefix}_{interim_id}_data_dp1.csv"),
    )
    n_pid = xi['pid'].nunique()
    interim_m = N_FULL - n_pid
    model = PartialCreditModel(dit=dit_df, dcati=xi, x_formula="~ time - 1",
                               seed=123)
    zi = model.get_interim_z_from_ypredi(
        os.path.join(DIR_RGE, f"{file_prefix}_{interim_id}_draws.zarr"),
        interim_m, pps_z_total=D, seed=123, keep_order=True,
    )
    ypred_cols = [f'ypred_{s}' for s in range(D)]
    out = np.zeros((D, len(cat_labels), 2), dtype=np.float32)
    for jj, lbl in enumerate(cat_labels):
        kmax = item_klevels[cat_j[jj]] - 1.0
        for t in (0, 1):
            sub = zi[(zi['item_label'] == lbl) & (zi['time'] == t)]
            vals = sub[ypred_cols].to_numpy(dtype=np.float32)   # (m, D), 1..K
            out[:, jj, t] = (vals - 1.0).mean(axis=0) / kmax
    return out


def _tensorise_interim(wa: pd.DataFrame, interim_id: int):
    """Return (tokens_base [D, J, 7], target [D, J], khat [J, J],
    aux_row [2]) for one interim.

    Token features (F = 8; the 8th, khat[j*, j], is query-dependent and
    appended at batch-construction time):
      0: w_x_baseline_j   -- OBSERVED cohort, s-invariant
      1: w_x_endline_j    -- OBSERVED cohort, s-invariant
      2: w_z_baseline_j^(s) -- future block z^(s)
      3: w_z_endline_j^(s)  -- future block z^(s)
      4: w_z_ratio_j^(s)    -- future block z^(s)
      5: c_j (item type), 6: d_j (direction)
    aux = (n / N_FULL, m / N_FULL).
    """
    w_x, khat, n_pid = _interim_x_features(interim_id)
    wa = wa.copy()
    wa['item_pos'] = wa['item_label'].map(item_pos)
    piv = wa.pivot_table(
        index='draw', columns='item_pos',
        values=['w_baseline', 'w_endline', 'w_ratio', 'pps_ratio_x'],
    )
    D = piv.index.size
    w_base = _rescale_w(
        piv['w_baseline'].reindex(columns=range(J)).to_numpy(dtype=np.float32))
    w_end = _rescale_w(
        piv['w_endline'].reindex(columns=range(J)).to_numpy(dtype=np.float32))
    w_rat = np.array(
        piv['w_ratio'].reindex(columns=range(J)).to_numpy(dtype=np.float32),
        copy=True,
    )
    # Categorical items: replace the thresholded-proportion z-side
    # summaries from the wa frames with mean-response summaries rebuilt
    # from the cached SVI draws (no input-side thresholding). The change
    # ratio is recomputed on the means (all CG-MH are lower_is_better).
    cat_pos = np.flatnonzero(item_type_idx == 1.0)
    if cat_pos.size:
        zc = _interim_z_cat_stats(interim_id, D)             # (D, n_cat, 2)
        w_base[:, cat_pos] = zc[:, :, 0]
        w_end[:, cat_pos] = zc[:, :, 1]
        wb = np.maximum(zc[:, :, 0], 1e-3)
        we = zc[:, :, 1]
        high = (item_high_idx[cat_pos] > 0.5)[None, :]
        w_rat[:, cat_pos] = np.where(high, we / wb - 1.0, 1.0 - we / wb)
    w_x_b = np.broadcast_to(w_x[None, :, 0], (D, J))
    w_x_e = np.broadcast_to(w_x[None, :, 1], (D, J))
    feats = np.stack([w_x_b, w_x_e, w_base, w_end, w_rat], axis=-1)  # (D, J, 5)
    target = piv['pps_ratio_x'].reindex(
        columns=range(J),
    ).to_numpy(dtype=np.float32)                            # (D, J)
    feats = np.nan_to_num(feats, nan=0.0, posinf=0.0, neginf=0.0)
    target = np.nan_to_num(target, nan=0.0, posinf=TARGET_CLIP, neginf=-TARGET_CLIP)
    target = np.clip(target, -TARGET_CLIP, TARGET_CLIP)
    metadata = np.stack([item_type_idx, item_high_idx], axis=-1)  # (J, 2)
    metadata_b = np.broadcast_to(
        metadata[None, :, :], (D, J, 2),
    ).astype(np.float32)
    tokens = np.concatenate([feats, metadata_b], axis=-1)   # (D, J, 7)
    aux_row = np.array([n_pid / N_FULL, (N_FULL - n_pid) / N_FULL],
                       dtype=np.float32)
    return tokens.astype(np.float32), target, khat, aux_row


def _fanout_with_khat(tok, khat, aux_row):
    """tok (D, J, 7), khat (J, J), aux_row (2,) -> per-query batch
    tensors: tokens (D*J, J, 8) with khat[q, j] as feature 8, mask,
    query_idx, aux (D*J, 2)."""
    D = tok.shape[0]
    B = D * J
    tokens_rep = np.empty((D, J, J, F_TOK), dtype=np.float32)
    tokens_rep[..., :7] = tok[:, None, :, :]
    tokens_rep[..., 7] = np.broadcast_to(khat[None, :, :], (D, J, J))
    tokens_rep = tokens_rep.reshape(B, J, F_TOK)
    mask_rep = np.ones((B, J), dtype=np.float32)
    query_idx = np.broadcast_to(
        np.arange(J)[None, :], (D, J),
    ).reshape(B).astype(np.int32)
    aux_rep = np.broadcast_to(aux_row[None, :], (B, 2)).astype(np.float32)
    return tokens_rep, mask_rep, query_idx, aux_rep


interim_tensors = {
    i: _tensorise_interim(wa, i) for i, wa in interim_wa.items()
}
for i, (t, y, kh, ax) in interim_tensors.items():
    print(f"  interim {i}: tokens {t.shape}, target {y.shape}, "
          f"khat mean|offdiag| "
          f"{np.abs(kh[~np.eye(J, dtype=bool)]).mean():.3f}, aux {ax}")

# %%

# =============================================================================
# Build training pool: concatenate all interims. At training we sample
# uniformly across (interim, draw); each sample fans out to J = 20
# per-item rows (one query per row).
# =============================================================================

train_tokens = np.concatenate(
    [interim_tensors[i][0] for i in INTERIM_IDS], axis=0,
)                                                           # (D_total, J, 7)
train_target = np.concatenate(
    [interim_tensors[i][1] for i in INTERIM_IDS], axis=0,
)                                                           # (D_total, J)
train_interim_ii = np.concatenate([
    np.full(interim_tensors[i][0].shape[0], k, dtype=np.int32)
    for k, i in enumerate(INTERIM_IDS)
])                                                          # (D_total,)
khat_arr = np.stack(
    [interim_tensors[i][2] for i in INTERIM_IDS], axis=0,
)                                                           # (C, J, J)
aux_arr = np.stack(
    [interim_tensors[i][3] for i in INTERIM_IDS], axis=0,
)                                                           # (C, 2)
D_total = train_tokens.shape[0]
print(f"\nTraining pool: {D_total} draws across {len(INTERIM_IDS)} interims.")

# Per-item target standardisation (§14.1 CG-MH resolution fix).
ITEM_STD_PATH = os.path.join(dir_out, f"{file_prefix}_item_std.npy")
if os.path.exists(ITEM_STD_PATH):
    item_std = np.load(ITEM_STD_PATH).astype(np.float32)
    assert item_std.shape == (J,)
    print(f"Loaded cached per-item target std: {ITEM_STD_PATH}")
else:
    item_std = np.nanstd(train_target, axis=0).astype(np.float32)
    item_std = np.where(item_std > 1e-6, item_std, 1.0).astype(np.float32)
    np.save(ITEM_STD_PATH, item_std)
    print(f"Computed + saved per-item target std to {ITEM_STD_PATH}")
print(f"item_std range [{item_std.min():.3g}, {item_std.max():.3g}], "
      f"mean {item_std.mean():.3g}")


def sample_fn(rng, S):
    """Sample S draws from the pooled tensors. Fan out to B = S * J
    per-item rows -- one query index per row, khat[q, :] appended as
    the query-dependent 8th token feature. Target is standardised per
    item by item_std."""
    draw_idx = rng.integers(0, D_total, size=S)              # (S,)
    tok = train_tokens[draw_idx]                            # (S, J, 7)
    kk = khat_arr[train_interim_ii[draw_idx]]               # (S, J, J)
    ax = aux_arr[train_interim_ii[draw_idx]]                # (S, 2)
    tgt = train_target[draw_idx] / item_std[None, :]        # (S, J) standardised
    B = S * J
    tokens_rep = np.empty((S, J, J, F_TOK), dtype=np.float32)
    tokens_rep[..., :7] = tok[:, None, :, :]
    tokens_rep[..., 7] = kk
    tokens_rep = tokens_rep.reshape(B, J, F_TOK)
    mask_rep = np.ones((B, J), dtype=np.float32)
    query_idx = np.broadcast_to(
        np.arange(J)[None, :], (S, J),
    ).reshape(B).astype(np.int32)
    aux_rep = np.broadcast_to(
        ax[:, None, :], (S, J, 2),
    ).reshape(B, 2).astype(np.float32)
    target_rho = tgt.reshape(B).astype(np.float32)
    return (
        {
            'tokens': tokens_rep,
            'mask': mask_rep,
            'query_idx': query_idx,
            'aux': aux_rep,
        },
        target_rho,
    )


def _predict_with_unstandardise(fit, batch, thresh):
    """Forward-pass batch, multiply predicted quantiles by
    item_std[query_idx] (undo per-item target standardisation), then
    apply monotone + interp to get p_h1_xz on the true target scale."""
    apply_fn = fit['apply_fn']
    params   = fit['params']
    taus     = np.asarray(fit['pps_ProbH1_lwr_quantiles_mesh'])
    from amortiser_common import _jax_pytree
    preds = np.asarray(apply_fn(params, _jax_pytree(batch)))    # (B, K_tau) standardised
    scale = item_std[np.asarray(batch['query_idx']).astype(int)]
    preds = preds * scale[:, None]                              # unstandardise
    preds = np.maximum.accumulate(preds, axis=1)
    F = np.array([
        float(np.interp(thresh, row, taus, left=0.0, right=1.0))
        for row in preds
    ])
    p_h1_xz = 1.0 - F
    med_idx = int(np.argmin(np.abs(taus - 0.5)))
    q_hat = preds[:, med_idx]
    return p_h1_xz, q_hat, preds

# %%

# =============================================================================
# Train amortised net once (or load).
# =============================================================================

NET_CKPT = os.path.join(dir_out, f"{file_prefix}_amortised_pps_net.pkl")

if os.path.exists(NET_CKPT):
    print(f"\nLoading cached amortised net from {NET_CKPT}")
    fit = load_fitted_model(NET_CKPT)
else:
    print(f"\nTraining Ukraine xcompAtt amortised net "
          f"({NET_STEPS} steps, batch S={NET_BATCH} -> B={NET_BATCH * J}) ...")
    fit = amortiser_train(
        sample_fn,
        Amortiser_PPS_features_itemScompAtt_qpsi_MLP_loss_multiquantilehead,
        net_kwargs={
            'q_tok_hidden': NET_TOK_HIDDEN,
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
# Per-interim deployment. For each interim, load its wa tensors, run one
# forward pass per draw with J = 20 query rows, aggregate into p_h1_xz +
# PPS. Same schema as the mquantile / gauss-approx regression scripts.
# =============================================================================

p_h1_xz_rows = []
pps_timing_rows = []
perf_rows = []
t_all0 = time.time()
for interim_id in INTERIM_IDS:
    t0 = time.time()
    tokens, _target, khat_i, aux_i = interim_tensors[interim_id]
    D = tokens.shape[0]
    S = D
    tokens_rep, mask_rep, query_idx, aux_rep = _fanout_with_khat(
        tokens, khat_i, aux_i,
    )
    batch = {
        'tokens':   tokens_rep,
        'mask':     mask_rep,
        'query_idx': query_idx,
        'aux':      aux_rep,
    }
    p_h1_xz, q_hat, _ = _predict_with_unstandardise(
        fit, batch, pps_H1_min_effect_size_thresh,
    )
    p_h1_xz = p_h1_xz.reshape(D, J)
    q_hat = q_hat.reshape(D, J)

    # Attach the observed pps_ratio_x target (from the source wa frame,
    # column already carries the per-(draw, item) ground truth) so we
    # can compute per-(interim, item) rho + R^2 in the aggregate step.
    wa_this = interim_wa[interim_id][
        ['draw', 'item_label', 'pps_ratio_x']
    ].copy()
    wa_this = wa_this.merge(
        items_df[['item_label', 'item_pos']], on='item_label', how='left',
    )
    target_per_pos = (
        wa_this.pivot_table(
            index='draw', columns='item_pos', values='pps_ratio_x',
        ).reindex(columns=range(J)).to_numpy(dtype=np.float64)
    )
    # Build long form (adds pps_ratio_x column for downstream rho/R^2).
    rows = []
    for j in range(J):
        rows.append(pd.DataFrame({
            'item_label':      [items_df.iloc[j]['item_label']] * D,
            'item_type':       [items_df.iloc[j]['item_type']] * D,
            'item_high_label': [items_df.iloc[j]['item_high_label']] * D,
            's':               np.arange(1, D + 1, dtype=int),
            'p_h1_xz':         p_h1_xz[:, j].astype(np.float64),
            'q_hat':           q_hat[:, j].astype(np.float64),
            'pps_ratio_x':     target_per_pos[:, j],
        }))
    df = pd.concat(rows, ignore_index=True)
    df['interim_id'] = interim_id
    p_h1_xz_rows.append(df)

    mins_interim = (time.time() - t0) / 60.0
    pps_timing_rows.append({
        'interim_id':      interim_id,
        'mins_interim_id': round(mins_interim, 3),
    })
    print(f"  interim {interim_id} xcompAtt amortised done in "
          f"{mins_interim:.3f} min (D={D} draws)")

# %%

# =============================================================================
# Aggregate + save.
# =============================================================================

dp_h1_xz = pd.concat(p_h1_xz_rows, ignore_index=True)
dp_h1_xz['pps_H1_min_effect_size_thresh'] = pps_H1_min_effect_size_thresh
dp_h1_xz['pps_ProbH1_target_lwr_quantile'] = pps_ProbH1_target_lwr_quantile
pps_timing = pd.DataFrame(pps_timing_rows)
pps_timing['mins_total'] = round((time.time() - t_all0) / 60.0, 3)

pkl_path = os.path.join(dir_out, f"{file_prefix}_pps_RGEE_p_h1_xz.pkl")
dp_h1_xz.to_pickle(pkl_path)
print(f"Saved RGEE P(H_1 | x, z) to: {pkl_path}")
pps_timing.to_csv(
    os.path.join(dir_out, f"{file_prefix}_pps_RGEE_timing.csv"), index=False,
)

pps_df = (
    dp_h1_xz.groupby(['interim_id', 'item_label', 'item_type',
                      'item_high_label'], observed=True)['p_h1_xz']
    .apply(lambda p: float(
        (p > pps_ProbH1_target_lwr_quantile).mean()
    ))
    .reset_index(name='pps')
)
pps_df['eta'] = pps_ProbH1_target_lwr_quantile
pps_df.to_csv(
    os.path.join(dir_out, f"{file_prefix}_pps_RGEE.csv"), index=False,
)
print(f"\nPPS table head:\n{pps_df.head(20).to_string(index=False)}")

# =============================================================================
# Per-(interim, item) rho + R^2 of predicted median quantile q_hat vs the
# observed pps_ratio_x. Mirrors the RGEM ``pps_RGEM_perf.csv`` schema so
# the compare-methods loader treats xcompAtt uniformly.
# =============================================================================

perf_rows = []
for (interim_id, item_label), g in dp_h1_xz.groupby(
    ['interim_id', 'item_label'], observed=True,
):
    y = g['pps_ratio_x'].to_numpy(dtype=np.float64)
    yh = g['q_hat'].to_numpy(dtype=np.float64)
    ok = np.isfinite(y) & np.isfinite(yh)
    if ok.sum() < 3:
        rho, r2 = float('nan'), float('nan')
    else:
        y, yh = y[ok], yh[ok]
        # Pearson rho + R^2 = 1 - SS_res / SS_tot.
        if y.std(ddof=0) == 0 or yh.std(ddof=0) == 0:
            rho = float('nan')
        else:
            rho = float(np.corrcoef(y, yh)[0, 1])
        ss_res = float(np.sum((y - yh) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float('nan')
    perf_rows.append({
        'item_label':  item_label,
        'item_type':   g['item_type'].iloc[0],
        'rho':         rho,
        'r2':          r2,
        'interim_id':  int(interim_id),
    })
perf = pd.DataFrame(perf_rows)

# Reuse the (interim_id -> interim_month_year, interim_date) mapping from
# the RGE perf-long pkl.
rge_perf_long = pd.read_pickle(
    os.path.join(DIR_RGE, f"{file_prefix}_pps_RGE_perf_long.pkl"),
)
_interim_map = (
    rge_perf_long[['interim_id', 'interim_month_year']]
    .drop_duplicates()
    .reset_index(drop=True)
)
# Grab interim_date from any regression-training pkl (all interims stamped
# on the same schema upstream).
_date_map = pd.DataFrame({
    'interim_id':   list(interim_wa.keys()),
    'interim_date': [pd.NaT] * len(interim_wa),
})
# Try to enrich with actual dates via RGE p_h1_xz pkl.
try:
    _rge_p = pd.read_pickle(
        os.path.join(DIR_RGE, f"{file_prefix}_pps_RGE_p_h1_xz.pkl"),
    )[['interim_id', 'interim_date']].drop_duplicates()
    _date_map = _rge_p
except Exception:
    pass
perf = perf.merge(_date_map, on='interim_id', how='left')
perf_csv = os.path.join(dir_out, f"{file_prefix}_pps_RGEE_perf.csv")
perf.to_csv(perf_csv, index=False)
print(f"Saved RGEE perf CSV to: {perf_csv}")

# =============================================================================
# Perf long-form: rho + R^2 (per item x interim) + time (per interim, broadcast
# across items). Matches the RGEM perf_long schema.
# =============================================================================

perf = perf.merge(_interim_map, on='interim_id', how='left')
perf = perf.merge(
    pps_timing[['interim_id', 'mins_interim_id']], on='interim_id', how='left',
)
metric_labels = {'rho': 'rho', 'r2': 'R^2', 'mins_interim_id': 'time (min)'}
perf_long = perf.melt(
    id_vars=['interim_id', 'interim_month_year', 'item_label'],
    value_vars=list(metric_labels.keys()),
    var_name='metric', value_name='value',
)
perf_long['metric'] = pd.Categorical(
    perf_long['metric'].map(metric_labels),
    categories=list(metric_labels.values()), ordered=True,
)
perf_long['method'] = 'Regression endpt amortised xcompAtt (Strong-Oakley)'
perf_long.to_pickle(
    os.path.join(dir_out, f"{file_prefix}_pps_RGEE_perf_long.pkl"),
)
print(f"Saved RGEE perf long-form pkl.")

# =============================================================================
# Box-stats + boxplot PDF of p(H_1 | x, z) per (item, interim). Mirrors
# ``pps_RGEM_p_h1_xz_boxplot.pkl / .pdf`` layout.
# =============================================================================

# Load dit metadata (group_label_long + item_label_short) from the RGE dir's
# cached data csv for interim 1 (schema stable across interims).
_dit_path = os.path.join(DIR_RGE, f'{file_prefix}_1_data_dit.csv')
try:
    dit_meta = pd.read_csv(_dit_path)[
        ['item_label', 'group_label_long', 'item_label_short']
    ].drop_duplicates()
except Exception as e:
    print(f"[warn] could not read dit meta at {_dit_path}: {e};"
          f" falling back to item_label as item_label_long.")
    dit_meta = pd.DataFrame({
        'item_label':       items_df['item_label'],
        'group_label_long': items_df['item_label'],
        'item_label_short': [''] * len(items_df),
    })

tmp = dp_h1_xz.merge(dit_meta, on='item_label', how='left')
tmp['item_label_long'] = tmp['group_label_long'].fillna(tmp['item_label']) + (
    np.where(
        tmp['item_label_short'].fillna('').astype(str).str.len() > 0,
        '\n' + tmp['item_label_short'].fillna('').astype(str),
        '',
    )
)
tmp = tmp.merge(_interim_map, on='interim_id', how='left')
tmp = tmp.merge(
    rge_perf_long[['interim_id', 'interim_month_year']].drop_duplicates(),
    on=['interim_id', 'interim_month_year'], how='left',
)
_interim_order = (
    rge_perf_long[['interim_id', 'interim_month_year']]
    .drop_duplicates()
    .sort_values('interim_id')['interim_month_year'].tolist()
)
tmp['interim_month_year'] = pd.Categorical(
    tmp['interim_month_year'], categories=_interim_order, ordered=True,
)

box_stats = (
    tmp.groupby(['item_label_long', 'interim_month_year'],
                observed=True)['p_h1_xz']
    .agg(
        min='min',
        q025=lambda s: s.quantile(0.025),
        q25 =lambda s: s.quantile(0.25),
        q50 =lambda s: s.quantile(0.50),
        q75 =lambda s: s.quantile(0.75),
        q975=lambda s: s.quantile(0.975),
        max='max',
    )
    .reset_index()
)
box_pkl = os.path.join(dir_out, f"{file_prefix}_pps_RGEE_p_h1_xz_boxplot.pkl")
box_stats.to_pickle(box_pkl)
print(f"Saved RGEE p(H_1 | x, z) box-stats to: {box_pkl}")

all_items = list(pd.unique(tmp['item_label_long']))
color_dict = dict(zip(all_items, _futurama_palette(len(all_items))))

p = (
    ggplot(box_stats, aes(x='interim_month_year', fill='item_label_long'))
    + geom_boxplot(
        aes(ymin='q025', lower='q25', middle='q50', upper='q75', ymax='q975',
            group='interim_month_year'),
        stat='identity',
    )
    + geom_hline(yintercept=pps_ProbH1_target_lwr_quantile,
                 colour='black', size=1.5)
    + facet_wrap('~ item_label_long', ncol=4)
    + scale_fill_manual(values=color_dict)
    + scale_y_continuous(
        limits=[0, 1],
        breaks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        labels=['0%', '20%', '40%', '60%', '80%', '100%'],
    )
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        legend_position='none',
        figure_size=(15, 15),
    )
    + labs(x='Interim',
           y='p(H_1 | x, z)  [Regression endpt amortised xcompAtt]')
)
box_pdf = os.path.join(dir_out, f"{file_prefix}_pps_RGEE_p_h1_xz_boxplot.pdf")
p.save(box_pdf, verbose=False, limitsize=False)
print(f"Saved RGEE p(H_1 | x, z) boxplot to: {box_pdf}")

# =============================================================================
# Per-interim diagnostic scatter: SVI posterior draw of rho (rendered as
# pps_ratio_x from p(theta | x) SVI draws) vs amortiser predicted median
# q_psi(rho | x), faceted per item, with per-item Pearson rho + R^2
# annotation. This is the direct "SVI vs amortiser" diagnostic: y-axis
# is truth (the SVI posterior sample of the endpoint ratio for that
# (interim, item, draw)); x-axis is the amortiser's median-quantile
# prediction of the same quantity conditional on the n interim
# observations x only.
# =============================================================================

# Single boxplot summary across all interims + items, mirroring the
# `mvn_J{J}_pps_RGEM_p_h1_xz_all.pdf` layout: x = interim_month_year,
# y = rho, one box per (interim, source) side by side, facet per item.
_source_svi = 'SVI  p(rho|x)'
_source_amo = 'amortiser  q_psi(rho|x)'
_dp_with_month = dp_h1_xz.merge(_interim_map, on='interim_id', how='left')
svi = _dp_with_month[[
    'interim_id', 'interim_month_year', 'item_label', 'pps_ratio_x',
]].rename(columns={'pps_ratio_x': 'rho'})
svi['source'] = _source_svi
amo = _dp_with_month[[
    'interim_id', 'interim_month_year', 'item_label', 'q_hat',
]].rename(columns={'q_hat': 'rho'})
amo['source'] = _source_amo
box_all = pd.concat([svi, amo], ignore_index=True)
_month_order = (
    _dp_with_month.sort_values('interim_id')['interim_month_year'].drop_duplicates().tolist()
)
box_all['interim_month_year'] = pd.Categorical(
    box_all['interim_month_year'], categories=_month_order, ordered=True,
)
box_all['source'] = pd.Categorical(
    box_all['source'], categories=[_source_svi, _source_amo], ordered=True,
)
# Precomputed quantiles: box = 25-75%, whiskers = 2.5-97.5%, no outliers.
box_stats = (
    box_all.groupby(['interim_month_year', 'item_label', 'source'], observed=True)['rho']
    .quantile([0.025, 0.25, 0.5, 0.75, 0.975])
    .unstack().reset_index()
    .rename(columns={0.025: 'q025', 0.25: 'q25', 0.5: 'q50',
                     0.75: 'q75', 0.975: 'q975'})
)
box_stats['grp'] = (
    box_stats['interim_month_year'].astype(str) + '|' + box_stats['source'].astype(str)
)
_items_sorted = sorted(dp_h1_xz['item_label'].unique())
_n_col = 4
_n_row = int(np.ceil(len(_items_sorted) / _n_col))
p = (
    ggplot(box_stats, aes(x='interim_month_year',
                          ymin='q025', lower='q25', middle='q50',
                          upper='q75', ymax='q975',
                          fill='source', group='grp'))
    + geom_boxplot(stat='identity',
                   position=position_dodge(width=0.8), width=0.7,
                   colour='#404040', size=0.3)
    + scale_fill_manual(values={_source_svi: '#1f77b4',
                                 _source_amo: '#d62728'})
    + facet_wrap('~ item_label', ncol=_n_col, scales='free_y')
    + theme_bw()
    + theme(axis_text_x=element_text(angle=45, vjust=1, hjust=1),
            figure_size=(3.0 * _n_col, 2.5 * _n_row),
            legend_position='top',
            strip_background=element_blank(),
            strip_text=element_text(face='bold'))
    + labs(x='Interim', y='rho  (box 25-75%, whiskers 2.5-97.5%)',
           title='SVI posterior of rho|x vs amortiser q_psi(rho|x) across interims')
)
pdf_path = os.path.join(
    dir_out, f"{file_prefix}_pps_svi_vs_amortiser_rho_all.pdf",
)
p.save(pdf_path, verbose=False, limitsize=False)
print(f"  saved combined boxplot: {pdf_path}")

# =============================================================================
# Attention-picked per-query summary scatter (§13.6 diagnostic).
# For each (interim, posterior draw s, item j*), compute the query-
# conditional attention output that q_psi actually consumes:
#     bar_h^(j*) = self_attention(h)[j*] in R^E   (self-attn + gather + B)
# then project all (interim, draw, item) bar_h vectors to R via a single
# global PCA (first component). Scatter y = SVI posterior draw of rho
# for that (interim, draw, item) vs x = PC1(bar_h^(j*)); facet per item.
# Panel Pearson rho then measures how well the per-query learned summary
# tracks the item's SVI posterior across draws.
# =============================================================================

import importlib

_net_module  = importlib.import_module(fit['net_class_module'])
_net_class   = getattr(_net_module, fit['net_class_name'])
_net_kwargs  = dict(fit['net_kwargs'])
_net_kwargs.setdefault('num_quantiles', len(fit['pps_ProbH1_lwr_quantiles_mesh']))
_net         = _net_class(**_net_kwargs)

_summary_rows = []
for interim_id in INTERIM_IDS:
    tokens_c, _t, khat_c, aux_c = interim_tensors[interim_id]  # (D, J, 7)
    D_i = tokens_c.shape[0]
    tokens_rep, mask_rep, qidx_rep, _aux_rep = _fanout_with_khat(
        tokens_c, khat_c, aux_c,
    )
    bar_h = _net.apply(
        fit['params'], tokens_rep, mask_rep, qidx_rep,
        method=lambda mdl, t, m, q: mdl.summary_for_query(t, m, q),
    )                                                        # (D*J, E)
    bar_h = np.asarray(bar_h).reshape(D_i, J, -1)            # (D, J, E)
    for d in range(D_i):
        for j in range(J):
            _summary_rows.append({
                'interim_id': interim_id,
                's':          d + 1,
                'item_label': items_df.iloc[j]['item_label'],
                'bar_h':      bar_h[d, j],
            })
_bar_all     = np.stack([r['bar_h'] for r in _summary_rows], axis=0)  # (N, E)
_bar_center  = _bar_all - _bar_all.mean(axis=0, keepdims=True)
_U, _s, _Vt  = np.linalg.svd(_bar_center, full_matrices=False)
_pc1_all     = (_U[:, 0] * _s[0]).astype(np.float64)                  # (N,)
_pc1_df  = pd.DataFrame({
    'interim_id': [r['interim_id'] for r in _summary_rows],
    's':          [r['s']          for r in _summary_rows],
    'item_label': [r['item_label'] for r in _summary_rows],
    'pc1':        _pc1_all,
})
_scatter = dp_h1_xz.merge(
    _pc1_df, on=['interim_id', 's', 'item_label'], how='inner',
)

for interim_id in sorted(_scatter['interim_id'].unique()):
    sub = _scatter[_scatter['interim_id'] == interim_id][
        ['item_label', 'pps_ratio_x', 'pc1']
    ].copy()
    if sub.empty:
        continue
    ann_rows = []
    for item in sub['item_label'].unique():
        g = sub[sub['item_label'] == item]
        rho_v = float(g['pc1'].corr(g['pps_ratio_x'])) if len(g) >= 2 else float('nan')
        ann_rows.append({
            'item_label': item,
            'x_left':  float(g['pc1'].min()),
            'y_top':   float(g['pps_ratio_x'].max()),
            'rho_label': f"rho(PC1,y)={rho_v:.2f}",
        })
    ann = pd.DataFrame(ann_rows)

    _items_sorted = sorted(sub['item_label'].unique())
    _item_colors  = dict(
        zip(_items_sorted, _futurama_palette(len(_items_sorted))),
    )
    p = (
        ggplot(sub, aes(x='pc1', y='pps_ratio_x', colour='item_label'))
        + geom_point(alpha=0.3, size=0.5)
        + geom_smooth(method='glm', method_args={'family': 'gaussian'},
                      se=True, colour='black', size=0.8)
        + geom_text(ann, aes(x='x_left', y='y_top', label='rho_label'),
                    ha='left', va='top', size=8, inherit_aes=False)
        + facet_wrap('~ item_label', ncol=4, scales='free')
        + scale_colour_manual(values=_item_colors)
        + theme_bw()
        + theme(
            figure_size=(15, 14),
            strip_background=element_rect(fill='none', colour='none'),
            panel_background=element_rect(fill='none', colour='none'),
            legend_position='none',
        )
        + labs(x=f'PC1 of attention-picked per-query summary bar_h^(j*) at interim {interim_id}',
               y=f'SVI posterior draw of rho from p(theta | x) at interim {interim_id}')
    )
    pdf_path = os.path.join(
        dir_out, f"{file_prefix}_i{interim_id}_svi_rho_vs_pc1_attn_summary.pdf",
    )
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"  saved PC1 diagnostic plot: {pdf_path}")

print("\nAll interims done -- Ukraine xcompAtt amortised PPS complete.")
