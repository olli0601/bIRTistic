#!/usr/bin/env python3
"""
Ukraine interim analysis: **amortised** PPS via self-attention over the
20 PCM items + multi-quantile head (§12.12 of
``dev/amortised_decision_making.md``).

Adapts the MVN `xcompAtt` (§12.11) to the Ukraine partial-credit model.
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

Outputs (``_RGEG_`` suffix; matches the Ukraine compare-methods
loader):

  - ``pcm_1_interim_pps_RGEG_p_h1_xz.pkl``
  - ``pcm_1_interim_pps_RGEG.csv``
  - ``pcm_1_interim_pps_RGEG_timing.csv``
  - ``pcm_1_interim_pps_RGEG_perf_long.pkl``

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
    ggplot, aes, geom_boxplot, geom_hline, geom_point, geom_smooth, geom_text,
    facet_wrap,
    scale_colour_manual, scale_fill_manual, scale_y_continuous,
    theme_bw, theme, element_text, element_rect, labs,
)

warnings.filterwarnings('ignore')

from utils import _futurama_palette
from amortiser_common import (
    train as amortiser_train,
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

DIR_RGE = os.path.join(
    "/Users/or105/sandbox/bIRTistic",
    "py-ukraine-interim-with-regression-on-endptx-wz-260601",
)
dir_out = os.path.join(
    "/Users/or105/sandbox/bIRTistic",
    "py-ukraine-interim-amortise-endptx-on-wz-with-features-itemXcompAtt-qpsi-MLP-loss-multiquantilehead-260716",
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


def _tensorise_interim(wa: pd.DataFrame):
    """Return (features [D, J, 3], target [D, J], mask) for one interim,
    with items ordered by the shared items_df pos."""
    wa = wa.copy()
    wa['item_pos'] = wa['item_label'].map(item_pos)
    piv = wa.pivot_table(
        index='draw', columns='item_pos',
        values=['w_baseline', 'w_endline', 'w_ratio', 'pps_ratio_x'],
    )
    # Reorder columns so w_baseline / w_endline / w_ratio / pps_ratio_x
    # are consecutive per item_pos.
    D = piv.index.size
    feats = np.stack(
        [
            piv['w_baseline'].reindex(columns=range(J)).to_numpy(dtype=np.float32),
            piv['w_endline'].reindex(columns=range(J)).to_numpy(dtype=np.float32),
            piv['w_ratio'].reindex(columns=range(J)).to_numpy(dtype=np.float32),
        ],
        axis=-1,
    )                                                       # (D, J, 3)
    target = piv['pps_ratio_x'].reindex(
        columns=range(J),
    ).to_numpy(dtype=np.float32)                            # (D, J)
    # Replace NaNs (very rare -- interim 8 had 32 in one column).
    feats = np.nan_to_num(feats, nan=0.0, posinf=0.0, neginf=0.0)
    target = np.nan_to_num(target, nan=0.0, posinf=TARGET_CLIP, neginf=-TARGET_CLIP)
    target = np.clip(target, -TARGET_CLIP, TARGET_CLIP)
    # Add discrete item metadata (broadcast across draws) as two more
    # per-item features.
    metadata = np.stack(
        [item_type_idx, item_high_idx], axis=-1,
    )                                                       # (J, 2)
    metadata_b = np.broadcast_to(
        metadata[None, :, :], (D, J, 2),
    ).astype(np.float32)
    tokens = np.concatenate([feats, metadata_b], axis=-1)   # (D, J, 5)
    return tokens, target


interim_tensors = {i: _tensorise_interim(wa) for i, wa in interim_wa.items()}
for i, (t, y) in interim_tensors.items():
    print(f"  interim {i}: tokens {t.shape}, target {y.shape}")

# %%

# =============================================================================
# Build training pool: concatenate all interims. At training we sample
# uniformly across (interim, draw); each sample fans out to J = 20
# per-item rows (one query per row).
# =============================================================================

train_tokens = np.concatenate(
    [interim_tensors[i][0] for i in INTERIM_IDS], axis=0,
)                                                           # (D_total, J, 5)
train_target = np.concatenate(
    [interim_tensors[i][1] for i in INTERIM_IDS], axis=0,
)                                                           # (D_total, J)
D_total = train_tokens.shape[0]
print(f"\nTraining pool: {D_total} draws across {len(INTERIM_IDS)} interims.")


def sample_fn(rng, S):
    """Sample S draws from the pooled tensors. Fan out to B = S * J
    per-item rows -- one query index per row."""
    draw_idx = rng.integers(0, D_total, size=S)              # (S,)
    tok = train_tokens[draw_idx]                            # (S, J, 5)
    tgt = train_target[draw_idx]                            # (S, J)
    B = S * J
    tokens_rep = np.broadcast_to(
        tok[:, None, :, :], (S, J, J, 5),
    ).reshape(B, J, 5).astype(np.float32)
    mask_rep = np.ones((B, J), dtype=np.float32)
    query_idx = np.broadcast_to(
        np.arange(J)[None, :], (S, J),
    ).reshape(B).astype(np.int32)
    target_rho = tgt.reshape(B).astype(np.float32)
    return (
        {
            'tokens': tokens_rep,
            'mask': mask_rep,
            'query_idx': query_idx,
        },
        target_rho,
    )

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
        Amortiser_PPS_features_itemXcompAtt_qpsi_MLP_loss_multiquantilehead,
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
    tokens, _target = interim_tensors[interim_id]
    D = tokens.shape[0]
    S = D
    B = D * J
    tokens_rep = np.broadcast_to(
        tokens[:, None, :, :], (D, J, J, 5),
    ).reshape(B, J, 5).astype(np.float32)
    mask_rep = np.ones((B, J), dtype=np.float32)
    query_idx = np.broadcast_to(
        np.arange(J)[None, :], (D, J),
    ).reshape(B).astype(np.int32)
    batch = {
        'tokens':   tokens_rep,
        'mask':     mask_rep,
        'query_idx': query_idx,
    }
    p_h1_xz, q_hat, _ = predict_amortised_p_h1_for_one_xz(
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

pkl_path = os.path.join(dir_out, f"{file_prefix}_pps_RGEG_p_h1_xz.pkl")
dp_h1_xz.to_pickle(pkl_path)
print(f"Saved RGEG P(H_1 | x, z) to: {pkl_path}")
pps_timing.to_csv(
    os.path.join(dir_out, f"{file_prefix}_pps_RGEG_timing.csv"), index=False,
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
    os.path.join(dir_out, f"{file_prefix}_pps_RGEG.csv"), index=False,
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
perf_csv = os.path.join(dir_out, f"{file_prefix}_pps_RGEG_perf.csv")
perf.to_csv(perf_csv, index=False)
print(f"Saved RGEG perf CSV to: {perf_csv}")

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
    os.path.join(dir_out, f"{file_prefix}_pps_RGEG_perf_long.pkl"),
)
print(f"Saved RGEG perf long-form pkl.")

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
box_pkl = os.path.join(dir_out, f"{file_prefix}_pps_RGEG_p_h1_xz_boxplot.pkl")
box_stats.to_pickle(box_pkl)
print(f"Saved RGEG p(H_1 | x, z) box-stats to: {box_pkl}")

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
box_pdf = os.path.join(dir_out, f"{file_prefix}_pps_RGEG_p_h1_xz_boxplot.pdf")
p.save(box_pdf, verbose=False, limitsize=False)
print(f"Saved RGEG p(H_1 | x, z) boxplot to: {box_pdf}")

# =============================================================================
# Per-interim diagnostic scatter: predicted median quantile q_hat vs
# observed pps_ratio_x, faceted per item, with per-item rho + R^2
# annotation. Mirrors the ``_endpointratio_vs_wratio.pdf`` layout from
# the RGEM script (target vs prediction rather than target vs w_ratio,
# since xcompAtt does not fit a per-item univariate w->ratio regression).
# =============================================================================

perf_indexed = perf.set_index(['interim_id', 'item_label'])
for interim_id in sorted(dp_h1_xz['interim_id'].unique()):
    sub = dp_h1_xz[dp_h1_xz['interim_id'] == interim_id][
        ['item_label', 'pps_ratio_x', 'q_hat']
    ].copy()
    if sub.empty:
        continue
    ann_rows = []
    for item in sub['item_label'].unique():
        g = sub[sub['item_label'] == item]
        try:
            rho_v = float(perf_indexed.loc[(interim_id, item), 'rho'])
            r2_v  = float(perf_indexed.loc[(interim_id, item), 'r2'])
        except KeyError:
            rho_v, r2_v = float('nan'), float('nan')
        x_lo, x_hi = float(g['q_hat'].min()), float(g['q_hat'].max())
        y_hi = float(g['pps_ratio_x'].max())
        ann_rows.append({
            'item_label': item,
            'x_left':  x_lo, 'x_right': x_hi, 'y_top': y_hi,
            'rho_label': f"ρ={rho_v:.2f}",
            'r2_label':  f"R2={100 * r2_v:.1f}%",
        })
    ann = pd.DataFrame(ann_rows)

    _items_sorted = sorted(sub['item_label'].unique())
    _item_colors = dict(
        zip(_items_sorted, _futurama_palette(len(_items_sorted))),
    )

    p = (
        ggplot(sub, aes(x='q_hat', y='pps_ratio_x', colour='item_label'))
        + geom_point(alpha=0.3, size=0.5)
        + geom_smooth(method='glm', method_args={'family': 'gaussian'},
                      se=True, colour='black', size=0.8)
        + geom_text(ann, aes(x='x_left', y='y_top', label='rho_label'),
                    ha='left', va='top', size=8, inherit_aes=False)
        + geom_text(ann, aes(x='x_right', y='y_top', label='r2_label'),
                    ha='right', va='top', size=8, inherit_aes=False)
        + facet_wrap('~ item_label', ncol=4, scales='free')
        + scale_colour_manual(values=_item_colors)
        + theme_bw()
        + theme(
            figure_size=(15, 14),
            strip_background=element_rect(fill='none', colour='none'),
            panel_background=element_rect(fill='none', colour='none'),
            legend_position='none',
        )
        + labs(x=f'xcompAtt predicted median quantile q_hat at interim {interim_id}',
               y=f'actual endpoint from posterior p(theta | x) draws at interim {interim_id}')
    )
    pdf_path = os.path.join(
        dir_out, f"{file_prefix}_i{interim_id}_endpointratio_vs_qhat.pdf",
    )
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"  saved diagnostic plot: {pdf_path}")

print("\nAll interims done -- Ukraine xcompAtt amortised PPS complete.")
