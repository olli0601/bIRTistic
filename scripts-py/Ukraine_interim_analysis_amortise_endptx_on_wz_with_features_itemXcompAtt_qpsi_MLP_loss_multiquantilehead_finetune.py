#!/usr/bin/env python3
"""
Deployment-time fine-tune of the prior-predictive-trained `itemXcompAtt`
amortiser on cached SVI interim data (remedy workstreams B2-alt and C of
the §14.1.12 calibration plan).

Loads the prior-trained net (encoders + head) and fine-tunes on the
(z^(s), rho_SVI-x) pairs of a chosen set of interims, then writes the
fine-tuned net + item_std into FT_DIR. The standard deploy + calibration
tests are then run against FT_DIR (via UKR_DIR) to score calibration.

Env vars:
  FT_BASE_DIR   prior-trained dir to load net + item_std from
  FT_DIR        output dir for the fine-tuned net (+ item_std copy)
  FT_FREEZE     'encoders' -> freeze q_tok, q_query, train only q_psi
                'none'      -> train all params
  FT_FIT        comma-sep interims to FIT on (e.g. '1' for B2-alt;
                '2,3,4,5,6,7,8' for a C leave-one-out fold that holds
                out interim 1)
  FT_STEPS      fine-tune steps (default 800)
  FT_LR         learning rate (default 3e-4)
  FT_BATCH      rows per step (default 2048)

Usage (B2-alt, head-only, fit interim 1):
    FT_BASE_DIR=<prior_dir> FT_DIR=<out_dir> FT_FREEZE=encoders FT_FIT=1 \
    pixi run python -u scripts-py/..._finetune.py
"""

import os
import sys
import shutil
from functools import partial
from pathlib import Path

try:
    project_root = Path(__file__).resolve().parent.parent
except NameError:
    project_root = Path.cwd()
sys.path.insert(0, str(project_root / 'python'))

import numpy as np
import pandas as pd
import jax
import jax.numpy as jnp
import optax
import importlib

from amortiser_common import load_fitted_model, save_trained_model
from model_pcm import PartialCreditModel

# ---- config ----
DIR_RGE = os.path.join(
    "/Users/or105/sandbox/bIRTistic",
    "py-ukraine-interim-with-regression-on-endptx-wz-260601",
)
FT_BASE_DIR = os.environ['FT_BASE_DIR']
FT_DIR = os.environ['FT_DIR']
FT_FREEZE = os.environ.get('FT_FREEZE', 'encoders')
FT_FIT = [int(x) for x in os.environ.get('FT_FIT', '1').split(',')]
FT_STEPS = int(os.environ.get('FT_STEPS', 800))
FT_LR = float(os.environ.get('FT_LR', 3e-4))
FT_BATCH = int(os.environ.get('FT_BATCH', 2048))
os.makedirs(FT_DIR, exist_ok=True)
file_prefix = "pcm_1_interim"
N_FULL = 503
KHAT_SHRINK_N0 = 50.0
F_TOK = 8
TARGET_CLIP = 20.0
NET_TAUS = np.array([0.05, 0.25, 0.5, 0.75, 0.95], dtype=np.float32)
print(f"Fine-tune: freeze={FT_FREEZE} fit_interims={FT_FIT} "
      f"steps={FT_STEPS} lr={FT_LR}\n base={FT_BASE_DIR}\n out={FT_DIR}")

# ---- item metadata (mirror the main script) ----
interim_wa = {i: pd.read_pickle(
    os.path.join(DIR_RGE, f"{file_prefix}_i{i}_regression_training.pkl"))
    for i in range(1, 9)
    if os.path.exists(os.path.join(
        DIR_RGE, f"{file_prefix}_i{i}_regression_training.pkl"))}
items_df = (interim_wa[1][['item_label', 'item_type', 'item_high_label']]
            .drop_duplicates()
            .sort_values(['item_type', 'item_label']).reset_index(drop=True))
items_df['item_pos'] = np.arange(len(items_df))
J = len(items_df)
item_pos = dict(zip(items_df['item_label'], items_df['item_pos']))
item_type_idx = items_df['item_type'].map(
    {'out-of-7': 0.0, 'categorical': 1.0}).to_numpy(np.float32)
item_high_idx = items_df['item_high_label'].map(
    {'higher_is_better': 1.0, 'lower_is_better': 0.0}).to_numpy(np.float32)
dit_df = pd.read_csv(os.path.join(DIR_RGE, f"{file_prefix}_1_data_dit.csv"))
items_df = items_df.merge(
    dit_df[['item_label', 'cat_length']].drop_duplicates(),
    on='item_label', how='left')
item_klevels = items_df['cat_length'].to_numpy(np.float32)

# ---- net + item_std from prior-trained dir ----
fit = load_fitted_model(os.path.join(
    FT_BASE_DIR, f"{file_prefix}_amortised_pps_net.pkl"))
item_std = np.load(os.path.join(
    FT_BASE_DIR, f"{file_prefix}_item_std.npy")).astype(np.float32)
net_mod = importlib.import_module(fit['net_class_module'])
Net = getattr(net_mod, fit['net_class_name'])
net_kwargs = dict(fit['net_kwargs'])
net_kwargs.setdefault('num_quantiles', len(fit['pps_ProbH1_lwr_quantiles_mesh']))
net = Net(**net_kwargs)


# ---- token construction (posterior-predictive z^(s); mirrors main) ----
def _x_feats(interim_id):
    xi = pd.read_csv(os.path.join(
        DIR_RGE, f"{file_prefix}_{interim_id}_data_dp1.csv"))
    piv = xi.pivot_table(index='pid', columns=['item_label', 'time'],
                         values='y')
    n_pid = piv.index.size
    w_x = np.zeros((J, 2), np.float32)
    chg = np.zeros((n_pid, J))
    for j, lbl in enumerate(items_df['item_label']):
        y0 = piv[(lbl, 0)].to_numpy(float); y1 = piv[(lbl, 1)].to_numpy(float)
        km = item_klevels[j] - 1.0
        w_x[j, 0] = np.nanmean(y0) / km; w_x[j, 1] = np.nanmean(y1) / km
        chg[:, j] = y1 - y0
    r = np.nan_to_num(pd.DataFrame(chg).corr(method='spearman').to_numpy())
    np.fill_diagonal(r, 1.0)
    lam = n_pid / (n_pid + KHAT_SHRINK_N0)
    return w_x, (lam * r + (1 - lam) * np.eye(J)).astype(np.float32), n_pid


def _z_means(interim_id, D):
    xi = pd.read_csv(os.path.join(
        DIR_RGE, f"{file_prefix}_{interim_id}_data_dp1.csv"))
    n_pid = xi['pid'].nunique(); m = N_FULL - n_pid
    model = PartialCreditModel(dit=dit_df, dcati=xi, x_formula="~ time - 1",
                               seed=123)
    zi = model.get_interim_z_from_ypredi(
        os.path.join(DIR_RGE, f"{file_prefix}_{interim_id}_draws.zarr"),
        m, pps_z_total=D, seed=123, keep_order=True)
    cols = [f'ypred_{s}' for s in range(D)]
    out = np.zeros((D, J, 2), np.float32)
    for j, lbl in enumerate(items_df['item_label']):
        km = item_klevels[j] - 1.0
        for tt in (0, 1):
            v = zi[(zi['item_label'] == lbl) & (zi['time'] == tt)][cols].to_numpy(np.float32)
            out[:, j, tt] = (v - 1.0).mean(0) / km
    return out


def _rows(interim_id):
    """Return (tokens_fanned (D*J,J,8), aux (D*J,2), qidx (D*J,),
    target_std (D*J,)) for one interim."""
    w_x, khat, n_pid = _x_feats(interim_id)
    wa = interim_wa[interim_id].copy()
    wa['item_pos'] = wa['item_label'].map(item_pos)
    piv = wa.pivot_table(index='draw', columns='item_pos', values='pps_ratio_x')
    D = piv.index.size
    tgt = np.clip(np.nan_to_num(
        piv.reindex(columns=range(J)).to_numpy(np.float32),
        posinf=TARGET_CLIP, neginf=-TARGET_CLIP), -TARGET_CLIP, TARGET_CLIP)
    zm = _z_means(interim_id, D)
    wb = np.maximum(zm[:, :, 0], 1e-3); we = zm[:, :, 1]
    hi = (item_high_idx > 0.5)[None, :]
    w_rat = np.clip(np.where(hi, we / wb - 1, 1 - we / wb),
                    -TARGET_CLIP, TARGET_CLIP)
    feats = np.stack([np.broadcast_to(w_x[None, :, 0], (D, J)),
                      np.broadcast_to(w_x[None, :, 1], (D, J)),
                      zm[:, :, 0], zm[:, :, 1], w_rat], -1)
    feats = np.nan_to_num(feats)
    meta = np.broadcast_to(
        np.stack([item_type_idx, item_high_idx], -1)[None], (D, J, 2))
    tok = np.concatenate([feats, meta], -1).astype(np.float32)   # (D,J,7)
    tok_f = np.empty((D, J, J, F_TOK), np.float32)
    tok_f[..., :7] = tok[:, None]; tok_f[..., 7] = khat[None]
    tok_f = tok_f.reshape(D * J, J, F_TOK)
    qidx = np.broadcast_to(np.arange(J)[None], (D, J)).reshape(D * J).astype(np.int32)
    aux = np.broadcast_to(
        np.array([n_pid / N_FULL, (N_FULL - n_pid) / N_FULL], np.float32)[None],
        (D * J, 2))
    tgt_std = (tgt / item_std[None, :]).reshape(D * J).astype(np.float32)
    return tok_f, aux.astype(np.float32), qidx, tgt_std


print("Building fine-tune rows ...")
TOK = []; AUX = []; QID = []; TGT = []
for i in FT_FIT:
    tk, ax, qi, tg = _rows(i)
    TOK.append(tk); AUX.append(ax); QID.append(qi); TGT.append(tg)
    print(f"  interim {i}: {tk.shape[0]} rows")
TOK = np.concatenate(TOK); AUX = np.concatenate(AUX)
QID = np.concatenate(QID); TGT = np.concatenate(TGT)
Nrow = TOK.shape[0]

# ---- optax with optional encoder freeze ----
taus = jnp.asarray(NET_TAUS)


def pinball(preds, y):
    preds = jnp.maximum.accumulate(preds, axis=1)          # monotone
    e = y[:, None] - preds
    return jnp.mean(jnp.maximum(taus[None, :] * e, (taus[None, :] - 1.0) * e))


def loss_fn(params, batch):
    preds = net.apply(params, batch)
    return pinball(preds, batch['target'])


def _label(params):
    # label top-level submodules: 'frozen' for encoders, 'train' for head
    def lab(path, _):
        name = path[1].key if len(path) > 1 else ''
        if FT_FREEZE == 'encoders' and name in ('q_tok', 'q_query'):
            return 'frozen'
        return 'train'
    return jax.tree_util.tree_map_with_path(lab, params)


tx = optax.multi_transform(
    {'train': optax.adam(FT_LR), 'frozen': optax.set_to_zero()},
    _label(fit['params']),
)
opt_state = tx.init(fit['params'])
params = fit['params']


@jax.jit
def step(params, opt_state, batch):
    l, g = jax.value_and_grad(loss_fn)(params, batch)
    updates, opt_state = tx.update(g, opt_state, params)
    params = optax.apply_updates(params, updates)
    return params, opt_state, l


rng = np.random.default_rng(0)
print(f"Fine-tuning ({FT_STEPS} steps, {Nrow} rows, batch {FT_BATCH}) ...")
for s in range(FT_STEPS):
    idx = rng.integers(0, Nrow, size=min(FT_BATCH, Nrow))
    batch = {'tokens': jnp.asarray(TOK[idx]), 'mask': jnp.ones((len(idx), J)),
             'query_idx': jnp.asarray(QID[idx]), 'aux': jnp.asarray(AUX[idx]),
             'target': jnp.asarray(TGT[idx])}
    params, opt_state, l = step(params, opt_state, batch)
    if s % 100 == 0 or s == FT_STEPS - 1:
        print(f"  step {s}/{FT_STEPS}  loss={float(l):.5f}")

# ---- save fine-tuned net + item_std copy ----
fit_ft = dict(fit)
fit_ft['params'] = jax.device_get(params)
save_trained_model(fit_ft, os.path.join(
    FT_DIR, f"{file_prefix}_amortised_pps_net.pkl"))
np.save(os.path.join(FT_DIR, f"{file_prefix}_item_std.npy"), item_std)
print(f"Saved fine-tuned net + item_std to {FT_DIR}")
print("Fine-tune complete.")
