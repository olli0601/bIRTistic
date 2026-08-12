#!/usr/bin/env python3
"""
Does contraction pressure make an INTERIM-1-ONLY head fine-tune enough
(§14.4.9 next-step test)? Fit q_psi on interim-1 SVI only, apply to ALL
8 monthly interims, measure coverage/PIT. Compare:

  baseline(new net)   contraction-retrained net, no recalibration
  head-ft-i1          q_psi refit on interim-1 embeddings only

Run on the contraction-retrained net (aux-dim-4, AUX_SQRT=1). Success =
near-nominal coverage at ALL interims from the single interim-1 fit
(previously, on the base net, i1-only fixed i1 but drifted elsewhere).
"""
import os
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'python'))
import numpy as np
import pandas as pd
import jax
import jax.numpy as jnp
import optax

warnings.filterwarnings('ignore')
from amortiser_common import load_fitted_model
from amortiser_pps_features_deepsetXcompAtt_qpsi_MLP_loss_multiquantilehead import (
    Amortiser_PPS_features_deepsetXcompAtt_qpsi_MLP_loss_multiquantilehead as DSNet, _MLP)
from model_pcm import PartialCreditModel

SB = "/Users/or105/sandbox/bIRTistic"
DIR_RGE = f"{SB}/py-ukraine-interim-with-regression-on-endptx-wz-260601"
BASE = os.environ.get('RGDX_BASE',
    f"{SB}/py-ukraine-interim-amortise-endptx-on-wz-with-features-deepsetXcompAtt-qpsi-MLP-loss-multiquantilehead-contraction-260811")
AUX_SQRT = os.environ.get('AUX_SQRT', '1') == '1'
CTAG = os.environ.get('RGDX_CTAG', 'contract')
file_prefix = "pcm_1_interim"
N_FULL, CLIP = 503, 20.0
S = int(os.environ.get('UKR_DEEPSET_S', 200))
TAUS = np.array([0.05, 0.25, 0.5, 0.75, 0.95], np.float32)
INTERIMS = list(range(1, 9))
CACHE = f"/private/tmp/claude-501/-Users-or105-git-bIRTistic/49372e22-d11b-4882-a442-d1c60bcbdfb0/scratchpad/headft_i1_cells_{CTAG}"
os.makedirs(CACHE, exist_ok=True)

dit = pd.read_csv(f"{DIR_RGE}/{file_prefix}_1_data_dit.csv")
_wa1 = pd.read_pickle(f"{DIR_RGE}/{file_prefix}_i1_regression_training.pkl")
items = (_wa1[['item_label', 'item_type', 'item_high_label']].drop_duplicates()
         .sort_values(['item_type', 'item_label']).reset_index(drop=True))
J = len(items); labels = items.item_label.tolist()
itype = items.item_type.map({'out-of-7': 0., 'categorical': 1.}).to_numpy(np.float32)
ihigh = items.item_high_label.map({'higher_is_better': 1., 'lower_is_better': 0.}).to_numpy(np.float32)
items = items.merge(dit[['item_label', 'cat_length']].drop_duplicates(), on='item_label', how='left')
kmax = items.cat_length.to_numpy(np.float32)
fit = load_fitted_model(f"{BASE}/{file_prefix}_amortised_pps_net.pkl")
nkw = dict(fit['net_kwargs']); nkw.setdefault('num_quantiles', 5); net = DSNet(**nkw)


def _pivot(df, col):
    piv = df.pivot_table(index='pid', columns=['item_label', 'time'], values=col)
    out = np.zeros((piv.index.size, J, 2), np.float32)
    for j, l in enumerate(labels):
        for t in (0, 1):
            out[:, j, t] = piv[(l, t)].to_numpy(np.float32)
    return np.nan_to_num((out - 1.) / (kmax[None, :, None] - 1.)), piv.index.to_numpy()


@jax.jit
def _fwd(params, batch):
    out, state = net.apply(params, batch, mutable=['intermediates'])
    return out, state['intermediates']['head_in'][0]


def build(k):
    xi = pd.read_csv(f"{DIR_RGE}/{file_prefix}_{k}_data_dp1.csv")
    x_raw, _ = _pivot(xi, 'y_stan'); n = x_raw.shape[0]; m = N_FULL - n
    model = PartialCreditModel(dit=dit, dcati=xi, x_formula="~ time - 1", seed=123)
    zi = model.get_interim_z_from_ypredi(f"{DIR_RGE}/{file_prefix}_{k}_draws.zarr",
                                         m, pps_z_total=S, seed=123, keep_order=True)
    x_pad = np.zeros((N_FULL, J, 2), np.float32); x_pad[:n] = x_raw
    mx = np.zeros(N_FULL, np.float32); mx[:n] = 1.
    mz = np.zeros(N_FULL, np.float32); mz[:m] = 1.
    xb = np.broadcast_to(x_pad[None], (J, N_FULL, J, 2)).astype(np.float32)
    mxb = np.broadcast_to(mx[None], (J, N_FULL)).astype(np.float32)
    mzb = np.broadcast_to(mz[None], (J, N_FULL)).astype(np.float32)
    metab = np.broadcast_to(np.stack([itype, ihigh], -1)[None], (J, J, 2)).astype(np.float32)
    av = [n/N_FULL, m/N_FULL]
    if AUX_SQRT:
        av += [1.0/np.sqrt(max(n, 1)), 1.0/np.sqrt(max(m, 1))]
    auxb = np.broadcast_to(np.array(av, np.float32)[None], (J, len(av))).astype(np.float32)
    qidb = np.arange(J, dtype=np.int32)
    qs = np.empty((S, J, 5), np.float64); hd = None
    for s in range(S):
        z_raw, _ = _pivot(zi.rename(columns={f'ypred_{s}': '_y'}), '_y')
        z_pad = np.zeros((N_FULL, J, 2), np.float32); z_pad[:m] = z_raw
        zb = np.broadcast_to(z_pad[None], (J, N_FULL, J, 2)).astype(np.float32)
        batch = {'x_responses': jnp.asarray(xb), 'mask_x': jnp.asarray(mxb),
                 'z_responses': jnp.asarray(zb), 'mask_z': jnp.asarray(mzb),
                 'item_metadata': jnp.asarray(metab), 'query_idx': jnp.asarray(qidb),
                 'aux': jnp.asarray(auxb)}
        out, head = _fwd(fit['params'], batch)
        qs[s] = np.maximum.accumulate(np.asarray(out), 1)
        if hd is None:
            hd = np.empty((S, J, head.shape[-1]), np.float32)
        hd[s] = np.asarray(head)
    wa = pd.read_pickle(f"{DIR_RGE}/{file_prefix}_i{k}_regression_training.pkl")
    piv = wa.pivot_table(index='draw', columns='item_label', values='pps_ratio_x').reindex(columns=labels)
    tgt = np.full((S, J), np.nan)
    idx = piv.index[piv.index < S].to_numpy()
    tgt[idx] = np.clip(piv.loc[idx].to_numpy(np.float64), -CLIP, CLIP)
    return qs, hd.astype(np.float32), tgt


print(f"Phase 0: forward monthly interims (AUX_SQRT={AUX_SQRT}) ...")
QS, HD, TGT = {}, {}, {}
for k in INTERIMS:
    cf = f"{CACHE}/cell_{k}.npz"
    if os.path.exists(cf):
        z = np.load(cf); QS[k], HD[k], TGT[k] = z['qs'], z['hd'], z['tgt']
    else:
        QS[k], HD[k], TGT[k] = build(k); np.savez(cf, qs=QS[k], hd=HD[k], tgt=TGT[k])
    print(f"  interim {k}: qs{QS[k].shape} head{HD[k].shape}")


def pit(qs, y):
    return np.array([np.interp(y[s], qs[s], TAUS, 0., 1.) for s in range(len(y))])


def cov_pit(get_qs):
    rows = []
    for k in INTERIMS:
        us = []
        for j in range(J):
            y = TGT[k][:, j]; ok = np.isfinite(y)
            u = pit(get_qs(k)[ok, j, :], y[ok]); us.append(u)
        u = np.concatenate(us)
        uu = np.sort(u); ec = np.arange(1, len(u)+1)/len(u)
        rows.append(dict(interim_id=k, cov5=float((u <= .5).mean()),
                         cov95=float((u <= .95).mean()),
                         pit_ks=float(np.max(np.abs(ec-uu)))))
    return pd.DataFrame(rows)


# ---- baseline (new net) ----
base = cov_pit(lambda k: QS[k]).assign(config='baseline(new net)')

# ---- head-ft on interim-1 ONLY, applied to all ----
taus_j = jnp.asarray(TAUS)
qpsi = _MLP(dims=(*nkw.get('hidden_dims', (64, 64)), 5))
qpsi0 = {'params': fit['params']['params']['q_psi']}


def pinball(params, hd, y):
    pr = jnp.maximum.accumulate(qpsi.apply(params, hd), 1)
    e = y[:, None] - pr
    return jnp.mean(jnp.maximum(taus_j[None]*e, (taus_j[None]-1)*e))


HD1 = HD[1].reshape(-1, HD[1].shape[-1]); TG1 = TGT[1].reshape(-1)
ok = np.isfinite(TG1); HD1, TG1 = HD1[ok].astype(np.float32), TG1[ok].astype(np.float32)
params = qpsi0; opt = optax.adam(3e-4); ost = opt.init(params); rng = np.random.default_rng(1)


@jax.jit
def step(params, ost, h, y):
    l, g = jax.value_and_grad(pinball)(params, h, y)
    up, ost = opt.update(g, ost, params)
    return optax.apply_updates(params, up), ost, l


bs = min(2048, HD1.shape[0])
for it in range(800):
    ix = rng.integers(0, HD1.shape[0], bs)
    params, ost, l = step(params, ost, jnp.asarray(HD1[ix]), jnp.asarray(TG1[ix]))
print(f"head-ft-i1 fit done (loss {float(l):.4f})")

HK = {}
for k in INTERIMS:
    pr = np.maximum.accumulate(np.asarray(qpsi.apply(params, jnp.asarray(HD[k].reshape(-1, HD[k].shape[-1])))), 1)
    HK[k] = pr.reshape(S, J, 5)
h1 = cov_pit(lambda k: HK[k]).assign(config='head-ft-i1(new net)')

out = pd.concat([base, h1], ignore_index=True)
out.to_csv(f"{BASE}/{file_prefix}_pps_RGDX_headft_i1_coverage.csv", index=False)
print("\n=== per-interim (nominal cov5=0.50 cov95=0.95; low pit_ks good) ===")
for cfg in ['baseline(new net)', 'head-ft-i1(new net)']:
    s = out[out.config == cfg]
    print(f"\n{cfg}:")
    print(s[['interim_id', 'cov5', 'cov95', 'pit_ks']].to_string(index=False))
print("\n=== MEAN over interims ===")
print(out.groupby('config')[['cov5', 'cov95', 'pit_ks']].mean().round(3).to_string())
print(f"\n-> {BASE}")
