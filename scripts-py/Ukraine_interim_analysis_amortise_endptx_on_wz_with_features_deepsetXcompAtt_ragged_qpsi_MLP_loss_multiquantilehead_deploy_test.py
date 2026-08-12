#!/usr/bin/env python3
"""
Deployment test of the RAGGED focused-proposal amortiser (§14.4.11) on
the weekly SVI grid: (a) posterior-contraction slope (predicted marginal
width vs n, against the SVI truth), and (b) interim-1-only head fine-tune
coverage at all interims. Success = W_marg slope -> ~ -0.59 and i1-only
head-ft near-nominal everywhere.
"""
import os
import sys
import glob
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
from amortiser_pps_features_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead import (
    Amortiser_PPS_features_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead as Net, _MLP)
from model_pcm import PartialCreditModel

SB = "/Users/or105/sandbox/bIRTistic"
WK = f"{SB}/py-ukraine-interim-weekly-svi-260811"
BASE = f"{SB}/py-ukraine-interim-amortise-endptx-on-wz-with-features-deepsetXcompAtt-ragged-qpsi-MLP-loss-multiquantilehead-260812"
file_prefix = "pcm_1_interim"
ANCHOR = 4                 # weekly interim used for i1-only head-ft
N_REF, CLIP = 503, 20.0
S = int(os.environ.get('RAG_S', 200))
TAUS = np.array([0.05, 0.25, 0.5, 0.75, 0.95], np.float32)
CACHE = "/private/tmp/claude-501/-Users-or105-git-bIRTistic/49372e22-d11b-4882-a442-d1c60bcbdfb0/scratchpad/ragged_deploy_cells"
os.makedirs(CACHE, exist_ok=True)

dit = pd.read_csv(f"{WK}/{file_prefix}_1_data_dit.csv")
_wa = sorted(glob.glob(f"{WK}/{file_prefix}_i*_regression_training.pkl"),
             key=lambda p: int(p.split('_i')[-1].split('_')[0]))
INTERIMS = [int(p.split('_i')[-1].split('_')[0]) for p in _wa]
items = (pd.read_pickle(_wa[0])[['item_label', 'item_type', 'item_high_label']].drop_duplicates()
         .sort_values(['item_type', 'item_label']).reset_index(drop=True))
J = len(items); labels = items.item_label.tolist()
itype = items.item_type.map({'out-of-7': 0., 'categorical': 1.}).to_numpy(np.float32)
ihigh = items.item_high_label.map({'higher_is_better': 1., 'lower_is_better': 0.}).to_numpy(np.float32)
items = items.merge(dit[['item_label', 'cat_length']].drop_duplicates(), on='item_label', how='left')
kmax = items.cat_length.to_numpy(np.float32)
META = np.stack([itype, ihigh], -1).astype(np.float32)
fit = load_fitted_model(f"{BASE}/{file_prefix}_amortised_pps_net.pkl")
sig = np.load(f"{BASE}/{file_prefix}_item_std.npy").astype(np.float32)
net = Net(num_quantiles=5)


def _pivot_pids(df, col, pids):
    piv = df.pivot_table(index='pid', columns=['item_label', 'time'], values=col).reindex(pids)
    out = np.zeros((len(pids), J, 2), np.float32)
    for j, l in enumerate(labels):
        for t in (0, 1):
            if (l, t) in piv.columns:
                out[:, j, t] = piv[(l, t)].to_numpy(np.float32)
    return np.nan_to_num((out - 1.) / (kmax[None, :, None] - 1.))


@jax.jit
def _fwd(params, batch):
    out, st = net.apply(params, batch, mutable=['intermediates'])
    return out, st['intermediates']['head_in'][0]


def build(k):
    xi = pd.read_csv(f"{WK}/{file_prefix}_{k}_data_dp1.csv")
    xpids = np.sort(xi.pid.unique()); n = len(xpids)
    x_raw = _pivot_pids(xi, 'y_stan', xpids)                # (n, J, 2)
    m = N_REF - n
    model = PartialCreditModel(dit=dit, dcati=xi, x_formula="~ time - 1", seed=123)
    zi = model.get_interim_z_from_ypredi(f"{WK}/{file_prefix}_{k}_draws.zarr",
                                         m, pps_z_total=S, seed=123, keep_order=True)
    zpids = np.sort(zi.pid.unique()); mm = len(zpids)
    aux = np.array([[1/np.sqrt(n), 1/np.sqrt(mm), n/N_REF, mm/N_REF]], np.float32)
    xb = jnp.asarray(x_raw); xseg = jnp.zeros(n, jnp.int32)
    metab = jnp.asarray(META[None]); qidb = jnp.arange(J)[None]; auxb = jnp.asarray(aux)
    qs = np.empty((S, J, 5), np.float64); hd = None
    for s in range(S):
        z_raw = _pivot_pids(zi.rename(columns={f'ypred_{s}': '_y'}), '_y', zpids)
        batch = dict(x_flat=xb, x_seg=xseg, z_flat=jnp.asarray(z_raw),
                     z_seg=jnp.zeros(mm, jnp.int32), item_metadata=metab,
                     query_idx=qidb, aux=auxb)
        out, head = _fwd(fit['params'], batch)
        qs[s] = np.maximum.accumulate(np.asarray(out)[0], 1) * sig[:, None]
        if hd is None:
            hd = np.empty((S, J, head.shape[-1]), np.float32)
        hd[s] = np.asarray(head)[0]
    wa = pd.read_pickle(f"{WK}/{file_prefix}_i{k}_regression_training.pkl")
    piv = wa.pivot_table(index='draw', columns='item_label', values='pps_ratio_x').reindex(columns=labels)
    tgt = np.full((S, J), np.nan)
    idx = piv.index[piv.index < S].to_numpy()
    tgt[idx] = np.clip(piv.loc[idx].to_numpy(np.float64), -CLIP, CLIP)
    return qs, hd.astype(np.float32), tgt, n


print(f"Phase 0: forward ragged net on {len(INTERIMS)} weekly interims (S={S}) ...")
QS, HD, TGT, NOBS = {}, {}, {}, {}
for k in INTERIMS:
    cf = f"{CACHE}/cell_{k}.npz"
    if os.path.exists(cf):
        z = np.load(cf); QS[k], HD[k], TGT[k], NOBS[k] = z['qs'], z['hd'], z['tgt'], int(z['n'])
    else:
        QS[k], HD[k], TGT[k], NOBS[k] = build(k)
        np.savez(cf, qs=QS[k], hd=HD[k], tgt=TGT[k], n=NOBS[k])
    print(f"  interim {k} (n={NOBS[k]}) done")


def marg_width(qs):
    lo = float(qs[:, 0].min()); hi = float(qs[:, 4].max())
    if hi <= lo:
        return np.nan
    gr = np.linspace(lo, hi, 400); F = np.zeros_like(gr)
    for s in range(qs.shape[0]):
        F += np.interp(gr, qs[s], TAUS, 0., 1.)
    F /= qs.shape[0]
    return float(np.interp(0.95, F, gr) - np.interp(0.05, F, gr))


def pit(qs, y):
    return np.array([np.interp(y[s], qs[s], TAUS, 0., 1.) for s in range(len(y))])


# ---- (a) contraction slopes on the weekly grid ----
BLOW = "CG-VIO_ph-punish"
rows = []
for k in INTERIMS:
    for j in range(J):
        y = TGT[k][:, j]; ok = np.isfinite(y); yv = y[ok]
        if yv.size < 10 or labels[j] == BLOW:
            continue
        rows.append(dict(interim_id=k, n=NOBS[k], item_label=labels[j],
                         W_svi=float(np.quantile(yv, .95) - np.quantile(yv, .05)),
                         W_marg=marg_width(QS[k][ok, j, :])))
cw = pd.DataFrame(rows)
cw.to_csv(f"{BASE}/{file_prefix}_pps_RAGD_contraction.csv", index=False)
def slope(col):
    d = cw[(cw[col] > 0) & (cw.n > 0)]
    return float(np.polyfit(np.log(d.n), np.log(d[col]), 1)[0])
print(f"\n(a) CONTRACTION slopes (weekly grid): W_svi={slope('W_svi'):.3f}  W_marg={slope('W_marg'):.3f}")
lo, hi = cw.n.min(), cw.n.max()
for c in ['W_svi', 'W_marg']:
    print(f"    {c}: n={lo}:{cw[cw.n==lo][c].median():.3f} -> n={hi}:{cw[cw.n==hi][c].median():.3f}")

# ---- (b) interim-1-only head fine-tune ----
taus_j = jnp.asarray(TAUS)
qpsi = _MLP(dims=(64, 64, 5)); qpsi0 = {'params': fit['params']['params']['q_psi']}


def pinball(params, hd, y):
    pr = jnp.maximum.accumulate(qpsi.apply(params, hd), 1)
    e = y[:, None] - pr
    return jnp.mean(jnp.maximum(taus_j[None]*e, (taus_j[None]-1)*e))


HD1 = HD[ANCHOR].reshape(-1, HD[ANCHOR].shape[-1])
TG1 = (TGT[ANCHOR] / sig[None]).reshape(-1)
okm = np.isfinite(TG1); HD1, TG1 = HD1[okm].astype(np.float32), TG1[okm].astype(np.float32)
params = qpsi0; opt = optax.adam(3e-4); ost = opt.init(params); rng = np.random.default_rng(1)


@jax.jit
def hstep(params, ost, h, y):
    l, g = jax.value_and_grad(pinball)(params, h, y)
    up, ost = opt.update(g, ost, params)
    return optax.apply_updates(params, up), ost, l


for it in range(800):
    ix = rng.integers(0, HD1.shape[0], min(2048, HD1.shape[0]))
    params, ost, l = hstep(params, ost, jnp.asarray(HD1[ix]), jnp.asarray(TG1[ix]))
print(f"\n(b) i1-only head-ft fit done (loss {float(l):.4f})")


def cov_pit(getq):
    out = []
    for k in INTERIMS:
        us = []
        for j in range(J):
            y = TGT[k][:, j]; ok = np.isfinite(y)
            us.append(pit(getq(k)[ok, j, :], y[ok]))
        u = np.concatenate(us); uu = np.sort(u); ec = np.arange(1, len(u)+1)/len(u)
        out.append(dict(interim_id=k, n=NOBS[k], cov5=float((u <= .5).mean()),
                        cov95=float((u <= .95).mean()), pit_ks=float(np.max(np.abs(ec-uu)))))
    return pd.DataFrame(out)


HK = {}
for k in INTERIMS:
    pr = np.maximum.accumulate(np.asarray(qpsi.apply(params, jnp.asarray(HD[k].reshape(-1, HD[k].shape[-1])))), 1)
    HK[k] = (pr.reshape(HD[k].shape[0], J, 5)) * sig[None, :, None]
base = cov_pit(lambda k: QS[k]).assign(config='baseline(ragged)')
h1 = cov_pit(lambda k: HK[k]).assign(config='head-ft-i1(ragged)')
res = pd.concat([base, h1], ignore_index=True)
res.to_csv(f"{BASE}/{file_prefix}_pps_RAGD_headft_i1_coverage.csv", index=False)
print("\n=== MEAN over weekly interims (nominal cov5=.50 cov95=.95) ===")
print(res.groupby('config')[['cov5', 'cov95', 'pit_ks']].mean().round(3).to_string())
print("\nhead-ft-i1 per interim:")
print(h1[['interim_id', 'n', 'cov5', 'cov95', 'pit_ks']].to_string(index=False))
print(f"\n-> {BASE}")
