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
from plotnine import (ggplot, aes, geom_abline, geom_line, geom_point, geom_smooth,
                      geom_histogram, geom_hline, geom_text, geom_boxplot, position_dodge,
                      facet_wrap, facet_grid,
                      theme_bw, theme, element_blank, element_text, labs,
                      scale_color_cmap, scale_linetype_manual, scale_color_manual,
                      scale_fill_manual, scale_fill_gradient2)
from amortiser_common import load_fitted_model
from amortiser_pps_features_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead import (
    Amortiser_PPS_features_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead as Net, _MLP)
from model_pcm import PartialCreditModel

SB = "/Users/or105/sandbox/bIRTistic"
# Source of interim artifacts: weekly grid by default; point RAGD_RGE at the
# monthly RGE dir for the 8-interim diagnostics that match the Xcomp/Scomp PDFs.
WK = os.environ.get('RAGD_RGE', f"{SB}/py-ukraine-interim-weekly-svi-260811")
BASE = os.environ.get('RAGD_BASE',
    f"{SB}/py-ukraine-interim-amortise-deepsetXcompAtt-net-260812")
file_prefix = "pcm_1_interim"
ANCHOR = int(os.environ.get('RAGD_ANCHOR', 4))   # interim used for i1-only head-ft
MAKE_PDFS = os.environ.get('RAGD_PDFS', '0') == '1'
OUT = os.environ.get('RAGD_OUT', BASE)           # write outputs here (BASE loads the net)
HEADFT = os.environ.get('RAGD_HEADFT', 'i1')     # i1 | expand (as §14.4.8)
os.makedirs(OUT, exist_ok=True)
SUF = os.environ.get('RAGD_SUF', 'RAGD')
N_REF = int(os.environ.get('RAGD_NREF', 503)); CLIP = 20.0   # trial total N (per application)
S = int(os.environ.get('RAG_S', 200))
TAUS = np.array([0.05, 0.25, 0.5, 0.75, 0.95], np.float32)
CACHE = f"/private/tmp/claude-501/-Users-or105-git-bIRTistic/49372e22-d11b-4882-a442-d1c60bcbdfb0/scratchpad/ragged_deploy_cells_{os.environ.get('RAGD_CTAG', 'focused')}"
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
net = Net(**dict(fit['net_kwargs']))             # respect head_mode / precision_pool
_mdf = f"{BASE}/{file_prefix}_meta_dim.npy"       # §14.4.9 item-amortised-family net
if os.path.exists(_mdf) or sig.size != J:
    _KM = 10; _Md = int(np.load(_mdf)[0]) if os.path.exists(_mdf) else 3
    sig = np.full(J, float(np.reshape(sig, -1)[0]), np.float32)   # scalar global sigma -> per-item
    _c = np.where(itype > .5, 2., 0.)             # caseness threshold: 2 categorical, 0 out-of-7
    _cols = [itype, ihigh, kmax / _KM, _c / _KM]  # canonical metadata order (type,dir,K/Kmax,c/Kmax)
    META = np.stack(_cols[:_Md], -1).astype(np.float32)
    print(f"  item-amortised net: global sigma={sig[0]:.3f}, META M={_Md}")


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
    # §14.4.25 Option A: empty future cohort -> p(rho|x) directly (single forward, no z^s mixing)
    aux_e = jnp.asarray(np.array([[1/np.sqrt(n), 0.0, n/N_REF, 0.0]], np.float32))
    be = dict(x_flat=xb, x_seg=xseg, z_flat=jnp.zeros((0, J, 2), np.float32),
              z_seg=jnp.zeros(0, jnp.int32), item_metadata=metab, query_idx=qidb, aux=aux_e)
    oute, heade = _fwd(fit['params'], be)
    qse = (np.maximum.accumulate(np.asarray(oute)[0], 1) * sig[:, None]).astype(np.float32)  # (J,5)
    hde = np.asarray(heade)[0].astype(np.float32)                                            # (J, A)
    wa = pd.read_pickle(f"{WK}/{file_prefix}_i{k}_regression_training.pkl")
    piv = wa.pivot_table(index='draw', columns='item_label', values='pps_ratio_x').reindex(columns=labels)
    tgt = np.full((S, J), np.nan)
    idx = piv.index[piv.index < S].to_numpy()
    tgt[idx] = np.clip(piv.loc[idx].to_numpy(np.float64), -CLIP, CLIP)
    return qs, hd.astype(np.float32), tgt, n, qse, hde


def build_splits(k, R=100):
    """§14.4.25 Option B (LOTP diagnostic): estimate p(rho|x) with NO simulation. Partition the
    REAL observed cohort x -> (x1, x2) many ways (balanced, R random permutations), feed x1 as
    current / x2 as future. Each forward conditions on x1 u x2 = x, so it IS p(rho|x); mixing
    over partitions averages the arbitrary split. n1+n2=n always -> no small-m inflation, and it
    uses the net's normal both-cohorts-present representation (unlike the empty-z bypass)."""
    xi = pd.read_csv(f"{WK}/{file_prefix}_{k}_data_dp1.csv")
    xpids = np.sort(xi.pid.unique()); n = len(xpids)
    x_raw = _pivot_pids(xi, 'y_stan', xpids)
    metab = jnp.asarray(META[None]); qidb = jnp.arange(J)[None]
    n1 = max(2, n // 2); n2 = n - n1
    aux = jnp.asarray(np.array([[1/np.sqrt(n1), 1/np.sqrt(max(n2, 1)), n1/N_REF, n2/N_REF]], np.float32))
    rngb = np.random.default_rng(1000 + k)
    qsb = np.empty((R, J, 5), np.float64); hdb = None
    for r in range(R):
        perm = rngb.permutation(n)
        bt = dict(x_flat=jnp.asarray(x_raw[perm[:n1]]), x_seg=jnp.zeros(n1, jnp.int32),
                  z_flat=jnp.asarray(x_raw[perm[n1:]]), z_seg=jnp.zeros(n2, jnp.int32),
                  item_metadata=metab, query_idx=qidb, aux=aux)
        o, h = _fwd(fit['params'], bt)
        qsb[r] = np.maximum.accumulate(np.asarray(o)[0], 1) * sig[:, None]
        if hdb is None:
            hdb = np.empty((R, J, h.shape[-1]), np.float32)
        hdb[r] = np.asarray(h)[0]
    return qsb.astype(np.float32), hdb


print(f"Phase 0: forward ragged net on {len(INTERIMS)} weekly interims (S={S}) ...")
QS, HD, TGT, NOBS, QSE, HDE = {}, {}, {}, {}, {}, {}
for k in INTERIMS:
    cf = f"{CACHE}/cell_{k}.npz"
    if os.path.exists(cf) and 'hde' in np.load(cf):
        z = np.load(cf)
        QS[k], HD[k], TGT[k], NOBS[k], QSE[k], HDE[k] = (
            z['qs'], z['hd'], z['tgt'], int(z['n']), z['qse'], z['hde'])
    else:
        QS[k], HD[k], TGT[k], NOBS[k], QSE[k], HDE[k] = build(k)
        np.savez(cf, qs=QS[k], hd=HD[k], tgt=TGT[k], n=NOBS[k], qse=QSE[k], hde=HDE[k])
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
cw.to_csv(f"{OUT}/{file_prefix}_pps_RAGD_contraction.csv", index=False)
def slope(col):
    d = cw[(cw[col] > 0) & (cw.n > 0)]
    return float(np.polyfit(np.log(d.n), np.log(d[col]), 1)[0])
print(f"\n(a) CONTRACTION slopes (weekly grid): W_svi={slope('W_svi'):.3f}  W_marg={slope('W_marg'):.3f}")
lo, hi = cw.n.min(), cw.n.max()
for c in ['W_svi', 'W_marg']:
    print(f"    {c}: n={lo}:{cw[cw.n==lo][c].median():.3f} -> n={hi}:{cw[cw.n==hi][c].median():.3f}")

# ---- (b) interim-1-only head fine-tune (mode-aware: plain | semiparam) ----
taus_j = jnp.asarray(TAUS)
HEAD_MODE = dict(fit['net_kwargs']).get('head_mode', 'plain')
HIDDEN = tuple(dict(fit['net_kwargs']).get('hidden_dims', (64, 64)))
EMB = int(dict(fit['net_kwargs']).get('embed_dim', 32))
A_AUX = HD[ANCHOR].shape[-1] - 2 * EMB            # aux tail width in head_in
ZQ = jnp.asarray([-1.6449, -0.6745, 0.0, 0.6745, 1.6449])[:len(TAUS)]
qpsi = _MLP(dims=(*HIDDEN, 5)); qpsi0 = {'params': fit['params']['params']['q_psi']}
g_head = _MLP(dims=(16, 1))                       # factored width-scaler (frozen)
g_params = {'params': fit['params']['params'].get('g_head', {})}
prec_anchor = 1.0 / np.sqrt(NOBS[ANCHOR])
print(f"head-ft mode={HEAD_MODE}  aux_dim={A_AUX}")


def _transform(raw, hd, prec):
    """standardised quantiles from q_psi raw output, per head_mode. hd is the
    head_in (its aux tail feeds the factored width-scaler)."""
    if HEAD_MODE == 'semiparam':                 # m + zq*softplus(s)*(1/sqrt n)
        return raw[..., 0:1] + ZQ[None, :] * jax.nn.softplus(raw[..., 1:2]) * prec
    if HEAD_MODE == 'factored':                  # med + (raw-med)*g(aux), g frozen
        med = raw[..., len(TAUS) // 2:len(TAUS) // 2 + 1]
        g = jax.nn.softplus(g_head.apply(g_params, hd[..., -A_AUX:]))
        return med + (raw - med) * g
    if HEAD_MODE == 'powerlaw':                  # m + zq*C*n^{-p}, learned C,p
        m = raw[..., 0:1]; C = jax.nn.softplus(raw[..., 1:2]); p = jax.nn.sigmoid(raw[..., 2:3])
        return m + ZQ[None, :] * (C * prec ** (2.0 * p))
    if HEAD_MODE == 'floor':                     # m + zq*sqrt(a^2 + b^2/n)
        m = raw[..., 0:1]; a = jax.nn.softplus(raw[..., 1:2]); b = jax.nn.softplus(raw[..., 2:3])
        return m + ZQ[None, :] * jnp.sqrt(a ** 2 + (b * prec) ** 2)
    return jnp.maximum.accumulate(raw, -1)       # plain: monotone cummax


def pinball(params, hd, y, prec):                # prec: (N,1) per-row 1/sqrt(n)
    pr = _transform(qpsi.apply(params, hd), hd, prec)
    e = y[:, None] - pr
    return jnp.mean(jnp.maximum(taus_j[None]*e, (taus_j[None]-1)*e))


opt0 = optax.adam(3e-4)
MARGW = float(os.environ.get('RAGD_MARGW', '0'))   # weight of marginal-CDF (mixture) loss
GG = 40                                            # grid points for the mixture CDF
SMARG = 64                                         # draws used to estimate the mixture (Fmix)


def _pool(ks):
    """pooled (head_in, standardised target, per-row 1/sqrt(n)) over interims ks."""
    HDp = np.concatenate([HD[i].reshape(-1, HD[i].shape[-1]) for i in ks])
    TGp = np.concatenate([(TGT[i] / sig[None]).reshape(-1) for i in ks])
    PRp = np.concatenate([np.full(HD[i].shape[0] * J, 1.0 / np.sqrt(NOBS[i]), np.float32) for i in ks])
    ok = np.isfinite(TGp)
    return HDp[ok].astype(np.float32), TGp[ok].astype(np.float32), PRp[ok, None].astype(np.float32)


# per-interim standardised marginal targets: grid + SVI empirical CDF (Gsvi) per item
MARG = {}
for _k in INTERIMS:
    GRk = np.zeros((J, GG), np.float32); Gs = np.zeros((J, GG), np.float32); drk = np.zeros(J, np.float32)
    for j in range(J):
        yv = (TGT[_k][:, j] / sig[j]); yv = yv[np.isfinite(yv)]
        if yv.size < 10:
            GRk[j] = np.linspace(-1, 1, GG); continue
        lo, hi = np.percentile(yv, [1, 99]); hi = max(hi, lo + 1e-3)
        GRk[j] = np.linspace(lo, hi, GG); drk[j] = (hi - lo) / GG
        Gs[j] = (yv[:, None] <= GRk[j][None]).mean(0)
    MARG[_k] = (jnp.asarray(GRk), jnp.asarray(Gs), jnp.asarray(drk))


HDM = {k: jnp.asarray(HD[k][:SMARG].reshape(-1, HD[k].shape[-1])) for k in INTERIMS}  # subsampled


def _marg_loss(params, k):
    """L2 distance between the amortiser mixture CDF and the SVI empirical CDF
    (standardised units), summed over items. The F-dependent part of the
    mixture CRPS; minimised when the marginal p(rho|x) matches SVI. Mixture
    estimated from SMARG draws for speed; the SVI target CDF uses all draws."""
    s = min(SMARG, HD[k].shape[0])
    Qk = _transform(qpsi.apply(params, HDM[k]), HDM[k], 1.0 / np.sqrt(NOBS[k])).reshape(s, J, 5)
    GR, Gs, dr = MARG[k]

    def perj(qj, grj):                              # qj (s,5), grj (GG,) -> mixture CDF (GG,)
        return jax.vmap(lambda q: jnp.interp(grj, q, taus_j, 0., 1.))(qj).mean(0)
    Fmix = jax.vmap(perj)(jnp.transpose(Qk, (1, 0, 2)), GR)   # (J, GG)
    return jnp.mean(jnp.sum((Fmix - Gs) ** 2 * dr[:, None], axis=1))


def _pool_e(ks):
    """empty-z pooled fit data: the single empty-z head_in per item, broadcast over the SVI
    marginal draws -> fits q_psi so p(rho|x) matches the SVI MARGINAL (spread + location)."""
    HDp = np.concatenate([np.broadcast_to(
        HDE[i][None], (TGT[i].shape[0], J, HDE[i].shape[-1])).reshape(-1, HDE[i].shape[-1]) for i in ks])
    TGp = np.concatenate([(TGT[i] / sig[None]).reshape(-1) for i in ks])
    PRp = np.concatenate([np.full(TGT[i].shape[0] * J, 1.0 / np.sqrt(NOBS[i]), np.float32) for i in ks])
    ok = np.isfinite(TGp)
    return HDp[ok].astype(np.float32), TGp[ok].astype(np.float32), PRp[ok, None].astype(np.float32)


def _fit_qpsi(ks, steps=800, poolf=None):
    HDp, TGp, PRp = (poolf or _pool)(ks)
    mks = [ks[i] for i in np.unique(np.linspace(0, len(ks) - 1, min(4, len(ks))).astype(int))]
    nstep = steps if MARGW == 0 else 300

    @jax.jit
    def step(params, ost, h, y, pr):
        def loss(p):
            lc = pinball(p, h, y, pr)
            if MARGW > 0:
                lm = jnp.mean(jnp.stack([_marg_loss(p, k) for k in mks]))
                return lc + MARGW * lm
            return lc
        l, g = jax.value_and_grad(loss)(params)
        up, ost = opt0.update(g, ost, params)
        return optax.apply_updates(params, up), ost, l

    params = qpsi0; ost = opt0.init(params); rng = np.random.default_rng(1); l = 0.0
    for it in range(nstep):
        ix = rng.integers(0, HDp.shape[0], min(2048, HDp.shape[0]))
        params, ost, l = step(params, ost, jnp.asarray(HDp[ix]), jnp.asarray(TGp[ix]), jnp.asarray(PRp[ix]))
    return params, float(l)


def _predict_hd(params, hd_k, k):
    """apply the fitted q_psi to any (M,J,A) head_in block -> (M,J,5) real scale."""
    flat = jnp.asarray(hd_k.reshape(-1, hd_k.shape[-1]))
    pr = np.asarray(_transform(qpsi.apply(params, flat), flat, 1.0 / np.sqrt(NOBS[k])))
    return pr.reshape(hd_k.shape[0], J, 5) * sig[None, :, None]


def _predict(params, k):
    return _predict_hd(params, HD[k], k)


HK, PKS = {}, {}                                 # PKS: expanding head-ft params per interim
if HEADFT == 'expand':                           # §14.4.8: refit q_psi on 1..k at each k
    for k in INTERIMS:
        pk, l = _fit_qpsi([i for i in INTERIMS if i <= k])
        HK[k] = _predict(pk, k); PKS[k] = pk
        print(f"  expand head-ft on 1..{k} done (loss {l:.4f})")
    HLABEL = 'head-ft-expand(ragged)'
else:                                            # interim-1(anchor)-only, reused for all k
    params, l = _fit_qpsi([ANCHOR]); print(f"\n(b) i1-only head-ft fit done (loss {l:.4f})")
    for k in INTERIMS:
        HK[k] = _predict(params, k); PKS[k] = params
    HLABEL = 'head-ft-i1(ragged)'

# §14.4.23 per-item EXPANDING median-shift (DEFAULT part of the head-ft): the prior-
# trained encoder systematically under-predicts every item (global shrinkage-to-prior
# bias); the shared head-ft fixes spread (PIT) but leaves a per-item location residual.
# Correct it with one scalar per item, fit on the SAME SVI 1..k the head-ft uses:
# Delta_j(k) = mean_{i<=k}(SVI marginal median - amortiser marginal median). Shifting ALL
# quantiles by Delta_j is a pure location move -> removes the offset without touching the
# (already-good) spread. On by default; RAGD_MEDSHIFT=0 disables it (ablation).
if os.environ.get('RAGD_MEDSHIFT', '1') == '1':
    def _margmed(qs):                            # marginal (mixture) median, real scale
        lo, hi = float(qs.min()), float(qs.max())
        if hi <= lo:
            return float(np.median(qs))
        gr = np.linspace(lo, hi, 200); F = np.zeros_like(gr)
        for s in range(qs.shape[0]):
            F += np.interp(gr, qs[s], TAUS, 0., 1.)
        return float(np.interp(0.5, F / qs.shape[0], gr))
    resid = {}                                   # residual per (interim, item), unshifted HK
    for k in INTERIMS:
        rj = np.full(J, np.nan)
        for j in range(J):
            y = TGT[k][:, j]; ok = np.isfinite(y); yv = y[ok]
            if yv.size >= 10:
                rj[j] = np.median(yv) - _margmed(HK[k][ok, j, :])
        resid[k] = rj
    ks_sorted = sorted(INTERIMS)
    HK0 = {k: HK[k] for k in INTERIMS}           # freeze the unshifted preds for residuals
    for k in ks_sorted:
        dj = np.nan_to_num(np.nanmean(
            np.stack([resid[i] for i in ks_sorted if i <= k]), axis=0))   # (J,)
        HK[k] = HK0[k] + dj[None, :, None]
    print("  applied per-item expanding median-shift (default; RAGD_MEDSHIFT=0 to disable)")
else:
    HLABEL = HLABEL.replace('(ragged)', '-noshift(ragged)')

# §14.4.25 Option A: effect-size marginal p(rho|x) from the EMPTY-z forward (no z^s mixing,
# so no finite-m Monte-Carlo inflation). Per-item expanding median-shift, then broadcast to
# (S,J,5) so the existing marg-KS / plot code sees a degenerate "mixture" = the single p(rho|x).
EMPTYZ = os.environ.get('RAGD_EMPTYZ', '0') == '1'
MARGSRC = None
if EMPTYZ:
    def _predict_e(params, k):                            # spread+location-calibrated empty-z (J,5)
        raw = qpsi.apply(params, jnp.asarray(HDE[k]))
        pr = np.asarray(_transform(raw, jnp.asarray(HDE[k]), 1.0 / np.sqrt(NOBS[k])))
        return np.maximum.accumulate(pr, 1) * sig[:, None]
    ks_sorted = sorted(INTERIMS); HKE = {}
    for k in ks_sorted:                                   # expanding head-ft on the EMPTY-z pathway
        pk, _ = _fit_qpsi([i for i in INTERIMS if i <= k], poolf=_pool_e)
        HKE[k] = np.broadcast_to(_predict_e(pk, k)[None], (S, J, 5)).copy()
    MARGSRC = HKE
    HLABEL = HLABEL.replace('(ragged)', '+emptyz(ragged)')
    print("  §14.4.25 Option A: effect-size marginal = spread-calibrated EMPTY-z forward (no z^s mixing)")

# §14.4.25 Option B: effect-size marginal from real-data SPLITS of x (no simulation). Forward the
# net on R balanced partitions (x1,x2) of the observed cohort, head-ft with the SAME conditional
# q_psi, mix over partitions; per-item affine shift as usual. n1+n2=n -> no small-m inflation.
if os.environ.get('RAGD_SPLITB', '0') == '1':
    print("  §14.4.25 Option B: building real-data split marginals (no simulation) ...")

    def _mixmed(qs):                                     # mixture median of (R,5) quantiles
        lo, hi = float(qs.min()), float(qs.max())
        if hi <= lo:
            return float(np.median(qs))
        gr = np.linspace(lo, hi, 200); F = np.zeros_like(gr)
        for s in range(qs.shape[0]):
            F += np.interp(gr, qs[s], TAUS, 0., 1.)
        return float(np.interp(0.5, F / qs.shape[0], gr))
    HKB0 = {}
    for k in INTERIMS:
        HKB0[k] = _predict_hd(PKS[k], build_splits(k)[1], k)   # (R,J,5) head-ft'd split conditionals
        print(f"  split-B interim {k} (n={NOBS[k]}) done")
    resid = {}
    for k in INTERIMS:
        rj = np.full(J, np.nan)
        for j in range(J):
            y = TGT[k][:, j]; ok = np.isfinite(y); yv = y[ok]
            if yv.size >= 10:
                rj[j] = np.median(yv) - _mixmed(HKB0[k][:, j, :])
        resid[k] = rj
    ks_sorted = sorted(INTERIMS); HKB = {}
    for k in ks_sorted:
        dj = np.nan_to_num(np.nanmean(np.stack([resid[i] for i in ks_sorted if i <= k]), 0))
        HKB[k] = HKB0[k] + dj[None, :, None]
    MARGSRC = HKB
    HLABEL = HLABEL.replace('(ragged)', '+splitB(ragged)')

# §14.4.26 BvM between-first self-consistency correction (post-hoc on the affine-shifted HK).
# Per item, fit the marginal power law SD=C n^-p to the SVI SD; then at each interim the BvM
# targets are T=(C n^-p)^2 (marginal) and W_true=(C (n+m)^-p)^2 (conditional at n+m=N_REF), so the
# non-spurious between is B*=T-W_true. Scale each conditional width by alpha=sqrt(W_true/W) and
# shrink the per-draw medians toward their mean by beta=sqrt(B*/B); then the mixture variance =
# W_true + B* = T exactly. Intrinsically justified for the power-law/floor heads. RAGD_BVM=1.
if os.environ.get('RAGD_BVM', '0') == '1':
    print("  §14.4.26 BvM between-first correction (per-item SVI power-law targets) ...")
    Cj = np.full(J, np.nan); pj = np.full(J, np.nan)          # per-item marginal power law
    for j in range(J):
        ns, sds = [], []
        for k in INTERIMS:
            yv = TGT[k][:, j]; yv = yv[np.isfinite(yv)]
            if yv.size >= 10 and np.std(yv) > 1e-9:
                ns.append(NOBS[k]); sds.append(np.std(yv))
        if len(ns) >= 3:
            b = np.polyfit(np.log(ns), np.log(sds), 1); pj[j] = -b[0]; Cj[j] = float(np.exp(b[1]))
    HKV = {}
    for k in INTERIMS:
        n = NOBS[k]; Hk = HK[k].copy()
        for j in range(J):
            if not np.isfinite(Cj[j]):
                continue
            C, p = Cj[j], pj[j]
            T = (C * n ** (-p)) ** 2                          # marginal var target
            Wt = (C * N_REF ** (-p)) ** 2                     # conditional var (n+m = N_REF)
            Bs = max(T - Wt, 0.0)                             # non-spurious between (>=0)
            qs = HK[k][:, j, :]; mu = qs[:, 2]; mub = float(mu.mean())
            sg = (qs[:, 4] - qs[:, 0]) / 3.2897
            W = float(np.mean(sg ** 2)); B = float(np.var(mu))
            a = np.sqrt(Wt / W) if W > 1e-12 else 1.0
            be = np.sqrt(Bs / B) if B > 1e-12 else 0.0
            Hk[:, j, :] = (mub + be * (mu - mub))[:, None] + a * (qs - mu[:, None])
        HKV[k] = Hk
    MARGSRC = HKV
    HLABEL = HLABEL.replace('(ragged)', '+bvm(ragged)')


# --- PPS across items/interims (§14.1.6): P(H1|x,z^s) = P(rho>eta0|x,z^s) via CDF interp;
#     PPS = mean over future draws s of 1{ P(H1|x,z^s) > etaH }. Uses the final per-draw
#     conditional quantiles (BvM-corrected when RAGD_BVM=1, else the affine-shifted HK). ---
ETA0 = float(os.environ.get('RAGD_ETA0', '0.5'))
ETAH = float(os.environ.get('RAGD_ETAH', '0.89'))
PPSSRC = HKV if os.environ.get('RAGD_BVM', '0') == '1' else HK
_ppsrows = []
for k in INTERIMS:
    for j in range(J):
        q = PPSSRC[k][:, j, :]                                    # (S,5) per-draw conditional quantiles
        ph1 = 1.0 - np.array([np.interp(ETA0, q[s], TAUS, 0.0, 1.0) for s in range(q.shape[0])])
        _ppsrows.append(dict(interim_id=k, n=NOBS[k], item_label=labels[j],
                             pps=float(np.mean(ph1 > ETAH)), eta0=ETA0, etaH=ETAH, method=HLABEL))
pd.DataFrame(_ppsrows).to_csv(f"{OUT}/{file_prefix}_pps_RAGD_pps_by_item.csv", index=False)
print(f"  wrote PPS by item ({len(INTERIMS)} interims x {J} items), eta0={ETA0} etaH={ETAH} -> {HLABEL}")

# --- p(H1|x,z) boxplot pkl (for cross-method comparison, §14.5): per (interim,item),
#     quantiles of p(H1|x,z^s)=P(rho>eta0|x,z^s) over the future draws s. Keyed on
#     item_label + interim_month_year so the compare script can dodge it beside SVI/HMC/IS. ---
_MONTHYR = {}
for k in INTERIMS:
    try:
        _dt = pd.to_datetime(pd.read_csv(f"{WK}/{file_prefix}_{k}_data_dp1.csv",
                                         usecols=['submission_date'])['submission_date'],
                             errors='coerce', utc=True).max()
        _MONTHYR[k] = _dt.strftime('%Y-%b') if pd.notna(_dt) else f"interim {k}"
    except Exception:
        _MONTHYR[k] = f"interim {k}"
_boxrows = []
for k in INTERIMS:
    for j in range(J):
        q = PPSSRC[k][:, j, :]
        ph1 = 1.0 - np.array([np.interp(ETA0, q[s], TAUS, 0.0, 1.0) for s in range(q.shape[0])])
        ph1 = ph1[np.isfinite(ph1)]
        if ph1.size < 5:
            continue
        qq = np.percentile(ph1, [2.5, 25, 50, 75, 97.5])
        _boxrows.append(dict(item_label=labels[j], interim_id=k, interim_month_year=_MONTHYR[k],
                             n=NOBS[k], q025=qq[0], q25=qq[1], q50=qq[2], q75=qq[3], q975=qq[4],
                             eta0=ETA0, method=HLABEL))
pd.DataFrame(_boxrows).to_pickle(f"{OUT}/{file_prefix}_pps_RAGD_p_h1_xz_boxplot.pkl")
print(f"  wrote p(H1|x,z) boxplot pkl ({len(_boxrows)} rows, eta0={ETA0})")


def cov_pit(getq, margq=None):
    out = []
    for k in INTERIMS:
        us = []; mks = []
        for j in range(J):
            y = TGT[k][:, j]; ok = np.isfinite(y); yv = y[ok]
            if yv.size < 10:
                continue
            q = getq(k)[ok, j, :]
            us.append(pit(q, yv))                          # conditional PIT (draw-aligned)
            qm = (margq(k) if margq is not None else getq(k))[:, j, :]    # marg source: ALL mixture rows
            qm = qm[np.isfinite(qm).all(1)]                # (not SVI-draw-aligned; R may != S)
            # marg-KS: SVI empirical CDF vs amortiser mixture CDF (the CDF-plot area)
            lo = min(yv.min(), qm[:, 0].min()); hi = max(yv.max(), qm[:, 4].max())
            gr = np.linspace(lo, hi, 200)
            Fs = (yv[:, None] <= gr[None]).mean(0)
            Fm = np.mean([np.interp(gr, qm[s], TAUS, 0., 1.) for s in range(qm.shape[0])], 0)
            mks.append(float(np.max(np.abs(Fs - Fm))))
        u = np.concatenate(us); uu = np.sort(u); ec = np.arange(1, len(u)+1)/len(u)
        out.append(dict(interim_id=k, n=NOBS[k], cov5=float((u <= .5).mean()),
                        cov95=float((u <= .95).mean()), pit_ks=float(np.max(np.abs(ec-uu))),
                        marg_ks=float(np.mean(mks))))
    return pd.DataFrame(out)


_MARGF = (lambda k: MARGSRC[k]) if MARGSRC is not None else None
base = cov_pit(lambda k: QS[k]).assign(config='baseline(ragged)')
h1 = cov_pit(lambda k: HK[k], _MARGF).assign(config=HLABEL)   # PIT=conditional, marg=empty-z if EMPTYZ
res = pd.concat([base, h1], ignore_index=True)
res.to_csv(f"{OUT}/{file_prefix}_pps_RAGD_headft_coverage.csv", index=False)
print(f"\n=== MEAN over interims (nominal cov5=.50 cov95=.95; head-ft={HEADFT}) ===")
print(res.groupby('config')[['cov5', 'cov95', 'pit_ks', 'marg_ks']].mean().round(3).to_string())
print(f"\n{HLABEL} per interim:")
print(h1[['interim_id', 'n', 'cov5', 'cov95', 'pit_ks', 'marg_ks']].to_string(index=False))


# ---- (c) diagnostic PDF suite (matches the Xcomp/Scomp format) ----
def _marg_cdf(qs, grid):
    F = np.zeros_like(grid)
    for s in range(qs.shape[0]):
        F += np.interp(grid, qs[s], TAUS, 0., 1.)
    return F / qs.shape[0]


def save_pdfs(getq, suf):
    covr, pk, mc, sc = [], [], [], []
    for k in INTERIMS:
        for j in range(J):
            y = TGT[k][:, j]; ok = np.isfinite(y); yv = y[ok]
            if yv.size < 10:
                continue
            qs = getq(k)[ok, j, :]; u = pit(qs, yv)
            for e in TAUS:
                covr.append(dict(interim_id=k, item_label=labels[j], eta_inpol=float(e),
                                 coverage=float((u <= e).mean())))
            uu = np.sort(u); ec = np.arange(1, len(u)+1)/len(u)
            pk.append(dict(item_label=labels[j], ks=float(np.max(np.abs(ec-uu))), pit=u))
            lo, hi = float(min(yv.min(), qs.min())), float(max(yv.max(), qs.max()))
            gr = np.linspace(lo, hi, 200); Fm = _marg_cdf(qs, gr); Fs = (yv[:, None] <= gr[None]).mean(0)
            idx = np.linspace(0, 199, 60).astype(int)
            for src, F in (('SVI  p(rho|x)', Fs), ('amortiser  p_hat(rho|x)', Fm)):
                for rv, cv in zip(gr[idx], F[idx]):
                    mc.append(dict(interim_id=k, item_label=labels[j], source=src, rho=float(rv), cdf=float(cv)))
            si = np.linspace(0, len(yv)-1, min(400, len(yv))).astype(int)
            for a, b2 in zip(qs[si, 2], yv[si]):
                sc.append(dict(interim_id=k, item_label=labels[j], med=float(a), y=float(b2)))
    d = OUT
    (ggplot(pd.DataFrame(covr), aes('eta_inpol', 'coverage', colour='factor(interim_id)', group='interim_id'))
     + geom_abline(intercept=0, slope=1, linetype='dashed', colour='black')
     + geom_line(size=.4, alpha=.8) + geom_point(size=.8) + facet_wrap('~ item_label', ncol=4)
     + theme_bw() + theme(figure_size=(12, 12.5), legend_position='top', strip_background=element_blank(),
       strip_text=element_text(face='bold')) + labs(x='nominal', y='empirical coverage', colour='interim',
       title=f'Coverage calibration — {suf}')).save(
        f"{d}/{file_prefix}_pps_{suf}_tests_coverage.pdf", verbose=False, limitsize=False)
    pl = pd.concat([pd.DataFrame({'item_label': r['item_label'], 'u': r['pit']}) for r in pk])
    (ggplot(pl, aes('u')) + geom_histogram(aes(y='..density..'), bins=20, fill='#1f77b4', colour='white', size=.2)
     + geom_hline(yintercept=1, linetype='dashed', colour='black') + facet_wrap('~ item_label', ncol=4)
     + theme_bw() + theme(figure_size=(12, 12.5), strip_background=element_blank(),
       strip_text=element_text(face='bold')) + labs(x='PIT u', y='density', title=f'PIT uniformity — {suf}')).save(
        f"{d}/{file_prefix}_pps_{suf}_tests_pit.pdf", verbose=False, limitsize=False)
    mcdf = pd.DataFrame(mc); scdf = pd.DataFrame(sc)
    for k in INTERIMS:
        s = mcdf[mcdf.interim_id == k]
        if len(s):
            (ggplot(s, aes('rho', 'cdf', colour='source', group='source')) + geom_line(size=.6)
             + facet_wrap('~ item_label', ncol=4, scales='free_x') + theme_bw()
             + theme(figure_size=(15, 14), legend_position='top', strip_background=element_blank(),
               strip_text=element_text(face='bold')) + labs(x=f'rho at interim {k}', y='CDF', colour='',
               title=f'Marginal p(rho|x) — {suf}, interim {k}')).save(
                f"{d}/{file_prefix}_i{k}_svi_vs_amortiser_marginal_cdf.pdf", verbose=False, limitsize=False)
        s2 = scdf[scdf.interim_id == k]
        if len(s2):
            (ggplot(s2, aes('med', 'y')) + geom_point(alpha=.25, size=.5, colour='#1f77b4')
             + geom_smooth(method='glm', method_args={'family': 'gaussian'}, se=False, colour='black', size=.6)
             + facet_wrap('~ item_label', ncol=4, scales='free') + theme_bw()
             + theme(figure_size=(15, 14), strip_background=element_blank(), strip_text=element_text(face='bold'))
             + labs(x=f'amortiser median rho_hat_3 at interim {k}', y=f'SVI rho at interim {k}',
               title=f'median association — {suf}, interim {k}')).save(
                f"{d}/{file_prefix}_i{k}_svi_rho_vs_amortiser_median.pdf", verbose=False, limitsize=False)
    print(f"  saved PDF suite ({suf})")


def save_contraction_cdf(getq, suf):
    """Per-item quantile boxplots of the marginal p(rho|x) across interims.
    ALL box statistics are PRECOMPUTED (geom_boxplot(stat='identity'), NOT the
    default Tukey whisker rule): whiskers = 10/90% quantiles, box = 25/75%,
    midline = median. x = interim id; y = p(rho|x) estimate; SVI (blue) vs
    amortiser (red) dodged side by side; one facet per item. Replaces the earlier
    unreadable CDF-facet plot. SVI quantiles are of the SVI draws; amortiser
    quantiles are read off the inverted mixture CDF Fmix (mean over z^(s) draws)."""
    QL = np.array([0.10, 0.25, 0.50, 0.75, 0.90])
    rows = []
    for j in range(J):
        for k in INTERIMS:
            y = TGT[k][:, j]; ok = np.isfinite(y); yv = y[ok]
            if yv.size < 10:
                continue
            qs = getq(k)[:, j, :]                                     # ALL marg-source rows (R may != S)
            qs = qs[np.isfinite(qs).all(1)]
            sq = np.percentile(yv, QL * 100.0)                        # SVI marginal quantiles
            lo, hi = float(min(yv.min(), qs.min())), float(max(yv.max(), qs.max()))
            gr = np.linspace(lo, hi, 400); Fm = _marg_cdf(qs, gr)
            aq = np.interp(QL, Fm, gr)                                # amortiser mixture quantiles
            for src, q in (('SVI', sq), ('amortiser', aq)):
                rows.append(dict(interim_id=int(k), item_label=labels[j], source=src,
                                 ymin=float(q[0]), lower=float(q[1]), middle=float(q[2]),
                                 upper=float(q[3]), ymax=float(q[4])))
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/{file_prefix}_pps_{suf}_contraction_cdf_by_item.csv", index=False)  # box data (calibrated marginal)
    df['interim_id'] = pd.Categorical(df.interim_id, categories=sorted(df.interim_id.unique()), ordered=True)
    (ggplot(df, aes('interim_id', ymin='ymin', lower='lower', middle='middle',
                    upper='upper', ymax='ymax', fill='source'))
     + geom_boxplot(stat='identity', position=position_dodge(width=0.78),
                    size=.3, width=.7, alpha=.85)
     + facet_wrap('~ item_label', ncol=4, scales='free_y')
     + scale_fill_manual(values={'SVI': CBLUE, 'amortiser': CRED})
     + theme_bw() + theme(figure_size=(16, 20), legend_position='top', panel_spacing=0.02,
       strip_background=element_blank(), strip_text=element_text(face='bold', size=8),
       axis_text_x=element_text(size=6))
     + labs(x='interim id', y='p(rho | x) estimate', fill='',
       title=(f'Marginal p(rho|x): SVI vs amortiser — {suf}  '
              + ('[RAW network: no head-ft, no power-law recalibration]'
                 if suf.endswith('base') else '[deployed: head-ft + affine shift]')
              + '  (box 25-75%, midline median, whiskers 10-90%)'))).save(
        f"{OUT}/{file_prefix}_pps_{suf}_contraction_cdf_by_item.pdf", verbose=False, limitsize=False)
    print(f"  saved contraction quantile-box ({suf})")


CBLUE, CRED = '#1f77b4', '#d62728'


def save_pit_box(getq, suf):
    """CONDITIONAL calibration (PIT) as quantile boxes — the PIT-KS analogue of the marginal
    box. Per item, x=interim; box the PIT values u = F_amortiser(rho_SVI | x, z^(s)) (whiskers
    10/90, box 25/75, midline median, ALL precomputed). A calibrated conditional has u~Uniform,
    so a perfectly calibrated box sits exactly on the dashed grey references (median 0.5, box
    0.25-0.75, whiskers 0.1-0.9). Deviations read directly: midline above 0.5 => amortiser
    median too low (SVI draws sit high); box WIDER than 0.25-0.75 => intervals too NARROW
    (over-confident); box narrower => under-confident. The KS gap to uniform is PIT-KS."""
    QL = np.array([0.10, 0.25, 0.50, 0.75, 0.90]); REF = QL
    rows = []
    for j in range(J):
        for k in INTERIMS:
            y = TGT[k][:, j]; ok = np.isfinite(y); yv = y[ok]
            if yv.size < 10:
                continue
            u = pit(getq(k)[ok, j, :], yv)               # PIT values (conditional)
            q = np.quantile(u, QL)
            rows.append(dict(interim_id=int(k), item_label=labels[j], ymin=float(q[0]),
                             lower=float(q[1]), middle=float(q[2]), upper=float(q[3]), ymax=float(q[4])))
    df = pd.DataFrame(rows)
    df['interim_id'] = pd.Categorical(df.interim_id, categories=sorted(df.interim_id.unique()), ordered=True)
    df['dev'] = df['middle'] - 0.5                        # median PIT deviation from calibrated 0.5
    (ggplot(df, aes('interim_id', ymin='ymin', lower='lower', middle='middle',
                    upper='upper', ymax='ymax', fill='dev'))
     + geom_hline(yintercept=list(REF), linetype='dashed', colour='#9e9e9e', size=.3)
     + geom_boxplot(stat='identity', alpha=.9, size=.3, width=.7)
     + scale_fill_gradient2(low='#2166ac', mid='#f7f7f7', high='#b2182b', midpoint=0.0,
                            name='PIT median − 0.5\n(blue: amortiser median too high;\nred: too low)')
     + facet_wrap('~ item_label', ncol=4)
     + theme_bw() + theme(figure_size=(16, 20), panel_spacing=0.02, legend_position='top',
       strip_background=element_blank(), strip_text=element_text(face='bold', size=8),
       axis_text_x=element_text(size=6))
     + labs(x='interim id', y='PIT   u = F_amortiser(rho_SVI | x, z)',
       title=f'Conditional calibration (PIT) — {suf}  '
             f'(box on dashed refs .1/.25/.5/.75/.9 = calibrated; fill = median bias, '
             f'box wider than .25-.75 = over-confident)')).save(
        f"{OUT}/{file_prefix}_pps_{suf}_pit_box_by_item.pdf", verbose=False, limitsize=False)
    print(f"  saved PIT quantile-box ({suf})")


def save_contraction_factor(getq, suf):
    """Posterior-contraction rate: SD of rho vs sqrt(n), one facet per item.
    SVI posterior SD (blue) vs amortiser marginal-predictive SD (red), each with
    a linear fit; per item R^2 and slope annotated (SVI top, amortiser bottom).
    Under BvM SD ~ 1/sqrt(n); a shallower amortiser slope = under-contraction."""
    rows = []
    for j in range(J):
        for k in INTERIMS:
            y = TGT[k][:, j]; ok = np.isfinite(y); yv = y[ok]
            if yv.size < 10:
                continue
            sn = float(np.sqrt(NOBS[k]))
            rows.append(dict(item_label=labels[j], sqrt_n=sn, sd=float(np.std(yv)), source='SVI'))
            qs = getq(k)[:, j, :]; qs = qs[np.isfinite(qs).all(1)]   # ALL marg rows (R may != S)
            mu = qs[:, 2]; sdw = (qs[:, 4] - qs[:, 0]) / 3.2897   # per-draw sd (90%/z)
            asd = float(np.sqrt(np.mean(sdw ** 2) + np.var(mu)))  # marginal-predictive sd
            rows.append(dict(item_label=labels[j], sqrt_n=sn, sd=asd, source='amortiser'))
    df = pd.DataFrame(rows)

    def _fit(g):
        b = np.polyfit(g.sqrt_n, g.sd, 1); pred = np.polyval(b, g.sqrt_n)
        ss = np.sum((g.sd - pred) ** 2); st = np.sum((g.sd - g.sd.mean()) ** 2)
        return float(b[0]), (1 - ss / st if st > 0 else np.nan)

    tl, tr, bl, br = [], [], [], []      # svi R2 / svi slope / amo R2 / amo slope
    for j, gi in df.groupby('item_label'):
        xlo, xhi = gi.sqrt_n.min(), gi.sqrt_n.max(); ylo, yhi = gi.sd.min(), gi.sd.max()
        for src in ('SVI', 'amortiser'):
            g = gi[gi.source == src]
            if len(g) < 3 or g.sd.std() == 0:
                continue
            sl, r2 = _fit(g)
            if src == 'SVI':
                tl.append(dict(item_label=j, sqrt_n=xlo, sd=yhi, label=f'R2 = {r2:.2f}'))
                tr.append(dict(item_label=j, sqrt_n=xhi, sd=yhi, label=f'slope = {sl:.3f}'))
            else:
                bl.append(dict(item_label=j, sqrt_n=xlo, sd=ylo, label=f'R2 = {r2:.2f}'))
                br.append(dict(item_label=j, sqrt_n=xhi, sd=ylo, label=f'slope = {sl:.3f}'))
    (ggplot(df, aes('sqrt_n', 'sd', colour='source'))
     + geom_point(size=1.5) + geom_smooth(method='lm', se=False, size=.5)
     + scale_color_manual(values={'SVI': CBLUE, 'amortiser': CRED})
     + geom_text(pd.DataFrame(tl), aes('sqrt_n', 'sd', label='label'), inherit_aes=False,
                 ha='left', va='top', size=8, colour=CBLUE)
     + geom_text(pd.DataFrame(tr), aes('sqrt_n', 'sd', label='label'), inherit_aes=False,
                 ha='right', va='top', size=8, colour=CBLUE)
     + geom_text(pd.DataFrame(bl), aes('sqrt_n', 'sd', label='label'), inherit_aes=False,
                 ha='left', va='bottom', size=8, colour=CRED)
     + geom_text(pd.DataFrame(br), aes('sqrt_n', 'sd', label='label'), inherit_aes=False,
                 ha='right', va='bottom', size=8, colour=CRED)
     + facet_wrap('~ item_label', ncol=4, scales='free') + theme_bw()
     + theme(figure_size=(15, 14), legend_position='top', strip_background=element_blank(),
       strip_text=element_text(face='bold'))
     + labs(x='sqrt(number of participants)  sqrt(n)', y='posterior SD of rho', colour='',
       title=f'Posterior contraction factor — {suf}  (SVI blue, amortiser red)')).save(
        f"{OUT}/{file_prefix}_pps_{suf}_posterior-contraction-factor.pdf",
        verbose=False, limitsize=False)
    print(f"  saved contraction-factor ({suf})")


def save_contraction_law(getq, suf):
    """Trained (amortiser, red) vs actual (SVI, blue) contraction with a power-
    law fit SD=C n^-p per item/series; exponent p annotated (SVI top, amortiser
    bottom). Shows whether the net learned the item-specific rate."""
    pts = []
    for j in range(J):
        for k in INTERIMS:
            y = TGT[k][:, j]; ok = np.isfinite(y); yv = y[ok]
            if yv.size < 10 or labels[j] == BLOW:
                continue
            q = getq(k)[:, j, :]; q = q[np.isfinite(q).all(1)]      # ALL marg rows (R may != S)
            mu = q[:, 2]; sdw = (q[:, 4] - q[:, 0]) / 3.2897
            pts.append(dict(item_label=labels[j], n=NOBS[k], sqrt_n=float(np.sqrt(NOBS[k])),
                            sd=float(np.std(yv)), source='SVI'))
            pts.append(dict(item_label=labels[j], n=NOBS[k], sqrt_n=float(np.sqrt(NOBS[k])),
                            sd=float(np.sqrt(np.mean(sdw ** 2) + np.var(mu))), source='amortiser'))
    df = pd.DataFrame(pts)

    def _pow(g):
        n = g.n.values.astype(float); y = np.maximum(g.sd.values, 1e-6)
        b = np.polyfit(np.log(n), np.log(y), 1); p = -b[0]; C = np.exp(b[1])
        pr = C * n ** (-p); ss = np.sum((y - pr) ** 2); st = np.sum((y - y.mean()) ** 2)
        return p, C, (1 - ss / st if st > 0 else np.nan)

    cur, atop, abot = [], [], []
    for j, gi in df.groupby('item_label'):
        ylo, yhi = gi.sd.min(), gi.sd.max(); xlo = gi.sqrt_n.min()
        for src, isS in (('SVI', True), ('amortiser', False)):
            g = gi[gi.source == src]
            if len(g) < 4:
                continue
            p, C, r2 = _pow(g); ng = np.linspace(g.n.min(), g.n.max(), 60)
            for xx, yy in zip(np.sqrt(ng), C * ng ** (-p)):
                cur.append(dict(item_label=j, source=src, sqrt_n=xx, sd=yy))
            (atop if isS else abot).append(dict(item_label=j, sqrt_n=xlo, sd=(yhi if isS else ylo),
                                                label=f'p={p:.2f} (R2={r2:.2f})'))
    (ggplot(df, aes('sqrt_n', 'sd', colour='source')) + geom_point(size=1.4)
     + geom_line(pd.DataFrame(cur), aes('sqrt_n', 'sd', colour='source'), size=.6)
     + scale_color_manual(values={'SVI': CBLUE, 'amortiser': CRED})
     + geom_text(pd.DataFrame(atop), aes('sqrt_n', 'sd', label='label'), inherit_aes=False, ha='left', va='top', size=8, colour=CBLUE)
     + geom_text(pd.DataFrame(abot), aes('sqrt_n', 'sd', label='label'), inherit_aes=False, ha='left', va='bottom', size=8, colour=CRED)
     + facet_wrap('~ item_label', ncol=4, scales='free') + theme_bw()
     + theme(figure_size=(15, 14), legend_position='top', strip_background=element_blank(),
       strip_text=element_text(face='bold'))
     + labs(x='sqrt(number of participants)  sqrt(n)', y='posterior SD of rho', colour='',
       title=f'Trained (amortiser) vs actual (SVI) contraction — {suf}  (power-law SD=C n^-p)')).save(
        f"{OUT}/{file_prefix}_pps_{suf}_contraction-law_trained-vs-svi.pdf", verbose=False, limitsize=False)
    print(f"  saved contraction-law ({suf})")


BLOW = "CG-VIO_ph-punish"
if MAKE_PDFS:
    print("\n(c) diagnostic PDF suite ...")
    _factoronly = os.environ.get('RAGD_FACTORONLY', '0') == '1'
    if not _factoronly and os.environ.get('RAGD_CDFONLY', '0') != '1':
        save_pdfs(lambda k: HK[k], SUF)          # head-ft-i1 (the deployable result)
        save_pdfs(lambda k: QS[k], f"{SUF}base")  # baseline (no recalibration)
    _mplot = (lambda k: MARGSRC[k]) if MARGSRC is not None else (lambda k: HK[k])  # marginal source
    if not _factoronly:
        save_contraction_cdf(_mplot, SUF)            # effect-size marginal (empty-z if EMPTYZ)
        save_contraction_cdf(lambda k: QS[k], f"{SUF}base")
        save_pit_box(lambda k: HK[k], SUF)           # CONDITIONAL calibration (PIT-KS) quantile-box
    save_contraction_factor(_mplot, SUF)             # SD-vs-n: empty-z marginal has no MC inflation
    save_contraction_law(_mplot, SUF)
if os.environ.get('RAGD_ETA0GRID', '0') == '1':      # eta_0 deployment sweep (§14.4.8 toolkit)
    from amortiser_diag_plots import eta0_sweep
    _eg = [float(x) for x in os.environ.get('RAGD_ETA0GRID_VALS', '0,0.25,0.5,0.75,1.0').split(',')]
    eta0_sweep(lambda k: PPSSRC[k], lambda k: TGT[k], lambda k: NOBS[k],
               labels, INTERIMS, np.asarray(TAUS), OUT, file_prefix, SUF, _eg,
               etaH=ETAH, date_of=lambda k: _MONTHYR.get(k, k))
print(f"\n-> {OUT}")
