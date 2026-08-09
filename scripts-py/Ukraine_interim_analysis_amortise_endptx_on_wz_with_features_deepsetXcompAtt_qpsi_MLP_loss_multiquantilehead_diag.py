#!/usr/bin/env python3
"""
Detailed calibration diagnostics + remedies for the prior-predictive-
trained nested-DeepSets amortiser (`deepsetXcompAtt`, §14.4 of
dev/amortised_decision_making.md).

Mirrors the item-level §14.1/§14.2 diagnostics (coverage curves, PIT
histograms, marginal-CDF overlays, median-association scatter + CSVs)
for the deepset architecture, plus the same three remedies:

  baseline  prior-predictive net as trained
  A1        conformal recalibration on interim-1 PIT, applied to all k
  A2        conformal recalibration EXPANDING: interim k recalibrated 1..k
  H         head-only fine-tune EXPANDING: interim k = q_psi refit on 1..k

Unlike the item-level script, the deepset net consumes the RAW cohorts
(x_obs, z^s), so Phase 0 forward-passes every draw once, capturing BOTH
the full 5-quantile predictions AND the frozen-encoder head input
`head_in` (sown in the net) so the head-only remedy can fine-tune q_psi
on cached embeddings without re-running the encoder.
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
from flax import linen as nn
from plotnine import (ggplot, aes, geom_abline, geom_line, geom_point,
                      geom_smooth, geom_histogram, geom_hline, facet_wrap,
                      theme_bw, theme, element_blank, element_text, labs)

warnings.filterwarnings('ignore')
from amortiser_common import load_fitted_model
from amortiser_pps_features_deepsetXcompAtt_qpsi_MLP_loss_multiquantilehead import (
    Amortiser_PPS_features_deepsetXcompAtt_qpsi_MLP_loss_multiquantilehead as DSNet,
    _MLP,
)
from model_pcm import PartialCreditModel

DIR_RGE = "/Users/or105/sandbox/bIRTistic/py-ukraine-interim-with-regression-on-endptx-wz-260601"
SB = "/Users/or105/sandbox/bIRTistic"
B = "py-ukraine-interim-amortise-endptx-on-wz-with-features-deepsetXcompAtt-qpsi-MLP-loss-multiquantilehead"
BASE = f"{SB}/{B}-260805"
file_prefix = "pcm_1_interim"
SUF = "RGDX"
N_FULL, CLIP = 503, 20.0
S = int(os.environ.get('UKR_DEEPSET_S', 200))
TAUS = np.array([0.05, 0.25, 0.5, 0.75, 0.95], np.float32)
INTERIMS = list(range(1, 9))
CACHE = "/private/tmp/claude-501/-Users-or105-git-bIRTistic/49372e22-d11b-4882-a442-d1c60bcbdfb0/scratchpad/deepset_cells"
os.makedirs(CACHE, exist_ok=True)

# ---- item metadata (shared ordering with the deploy script) ----
dit = pd.read_csv(f"{DIR_RGE}/{file_prefix}_1_data_dit.csv")
_wa1 = pd.read_pickle(f"{DIR_RGE}/{file_prefix}_i1_regression_training.pkl")
items = (_wa1[['item_label', 'item_type', 'item_high_label']].drop_duplicates()
         .sort_values(['item_type', 'item_label']).reset_index(drop=True))
J = len(items)
labels = items.item_label.tolist()
itype = items.item_type.map({'out-of-7': 0., 'categorical': 1.}).to_numpy(np.float32)
ihigh = items.item_high_label.map({'higher_is_better': 1., 'lower_is_better': 0.}).to_numpy(np.float32)
items = items.merge(dit[['item_label', 'cat_length']].drop_duplicates(), on='item_label', how='left')
kmax = items.cat_length.to_numpy(np.float32)

fit = load_fitted_model(f"{BASE}/{file_prefix}_amortised_pps_net.pkl")
nkw = dict(fit['net_kwargs']); nkw.setdefault('num_quantiles', 5)
net = DSNet(**nkw)


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
    """Forward-pass all S draws -> qs (S,J,5), head_in (S,J,H), tgt (S,J)."""
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
    auxb = np.broadcast_to(np.array([n/N_FULL, m/N_FULL], np.float32)[None], (J, 2)).astype(np.float32)
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
        out = np.maximum.accumulate(np.asarray(out), 1)      # monotone quantiles
        qs[s] = out
        if hd is None:
            hd = np.empty((S, J, head.shape[-1]), np.float32)
        hd[s] = np.asarray(head)
    # targets (draw-aligned) from the cached wa frame
    wa = pd.read_pickle(f"{DIR_RGE}/{file_prefix}_i{k}_regression_training.pkl")
    piv = wa.pivot_table(index='draw', columns='item_label', values='pps_ratio_x').reindex(columns=labels)
    tgt = np.full((S, J), np.nan)
    idx = piv.index[piv.index < S].to_numpy()
    tgt[idx] = np.clip(piv.loc[idx].to_numpy(np.float64), -CLIP, CLIP)
    return qs, hd.astype(np.float32), tgt


print(f"Phase 0: forward-pass all draws (S={S}) capturing qs + head_in ...")
QS, HD, TGT = {}, {}, {}
for k in INTERIMS:
    cf = f"{CACHE}/deepset_cell_{k}.npz"
    if os.path.exists(cf):
        z = np.load(cf); QS[k], HD[k], TGT[k] = z['qs'], z['hd'], z['tgt']
    else:
        QS[k], HD[k], TGT[k] = build(k)
        np.savez(cf, qs=QS[k], hd=HD[k], tgt=TGT[k])
    print(f"  interim {k}: qs{QS[k].shape} head{HD[k].shape}")


# ---- numpy diagnostic core (net-agnostic given qs + tgt) ----
def pit(qs, y):
    return np.array([np.interp(y[s], qs[s], TAUS, 0., 1.) for s in range(len(y))])


def marg_cdf(qs, grid):
    F = np.zeros_like(grid)
    for s in range(qs.shape[0]):
        F += np.interp(grid, qs[s], TAUS, 0., 1.)
    return F / qs.shape[0]


def med_at(qs, L):
    return np.array([np.interp(L, TAUS, qs[s]) for s in range(qs.shape[0])])


def ecdf_map(u):
    us = np.sort(u); n = len(us); yy = np.arange(1, n+1)/n
    return lambda x: np.interp(x, us, yy, left=0., right=1.)


def _fin(k, j):
    y = TGT[k][:, j]; ok = np.isfinite(y)
    return ok, y[ok]


# baseline PIT per (k,j) on finite targets
PB = {k: [pit(QS[k][:, j, :], TGT[k][:, j])[np.isfinite(TGT[k][:, j])] for j in range(J)]
      for k in INTERIMS}


def diagnose(pred_qs, get_u_med):
    """pred_qs(k,j)->(S,5) quantiles used for marginal/scatter.
    get_u_med(k,j)->(u, med, g_or_None). Returns row lists."""
    cov, pk, mk, mc, perf, sc = [], [], [], [], [], []
    for k in INTERIMS:
        for j in range(J):
            lbl = labels[j]; ok, y = _fin(k, j)
            qs = pred_qs(k, j)[ok]
            u, med, g = get_u_med(k, j)
            for e in TAUS:
                cov.append(dict(interim_id=k, item_label=lbl, eta_inpol=float(e),
                                coverage=float((u <= e).mean())))
            uu = np.sort(u); ec = np.arange(1, len(u)+1)/len(u)
            pk.append(dict(interim_id=k, item_label=lbl, ks=float(np.max(np.abs(ec-uu))), pit=u))
            lo, hi = float(min(y.min(), qs.min())), float(max(y.max(), qs.max()))
            gr = np.linspace(lo, hi, 200); Fm = marg_cdf(qs, gr)
            if g is not None:
                Fm = g(Fm)
            Fs = (y[:, None] <= gr[None]).mean(0)
            mk.append(dict(interim_id=k, item_label=lbl, ks_marginal=float(np.max(np.abs(Fm-Fs)))))
            idx = np.linspace(0, 199, 60).astype(int)
            for src, F in (('SVI  p(rho|x)', Fs), ('amortiser  p_hat(rho|x)', Fm)):
                for rv, cv in zip(gr[idx], F[idx]):
                    mc.append(dict(interim_id=k, item_label=lbl, source=src, rho=float(rv), cdf=float(cv)))
            perf.append(dict(interim_id=k, item_label=lbl,
                             rho=float(np.corrcoef(med, y)[0, 1]) if np.std(med) > 0 else np.nan))
            si = np.linspace(0, len(y)-1, min(400, len(y))).astype(int)
            for a, b2 in zip(med[si], y[si]):
                sc.append(dict(interim_id=k, item_label=lbl, med=float(a), y=float(b2)))
    return cov, pk, mk, mc, perf, sc


def save(dirn, cov, pk, mk, mc, perf, sc, title):
    d = f"{SB}/{dirn}"; os.makedirs(d, exist_ok=True)
    pd.DataFrame(cov).to_csv(f"{d}/{file_prefix}_pps_{SUF}_tests_coverage.csv", index=False)
    pd.DataFrame([{a: r[a] for a in ('interim_id', 'item_label', 'ks')} for r in pk]
                 ).to_csv(f"{d}/{file_prefix}_pps_{SUF}_tests_pit_ks.csv", index=False)
    pd.DataFrame(mk).to_csv(f"{d}/{file_prefix}_pps_{SUF}_tests_marginal_ks.csv", index=False)
    pd.DataFrame(perf).assign(item_type=lambda x: x.item_label.map(dict(zip(items.item_label, items.item_type)))
                 ).to_csv(f"{d}/{file_prefix}_pps_{SUF}_tests_perf.csv", index=False)
    covdf = pd.DataFrame(cov)
    (ggplot(covdf, aes('eta_inpol', 'coverage', colour='factor(interim_id)', group='interim_id'))
     + geom_abline(intercept=0, slope=1, linetype='dashed', colour='black')
     + geom_line(size=.4, alpha=.8) + geom_point(size=.8) + facet_wrap('~ item_label', ncol=4)
     + theme_bw() + theme(figure_size=(12, 12.5), legend_position='top',
       strip_background=element_blank(), strip_text=element_text(face='bold'))
     + labs(x='nominal', y='empirical coverage', colour='interim',
       title=f'Coverage calibration — {title}')).save(
        f"{d}/{file_prefix}_pps_{SUF}_tests_coverage.pdf", verbose=False, limitsize=False)
    pl = pd.concat([pd.DataFrame({'item_label': r['item_label'], 'u': r['pit']}) for r in pk])
    (ggplot(pl, aes('u')) + geom_histogram(aes(y='..density..'), bins=20, fill='#1f77b4',
       colour='white', size=.2) + geom_hline(yintercept=1, linetype='dashed', colour='black')
     + facet_wrap('~ item_label', ncol=4) + theme_bw()
     + theme(figure_size=(12, 12.5), strip_background=element_blank(), strip_text=element_text(face='bold'))
     + labs(x='PIT u', y='density', title=f'PIT uniformity — {title}')).save(
        f"{d}/{file_prefix}_pps_{SUF}_tests_pit.pdf", verbose=False, limitsize=False)
    mcdf = pd.DataFrame(mc)
    for k in INTERIMS:
        sub = mcdf[mcdf.interim_id == k]
        (ggplot(sub, aes('rho', 'cdf', colour='source', group='source')) + geom_line(size=.6)
         + facet_wrap('~ item_label', ncol=4, scales='free_x') + theme_bw()
         + theme(figure_size=(15, 14), legend_position='top', strip_background=element_blank(),
           strip_text=element_text(face='bold'))
         + labs(x=f'rho at interim {k}', y='CDF', colour='',
           title=f'Marginal p(rho|x) — {title}, interim {k}')).save(
            f"{d}/{file_prefix}_i{k}_svi_vs_amortiser_marginal_cdf.pdf", verbose=False, limitsize=False)
    scdf = pd.DataFrame(sc)
    for k in INTERIMS:
        s = scdf[scdf.interim_id == k]
        (ggplot(s, aes('med', 'y')) + geom_point(alpha=.25, size=.5, colour='#1f77b4')
         + geom_smooth(method='glm', method_args={'family': 'gaussian'}, se=False, colour='black', size=.6)
         + facet_wrap('~ item_label', ncol=4, scales='free') + theme_bw()
         + theme(figure_size=(15, 14), strip_background=element_blank(), strip_text=element_text(face='bold'))
         + labs(x=f'amortiser median rho_hat_3 at interim {k}', y=f'SVI rho at interim {k}',
                title=f'median association — {title}')).save(
            f"{d}/{file_prefix}_i{k}_svi_rho_vs_amortiser_median.pdf", verbose=False, limitsize=False)
    print(f"  saved {title} -> {dirn}")


# ---- baseline ----
print("baseline diagnostics ...")
def base_qs(k, j):
    return QS[k][:, j, :]
def base_u(k, j):
    ok, y = _fin(k, j); qs = QS[k][ok, j, :]
    return pit(qs, y), qs[:, 2], None
save(f"{B}-260805", *diagnose(base_qs, base_u), "deepset baseline (prior)")

# ---- A1: conformal on interim-1 PIT ----
print("A1 conformal-i1 ...")
g1 = {j: ecdf_map(PB[1][j]) for j in range(J)}
lvl1 = {j: float(np.median(PB[1][j])) for j in range(J)}
def a1_u(k, j):
    ok, y = _fin(k, j); qs = QS[k][ok, j, :]
    return g1[j](PB[k][j]), med_at(qs, lvl1[j]), g1[j]
save(f"{B}_conformal-i1_260809", *diagnose(base_qs, a1_u), "deepset conformal (interim-1)")

# ---- A2: conformal expanding 1..k ----
print("A2 conformal-expand ...")
gx, lvx = {}, {}
for k in INTERIMS:
    for j in range(J):
        pool = np.concatenate([PB[i][j] for i in range(1, k+1)])
        gx[(k, j)] = ecdf_map(pool); lvx[(k, j)] = float(np.median(pool))
def a2_u(k, j):
    ok, y = _fin(k, j); qs = QS[k][ok, j, :]
    return gx[(k, j)](PB[k][j]), med_at(qs, lvx[(k, j)]), gx[(k, j)]
save(f"{B}_conformal-expand_260809", *diagnose(base_qs, a2_u), "deepset conformal (expanding 1..k)")

# ---- H: head-only fine-tune expanding 1..k on cached head_in ----
print("H head-ft-expanding (q_psi on cached embeddings) ...")
taus_j = jnp.asarray(TAUS)
qpsi = _MLP(dims=(*nkw.get('hidden_dims', (64, 64)), 5))
qpsi0 = {'params': fit['params']['params']['q_psi']}


def pinball(params, hd, y):
    pr = jnp.maximum.accumulate(qpsi.apply(params, hd), 1)
    e = y[:, None] - pr
    return jnp.mean(jnp.maximum(taus_j[None]*e, (taus_j[None]-1)*e))


HK = {}
for k in INTERIMS:
    HDk = np.concatenate([HD[i].reshape(-1, HD[i].shape[-1]) for i in range(1, k+1)])
    TGk = np.concatenate([TGT[i].reshape(-1) for i in range(1, k+1)])
    ok = np.isfinite(TGk); HDk, TGk = HDk[ok].astype(np.float32), TGk[ok].astype(np.float32)
    params = qpsi0; opt = optax.adam(3e-4); ost = opt.init(params)
    rng = np.random.default_rng(k)

    @jax.jit
    def step(params, ost, h, y):
        l, g = jax.value_and_grad(pinball)(params, h, y)
        up, ost = opt.update(g, ost, params)
        return optax.apply_updates(params, up), ost, l
    bs = min(2048, HDk.shape[0]); l = 0.
    for it in range(800):
        ix = rng.integers(0, HDk.shape[0], bs)
        params, ost, l = step(params, ost, jnp.asarray(HDk[ix]), jnp.asarray(TGk[ix]))
    pr = np.maximum.accumulate(np.asarray(qpsi.apply(params, jnp.asarray(HD[k].reshape(-1, HD[k].shape[-1])))), 1)
    HK[k] = pr.reshape(S, J, 5)
    print(f"  interim {k}: head-ft on 1..{k} done (loss {float(l):.4f})")


def h_qs(k, j):
    return HK[k][:, j, :]
def h_u(k, j):
    ok, y = _fin(k, j); qs = HK[k][ok, j, :]
    return pit(qs, y), qs[:, 2], None
save(f"{B}_ftheadexpand_260809", *diagnose(h_qs, h_u), "deepset head-ft expanding 1..k")
print("Deepset diagnostics complete.")
