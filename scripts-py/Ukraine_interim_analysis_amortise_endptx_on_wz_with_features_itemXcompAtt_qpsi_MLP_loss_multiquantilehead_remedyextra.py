#!/usr/bin/env python3
"""
Extra calibration-remedy experiments on the prior-predictive-trained
`itemXcompAtt` amortiser (§14.2 of dev/amortised_decision_making.md):

  C1  conformal recalibration using ONLY interim-1 SVI, applied to all k
  C2  conformal recalibration EXPANDING: interim k recalibrated on 1..k
  H   head-only fine-tune EXPANDING: interim k = head fit on 1..k SVI

All three emit the SAME diagnostic outputs as the main tests script
(coverage curves, PIT histograms, marginal-CDF overlays + CSVs), one
result dir each, so collate_remedy.py picks them up.

Conformal (Kuleshov et al. 2018, distribution calibration): the PIT
u = F_hat(rho_SVI-x | x, z^s) carries everything (coverage@eta =
ecdf_u(eta)); a per-item monotone map g = ecdf(u_cal) recalibrates to
u' = g(u); recalibrated coverage/PIT/marginal follow from u' and the
per-s quantile CDFs.
"""
import os, sys, importlib
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'python'))
import numpy as np, pandas as pd, jax, jax.numpy as jnp, optax
from plotnine import (ggplot, aes, geom_abline, geom_line, geom_point, geom_smooth, geom_text,
                      geom_histogram, geom_hline, facet_wrap, theme_bw,
                      theme, element_blank, element_text, labs)
from amortiser_common import load_fitted_model, _jax_pytree
from model_pcm import PartialCreditModel

DIR_RGE = "/Users/or105/sandbox/bIRTistic/py-ukraine-interim-with-regression-on-endptx-wz-260601"
SB = "/Users/or105/sandbox/bIRTistic"
# B = xcomp default; OUTB names the output dirs, RX_BASE the net to load.
# Both parameterised so this script serves itemXcompAtt AND itemScompAtt
# (same F=8 token schema; only the trained net differs).
B = "py-ukraine-interim-amortise-endptx-on-wz-with-features-itemXcompAtt-qpsi-MLP-loss-multiquantilehead"
OUTB = os.environ.get('RX_OUTB', B)
BASE = os.environ.get('RX_BASE', f"{SB}/{B}-260716")
DO = set(os.environ.get('RX_DO', 'C1,C2,H').split(','))
# net-specific pred cache tag (tokens are net-independent, preds are not)
CTAG = os.environ.get('RX_CTAG', 'xcomp')
file_prefix = "pcm_1_interim"
N_FULL, KHAT_N0, F_TOK, CLIP = 503, 50.0, 8, 20.0
TAUS = np.array([0.05, 0.25, 0.5, 0.75, 0.95], np.float32)
INTERIMS = list(range(1, 9))
SUF = "RGEG"

# ---- metadata + net ----
iw = {i: pd.read_pickle(f"{DIR_RGE}/{file_prefix}_i{i}_regression_training.pkl")
      for i in INTERIMS}
items = (iw[1][['item_label', 'item_type', 'item_high_label']].drop_duplicates()
         .sort_values(['item_type', 'item_label']).reset_index(drop=True))
items['item_pos'] = np.arange(len(items)); J = len(items)
ipos = dict(zip(items.item_label, items.item_pos))
itype = items.item_type.map({'out-of-7': 0., 'categorical': 1.}).to_numpy(np.float32)
ihigh = items.item_high_label.map({'higher_is_better': 1., 'lower_is_better': 0.}).to_numpy(np.float32)
dit = pd.read_csv(f"{DIR_RGE}/{file_prefix}_1_data_dit.csv")
items = items.merge(dit[['item_label', 'cat_length']].drop_duplicates(), on='item_label', how='left')
klev = items.cat_length.to_numpy(np.float32)
labels = items.item_label.tolist()
fit = load_fitted_model(f"{BASE}/{file_prefix}_amortised_pps_net.pkl")
istd = np.load(f"{BASE}/{file_prefix}_item_std.npy").astype(np.float32)
netmod = importlib.import_module(fit['net_class_module'])
Net = getattr(netmod, fit['net_class_name'])
nkw = dict(fit['net_kwargs']); nkw.setdefault('num_quantiles', 5)
net = Net(**nkw)


def _xfeat(k):
    xi = pd.read_csv(f"{DIR_RGE}/{file_prefix}_{k}_data_dp1.csv")
    piv = xi.pivot_table(index='pid', columns=['item_label', 'time'], values='y')
    n = piv.index.size; wx = np.zeros((J, 2), np.float32); chg = np.zeros((n, J))
    for j, l in enumerate(labels):
        y0 = piv[(l, 0)].to_numpy(float); y1 = piv[(l, 1)].to_numpy(float)
        km = klev[j] - 1; wx[j, 0] = np.nanmean(y0) / km; wx[j, 1] = np.nanmean(y1) / km
        chg[:, j] = y1 - y0
    r = np.nan_to_num(pd.DataFrame(chg).corr('spearman').to_numpy()); np.fill_diagonal(r, 1)
    lam = n / (n + KHAT_N0)
    return wx, (lam * r + (1 - lam) * np.eye(J)).astype(np.float32), n


def _zmean(k, D):
    xi = pd.read_csv(f"{DIR_RGE}/{file_prefix}_{k}_data_dp1.csv")
    n = xi.pid.nunique()
    m = PartialCreditModel(dit=dit, dcati=xi, x_formula="~ time - 1", seed=123)
    zi = m.get_interim_z_from_ypredi(f"{DIR_RGE}/{file_prefix}_{k}_draws.zarr",
                                     N_FULL - n, pps_z_total=D, seed=123, keep_order=True)
    cols = [f'ypred_{s}' for s in range(D)]; out = np.zeros((D, J, 2), np.float32)
    for j, l in enumerate(labels):
        km = klev[j] - 1
        for tt in (0, 1):
            v = zi[(zi.item_label == l) & (zi.time == tt)][cols].to_numpy(np.float32)
            out[:, j, tt] = (v - 1).mean(0) / km
    return out


def build(k):
    wx, kh, n = _xfeat(k)
    wa = iw[k].copy(); wa['item_pos'] = wa.item_label.map(ipos)
    piv = wa.pivot_table(index='draw', columns='item_pos', values='pps_ratio_x')
    D = piv.index.size
    tgt = np.clip(np.nan_to_num(piv.reindex(columns=range(J)).to_numpy(np.float32),
                  posinf=CLIP, neginf=-CLIP), -CLIP, CLIP)
    zm = _zmean(k, D); wb = np.maximum(zm[:, :, 0], 1e-3); we = zm[:, :, 1]
    hi = (ihigh > .5)[None]; wr = np.clip(np.where(hi, we/wb-1, 1-we/wb), -CLIP, CLIP)
    feats = np.nan_to_num(np.stack([np.broadcast_to(wx[None, :, 0], (D, J)),
            np.broadcast_to(wx[None, :, 1], (D, J)), zm[:, :, 0], zm[:, :, 1], wr], -1))
    meta = np.broadcast_to(np.stack([itype, ihigh], -1)[None], (D, J, 2))
    tok = np.concatenate([feats, meta], -1).astype(np.float32)
    tf = np.empty((D, J, J, F_TOK), np.float32); tf[..., :7] = tok[:, None]; tf[..., 7] = kh[None]
    tf = tf.reshape(D*J, J, F_TOK)
    qid = np.broadcast_to(np.arange(J)[None], (D, J)).reshape(D*J).astype(np.int32)
    aux = np.broadcast_to(np.array([n/N_FULL, (N_FULL-n)/N_FULL], np.float32)[None], (D*J, 2)).astype(np.float32)
    return dict(D=D, tok=tf, aux=aux, qid=qid, tgt=tgt)


def preds_of(params, cell):
    p = np.asarray(net.apply(params, {'tokens': jnp.asarray(cell['tok']),
        'mask': jnp.ones((cell['tok'].shape[0], J)), 'query_idx': jnp.asarray(cell['qid']),
        'aux': jnp.asarray(cell['aux'])}))
    p = p * istd[cell['qid']][:, None]
    p = np.maximum.accumulate(p, 1)
    return p.reshape(cell['D'], J, 5)


print("Phase 0: build cells + baseline preds (cached) ...")
CACHE = "/private/tmp/claude-501/-Users-or105-git-bIRTistic/49372e22-d11b-4882-a442-d1c60bcbdfb0/scratchpad/remedy_cells"
os.makedirs(CACHE, exist_ok=True)
CELL = {}
for k in INTERIMS:
    cf = f"{CACHE}/cell_{k}.npz"                 # tokens: net-independent
    if os.path.exists(cf):
        z = np.load(cf)
        CELL[k] = dict(D=int(z['D']), tok=z['tok'], aux=z['aux'],
                       qid=z['qid'], tgt=z['tgt'])
    else:
        c = build(k)
        np.savez(cf, D=c['D'], tok=c['tok'], aux=c['aux'],
                 qid=c['qid'], tgt=c['tgt'])
        CELL[k] = c
    # preds ARE net-specific -> always recompute for the loaded net.
    CELL[k]['qs'] = preds_of(fit['params'], CELL[k])
    print(f"  interim {k} (D={CELL[k]['D']}, net={CTAG})")



def pit(qs, y):  # (D,Q),(D,) -> u (D,)
    return np.array([np.interp(y[s], qs[s], TAUS, 0., 1.) for s in range(len(y))])


def marg_cdf(qs, grid):  # amortiser marginal CDF on grid
    F = np.zeros_like(grid)
    for s in range(qs.shape[0]):
        F += np.interp(grid, qs[s], TAUS, 0., 1.)
    return F / qs.shape[0]
# Precompute baseline PIT per (k,j) + calibration levels (hoist the
# per-draw recompute that made C1 pathological).
PB = {k: [pit(CELL[k]['qs'][:, j, :], CELL[k]['tgt'][:, j]) for j in range(J)]
      for k in INTERIMS}
def med_at(qs, L):
    return np.array([np.interp(L, TAUS, qs[s]) for s in range(qs.shape[0])])


def diag_rows(get_u_qs_y):
    """get_u_qs_y(k,j)-> (u_recal, qs, y, g_or_None). Returns cov/pit/marg
    row lists over all (k,j)."""
    cov, pk, mk, mc, perf, sc = [], [], [], [], [], []
    for k in INTERIMS:
        for j in range(J):
            lbl = labels[j]; y = CELL[k]['tgt'][:, j]; qs = CELL[k]['qs'][:, j, :]
            u, qs_use, g = get_u_qs_y(k, j)
            for e in TAUS:
                cov.append(dict(interim_id=k, item_label=lbl, eta_inpol=float(e),
                                coverage=float((u <= e).mean())))
            uu = np.sort(u); ec = np.arange(1, len(u)+1)/len(u)
            pk.append(dict(interim_id=k, item_label=lbl, ks=float(np.max(np.abs(ec-uu))), pit=u))
            lo, hi = float(min(y.min(), qs.min())), float(max(y.max(), qs.max()))
            gr = np.linspace(lo, hi, 200)
            Fm = marg_cdf(qs, gr)
            if g is not None:      # conformal: recalibrate the marginal CDF too
                Fm = g(Fm)
            Fs = (y[:, None] <= gr[None]).mean(0)
            mk.append(dict(interim_id=k, item_label=lbl,
                           ks_marginal=float(np.max(np.abs(Fm-Fs)))))
            idx = np.linspace(0, 199, 60).astype(int)
            for src, F in (('SVI  p(rho|x)', Fs), ('amortiser  p_hat(rho|x)', Fm)):
                for rv, cv in zip(gr[idx], F[idx]):
                    mc.append(dict(interim_id=k, item_label=lbl, source=src,
                                   rho=float(rv), cdf=float(cv)))
            # correlation of recalibrated median vs target
            med = qs_use
            perf.append(dict(interim_id=k, item_label=lbl,
                             rho=float(np.corrcoef(med, y)[0, 1]) if np.std(med) > 0 else np.nan))
            si = np.linspace(0, len(y)-1, 400).astype(int)
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
                 ).to_csv(f"{d}/{file_prefix}_pps_{SUF}_perf.csv", index=False)
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


def ecdf_map(u_cal):
    us = np.sort(u_cal); n = len(us); yy = (np.arange(1, n+1))/n
    return lambda x: np.interp(x, us, yy, left=0., right=1.)


# ---- C1: conformal on interim-1 PIT, applied to all ----
if 'C1' in DO:
    print("C1 conformal-i1 ...")
    g1 = {j: ecdf_map(PB[1][j]) for j in range(J)}
    lvl1 = {j: float(np.median(PB[1][j])) for j in range(J)}
    def c1(k, j):
        qs = CELL[k]['qs'][:, j, :]
        return g1[j](PB[k][j]), med_at(qs, lvl1[j]), g1[j]
    save(f"{OUTB}_conformal-i1_260809", *diag_rows(c1), "conformal (interim-1)")

# ---- C2: conformal expanding 1..k ----
if 'C2' in DO:
    print("C2 conformal-expand ...")
    gexp, lvlexp = {}, {}
    for k in INTERIMS:
        for j in range(J):
            pool = np.concatenate([PB[i][j] for i in range(1, k+1)])
            gexp[(k, j)] = ecdf_map(pool); lvlexp[(k, j)] = float(np.median(pool))
    def c2(k, j):
        qs = CELL[k]['qs'][:, j, :]
        return gexp[(k, j)](PB[k][j]), med_at(qs, lvlexp[(k, j)]), gexp[(k, j)]
    save(f"{OUTB}_conformal-expand_260809", *diag_rows(c2), "conformal (expanding 1..k)")

# ---- H: head-only fine-tune expanding 1..k ----
if 'H' not in DO:
    print("Remedy-extra complete (H skipped).")
    sys.exit(0)
print("H head-ft-expanding ...")
taus_j = jnp.asarray(TAUS)
def pinball(pr, y):
    pr = jnp.maximum.accumulate(pr, 1); e = y[:, None] - pr
    return jnp.mean(jnp.maximum(taus_j[None]*e, (taus_j[None]-1)*e))
def lab(params):
    # Universal head-only rule: train ONLY q_psi, freeze every encoder
    # block (xcomp: q_tok,q_query; scomp: q_tok,query_bias_alpha).
    return jax.tree_util.tree_map_with_path(
        lambda p, _: 'train' if (len(p) > 1 and p[1].key == 'q_psi') else 'frozen', params)
HK = {}
for k in INTERIMS:
    TOK = np.concatenate([CELL[i]['tok'] for i in range(1, k+1)])
    AUX = np.concatenate([CELL[i]['aux'] for i in range(1, k+1)])
    QID = np.concatenate([CELL[i]['qid'] for i in range(1, k+1)])
    TGT = np.concatenate([(CELL[i]['tgt']/istd[None]).reshape(-1) for i in range(1, k+1)]).astype(np.float32)
    tx = optax.multi_transform({'train': optax.adam(3e-4), 'frozen': optax.set_to_zero()}, lab(fit['params']))
    params = fit['params']; ost = tx.init(params); rng = np.random.default_rng(k)
    @jax.jit
    def step(params, ost, b):
        l, g = jax.value_and_grad(lambda pp: pinball(net.apply(pp, b), b['target']))(params)
        up, ost = tx.update(g, ost, params); return optax.apply_updates(params, up), ost, l
    for s in range(800):
        ix = rng.integers(0, TOK.shape[0], 2048)
        b = {'tokens': jnp.asarray(TOK[ix]), 'mask': jnp.ones((2048, J)),
             'query_idx': jnp.asarray(QID[ix]), 'aux': jnp.asarray(AUX[ix]), 'target': jnp.asarray(TGT[ix])}
        params, ost, l = step(params, ost, b)
    HK[k] = preds_of(params, CELL[k])   # eval interim k with head-fit-on-1..k
    print(f"  interim {k}: head-ft on 1..{k} done (loss {float(l):.4f})")
def hexp(k, j):
    qs = HK[k][:, j, :]; y = CELL[k]['tgt'][:, j]
    return pit(qs, y), qs[:, 2], None
# override diag to use HK preds for qs/marginal
def diag_H():
    cov, pk, mk, mc, perf, sc = [], [], [], [], [], []
    for k in INTERIMS:
        for j in range(J):
            lbl = labels[j]; y = CELL[k]['tgt'][:, j]; qs = HK[k][:, j, :]
            u = pit(qs, y)
            for e in TAUS: cov.append(dict(interim_id=k, item_label=lbl, eta_inpol=float(e), coverage=float((u <= e).mean())))
            uu = np.sort(u); ec = np.arange(1, len(u)+1)/len(u)
            pk.append(dict(interim_id=k, item_label=lbl, ks=float(np.max(np.abs(ec-uu))), pit=u))
            lo, hi = float(min(y.min(), qs.min())), float(max(y.max(), qs.max())); gr = np.linspace(lo, hi, 200)
            Fm = marg_cdf(qs, gr); Fs = (y[:, None] <= gr[None]).mean(0)
            mk.append(dict(interim_id=k, item_label=lbl, ks_marginal=float(np.max(np.abs(Fm-Fs)))))
            idx = np.linspace(0, 199, 60).astype(int)
            for src, F in (('SVI  p(rho|x)', Fs), ('amortiser  p_hat(rho|x)', Fm)):
                for rv, cv in zip(gr[idx], F[idx]): mc.append(dict(interim_id=k, item_label=lbl, source=src, rho=float(rv), cdf=float(cv)))
            perf.append(dict(interim_id=k, item_label=lbl, rho=float(np.corrcoef(qs[:, 2], y)[0, 1]) if np.std(qs[:, 2]) > 0 else np.nan))
            si = np.linspace(0, len(y)-1, 400).astype(int)
            for a, b2 in zip(qs[si, 2], y[si]):
                sc.append(dict(interim_id=k, item_label=lbl, med=float(a), y=float(b2)))
    return cov, pk, mk, mc, perf, sc
save(f"{OUTB}_ftheadexpand_260809", *diag_H(), "head-ft expanding 1..k")
print("Remedy-extra complete.")
