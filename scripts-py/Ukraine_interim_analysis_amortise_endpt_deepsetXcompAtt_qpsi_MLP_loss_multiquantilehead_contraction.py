#!/usr/bin/env python3
"""
Posterior-contraction diagnostic for the deepset PPS amortiser
(§14.4.9 of dev/amortised_decision_making.md).

Question (memory note `deepset-contraction-calibration`): does the
amortiser's predicted spread for rho shrink with the observed cohort
size n at the SAME rate as the true SVI posterior p(rho | x)? If yes,
the interim-heteroscedastic calibration bias is a fixed function of n
that the encoder already represents, and a SINGLE calibration point
(interim-1 head fine-tune) could pin all interims.

Uses the WEEKLY-cadence SVI interims (Ukraine_interim_weekly_svi_for_contraction.py)
for many (n, width) points. Per (weekly interim k, item j) we compute:

  n_obs                       observed cohort size
  W_svi   = q95 - q05 of the SVI draws rho_SVI-x           (true p(rho|x) width)
  W_marg  = 5..95 width of the amortiser MARGINAL p(rho|x) (mixture over draws
            of the per-draw conditionals p(rho|x,z^s))
  W_cond  = mean_s (q95^s - q05^s)                          (amortiser per-slice width)

Log-log slope of W vs n estimates the contraction rate (Bernstein-von
Mises: ~ -1/2). Matching slopes (W_marg vs W_svi) => encoder learned
contraction. Outputs into the deepset baseline dir.
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
from plotnine import (ggplot, aes, geom_point, geom_line, geom_smooth,
                      facet_wrap, facet_grid, scale_x_log10, scale_y_log10, theme_bw,
                      theme, element_blank, element_text, labs,
                      scale_color_cmap, scale_linetype_manual, scale_color_manual)

warnings.filterwarnings('ignore')
from amortiser_common import load_fitted_model
from amortiser_pps_features_deepsetXcompAtt_qpsi_MLP_loss_multiquantilehead import (
    Amortiser_PPS_features_deepsetXcompAtt_qpsi_MLP_loss_multiquantilehead as DSNet,
)
from model_pcm import PartialCreditModel

DIR_RGE = os.environ.get(
    'UKR_WEEKLY_DIR',
    "/Users/or105/sandbox/bIRTistic/py-ukraine-interim-weekly-svi-260811")
SB = "/Users/or105/sandbox/bIRTistic"
B = "py-ukraine-interim-amortise-endptx-on-wz-with-features-deepsetXcompAtt-qpsi-MLP-loss-multiquantilehead"
BASE = os.environ.get('RGDX_BASE', f"{SB}/{B}-260805")
AUX_SQRT = os.environ.get('AUX_SQRT', '0') == '1'   # dim-4 aux for contraction net
CTAG = os.environ.get('RGDX_CTAG', 'base')
file_prefix = "pcm_1_interim"
N_FULL, CLIP = 503, 20.0
S = int(os.environ.get('UKR_DEEPSET_S', 200))
TAUS = np.array([0.05, 0.25, 0.5, 0.75, 0.95], np.float32)
CACHE = f"/private/tmp/claude-501/-Users-or105-git-bIRTistic/49372e22-d11b-4882-a442-d1c60bcbdfb0/scratchpad/weekly_contraction_cells_{CTAG}"
os.makedirs(CACHE, exist_ok=True)

# item metadata (shared ordering with deploy)
dit = pd.read_csv(f"{DIR_RGE}/{file_prefix}_1_data_dit.csv")
_wa_glob = sorted(glob.glob(f"{DIR_RGE}/{file_prefix}_i*_regression_training.pkl"),
                  key=lambda p: int(p.split('_i')[-1].split('_')[0]))
INTERIMS = [int(p.split('_i')[-1].split('_')[0]) for p in _wa_glob]
print(f"weekly interims found: {INTERIMS}")
_wa1 = pd.read_pickle(_wa_glob[0])
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
    return net.apply(params, batch)


def build(k):
    """Forward-pass S draws -> qs (S,J,5); return qs, tgt (S,J), n_obs."""
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
        av += [1.0/np.sqrt(max(n, 1)), 1.0/np.sqrt(max(m, 1))]   # BvM features
    auxb = np.broadcast_to(np.array(av, np.float32)[None], (J, len(av))).astype(np.float32)
    qidb = np.arange(J, dtype=np.int32)
    qs = np.empty((S, J, 5), np.float64)
    for s in range(S):
        z_raw, _ = _pivot(zi.rename(columns={f'ypred_{s}': '_y'}), '_y')
        z_pad = np.zeros((N_FULL, J, 2), np.float32); z_pad[:m] = z_raw
        zb = np.broadcast_to(z_pad[None], (J, N_FULL, J, 2)).astype(np.float32)
        batch = {'x_responses': jnp.asarray(xb), 'mask_x': jnp.asarray(mxb),
                 'z_responses': jnp.asarray(zb), 'mask_z': jnp.asarray(mzb),
                 'item_metadata': jnp.asarray(metab), 'query_idx': jnp.asarray(qidb),
                 'aux': jnp.asarray(auxb)}
        qs[s] = np.maximum.accumulate(np.asarray(_fwd(fit['params'], batch)), 1)
    wa = pd.read_pickle(f"{DIR_RGE}/{file_prefix}_i{k}_regression_training.pkl")
    piv = wa.pivot_table(index='draw', columns='item_label', values='pps_ratio_x').reindex(columns=labels)
    tgt = np.full((S, J), np.nan)
    idx = piv.index[piv.index < S].to_numpy()
    tgt[idx] = np.clip(piv.loc[idx].to_numpy(np.float64), -CLIP, CLIP)
    return qs, tgt, n


def marg_width(qs):
    """5-95 width of the amortiser marginal p(rho|x) (mixture over draws)."""
    lo = float(qs[:, 0].min()); hi = float(qs[:, 4].max())
    if not np.isfinite(lo) or hi <= lo:
        return np.nan
    gr = np.linspace(lo, hi, 400)
    F = np.zeros_like(gr)
    for s in range(qs.shape[0]):
        F += np.interp(gr, qs[s], TAUS, 0., 1.)
    F /= qs.shape[0]
    r05 = np.interp(0.05, F, gr); r95 = np.interp(0.95, F, gr)
    return float(r95 - r05)


print(f"Phase 0: forward-pass weekly interims (S={S}) ...")
rows = []
QS, TGT = {}, {}
for k in INTERIMS:
    cf = f"{CACHE}/wk_{k}.npz"
    if os.path.exists(cf):
        z = np.load(cf); qs, tgt, n = z['qs'], z['tgt'], int(z['n'])
    else:
        qs, tgt, n = build(k)
        np.savez(cf, qs=qs, tgt=tgt, n=n)
    QS[k], TGT[k] = qs, tgt
    for j in range(J):
        y = tgt[:, j]; ok = np.isfinite(y); yv = y[ok]
        if yv.size < 10:
            continue
        w_svi = float(np.quantile(yv, 0.95) - np.quantile(yv, 0.05))
        w_cond = float(np.mean(qs[:, j, 4] - qs[:, j, 0]))
        w_marg = marg_width(qs[:, j, :])
        rows.append(dict(interim_id=k, n=n, item_label=labels[j],
                         item_type=items.item_type.iloc[j],
                         W_svi=w_svi, W_marg=w_marg, W_cond=w_cond))
    print(f"  week {k} (n={n}) done")

df = pd.DataFrame(rows)
BLOW = "CG-VIO_ph-punish"
out = f"{BASE}"
df.to_csv(f"{out}/{file_prefix}_pps_RGDX_contraction.csv", index=False)

# ---- log-log contraction slope per (item, source) ----
def slope(sub, col):
    d = sub[(sub[col] > 0) & (sub.n > 0)]
    if len(d) < 3:
        return np.nan
    return float(np.polyfit(np.log(d.n), np.log(d[col]), 1)[0])


slopes = []
for j, lbl in enumerate(labels):
    sub = df[df.item_label == lbl]
    slopes.append(dict(item_label=lbl, item_type=items.item_type.iloc[j],
                       slope_svi=slope(sub, 'W_svi'),
                       slope_marg=slope(sub, 'W_marg'),
                       slope_cond=slope(sub, 'W_cond')))
sl = pd.DataFrame(slopes)
sl.to_csv(f"{out}/{file_prefix}_pps_RGDX_contraction_slopes.csv", index=False)
ok = sl[sl.item_label != BLOW]
print("\nmean log-log slope (contraction rate; BvM ~ -0.5):")
print(ok[['slope_svi', 'slope_marg', 'slope_cond']].mean().round(3).to_string())

# ---- plot: width vs n, faceted by item, 3 sources ----
long = df.melt(id_vars=['interim_id', 'n', 'item_label', 'item_type'],
               value_vars=['W_svi', 'W_marg', 'W_cond'],
               var_name='source', value_name='width')
long['source'] = long['source'].map({
    'W_svi': 'SVI p(rho|x)', 'W_marg': 'amortiser marginal p(rho|x)',
    'W_cond': 'amortiser per-draw p(rho|x,z^s)'})
long = long[(long.width > 0) & long.width.notna()]
(ggplot(long, aes('n', 'width', colour='source', group='source'))
 + geom_point(size=.8, alpha=.7) + geom_line(size=.4, alpha=.7)
 + facet_wrap('~ item_label', ncol=4, scales='free_y')
 + scale_x_log10() + scale_y_log10() + theme_bw()
 + theme(figure_size=(15, 15), legend_position='top', strip_background=element_blank(),
         strip_text=element_text(face='bold'))
 + labs(x='observed cohort size n (log)', y='90% interval width (log)', colour='',
        title='Posterior contraction: interval width vs n (weekly interims)')).save(
    f"{out}/{file_prefix}_pps_RGDX_contraction_by_item.pdf", verbose=False, limitsize=False)

# pooled panel: mean width across items (excl blow-up) vs n
pooled = (long[long.item_label != BLOW].groupby(['n', 'source'])['width']
          .median().reset_index())
(ggplot(pooled, aes('n', 'width', colour='source', group='source'))
 + geom_point(size=1.5) + geom_line(size=.6)
 + geom_smooth(method='lm', se=False, linetype='dashed', size=.4)
 + scale_x_log10() + scale_y_log10() + theme_bw()
 + theme(figure_size=(8, 6), legend_position='top')
 + labs(x='observed cohort size n (log)', y='median 90% interval width (log)', colour='',
        title='Posterior contraction (median over items) vs n')).save(
    f"{out}/{file_prefix}_pps_RGDX_contraction_pooled.pdf", verbose=False, limitsize=False)

# ---- contraction-CDF grid: items in rows, interims in columns (8 evenly-
# spaced from the weekly grid), SVI solid vs amortiser dotted ----
SEL = [INTERIMS[i] for i in np.linspace(0, len(INTERIMS) - 1, 8).astype(int)]
ccrows = []
for j in range(J):
    if labels[j] == BLOW:
        continue
    vals = []
    for k in SEL:
        y = TGT[k][:, j]; ok = np.isfinite(y)
        if ok.sum() >= 10:
            vals.append(y[ok]); vals.append(QS[k][ok, j, :].ravel())
    if not vals:
        continue
    lo, hi = np.percentile(np.concatenate(vals), [1, 99])
    if hi <= lo:
        continue
    gr = np.linspace(lo, hi, 120); grs = (gr - lo) / (hi - lo)
    for k in SEL:
        y = TGT[k][:, j]; ok = np.isfinite(y); yv = y[ok]
        if yv.size < 10:
            continue
        qs = QS[k][ok, j, :]
        Fs = (yv[:, None] <= gr[None]).mean(0)
        Fm = np.zeros_like(gr)
        for s in range(qs.shape[0]):
            Fm += np.interp(gr, qs[s], TAUS, 0., 1.)
        Fm /= qs.shape[0]
        for xs, fs, fm in zip(grs, Fs, Fm):
            ccrows.append(dict(interim_id=k, item_label=labels[j], x=float(xs), cdf=float(fs), source='SVI'))
            ccrows.append(dict(interim_id=k, item_label=labels[j], x=float(xs), cdf=float(fm), source='amortiser'))
ccdf = pd.DataFrame(ccrows)
(ggplot(ccdf, aes('x', 'cdf', colour='source', linetype='source', group='source'))
 + geom_line(size=.55) + facet_grid('item_label ~ interim_id')
 + scale_color_manual(values={'SVI': '#1f77b4', 'amortiser': '#d62728'})
 + scale_linetype_manual(values={'SVI': 'solid', 'amortiser': 'dotted'})
 + theme_bw() + theme(figure_size=(13, 24), legend_position='top', panel_spacing=0.015,
   strip_background=element_blank(), strip_text_x=element_text(face='bold', size=8),
   strip_text_y=element_text(angle=0, ha='left', size=6))
 + labs(x='rho (per-item min-max scaled)', y='CDF',
   title='Posterior contraction — contraction-pressure retrain  (rows=item, cols=weekly interim)')).save(
    f"{out}/{file_prefix}_pps_RGDX_contraction_cdf_by_item.pdf", verbose=False, limitsize=False)
print(f"\nContraction diagnostic complete -> {out}")
