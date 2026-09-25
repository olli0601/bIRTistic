"""
Shared amortiser calibration/contraction diagnostic plots (§14.4.22-25), so the
deepset (`deepsetXcompAtt`) and hand-token (`itemXcompAtt`) deploys produce the
*same* figures for a like-for-like comparison.

Each function takes callables over interim ids:
  qs_of(k)  -> (D, J, 5)  amortiser quantile predictions (the DEPLOYABLE ones)
  tgt_of(k) -> (D, J)     SVI draws of rho (the target marginal)
  n_of(k)   -> int        observed-cohort size at interim k
plus `labels` (J item names), `interims` (sorted ids), `taus` (the 5 quantile
levels), `out` (dir), `prefix` (file prefix), `suf` (title/file tag), and
`blow` (a degenerate item label to drop from SD/contraction fits).
"""

import numpy as np
import pandas as pd
from plotnine import (ggplot, aes, geom_boxplot, geom_hline, geom_point, geom_line,
                      geom_text, position_dodge, facet_wrap, theme_bw, theme, labs,
                      element_blank, element_text, scale_fill_manual, scale_fill_gradient2)

CBLUE, CRED = '#1f77b4', '#d62728'


def _marg_cdf(qs, gr, taus):
    F = np.zeros_like(gr)
    for s in range(qs.shape[0]):
        F += np.interp(gr, qs[s], taus, 0., 1.)
    return F / qs.shape[0]


def contraction_cdf(qs_of, tgt_of, labels, interims, taus, out, prefix, suf):
    """Marginal p(rho|x) quantile-boxes: x=interim, per item, SVI (blue) vs
    amortiser (red) dodged; whiskers 10/90, box 25/75, midline median (all
    precomputed). SVI from the draws; amortiser from the inverted mixture CDF."""
    J = len(labels); QL = np.array([0.10, 0.25, 0.50, 0.75, 0.90]); rows = []
    for j in range(J):
        for k in interims:
            y = tgt_of(k)[:, j]; ok = np.isfinite(y); yv = y[ok]
            if yv.size < 10:
                continue
            qs = qs_of(k)[:, j, :]; qs = qs[np.isfinite(qs).all(1)]
            sq = np.percentile(yv, QL * 100.0)
            lo, hi = float(min(yv.min(), qs.min())), float(max(yv.max(), qs.max()))
            gr = np.linspace(lo, hi, 400); aq = np.interp(QL, _marg_cdf(qs, gr, taus), gr)
            for src, q in (('SVI', sq), ('amortiser', aq)):
                rows.append(dict(interim_id=int(k), item_label=labels[j], source=src,
                                 ymin=float(q[0]), lower=float(q[1]), middle=float(q[2]),
                                 upper=float(q[3]), ymax=float(q[4])))
    df = pd.DataFrame(rows)
    df.to_csv(f"{out}/{prefix}_pps_{suf}_contraction_cdf_by_item.csv", index=False)  # box data (SVI vs amortiser marginal)
    df['interim_id'] = pd.Categorical(df.interim_id, categories=sorted(df.interim_id.unique()), ordered=True)
    (ggplot(df, aes('interim_id', ymin='ymin', lower='lower', middle='middle',
                    upper='upper', ymax='ymax', fill='source'))
     + geom_boxplot(stat='identity', position=position_dodge(width=0.78), size=.3, width=.7, alpha=.85)
     + facet_wrap('~ item_label', ncol=4, scales='free_y')
     + scale_fill_manual(values={'SVI': CBLUE, 'amortiser': CRED})
     + theme_bw() + theme(figure_size=(16, 20), legend_position='top', panel_spacing=0.02,
       strip_background=element_blank(), strip_text=element_text(face='bold', size=8),
       axis_text_x=element_text(size=6))
     + labs(x='interim id', y='p(rho | x) estimate', fill='',
       title=f'Marginal p(rho|x): SVI vs amortiser — {suf}  '
             f'(box 25-75%, midline median, whiskers 10-90%)')).save(
        f"{out}/{prefix}_pps_{suf}_contraction_cdf_by_item.pdf", verbose=False, limitsize=False)


def pit_box(qs_of, tgt_of, labels, interims, taus, out, prefix, suf):
    """Conditional calibration (PIT) quantile-boxes: x=interim, per item, box of
    u=F_amortiser(rho_SVI|x,z); dashed refs .1/.25/.5/.75/.9 = calibrated;
    diverging fill on median-0.5 (blue: median too high, red: too low); box wider
    than .25-.75 = over-confident."""
    J = len(labels); QL = np.array([0.10, 0.25, 0.50, 0.75, 0.90]); rows = []
    for j in range(J):
        for k in interims:
            y = tgt_of(k)[:, j]; ok = np.isfinite(y); yv = y[ok]
            if yv.size < 10:
                continue
            qs = qs_of(k)[ok, j, :]
            u = np.array([np.interp(yv[s], qs[s], taus, 0., 1.) for s in range(len(yv))])
            q = np.quantile(u, QL)
            rows.append(dict(interim_id=int(k), item_label=labels[j], ymin=float(q[0]),
                             lower=float(q[1]), middle=float(q[2]), upper=float(q[3]), ymax=float(q[4])))
    df = pd.DataFrame(rows)
    df['interim_id'] = pd.Categorical(df.interim_id, categories=sorted(df.interim_id.unique()), ordered=True)
    df['dev'] = df['middle'] - 0.5
    (ggplot(df, aes('interim_id', ymin='ymin', lower='lower', middle='middle',
                    upper='upper', ymax='ymax', fill='dev'))
     + geom_hline(yintercept=[0.1, 0.25, 0.5, 0.75, 0.9], linetype='dashed', colour='#9e9e9e', size=.3)
     + geom_boxplot(stat='identity', alpha=.9, size=.3, width=.7)
     + scale_fill_gradient2(low='#2166ac', mid='#f7f7f7', high='#b2182b', midpoint=0.0,
                            name='PIT median − 0.5\n(blue: median too high;\nred: too low)')
     + facet_wrap('~ item_label', ncol=4)
     + theme_bw() + theme(figure_size=(16, 20), panel_spacing=0.02, legend_position='top',
       strip_background=element_blank(), strip_text=element_text(face='bold', size=8),
       axis_text_x=element_text(size=6))
     + labs(x='interim id', y='PIT   u = F_amortiser(rho_SVI | x, z)',
       title=f'Conditional calibration (PIT) — {suf}  '
             f'(box on dashed refs .1/.25/.5/.75/.9 = calibrated; fill = median bias, '
             f'box wider than .25-.75 = over-confident)')).save(
        f"{out}/{prefix}_pps_{suf}_pit_box_by_item.pdf", verbose=False, limitsize=False)


def _sdrows(qs_of, tgt_of, n_of, labels, interims, blow):
    J = len(labels); rows = []
    for j in range(J):
        for k in interims:
            if labels[j] == blow:
                continue
            y = tgt_of(k)[:, j]; ok = np.isfinite(y); yv = y[ok]
            if yv.size < 10:
                continue
            sn = float(np.sqrt(n_of(k)))
            qs = qs_of(k)[:, j, :]; qs = qs[np.isfinite(qs).all(1)]
            mu = qs[:, 2]; sdw = (qs[:, 4] - qs[:, 0]) / 3.2897
            asd = float(np.sqrt(np.mean(sdw ** 2) + np.var(mu)))
            rows.append(dict(item_label=labels[j], n=n_of(k), sqrt_n=sn,
                             sd_svi=float(np.std(yv)), sd_amo=asd))
    return pd.DataFrame(rows)


def contraction_factor(qs_of, tgt_of, n_of, labels, interims, out, prefix, suf, blow):
    """SD of rho vs sqrt(n), per item: SVI (blue) vs amortiser marginal-predictive
    SD (red), each with a linear fit; slope/R2 annotated (SVI top, amortiser bottom)."""
    d = _sdrows(qs_of, tgt_of, n_of, labels, interims, blow)
    long = pd.concat([d[['item_label', 'sqrt_n', 'sd_svi']].rename(columns={'sd_svi': 'sd'}).assign(source='SVI'),
                      d[['item_label', 'sqrt_n', 'sd_amo']].rename(columns={'sd_amo': 'sd'}).assign(source='amortiser')])

    def _fit(g):
        b = np.polyfit(g.sqrt_n, g.sd, 1); pred = np.polyval(b, g.sqrt_n)
        ss = np.sum((g.sd - pred) ** 2); st = np.sum((g.sd - g.sd.mean()) ** 2)
        return float(b[0]), (1 - ss / st if st > 0 else np.nan)
    ann = []
    for j, gj in long.groupby('item_label'):
        for src, yy in (('SVI', 0.95), ('amortiser', 0.05)):
            s, r2 = _fit(gj[gj.source == src])
            ann.append(dict(item_label=j, source=src, txt=f"slope={s:+.3f} R2={r2:.2f}",
                            x=gj.sqrt_n.min(), y=gj.sd.min() + yy * (gj.sd.max() - gj.sd.min())))
    (ggplot(long, aes('sqrt_n', 'sd', colour='source'))
     + geom_point(size=.7, alpha=.6) + geom_line(stat='smooth', method='lm', se=False, size=.6)
     + geom_text(aes(x='x', y='y', label='txt', colour='source'), data=pd.DataFrame(ann),
                 ha='left', size=6, show_legend=False)
     + facet_wrap('~ item_label', ncol=4, scales='free')
     + scale_fill_manual(values={'SVI': CBLUE, 'amortiser': CRED})
     + theme_bw() + theme(figure_size=(16, 20), legend_position='top', panel_spacing=0.02,
       strip_background=element_blank(), strip_text=element_text(face='bold', size=8))
     + labs(x='sqrt(number of participants)  sqrt(n)', y='posterior SD of rho', colour='',
       title=f'Posterior-contraction rate — {suf}  (SVI blue vs amortiser red; under BvM SD ~ 1/sqrt n)')).save(
        f"{out}/{prefix}_pps_{suf}_posterior-contraction-factor.pdf", verbose=False, limitsize=False)


def contraction_law(qs_of, tgt_of, n_of, labels, interims, out, prefix, suf, blow):
    """Power-law fit SD = C n^-p: SVI (blue) vs amortiser (red) SD vs sqrt(n),
    with the fitted power law and per-item p / R2 annotated."""
    d = _sdrows(qs_of, tgt_of, n_of, labels, interims, blow)
    long = pd.concat([d[['item_label', 'n', 'sd_svi']].rename(columns={'sd_svi': 'sd'}).assign(source='SVI'),
                      d[['item_label', 'n', 'sd_amo']].rename(columns={'sd_amo': 'sd'}).assign(source='amortiser')])
    long['sqrt_n'] = np.sqrt(long.n)
    curves, ann = [], []
    for (j, src), g in long.groupby(['item_label', 'source']):
        gp = g[g.sd > 0]
        if len(gp) < 3:
            continue
        b = np.polyfit(np.log(gp.n), np.log(gp.sd), 1); p = -b[0]
        pred = np.exp(np.polyval(b, np.log(gp.n)))
        ss = np.sum((np.log(gp.sd) - np.log(pred)) ** 2); st = np.sum((np.log(gp.sd) - np.log(gp.sd).mean()) ** 2)
        r2 = 1 - ss / st if st > 0 else np.nan
        ng = np.linspace(gp.n.min(), gp.n.max(), 40)
        curves.append(pd.DataFrame(dict(item_label=j, source=src, sqrt_n=np.sqrt(ng),
                                        sd=np.exp(np.polyval(b, np.log(ng))))))
        ann.append(dict(item_label=j, source=src, txt=f"p={p:.2f} (R2={r2:.2f})",
                        x=np.sqrt(gp.n.min()), y=(0.95 if src == 'SVI' else 0.05)))
    cdf = pd.concat(curves) if curves else pd.DataFrame(columns=['item_label', 'source', 'sqrt_n', 'sd'])
    (ggplot(long, aes('sqrt_n', 'sd', colour='source'))
     + geom_point(size=.7, alpha=.6) + geom_line(data=cdf, size=.6)
     + facet_wrap('~ item_label', ncol=4, scales='free')
     + theme_bw() + theme(figure_size=(16, 20), legend_position='top', panel_spacing=0.02,
       strip_background=element_blank(), strip_text=element_text(face='bold', size=8))
     + labs(x='sqrt(number of participants)  sqrt(n)', y='posterior SD of rho', colour='',
       title=f'Trained (amortiser) vs actual (SVI) contraction — {suf}  (power-law SD=C n^-p)')).save(
        f"{out}/{prefix}_pps_{suf}_contraction-law_trained-vs-svi.pdf", verbose=False, limitsize=False)


def eta0_sweep(qs_of, tgt_of, n_of, labels, interims, taus, out, prefix, suf, eta_raw,
               etaH=0.89, date_of=None, blow=None):
    """eta_0 deployment sweep over a raw endpoint-change grid ``eta_raw`` (fractions).
    Post-hoc from the calibrated per-cohort quantiles ``qs_of(k)`` (S,J,5):
      PPS(eta0) = mean_s 1{ P(rho>eta0 | x,z^s) > etaH },  P(rho>eta0|x,z^s)=1-F_qs[s](eta0).
    Emits #6 PPS-by-item x eta0 grid, #7 SVI rho predictive with eta0 threshold lines, and a
    PPS-vs-eta0 decay curve. ``date_of(k)`` optionally labels the x-axis (else interim id)."""
    from plotnine import (geom_col, geom_boxplot, geom_hline, geom_line, geom_point,
                          facet_grid, scale_x_discrete, scale_colour_brewer, labeller)
    G = len(eta_raw); dlab = (lambda k: str(date_of(k))) if date_of is not None else (lambda k: str(k))
    pps_rows, rho_rows = [], []
    for k in interims:
        for j, lab in enumerate(labels):
            if lab == blow:
                continue
            q = qs_of(k)[:, j, :]; q = q[np.isfinite(q).all(1)]
            for e in eta_raw:
                ph1 = 1.0 - np.array([np.interp(e, q[s], taus, 0., 1.) for s in range(q.shape[0])])
                pps_rows.append(dict(interim=k, date=dlab(k), item=lab, eta0_pct=int(round(e * 100)),
                                     pps=float(np.mean(ph1 > etaH))))
            y = tgt_of(k)[:, j]; y = y[np.isfinite(y)]
            for r in y:
                rho_rows.append(dict(interim=k, date=dlab(k), item=lab, rho_pct=float(r) * 100.0))
    pdf = pd.DataFrame(pps_rows); rdf = pd.DataFrame(rho_rows)
    pdf.to_csv(f"{out}/{prefix}_pps_{suf}_eta0_pps.csv", index=False)
    rdf.to_csv(f"{out}/{prefix}_pps_{suf}_eta0_rho.csv", index=False)
    order_k = [f"{k:02d}" for k in interims]; date_lab = [dlab(k) for k in interims]
    for _d in (pdf, rdf):
        _d['xkey'] = pd.Categorical(_d['interim'].map(lambda k: f"{int(k):02d}"), categories=order_k, ordered=True)
    ldf = pd.DataFrame([dict(item=l, eta0_pct=int(round(e * 100)), rho_pct=float(e) * 100.0)
                        for l in rdf['item'].unique() for e in eta_raw])
    nit = pdf.item.nunique()
    (ggplot(pdf, aes('xkey', 'pps')) + geom_col(fill=CBLUE, width=.8)
     + facet_grid('item ~ eta0_pct', labeller=labeller(cols=lambda v: f'eta0={v}%'))
     + scale_x_discrete(breaks=order_k, labels=date_lab)
     + theme_bw() + theme(figure_size=(2 + 1.6 * G, 1.0 * nit), axis_text_x=element_text(rotation=90, size=4),
                          strip_text_y=element_text(angle=0, size=6), strip_text_x=element_text(size=7))
     + labs(x='interim', y='PPS', title=f'PPS by item x eta_0 (% endpoint change) — {suf}')).save(
        f"{out}/{prefix}_pps_{suf}_pps_by_item_eta0_grid.pdf", verbose=False, limitsize=False)
    (ggplot(rdf, aes('xkey', 'rho_pct')) + geom_boxplot(outlier_size=.2, fill='#d9d9d9', size=.3)
     + geom_hline(ldf, aes(yintercept='rho_pct', colour='factor(eta0_pct)'), size=.5)
     + facet_wrap('~ item', ncol=4, scales='free_y') + scale_x_discrete(breaks=order_k, labels=date_lab)
     + scale_colour_brewer(type='seq', palette='YlOrRd', name='eta_0 (%)')
     + theme_bw() + theme(figure_size=(16, 2.2 * ((nit + 3) // 4)), axis_text_x=element_text(rotation=90, size=4),
                          strip_text=element_text(size=6), legend_position='top')
     + labs(x='interim', y='rho = endpoint % change',
            title=f'SVI rho predictive vs eta_0 thresholds — {suf}')).save(
        f"{out}/{prefix}_pps_{suf}_rho_vs_eta0_lines_by_item.pdf", verbose=False, limitsize=False)
    sw = pdf.groupby(['eta0_pct', 'xkey'], observed=True)['pps'].mean().reset_index()
    (ggplot(sw, aes('eta0_pct', 'pps', colour='xkey', group='xkey'))
     + geom_line(size=.6) + geom_point(size=1.2)
     + scale_colour_brewer(type='seq', palette='Blues', name='interim')
     + theme_bw() + theme(figure_size=(8, 5), legend_position='right')
     + labs(x='eta_0 (% endpoint change)', y='mean PPS over items',
            title=f'PPS decay vs eta_0 — {suf}')).save(
        f"{out}/{prefix}_pps_{suf}_eta0_sweep_compare.pdf", verbose=False, limitsize=False)
    print(f"  saved eta0-sweep (3 plots, grid={[int(round(e*100)) for e in eta_raw]}%) ({suf}) -> {out}")


def all_plots(qs_of, tgt_of, n_of, labels, interims, taus, out, prefix, suf, blow):
    contraction_cdf(qs_of, tgt_of, labels, interims, taus, out, prefix, suf)
    pit_box(qs_of, tgt_of, labels, interims, taus, out, prefix, suf)
    contraction_factor(qs_of, tgt_of, n_of, labels, interims, out, prefix, suf, blow)
    contraction_law(qs_of, tgt_of, n_of, labels, interims, out, prefix, suf, blow)
    print(f"  saved 4 comparison plots ({suf}) -> {out}")


# =============================================================================
# FEDERATED combined diagnostics (§17.5 / §3.18-3.19): pool each endpoint's tidy plot data
# (RAGD_PLOTDATA) across the delegated amortisers of one federated parent and render, per
# diagnostic type, ONE combined figure faceted item (rows) x endpoint rho (columns). Shared by
# every app (SDY269 / SDY312 variants / CAVD ...); the per-app config lives in data_loading.
# Mixed-unit apps (columns in rate vs fold vs difference) stitch one own-y-scale column per
# endpoint (facet_grid free_y squashes across a row); single-/like-unit apps use one facet_grid.
# =============================================================================
import os as _os
import glob as _glob
import re as _re
import pickle as _pickle


def _fed_pal(keys):
    import ggsci, matplotlib.colors as mcolors
    pal = [mcolors.to_hex(c) for c in ggsci.pal_futurama("planetexpress")(12)]
    return {k: pal[i % len(pal)] for i, k in enumerate(keys)}


class FederatedDiagnostics:
    """Render the combined amortiser diagnostic suite for one federated parent dir. Construct from
    (sandbox root, app-config dict from data_loading.amortiser_app) and call .run()."""

    def __init__(self, sb, cfg, verbose=True):
        self.v = verbose; self.cfg = cfg
        self.FED = f"{sb}/{cfg['fed']}"; self.SRC = f"{sb}/{cfg['svi']}"; self.OUT = self.FED
        self.pfx = cfg.get('pfx', 'pcm_1_interim')
        self.title = cfg['title']; self.mixed = bool(cfg.get('mixed_units', False))
        self.eta_units = cfg.get('eta0_units', 'endpoint units x100')
        self.calib = cfg.get('calib_prefix', 'amortiser')
        self.item_word = cfg.get('item_kind', 'item')
        self.eps = self._resolve(cfg['endpoints'])
        self.RHO = [e['rho'] for e in self.eps]
        self.DIR = {e['rho']: f"{self.FED}/{e['rho']}" for e in self.eps}
        self.INST = {e['rho']: e['instance'] for e in self.eps}
        self.LONG = {e['rho']: e['long'] for e in self.eps}
        self.LONG_ORDER = [self.LONG[r] for r in self.RHO]
        self.ITEMS = list(pd.read_csv(f"{self.SRC}/{self.pfx}_1_data_dit.csv").item_label)
        self.ILAB = self._ilab(); self.IORDER = [self.ILAB[k] for k in sorted(self.ILAB)]

    # ---- config / label resolution -------------------------------------------------------
    def _resolve(self, eps):
        fs = sorted(_glob.glob(f"{self.SRC}/{self.pfx}_i*_regression_training.pkl"),
                    key=lambda f: int(_re.search(r'_i(\d+)_', f).group(1)))
        x = pd.read_pickle(fs[-1]) if fs else pd.DataFrame()
        # single-rho legacy pkls (Ukraine/covid) carry neither rho_id nor rho_label_long; guard on
        # column presence so the label lookup degrades to the endpoint's own `long`/`rho` (below).
        has = lambda cols: len(x) and set(cols).issubset(x.columns)
        by_lab = (x[['rho_label', 'rho_label_long']].drop_duplicates()
                  .set_index('rho_label')['rho_label_long'].to_dict()) if has(['rho_label', 'rho_label_long']) else {}
        by_id = (x[['rho_id', 'rho_label_long']].drop_duplicates()
                 .set_index('rho_id')['rho_label_long'].to_dict()) if has(['rho_id', 'rho_label_long']) else {}
        out = []
        for e in eps:
            e = dict(e)
            if 'long' not in e:
                e['long'] = by_id.get(e['rho_id']) if e.get('rho_id') is not None else by_lab.get(e['rho'], e['rho'])
            out.append(e)
        return out

    def _ilab(self):
        for r in self.RHO:
            p = f"{self.DIR[r]}/{self.pfx}_pps_RAGD_pps_by_item.csv"
            if _os.path.exists(p):
                m = pd.read_csv(p)[['interim_id', 'n']].drop_duplicates().set_index('interim_id')['n'].to_dict()
                return {int(k): f"interim {int(k)}\n(n={int(vv)})" for k, vv in m.items()}
        return {}

    # ---- small helpers -------------------------------------------------------------------
    def _catcols(self, df):
        df = df.copy()
        if 'item_label' in df:
            df['item_label'] = pd.Categorical(df['item_label'], categories=self.ITEMS, ordered=True)
        df['rho_label_long'] = pd.Categorical(
            df['rho_label_long'], ordered=True,
            categories=[c for c in self.LONG_ORDER if c in set(df['rho_label_long'])])
        return df

    def _xcat(self, df, col):
        df = df.copy()
        df['xlab'] = pd.Categorical(df[col].astype(int).map(self.ILAB), categories=self.IORDER, ordered=True)
        return df

    def _rot(self):
        return theme(axis_text_x=element_text(rotation=60, ha='right', size=6))

    def _grid(self, scales='fixed'):
        from plotnine import facet_grid
        return facet_grid('item_label ~ rho_label_long', scales=scales)

    def _theme(self, h=None, w=None):
        h = 1.95 * len(self.ITEMS) if h is None else h
        return (theme_bw() + theme(figure_size=(3.6 * len(self.RHO) if w is None else w, h),
                                   strip_background=element_blank(),
                                   strip_text=element_text(face='bold', size=7),
                                   strip_text_y=element_text(angle=270, size=7),
                                   legend_position='top', panel_spacing=0.02))

    def _rhos_in(self, df):
        return [r for r in self.LONG_ORDER if r in set(df['rho_label_long'].dropna())]

    # ---- data loading --------------------------------------------------------------------
    def load_plotdata(self):
        store = {}
        for r in self.RHO:
            p = f"{self.DIR[r]}/{self.pfx}_pps_RAGD_plotdata.pkl"
            if not _os.path.exists(p):
                continue
            with open(p, 'rb') as fh:
                frames = _pickle.load(fh)
            for key, df in frames.items():
                suf, name = key.split('::')
                df = df.copy(); df['rho_label'] = r; df['rho_label_long'] = self.LONG[r]
                store.setdefault(name, {}).setdefault(suf, []).append(df)
        return {name: {suf: pd.concat(v, ignore_index=True) for suf, v in bysuf.items()}
                for name, bysuf in store.items()}

    def _read_eta(self, name):
        rows = []
        for r in self.RHO:
            p = f"{self.DIR[r]}/{self.pfx}_pps_RAGD_{name}.csv"
            if not _os.path.exists(p):
                continue
            e = pd.read_csv(p).rename(columns={'item': 'item_label'})
            nmap = (pd.read_csv(f"{self.DIR[r]}/{self.pfx}_pps_RAGD_pps_by_item.csv")[['interim_id', 'n']]
                    .drop_duplicates().rename(columns={'interim_id': 'interim'}))
            e = e.merge(nmap, on='interim', how='left'); e['rho_label_long'] = self.LONG[r]; rows.append(e)
        return self._catcols(pd.concat(rows, ignore_index=True)) if rows else pd.DataFrame()

    def load_calibration(self):
        cal = []
        for r in self.RHO:
            p = f"{self.DIR[r]}/{self.pfx}_pps_RAGD_pit_by_item.csv"
            if _os.path.exists(p):
                df = pd.read_csv(p); m = df[['pit_ks', 'marg_ks', 'cov5', 'cov95']].mean()
                cal.append(dict(rho=r, instance=self.INST[r], n_item=len(df),
                                pit_ks=m.pit_ks, marg_ks=m.marg_ks, cov5=m.cov5, cov95=m.cov95))
        return pd.DataFrame(cal)

    # ---- stitch (mixed-units): one own-y-scale facet_grid column per endpoint ------------
    def _fill(self, d, rho_long):
        present = set(d['item_label'].dropna().astype(str))
        miss = [s for s in self.ITEMS if s not in present]
        if miss:
            pad = pd.DataFrame({'item_label': miss}); pad['rho_label_long'] = rho_long
            d = pd.concat([d, pad], ignore_index=True)
        if 'source' in d.columns:
            d['source'] = d['source'].fillna('SVI')
        d = self._catcols(d)
        if 'xlab' in d.columns and len(self.IORDER):
            d['xlab'] = pd.Categorical(d['xlab'].astype('object').where(d['xlab'].notna(), self.IORDER[0]),
                                       categories=self.IORDER, ordered=True)
        return d

    @staticmethod
    def _wrapcols(rhos, width=26):
        import textwrap
        w = {r: ('\n'.join(textwrap.wrap(str(r), width)) or str(r)) for r in rhos}
        m = max((s.count('\n') for s in w.values()), default=0)
        return {r: s + '\n' * (m - s.count('\n')) for r, s in w.items()}

    @staticmethod
    def _strip(i, n):
        ex = {}
        if i != n - 1:
            ex['strip_text_y'] = element_blank()
        if i != 0:
            ex['axis_title_y'] = element_blank()
        return theme(**ex) if ex else theme()

    def _stitch_save(self, subplots, out, h, title='', legend=None, legend_title='', w_each=4.7):
        import tempfile, matplotlib.pyplot as plt, matplotlib.image as mpimg
        from matplotlib.patches import Patch
        dpi = 150; tmp, imgs = [], []
        for p in subplots:
            f = tempfile.NamedTemporaryFile(suffix='.png', delete=False).name; tmp.append(f)
            p.save(f, width=w_each, height=h, dpi=dpi, verbose=False, limitsize=False)
            imgs.append(mpimg.imread(f))
        ws = [im.shape[1] for im in imgs]; hpx = max(im.shape[0] for im in imgs)
        band_t = 0.75 if title else 0.0; band_b = 0.65 if legend else 0.0
        Wax, Hax = sum(ws) / dpi, hpx / dpi; H = Hax + band_t + band_b
        fig, axes = plt.subplots(1, len(imgs), figsize=(Wax, H), gridspec_kw={'width_ratios': ws})
        axes = np.atleast_1d(axes)
        for ax, im in zip(axes, imgs):
            ax.imshow(im); ax.axis('off')
        fig.subplots_adjust(left=0, right=1, top=(Hax + band_b) / H, bottom=band_b / H, wspace=0.01)
        if title:
            fig.suptitle(title, y=1 - band_t / H * 0.45, fontsize=13)
        if legend:
            fig.legend(handles=[Patch(facecolor=c, label=l) for l, c in legend], title=legend_title,
                       loc='lower center', ncol=min(len(legend), 8), frameon=False, fontsize=9)
        fig.savefig(out, dpi=dpi); fig.savefig(out[:-4] + '.png', dpi=dpi); plt.close(fig)
        for f in tmp:
            _os.remove(f)

    def _grid_or_stitch(self, df, out, *, geoms, aes0, title, y_title, x_title='interim',
                        scales='free_y', rot=True, legend=None, legend_title='', hlines=None):
        """If mixed_units: stitch per-endpoint columns (own y). Else: one facet_grid. `geoms(p, d, r)`
        adds the layers to a base ggplot(d, aes0); `hlines(r)` optionally returns a per-endpoint frame."""
        rhos = self._rhos_in(df)
        if not self.mixed:
            d = self._xcat(self._catcols(df), 'interim_id') if 'interim_id' in df else self._catcols(df)
            p = ggplot(d, aes(**aes0)); p = geoms(p, d, None)
            p = p + self._grid(scales) + self._theme(w=3.6 * len(self.RHO)) + (self._rot() if rot else theme())
            p = p + labs(x=x_title, y=y_title, title=self.title + ' — ' + title)
            p.save(out, verbose=False, limitsize=False); p.save(out[:-4] + '.png', dpi=120, verbose=False, limitsize=False)
            return
        if 'interim_id' in df:                          # xlab on the REAL rows first; _fill pads them valid
            df = self._xcat(df, 'interim_id')
        wl = self._wrapcols(rhos); subs = []
        for i, r in enumerate(rhos):
            dd = self._fill(df[df.rho_label_long == r].copy(), r).assign(rho_label_long=wl[r])
            p = ggplot(dd, aes(**aes0)); p = geoms(p, dd, r)
            if hlines is not None:
                p = hlines(p, r, wl[r])
            p = (p + self._grid(scales) + self._theme(w=4.7) + (self._rot() if rot else theme())
                 + self._strip(i, len(rhos)) + theme(legend_position='none')
                 + labs(x=x_title, y=(y_title if i == 0 else '')))
            subs.append(p)
        self._stitch_save(subs, out, 2.6 * len(self.ITEMS), title=self.title + ' — ' + title,
                          legend=legend, legend_title=legend_title)

    # ---- the figures ---------------------------------------------------------------------
    def fig_tests_pit(self, store, suf):
        from plotnine import geom_histogram
        df = store.get('pit_u', {}).get(suf)
        if df is None:
            return
        p = (ggplot(self._catcols(df), aes('u')) + geom_histogram(aes(y='..density..'), bins=20,
             fill=CBLUE, colour='white', size=.2) + geom_hline(yintercept=1, linetype='dashed', colour='black')
             + self._grid() + self._theme() + labs(x='PIT  u = F_amortiser(rho_SVI | x, z)', y='density',
               title=f'{self.title} — PIT uniformity ({suf}); flat = calibrated'))
        p.save(f"{self.OUT}/{self.pfx}_pps_{suf}_tests_pit.pdf", verbose=False, limitsize=False)

    def fig_tests_coverage(self, store, suf):
        from plotnine import geom_abline, scale_color_manual
        df = store.get('coverage', {}).get(suf)
        if df is None:
            return
        df = self._catcols(df)
        df['ilab'] = pd.Categorical(df.interim_id.astype(int).map(self.ILAB), categories=self.IORDER, ordered=True)
        p = (ggplot(df, aes('eta_inpol', 'coverage', colour='ilab', group='ilab'))
             + geom_abline(intercept=0, slope=1, linetype='dashed', colour='black')
             + geom_line(size=.4, alpha=.85) + geom_point(size=.8)
             + scale_color_manual(values=_fed_pal(self.IORDER), name='interim') + self._grid('free_y') + self._theme()
             + labs(x='nominal coverage', y='empirical coverage',
               title=f'{self.title} — coverage calibration ({suf}); on the diagonal = calibrated'))
        p.save(f"{self.OUT}/{self.pfx}_pps_{suf}_tests_coverage.pdf", verbose=False, limitsize=False)

    def fig_contraction_cdf(self, store, suf):
        df = store.get('contbox', {}).get(suf)
        if df is None:
            return
        tag = '[deployed: head-ft + affine + BvM]' if suf == 'RAGD' else '[RAW net: no recalibration]'

        def geoms(p, d, r):
            return (p + geom_boxplot(stat='identity', position=position_dodge(width=0.78), size=.3,
                                     width=.7, alpha=.85)
                    + scale_fill_manual(values={'SVI': CBLUE, 'amortiser': CRED}, name=''))
        self._grid_or_stitch(
            df, f"{self.OUT}/{self.pfx}_pps_{suf}_contraction_cdf_by_item.pdf", geoms=geoms,
            aes0=dict(x='xlab', ymin='ymin', lower='lower', middle='middle', upper='upper',
                      ymax='ymax', fill='source'),
            title=f'marginal p(rho|x): SVI vs amortiser ({suf}) {tag}', y_title='p(rho | x) estimate',
            legend=[('SVI', CBLUE), ('amortiser', CRED)], legend_title='source')

    def fig_pit_box(self, store, suf):
        from plotnine import scale_fill_gradient2
        from mizani.bounds import squish
        df = store.get('pitbox', {}).get(suf)
        if df is None:
            return
        df = self._xcat(self._catcols(df), 'interim_id'); df['dev'] = df['middle'] - 0.5
        p = (ggplot(df, aes('xlab', ymin='ymin', lower='lower', middle='middle',
                            upper='upper', ymax='ymax', fill='dev'))
             + geom_hline(yintercept=[0.10, 0.25, 0.50, 0.75, 0.90], linetype='dashed', colour='#9e9e9e', size=.3)
             + geom_boxplot(stat='identity', alpha=.9, size=.3, width=.7)
             + scale_fill_gradient2(low='#2166ac', mid='#f7f7f7', high='#b2182b', midpoint=0.0,
                 limits=[-0.3, 0.3], oob=squish, name='PIT median − 0.5')
             + self._grid() + self._theme(w=3.6 * len(self.RHO)) + self._rot()
             + labs(x='interim', y='PIT  u = F_amortiser(rho_SVI | x, z)',
               title=f'{self.title} — conditional calibration (PIT box on dashed refs = calibrated)'))
        p.save(f"{self.OUT}/{self.pfx}_pps_{suf}_pit_box_by_item.pdf", verbose=False, limitsize=False)

    def fig_contraction_law(self, store, suf='RAGD'):
        from plotnine import geom_text, scale_color_manual
        df = store.get('contlaw', {}).get(suf)
        if df is None:
            return
        df = self._catcols(df)

        def _pow(g):
            n = g.n.values.astype(float); y = np.maximum(g.sd.values, 1e-6)
            b = np.polyfit(np.log(n), np.log(y), 1); pp = -b[0]; C = np.exp(b[1])
            pr = C * n ** (-pp); ss = np.sum((y - pr) ** 2); st = np.sum((y - y.mean()) ** 2)
            return pp, C, (1 - ss / st if st > 0 else np.nan)
        cur, ann = [], []
        for (it, r), gi in df.groupby(['item_label', 'rho_label_long'], observed=True):
            yhi, ylo = gi.sd.max(), gi.sd.min(); xlo = gi.sqrt_n.min()
            for src, isS in (('SVI', True), ('amortiser', False)):
                g = gi[gi.source == src]
                if len(g) < 4:
                    continue
                pp, C, r2 = _pow(g); ng = np.linspace(g.n.min(), g.n.max(), 60)
                for xx, yy in zip(np.sqrt(ng), C * ng ** (-pp)):
                    cur.append(dict(item_label=it, rho_label_long=r, source=src, sqrt_n=xx, sd=yy))
                ann.append(dict(item_label=it, rho_label_long=r, sqrt_n=xlo, sd=(yhi if isS else ylo),
                                src=src, label=f'{src[:3]}: p={pp:.2f} (R2={r2:.2f})'))
        cur = self._catcols(pd.DataFrame(cur)) if cur else pd.DataFrame()
        a = self._catcols(pd.DataFrame(ann)) if ann else pd.DataFrame()
        out = f"{self.OUT}/{self.pfx}_pps_{suf}_contraction-law_trained-vs-svi.pdf"
        rhos = self._rhos_in(df)
        if not self.mixed:
            p = (ggplot(df, aes('sqrt_n', 'sd', colour='source')) + geom_point(size=1.3)
                 + scale_color_manual(values={'SVI': CBLUE, 'amortiser': CRED}, name=''))
            if len(cur):
                p = p + geom_line(cur, aes('sqrt_n', 'sd', colour='source'), size=.6)
            for src, col in (('SVI', CBLUE), ('amortiser', CRED)):
                s = a[a.src == src] if len(a) else a
                if len(s):
                    p = p + geom_text(s, aes('sqrt_n', 'sd', label='label'), inherit_aes=False,
                                      ha='left', va='top', size=6, colour=col)
            p = (p + self._grid('free') + self._theme(w=3.6 * len(self.RHO))
                 + labs(x='sqrt(participants)  sqrt(n)', y='posterior SD of rho',
                   title=f'{self.title} — trained (amortiser) vs actual (SVI) contraction (SD = C n^-p)'))
            p.save(out, verbose=False, limitsize=False); p.save(out[:-4] + '.png', dpi=120, verbose=False, limitsize=False)
            return
        wl = self._wrapcols(rhos); subs = []
        for i, r in enumerate(rhos):
            rw = wl[r]
            d = self._fill(df[df.rho_label_long == r].copy(), r).assign(rho_label_long=rw)
            cr = cur[cur.rho_label_long == r].assign(rho_label_long=rw) if len(cur) else cur
            ar = a[a.rho_label_long == r].assign(rho_label_long=rw) if len(a) else a
            p = (ggplot(d, aes('sqrt_n', 'sd', colour='source')) + geom_point(size=1.3)
                 + scale_color_manual(values={'SVI': CBLUE, 'amortiser': CRED}, name=''))
            if len(cr):
                p = p + geom_line(cr, aes('sqrt_n', 'sd', colour='source'), size=.6)
            for src, col in (('SVI', CBLUE), ('amortiser', CRED)):
                s = ar[ar.src == src] if len(ar) else ar
                if len(s):
                    p = p + geom_text(s, aes('sqrt_n', 'sd', label='label'), inherit_aes=False,
                                      ha='left', va='top', size=6, colour=col)
            p = (p + self._grid('free') + self._theme(w=4.7) + self._strip(i, len(rhos)) + theme(legend_position='none')
                 + labs(x='sqrt(participants)  sqrt(n)', y=('posterior SD of rho' if i == 0 else '')))
            subs.append(p)
        self._stitch_save(subs, out, 2.6 * len(self.ITEMS),
                          title=f'{self.title} — trained (amortiser) vs actual (SVI) contraction (SD = C n^-p)',
                          legend=[('SVI', CBLUE), ('amortiser', CRED)], legend_title='source')

    def fig_rho_vs_eta0(self):
        from plotnine import scale_color_manual
        rdf = self._read_eta('eta0_rho')
        if rdf.empty:
            return
        pdf = self._read_eta('eta0_pps')
        ldf = (pdf[['item_label', 'rho_label_long', 'eta0_pct']].drop_duplicates()
               .assign(rho_pct=lambda d: d.eta0_pct.astype(float)))
        ldf['eta0'] = ldf.eta0_pct.astype(int).astype(str)
        thr = sorted(ldf.eta0_pct.unique()); cd = _fed_pal([str(int(t)) for t in thr])
        out = f"{self.OUT}/{self.pfx}_pps_RAGD_rho_vs_eta0_lines_by_item.pdf"
        rhos = self._rhos_in(rdf)
        if not self.mixed:
            d = self._xcat(rdf, 'interim')
            p = (ggplot(d, aes('xlab', 'rho_pct')) + geom_boxplot(outlier_size=.15, fill='#d9d9d9', size=.3)
                 + geom_hline(ldf, aes(yintercept='rho_pct', colour='eta0'), size=.5)
                 + scale_color_manual(values=cd, name=f'success threshold eta0 ({self.eta_units})')
                 + self._grid('free_y') + self._theme(w=3.6 * len(self.RHO)) + self._rot()
                 + labs(x='interim', y=f'rho  ({self.eta_units})',
                   title=f'{self.title} — SVI rho predictive vs eta0 thresholds'))
            p.save(out, verbose=False, limitsize=False); p.save(out[:-4] + '.png', dpi=120, verbose=False, limitsize=False)
            return
        wl = self._wrapcols(rhos); subs = []
        rdf = self._xcat(rdf, 'interim')
        for i, r in enumerate(rhos):
            rw = wl[r]
            d = self._fill(rdf[rdf.rho_label_long == r].copy(), r).assign(rho_label_long=rw)
            l = ldf[ldf.rho_label_long == r].assign(rho_label_long=rw)
            p = (ggplot(d, aes('xlab', 'rho_pct')) + geom_boxplot(outlier_size=.15, fill='#d9d9d9', size=.3)
                 + geom_hline(l, aes(yintercept='rho_pct', colour='eta0'), size=.5)
                 + scale_color_manual(values=cd, name=f'success threshold eta0 ({self.eta_units})')
                 + self._grid('free_y') + self._theme(w=4.7) + self._rot() + self._strip(i, len(rhos))
                 + theme(legend_position='none')
                 + labs(x='interim', y=(f'rho  ({self.eta_units})' if i == 0 else '')))
            subs.append(p)
        self._stitch_save(subs, out, 2.6 * len(self.ITEMS),
                          title=f'{self.title} — SVI rho predictive (grey box) vs eta0 thresholds',
                          legend=[('SVI rho', '#d9d9d9')] + [(t, cd[t]) for t in [str(int(x)) for x in thr]],
                          legend_title=f'eta0 threshold ({self.eta_units})')

    def fig_eta0_sweep(self):
        from plotnine import facet_grid, scale_color_manual
        pdf = self._read_eta('eta0_pps')
        if pdf.empty:
            return
        pdf['ilab'] = pd.Categorical(pdf.interim.astype(int).map(self.ILAB), categories=self.IORDER, ordered=True)
        p = (ggplot(pdf, aes('eta0_pct', 'pps', colour='ilab', group='ilab'))
             + geom_line(size=.5) + geom_point(size=1.1)
             + scale_color_manual(values=_fed_pal(self.IORDER), name='interim')
             + facet_grid('item_label ~ rho_label_long', scales='free_x') + self._theme()
             + labs(x=f'success threshold eta0 ({self.eta_units}, %)', y='PPS:  P( P(rho > eta0 | x) > 0.89 )',
               title=f'{self.title} — PPS decay vs eta0 threshold'))
        p.save(f"{self.OUT}/{self.pfx}_pps_RAGD_eta0_sweep_compare.pdf", verbose=False, limitsize=False)

    def fig_trajectory(self):
        from plotnine import facet_grid, scale_y_continuous, scale_color_manual
        import matplotlib.cm as cm, matplotlib.colors as mcolors
        df = self._read_eta('eta0_pps')
        if df.empty:
            return
        df = self._xcat(df, 'interim')
        thr = sorted(df.eta0_pct.unique()); cmap = cm.get_cmap('viridis'); n = max(1, len(thr) - 1)
        cd = {str(int(t)): mcolors.to_hex(cmap(1 - i / n)) for i, t in enumerate(thr)}
        df['eta0'] = pd.Categorical(df.eta0_pct.astype(int).astype(str),
                                    categories=[str(int(t)) for t in thr], ordered=True)
        p = (ggplot(df, aes('xlab', 'pps', colour='eta0', group='eta0'))
             + geom_hline(yintercept=[0.1, 0.9], linetype='dashed', colour='#999999')
             + geom_line() + geom_point(size=1.3)
             + scale_color_manual(values=cd, name=f'success threshold eta0 ({self.eta_units}; darker = harder)')
             + scale_y_continuous(limits=[0, 1], labels=lambda l: [f'{v:.0%}' for v in l])
             + facet_grid('item_label ~ rho_label_long') + self._theme() + self._rot()
             + labs(x='interim', y='amortised PPS',
               title=f'{self.title} — amortised PPS trajectory (dashed 10%/90% go/no-go guides)'))
        p.save(f"{self.OUT}/{self.pfx}_pps_RAGD_trajectory.pdf", verbose=False, limitsize=False)
        p.save(f"{self.OUT}/{self.pfx}_pps_RAGD_trajectory.png", dpi=110, verbose=False, limitsize=False)

    def fig_calibration(self, cal):
        from plotnine import geom_col, facet_wrap, coord_flip
        m = cal.melt(id_vars=['rho', 'instance'], value_vars=['pit_ks', 'marg_ks'],
                     var_name='metric', value_name='ks')
        m['rho_f'] = pd.Categorical(m['rho'], categories=self.RHO[::-1], ordered=True)
        p = (ggplot(m, aes(x='rho_f', y='ks', fill='instance')) + geom_col(show_legend=True) + coord_flip()
             + geom_hline(yintercept=0.10, linetype='dashed', color='#555555') + facet_wrap('metric', ncol=2)
             + scale_fill_manual(values=_fed_pal(cal['instance'].unique()), name='amortiser')
             + theme_bw() + theme(figure_size=(9, 3), strip_text=element_text(size=9))
             + labs(x='endpoint (federated amortiser)', y='KS vs SVI (dashed = 0.10)',
                    title=f'{self.title}: federated amortiser calibration (mean over {self.item_word}s)'))
        p.save(f"{self.OUT}/{self.calib}_calibration.pdf", verbose=False)
        p.save(f"{self.OUT}/{self.calib}_calibration.png", dpi=120, verbose=False)

    # ---- orchestration -------------------------------------------------------------------
    def run(self):
        cal = self.load_calibration()
        cal.to_csv(f"{self.OUT}/{self.calib}_calibration.csv", index=False)
        if self.v:
            pd.set_option('display.width', 160)
            print(f"=== {self.title} federated amortiser calibration (mean over {self.item_word}s) "
                  f"[{_os.path.basename(self.FED)}] ===")
            print(cal.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
        if not cal.empty:
            self.fig_calibration(cal)
        store = self.load_plotdata()
        for suf in ('RAGD', 'RAGDbase'):
            self.fig_tests_pit(store, suf); self.fig_tests_coverage(store, suf); self.fig_contraction_cdf(store, suf)
        self.fig_pit_box(store, 'RAGD'); self.fig_contraction_law(store, 'RAGD')
        self.fig_trajectory(); self.fig_rho_vs_eta0(); self.fig_eta0_sweep()
        for f in _glob.glob(f"{self.FED}/*/*.pdf"):          # subdirs are data-only
            _os.remove(f)
        if self.v:
            print(f"combined {self.item_word} x rho figures -> {self.OUT}/{self.pfx}_pps_*.pdf (+ {self.calib}_calibration)")
