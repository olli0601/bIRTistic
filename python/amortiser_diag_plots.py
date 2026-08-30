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


def all_plots(qs_of, tgt_of, n_of, labels, interims, taus, out, prefix, suf, blow):
    contraction_cdf(qs_of, tgt_of, labels, interims, taus, out, prefix, suf)
    pit_box(qs_of, tgt_of, labels, interims, taus, out, prefix, suf)
    contraction_factor(qs_of, tgt_of, n_of, labels, interims, out, prefix, suf, blow)
    contraction_law(qs_of, tgt_of, n_of, labels, interims, out, prefix, suf, blow)
    print(f"  saved 4 comparison plots ({suf}) -> {out}")
