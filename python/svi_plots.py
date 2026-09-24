"""Generic SVI-side plots shared by the applications (the amortiser diagnostics live in
amortiser_diag_plots; the per-category fit bars in utils._plot_prob_barplots). Currently:

  rho_contraction(svi_dir, ...) — the p(rho|x) contraction figure: for each endpoint rho declared in
  the per-interim i{k}_regression_training.pkl, box the posterior of rho across draws at each interim,
  faceted item/strain (row) x rho (column, free y). Dashed line = the H1 threshold inferred from the
  rho_label (SPR>0.70, GMFR>2.5, difference at 0). Writes pcm_..._p_rho_x_by_item.{csv,pdf,png}.

CLI (back-compat with the old IMMPORT_flu_rho_plots.py): default = SDY312/SDY314; override with
IMMPORT_RHO_DIRS='label:subdir[:file_prefix[:dp_prefix]],...'.
"""
import os
import glob
import re
import numpy as np
import pandas as pd

DEFAULT_PREFIX = "pcm_1_interim"


def _h1(rho_label):
    """H1 reference line inferred from the rho_label (short 'LAIV_spr'/'TIV_gmfr'/*_diff or long)."""
    rl = str(rho_label).lower()
    if 'diff' in rl:
        return 0.0
    if 'spr' in rl or 'seroprot' in rl:
        return 0.70
    if 'gmfr' in rl or 'gmt' in rl or 'fold' in rl:
        return 2.5
    return 0.0


def rho_contraction_summary(svi_dir, file_prefix=DEFAULT_PREFIX, dp_prefix=None):
    """Per (interim, n, item, rho) boxplot quantiles of the rho posterior across draws."""
    dp_prefix = dp_prefix or file_prefix
    fs = sorted(glob.glob(f"{svi_dir}/{file_prefix}_i*_regression_training.pkl"),
                key=lambda p: int(re.search(r'_i(\d+)_', p).group(1)))
    rows = []
    for f in fs:
        k = int(re.search(r'_i(\d+)_', f).group(1))
        n = pd.read_csv(f"{svi_dir}/{dp_prefix}_{k}_data_dp1.csv")['pid'].nunique()
        x = pd.read_pickle(f).copy()
        if 'rho_label_long' not in x.columns:
            x['rho_label_long'] = x['rho_label']
        for (item, rid, rlab, rll), g in x.groupby(['item_label', 'rho_id', 'rho_label', 'rho_label_long']):
            q = g['pps_rho_x'].quantile([0.025, 0.25, 0.5, 0.75, 0.975]).to_numpy()
            rows.append(dict(interim=k, n=n, item=item, rho_id=rid, rho_label=rlab, rho_label_long=rll,
                             ymin=q[0], lower=q[1], middle=q[2], upper=q[3], ymax=q[4]))
    return pd.DataFrame(rows)


def rho_contraction(svi_dir, study=None, file_prefix=DEFAULT_PREFIX, dp_prefix=None, height_mult=1.0):
    """SVI p(rho|x) contraction figure, faceted item (row) x rho (column). Saves into svi_dir."""
    import ggsci, matplotlib.colors as mcolors
    from plotnine import (ggplot, aes, geom_boxplot, geom_hline, facet_wrap, theme_bw, theme,
                          labs, scale_fill_manual, element_text)
    study = study or os.path.basename(svi_dir)
    df = rho_contraction_summary(svi_dir, file_prefix, dp_prefix)
    if df.empty:
        print(f"{study}: no endpoint pkls"); return
    df.to_csv(f"{svi_dir}/{file_prefix}_p_rho_x_by_item.csv", index=False)
    items = sorted(df['item'].unique())
    il = df[['interim', 'n']].drop_duplicates().sort_values('interim')
    xlab = {int(r.interim): f"interim {int(r.interim)}\n(n={int(r.n)})" for r in il.itertuples()}
    xorder = [xlab[int(i)] for i in il['interim']]
    rho_ord = df[['rho_id', 'rho_label', 'rho_label_long']].drop_duplicates().sort_values('rho_id')
    rhos = rho_ord['rho_label_long'].tolist(); nrho = len(rhos)
    df['n_f'] = pd.Categorical(df['interim'].astype(int).map(xlab), categories=xorder, ordered=True)
    df['item_f'] = pd.Categorical(df['item'], categories=items, ordered=True)
    df['rho_f'] = pd.Categorical(df['rho_label_long'], categories=rhos, ordered=True)
    lab_of = dict(zip(rho_ord['rho_label_long'], rho_ord['rho_label']))
    thr = pd.DataFrame([(it, rl, _h1(lab_of[rl])) for it in items for rl in rhos],
                       columns=['item_f', 'rho_f', 'yint'])
    thr['item_f'] = pd.Categorical(thr['item_f'], categories=items, ordered=True)
    thr['rho_f'] = pd.Categorical(thr['rho_f'], categories=rhos, ordered=True)
    pal = ggsci.pal_futurama("planetexpress")(12)
    cols = [mcolors.to_hex(c) for c in
            mcolors.LinearSegmentedColormap.from_list("f", pal)(np.linspace(0, 1, len(items)))]
    p = (ggplot(df, aes(x='n_f'))
         + geom_boxplot(aes(ymin='ymin', lower='lower', middle='middle', upper='upper',
                            ymax='ymax', fill='item_f'), stat='identity', alpha=0.8, width=0.7)
         + geom_hline(thr, aes(yintercept='yint'), linetype='dashed', color='#555555')
         + facet_wrap(['item_f', 'rho_f'], ncol=nrho, scales='free_y')
         + scale_fill_manual(values=dict(zip(items, cols)), guide=None)
         + theme_bw()
         + theme(axis_text_x=element_text(angle=60, vjust=1, hjust=1, size=6),
                 strip_text=element_text(size=7),
                 figure_size=(max(6, 3.2 * nrho), max(10, 1.9 * len(items)) * height_mult))
         + labs(x='interim', y='p(rho | x)',
                title=f'{study}: SVI p(rho|x) contraction by item (row) x rho (column) (dashed = H1 threshold)'))
    p.save(f"{svi_dir}/{file_prefix}_p_rho_x_by_item.pdf", verbose=False, limitsize=False)
    p.save(f"{svi_dir}/{file_prefix}_p_rho_x_by_item.png", dpi=110, verbose=False, limitsize=False)
    print(f"{study}: {df['n'].nunique()} interims x {len(items)} items x {nrho} rho [{file_prefix}] -> {svi_dir}")


if __name__ == "__main__":
    SB = os.environ.get('IMMPORT_SB', "/Users/or105/sandbox/bIRTistic")
    hm = float(os.environ.get('PROB_FIT_HEIGHT_MULT', '1.0'))
    if os.environ.get('IMMPORT_RHO_DIRS'):
        for x in os.environ['IMMPORT_RHO_DIRS'].split(','):
            parts = x.split(':'); lab, sub = parts[0], parts[1]
            fp = parts[2] if len(parts) > 2 else DEFAULT_PREFIX
            dpp = parts[3] if len(parts) > 3 else fp
            rho_contraction(f"{SB}/{sub}", study=lab, file_prefix=fp, dp_prefix=dpp, height_mult=hm)
    else:
        for study in ('SDY312', 'SDY314'):
            rho_contraction(f"{SB}/py-immport-{study}_260918", study=study, height_mult=hm)
