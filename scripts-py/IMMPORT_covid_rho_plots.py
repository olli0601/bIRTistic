"""p(rho|x) contraction plot for the ImmPort COVID subpopulation-contrast grids (§3.20):
the between-group endpoint rho (group2-vs-group1 relative shift in mean log2 ID50 titre) as it
evolves over accruing participants. x = cumulative participants, y = p(rho|x) posterior box
across draws (contracting with n), one panel per item. Dashed line at rho=0 (no group
difference). Writes pcm_1_interim_p_rho_x_by_item.{csv,pdf,png} in each COVID grid dir."""
# ---- boilerplate ----

import os, sys, glob, re
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
import numpy as np, pandas as pd

SB = os.environ.get('COVID_SB', "/Users/or105/sandbox/bIRTistic")
DIRS = sorted(glob.glob(f"{SB}/py-immport-covid-SDY*_260920"))
file_prefix = "pcm_1_interim"
_hmult = float(os.environ.get('PROB_FIT_HEIGHT_MULT', '1.0'))


def build_summary(d):
    fs = sorted(glob.glob(f"{d}/{file_prefix}_i*_regression_training.pkl"),
                key=lambda p: int(re.search(r'_i(\d+)_', p).group(1)))
    rows = []
    for f in fs:
        k = int(re.search(r'_i(\d+)_', f).group(1))
        n = pd.read_csv(f"{d}/{file_prefix}_{k}_data_dp1.csv")['pid'].nunique()
        x = pd.read_pickle(f)
        for item, g in x.groupby('item_label'):
            q = g['pps_ratio_x'].quantile([0.025, 0.25, 0.5, 0.75, 0.975]).to_numpy()
            rows.append(dict(n=n, item=item, ymin=q[0], lower=q[1], middle=q[2],
                             upper=q[3], ymax=q[4]))
    return pd.DataFrame(rows)


def plot_dir(d):
    import ggsci, matplotlib.colors as mcolors
    from plotnine import (ggplot, aes, geom_boxplot, geom_hline, facet_wrap, theme_bw, theme,
                          labs, scale_fill_manual, element_text)
    df = build_summary(d)
    if df.empty:
        return
    df.to_csv(f"{d}/{file_prefix}_p_rho_x_by_item.csv", index=False)
    ns = sorted(df['n'].unique()); items = sorted(df['item'].unique())
    df['n_f'] = pd.Categorical(df['n'].astype(str), categories=[str(x) for x in ns], ordered=True)
    pal = ggsci.pal_futurama("planetexpress")(12)
    cols = [mcolors.to_hex(c) for c in
            mcolors.LinearSegmentedColormap.from_list("f", pal)(np.linspace(0, 1, max(len(items), 2)))]
    cdict = dict(zip(items, cols))
    name = os.path.basename(d).replace('py-immport-covid-', '').replace('_260920', '')
    p = (ggplot(df, aes(x='n_f'))
         + geom_hline(yintercept=0, linetype='dashed', color='#555555')     # rho=0: no group difference
         + geom_boxplot(aes(ymin='ymin', lower='lower', middle='middle', upper='upper',
                            ymax='ymax', fill='item'), stat='identity', alpha=0.8, width=0.7)
         + facet_wrap('item', ncol=1, scales='free_y')
         + scale_fill_manual(values=cdict, guide=None)
         + theme_bw()
         + theme(axis_text_x=element_text(angle=45, vjust=1, hjust=1, size=7),
                 strip_text=element_text(size=8), figure_size=(9, 4.2 * len(items) * _hmult))
         + labs(x='interim (cumulative participants)',
                y='p(rho | x)  —  group2-vs-group1 relative shift in mean log2 titre',
                title=f'{name}: SVI p(rho|x) contraction (dashed = 0, no group difference)'))
    p.save(f"{d}/{file_prefix}_p_rho_x_by_item.pdf", verbose=False)
    p.save(f"{d}/{file_prefix}_p_rho_x_by_item.png", dpi=110, verbose=False)
    print(f"{name}: p(rho|x) plot ({df['n'].nunique()} interims x {len(items)} item) -> {d}")


if __name__ == "__main__":
    for d in DIRS:
        plot_dir(d)
