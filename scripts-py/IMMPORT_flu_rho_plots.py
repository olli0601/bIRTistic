"""p(rho|x) contraction plots for the ImmPort HAI applications (§3.18), like
py-icrc-dass-drc_260902/pcm_1_interim_p_rho_x_by_item.pdf but with the TWO HAI endpoints
(SPR, GMFR) as facet COLUMNS. Reads the per-interim i{k}_regression_training.pkl (draw x
item x rho), summarises the posterior of each rho as a boxplot across draws at each interim,
and facets strain (row) x rho measure (column, free y-scale: SPR is a probability, GMFR a
fold). Dashed line = the clinical H1 threshold per rho (SPR>0.70, GMFR>2.5). Emits
pcm_1_interim_p_rho_x_by_item.{csv,pdf,png} in each study dir."""
# ---- boilerplate ----

import os, sys, glob, re
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
import numpy as np, pandas as pd

SB = os.environ.get('IMMPORT_SB', "/Users/or105/sandbox/bIRTistic")
DEFAULT_PREFIX = "pcm_1_interim"
# default = §3.18 SDY312/314 (one dir, prefix pcm_1_interim). Override with
# IMMPORT_RHO_DIRS='label:subdir[:file_prefix],...'; an explicit prefix (e.g. the SDY269
# per-endpoint-arm series pcm_SPR_TIV_interim) makes a single-rho contraction named by that prefix.
STUDIES = {'SDY312': ('py-immport-SDY312_260918', DEFAULT_PREFIX),
           'SDY314': ('py-immport-SDY314_260918', DEFAULT_PREFIX)}
# entry: 'label:subdir[:file_prefix[:dp_prefix]]'. dp_prefix (where the _{k}_data_dp1.csv giving the
# interim n lives) defaults to file_prefix; for a joint fit read out under several endpoint prefixes
# it is the single fit prefix (e.g. SDY269: file_prefix=pcm_SPR_interim, dp_prefix=pcm_1_interim).
STUDIES = {k: v + (v[1],) for k, v in STUDIES.items()}            # -> (subdir, file_prefix, dp_prefix)
if os.environ.get('IMMPORT_RHO_DIRS'):
    STUDIES = {}
    for x in os.environ['IMMPORT_RHO_DIRS'].split(','):
        parts = x.split(':'); lab, sub = parts[0], parts[1]
        fp = parts[2] if len(parts) > 2 else DEFAULT_PREFIX
        STUDIES[lab] = (sub, fp, parts[3] if len(parts) > 3 else fp)
# H1 threshold line per rho_label; unknown labels (e.g. the head-to-head _diff endpoints) -> 0
H1 = {'SPR: P(titre>=1:40) at endline': 0.70, 'GMT fold-rise (endline/baseline)': 2.5}
_h1 = lambda rl: H1.get(rl, 0.0)
RHO_ORDER = ['SPR: P(titre>=1:40) at endline', 'GMT fold-rise (endline/baseline)']
_hmult = float(os.environ.get('PROB_FIT_HEIGHT_MULT', '1.0'))


def build_summary(d, file_prefix=DEFAULT_PREFIX, dp_prefix=None):
    """Per (interim n, strain, rho) boxplot quantiles of the rho posterior across draws."""
    dp_prefix = dp_prefix or file_prefix
    fs = sorted(glob.glob(f"{d}/{file_prefix}_i*_regression_training.pkl"),
                key=lambda p: int(re.search(r'_i(\d+)_', p).group(1)))
    rows = []
    for f in fs:
        k = int(re.search(r'_i(\d+)_', f).group(1))
        n = pd.read_csv(f"{d}/{dp_prefix}_{k}_data_dp1.csv")['pid'].nunique()
        x = pd.read_pickle(f).copy()
        # disambiguate an item measured under >1 arm (e.g. the head-to-head shared strain, one
        # item_label per arm) so its trajectories do not pool; single-arm items keep the bare label
        if 'arm' in x.columns and x.groupby('item_label')['arm'].nunique().max() > 1:
            x['item_label'] = x['item_label'] + ' [' + x['arm'].astype(str) + ']'
        for (item, rid, rlab), g in x.groupby(['item_label', 'rho_id', 'rho_label']):
            q = g['pps_rho_x'].quantile([0.025, 0.25, 0.5, 0.75, 0.975]).to_numpy()
            rows.append(dict(n=n, item=item, rho_id=rid, rho_label=rlab,
                             ymin=q[0], lower=q[1], middle=q[2], upper=q[3], ymax=q[4]))
    return pd.DataFrame(rows)


def plot_study(study, d, file_prefix=DEFAULT_PREFIX, dp_prefix=None):
    import ggsci, matplotlib.colors as mcolors
    from plotnine import (ggplot, aes, geom_boxplot, geom_hline, facet_wrap, theme_bw, theme,
                          labs, scale_fill_manual, element_text)
    df = build_summary(d, file_prefix, dp_prefix)
    if df.empty:
        print(f"{study}: no endpoint pkls"); return
    df.to_csv(f"{d}/{file_prefix}_p_rho_x_by_item.csv", index=False)
    ns = sorted(df['n'].unique()); items = sorted(df['item'].unique())
    present = set(df['rho_label'])                                 # keep known order, append any extras
    rhos = [r for r in RHO_ORDER if r in present] + sorted(present - set(RHO_ORDER))
    df['n_f'] = pd.Categorical(df['n'].astype(str), categories=[str(x) for x in ns], ordered=True)
    df['item_f'] = pd.Categorical(df['item'], categories=items, ordered=True)
    df['rho_f'] = pd.Categorical(df['rho_label'], categories=rhos, ordered=True)
    # per-panel H1 threshold line (0 for the head-to-head _diff endpoints)
    thr = pd.DataFrame([(it, rl, _h1(rl)) for it in items for rl in rhos],
                       columns=['item_f', 'rho_f', 'yint'])
    thr['item_f'] = pd.Categorical(thr['item_f'], categories=items, ordered=True)
    thr['rho_f'] = pd.Categorical(thr['rho_f'], categories=rhos, ordered=True)
    pal = ggsci.pal_futurama("planetexpress")(12)
    cols = [mcolors.to_hex(c) for c in
            mcolors.LinearSegmentedColormap.from_list("f", pal)(np.linspace(0, 1, len(items)))]
    cdict = dict(zip(items, cols))
    p = (ggplot(df, aes(x='n_f'))
         + geom_boxplot(aes(ymin='ymin', lower='lower', middle='middle', upper='upper',
                            ymax='ymax', fill='item_f'), stat='identity', alpha=0.8, width=0.7)
         + geom_hline(thr, aes(yintercept='yint'), linetype='dashed', color='#555555')
         + facet_wrap(['item_f', 'rho_f'], ncol=len(rhos), scales='free_y')   # rows=strain, cols=rho
         + scale_fill_manual(values=cdict, guide=None)
         + theme_bw()
         + theme(axis_text_x=element_text(angle=45, vjust=1, hjust=1, size=7),
                 strip_text=element_text(size=8),
                 figure_size=(5.5 * len(rhos), max(14, 2.4 * len(items)) * _hmult))
         + labs(x='interim (cumulative participants)',
                y='p(rho | x)  —  SPR (prob) / GMFR (fold)',
                title=f'{study} HAI: SVI p(rho|x) contraction by strain and endpoint '
                      f'(dashed = H1: SPR>0.70, GMFR>2.5)'))
    p.save(f"{d}/{file_prefix}_p_rho_x_by_item.pdf", verbose=False)
    p.save(f"{d}/{file_prefix}_p_rho_x_by_item.png", dpi=110, verbose=False)
    print(f"{study}: {df['n'].nunique()} interims x {len(items)} strains x {len(rhos)} rho "
          f"[{file_prefix}] -> {d}")


if __name__ == "__main__":
    for study, (sub, prefix, dpp) in STUDIES.items():
        plot_study(study, f"{SB}/{sub}", prefix, dpp)
