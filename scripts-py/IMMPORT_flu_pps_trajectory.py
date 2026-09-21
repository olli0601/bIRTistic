"""PPS trajectory summary for the ImmPort HAI amortiser deploys (§3.18): amortised
P(H1|x) vs interim n, one line per strain, faceted by success threshold eta0. BOTH endpoints
live in the single per-study amortiser dir, distinguished by prefix (pcm_gmfr_interim_* /
pcm_spr_interim_*). Writes pcm_<rho>_interim_pps_RAGD_trajectory.{pdf,png} in that dir."""
# ---- boilerplate ----

import os, sys, glob
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
import numpy as np, pandas as pd

SB = os.environ.get('IMMPORT_SB', "/Users/or105/sandbox/bIRTistic")
DEPLOYS = sorted(glob.glob(f"{SB}/py-immport-SDY*-amortise-*bvm*"))
ENDPOINTS = {'pcm_gmfr_interim': ('GMFR', 2.5), 'pcm_spr_interim': ('SPR', 0.70)}
_hmult = float(os.environ.get('PROB_FIT_HEIGHT_MULT', '1.0'))


def plot_trajectory(d, prefix, ep, chmp):
    import ggsci, matplotlib.colors as mcolors
    from plotnine import (ggplot, aes, geom_line, geom_point, geom_hline, facet_wrap,
                          theme_bw, theme, labs, scale_color_manual, element_text,
                          scale_y_continuous)
    ecsv = f"{d}/{prefix}_pps_RAGD_eta0_pps.csv"
    if not os.path.exists(ecsv):
        return
    e = pd.read_csv(ecsv)
    nmap = (pd.read_csv(f"{d}/{prefix}_pps_RAGD_pps_by_item.csv")
            [['interim_id', 'n']].drop_duplicates().rename(columns={'interim_id': 'interim'}))
    e = e.merge(nmap, on='interim', how='left')
    e['thr'] = e['eta0_pct'] / 100.0
    _lab = lambda f: f"{ep} > {f:g}" + ("  (CHMP)" if abs(f - chmp) < 1e-6 else "")
    order = sorted(e['thr'].unique())
    e['thr_f'] = pd.Categorical(e['thr'].map(_lab), categories=[_lab(f) for f in order], ordered=True)
    items = sorted(e['item'].unique())
    pal = ggsci.pal_futurama("planetexpress")(12)
    cols = [mcolors.to_hex(c) for c in
            mcolors.LinearSegmentedColormap.from_list("f", pal)(np.linspace(0, 1, len(items)))]
    cdict = dict(zip(items, cols))
    study = os.path.basename(d).split('-amortise')[0].replace('py-immport-', '')
    p = (ggplot(e, aes('n', 'pps', color='item'))
         + geom_hline(yintercept=[0.1, 0.9], linetype='dashed', color='#999999')
         + geom_line() + geom_point(size=1.4)
         + facet_wrap('thr_f', ncol=len(order))
         + scale_color_manual(values=cdict, name='strain')
         + scale_y_continuous(limits=[0, 1], labels=lambda l: [f'{v:.0%}' for v in l])
         + theme_bw()
         + theme(legend_position='bottom', strip_text=element_text(size=8),
                 axis_text_x=element_text(size=7), figure_size=(3.2 * len(order), 4.2 * _hmult))
         + labs(x='interim (cumulative participants)', y='amortised PPS = P(H1 | x)',
                title=f'{study} HAI: amortised PPS trajectory for {ep} '
                      f'(J64 deepsetXcompAtt + head BvM; dashed 10%/90% go/no-go guides)'))
    p.save(f"{d}/{prefix}_pps_RAGD_trajectory.pdf", verbose=False)
    p.save(f"{d}/{prefix}_pps_RAGD_trajectory.png", dpi=110, verbose=False)
    print(f"{study} {ep}: PPS trajectory -> {prefix}_pps_RAGD_trajectory")


if __name__ == "__main__":
    for d in DEPLOYS:
        for prefix, (ep, chmp) in ENDPOINTS.items():
            plot_trajectory(d, prefix, ep, chmp)
