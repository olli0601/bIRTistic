"""Head-to-head LAIV-vs-TIV comparison for SDY269 (§3.21). Reads the ONE joint per-draw rho frame
written by IMMPORT_flu_headhead_svi.py (dir `py-immport-SDY269`,
`pcm_1_interim_i{k}_regression_training.pkl`): six rho_labels — LAIV_spr, TIV_spr, LAIV_gmfr,
TIV_gmfr, TIV-LAIV_spr_diff, TIV-LAIV_gmfr_diff — all indexed on the same Monte-Carlo draw. The
strain x rho_label contraction plot is made by IMMPORT_flu_rho_plots.py; this script writes the
numeric summaries into the same dir:
  headhead_final_summary.csv       per (rho_label, strain) median + 95% CI + P(H1) at full accrual
  headhead_shared_strain_diff.csv / .{pdf,png}   the two _diff rhos (TIV-LAIV) over interims
and echoes the composite CHMP 'at least one of SPR/GMFR met' PPS (headhead_any_met.csv). Because
both arms share the joint fit, the _diff rhos are EXACT per-draw contrasts (no random pairing)."""
# ---- boilerplate ----

import os, sys, glob, re
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
import numpy as np, pandas as pd

SB = os.environ.get('IMMPORT_SB', "/Users/or105/sandbox/bIRTistic")
DIR = f"{SB}/py-immport-SDY269"; OUT = DIR; fit_prefix = "pcm_1_interim"


def _interims():
    return sorted(int(re.search(r'_i(\d+)_', f).group(1))
                  for f in glob.glob(f"{DIR}/{fit_prefix}_i*_regression_training.pkl"))


def _n_per_arm(k):
    dp = pd.read_csv(f"{DIR}/{fit_prefix}_{k}_data_dp1.csv")
    return {a: g['pid'].nunique() for a, g in dp.groupby('arm')}


def load():
    """The joint per-draw rho frame across interims (all six rho_labels, shared draw index)."""
    rows = []
    for k in _interims():
        x = pd.read_pickle(f"{DIR}/{fit_prefix}_i{k}_regression_training.pkl").copy()
        na = _n_per_arm(k); x['interim'] = k; x['n_arm'] = min(na.values())
        rows.append(x)
    return pd.concat(rows, ignore_index=True)


def summarise_final(df):
    k = df['interim'].max(); g = df[df.interim == k]
    rows = []
    for (rl, rll, item, arm), gg in g.groupby(['rho_label', 'rho_label_long', 'item_label', 'arm']):
        q = gg['pps_rho_x'].quantile([0.025, 0.5, 0.975]).to_numpy()
        rows.append(dict(rho_label=rl, rho_label_long=rll, strain=item, arm=arm,
                         median=q[1], lo=q[0], hi=q[2], p_H1=float(gg['pps_H1_x'].mean())))
    return pd.DataFrame(rows).sort_values(['rho_label', 'strain'])


def diffs_summary_plot(df):
    import ggsci, matplotlib.colors as mcolors
    from plotnine import (ggplot, aes, geom_boxplot, geom_hline, facet_wrap, theme_bw, theme,
                          labs, scale_fill_manual, element_text)
    d = df[df['rho_label'].str.contains('diff')].copy()
    rows_s, rows_b = [], []
    for (rl, rll, n), g in d.groupby(['rho_label', 'rho_label_long', 'n_arm']):
        q = g['pps_rho_x'].quantile([0.025, 0.25, 0.5, 0.75, 0.975]).to_numpy()
        rows_b.append(dict(rho_label=rl, rho_label_long=rll, n=n, ymin=q[0], lower=q[1],
                           middle=q[2], upper=q[3], ymax=q[4]))
        rows_s.append(dict(rho_label=rl, n_arm=n, LAIV_median=float(g['group1'].median()),
                           TIV_median=float(g['group2'].median()), diff_median=q[2],
                           diff_lo=q[0], diff_hi=q[4], P_TIV_gt_LAIV=float((g['pps_rho_x'] > 0).mean())))
    pd.DataFrame(rows_s).sort_values(['rho_label', 'n_arm']).to_csv(
        f"{OUT}/headhead_shared_strain_diff.csv", index=False)
    b = pd.DataFrame(rows_b); ns = sorted(b['n'].unique())
    b['n_f'] = pd.Categorical(b['n'].astype(str), categories=[str(x) for x in ns], ordered=True)
    pal = [mcolors.to_hex(c) for c in ggsci.pal_futurama("planetexpress")(12)]
    cdict = {p: pal[i] for i, p in enumerate(sorted(b['rho_label_long'].unique()))}
    p = (ggplot(b, aes(x='n_f'))
         + geom_hline(yintercept=0, linetype='dashed', color='#555555')     # 0 = arms equal
         + geom_boxplot(aes(ymin='ymin', lower='lower', middle='middle', upper='upper',
                            ymax='ymax', fill='rho_label_long'), stat='identity', alpha=0.85, width=0.6)
         + facet_wrap('rho_label_long', ncol=2, scales='free_y')
         + scale_fill_manual(values=cdict, guide=None)
         + theme_bw()
         + theme(axis_text_x=element_text(size=8), strip_text=element_text(size=9),
                 figure_size=(11, 4.2))
         + labs(x='interim (paired participants per arm)',
                y='rho = TIV - LAIV  (shared strain A/Uruguay H3N2, exact per-draw contrast)',
                title='SDY269 head-to-head cross-arm rho on the shared strain (JOINT fit; '
                      'dashed = 0: arms equal; box = median/IQR, whiskers 2.5-97.5%)'))
    p.save(f"{OUT}/headhead_shared_strain_diff.pdf", verbose=False)
    p.save(f"{OUT}/headhead_shared_strain_diff.png", dpi=120, verbose=False)


if __name__ == "__main__":
    df = load()
    s = summarise_final(df); s.to_csv(f"{OUT}/headhead_final_summary.csv", index=False)
    diffs_summary_plot(df)
    pd.set_option('display.width', 200, 'display.max_columns', 20)
    print("=== final per-rho seroresponse (median [95% CI], P(H1)) ===")
    print(s[['rho_label', 'strain', 'arm', 'median', 'lo', 'hi', 'p_H1']].to_string(
        index=False, float_format=lambda v: f"{v:.3f}"))
    am = pd.read_csv(f"{OUT}/headhead_any_met.csv")
    fin = am[am.interim == am.interim.max()]
    print("\n=== composite CHMP 'at least one of SPR/GMFR met' at full accrual (joint per-draw, per strain x arm) ===")
    print(fin[['arm', 'item_label', 'pps_any_met']].sort_values(['arm', 'item_label'])
          .to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(f"\n-> {OUT}")
