"""Head-to-head LAIV-vs-TIV seroresponse comparison for SDY269 (§3.21), JOINT-FIT version.
Reads the four endpoint posteriors written by IMMPORT_flu_headhead_svi.py off the ONE joint fit
in dir `py-immport-SDY269`:
  pcm_SPR_interim_ / pcm_GMFR_interim_            per (strain, arm) level endpoints
  pcm_SPR_diff_interim_ / pcm_GMFR_diff_interim_  TIV-minus-LAIV on the SHARED strain (per draw)
Produces (into the same dir):
  headhead_final_summary.csv     per (role, strain, arm, endpoint) median + 95% CI + P(H1) at full accrual
  headhead_compare.{pdf,png}     facet endpoint x role, x=arm, posterior boxplots + H1 line
  headhead_shared_strain_diff.csv / .{pdf,png}   the two _diff endpoints (TIV-LAIV) over interims
Because both arms share the joint fit, the _diff is an EXACT per-draw contrast (no random pairing)."""
# ---- boilerplate ----

import os, sys, glob, re
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
import numpy as np, pandas as pd

SB = os.environ.get('IMMPORT_SB', "/Users/or105/sandbox/bIRTistic")
DIR = f"{SB}/py-immport-SDY269"; OUT = DIR; fit_prefix = "pcm_1_interim"
H1 = {1: 0.70, 2: 2.5}                                  # SPR>0.70, GMFR>2.5 (levels); diffs use >0
RHO_ORDER = ['SPR: P(titre>=1:40) at endline', 'GMT fold-rise (endline/baseline)']
ROLE_ORDER = ['H1N1', 'H3N2 (shared)', 'B']


def _role(strain):
    s = str(strain)
    if 'H1N1' in s: return 'H1N1'
    if 'H3N2' in s: return 'H3N2 (shared)'
    if s.startswith('B/'): return 'B'
    return s


def _interims(tag):
    return sorted(int(re.search(r'_i(\d+)_', f).group(1))
                  for f in glob.glob(f"{DIR}/pcm_{tag}_interim_i*_regression_training.pkl"))


def _n_per_arm(k):
    dp = pd.read_csv(f"{DIR}/{fit_prefix}_{k}_data_dp1.csv")
    return {a: g['pid'].nunique() for a, g in dp.groupby('arm')}


def load(tags):
    """Concat the given endpoint tags across all interims -> one long frame with interim n."""
    rows = []
    for tag in tags:
        for k in _interims(tag):
            x = pd.read_pickle(f"{DIR}/pcm_{tag}_interim_i{k}_regression_training.pkl").copy()
            na = _n_per_arm(k); x['interim'] = k; x['n_arm'] = min(na.values())
            x['strain'] = x['item_label']                    # item_label IS the strain (pure item)
            x['role'] = x['strain'].map(_role); x['tag'] = tag
            rows.append(x)
    return pd.concat(rows, ignore_index=True)


def summarise_final(levels):
    k = levels['interim'].max(); g = levels[levels.interim == k]
    rows = []
    for (arm, role, strain, rid, rlab), gg in g.groupby(['arm', 'role', 'strain', 'rho_id', 'rho_label']):
        q = gg['pps_rho_x'].quantile([0.025, 0.5, 0.975]).to_numpy()
        rows.append(dict(arm=arm, role=role, strain=strain, rho_id=rid, rho_label=rlab,
                         median=q[1], lo=q[0], hi=q[2], p_H1=float((gg['pps_rho_x'] > H1[rid]).mean())))
    return pd.DataFrame(rows).sort_values(['role', 'rho_id', 'arm'])


def plot_levels(levels):
    import ggsci, matplotlib.colors as mcolors
    from plotnine import (ggplot, aes, geom_boxplot, geom_hline, facet_grid, theme_bw, theme,
                          labs, scale_fill_manual, element_text)
    k = levels['interim'].max(); g = levels[levels.interim == k]
    rows = []
    for (role, rid, rlab, arm, strain), gg in g.groupby(['role', 'rho_id', 'rho_label', 'arm', 'strain']):
        q = gg['pps_rho_x'].quantile([0.025, 0.25, 0.5, 0.75, 0.975]).to_numpy()
        rows.append(dict(role=role, rho_id=rid, rho_label=rlab, arm=arm, strain=strain,
                         ymin=q[0], lower=q[1], middle=q[2], upper=q[3], ymax=q[4]))
    b = pd.DataFrame(rows)
    b['x'] = b['arm'] + '\n' + b['strain'].map(lambda s: re.sub(r'\(.*\)', '', str(s)).strip())
    b['role_f'] = pd.Categorical(b['role'], categories=ROLE_ORDER, ordered=True)
    b['rho_f'] = pd.Categorical(b['rho_label'], categories=RHO_ORDER, ordered=True)
    thr = pd.DataFrame([(rl, H1[rid]) for rid, rl in zip([1, 2], RHO_ORDER)], columns=['rho_f', 'yint'])
    thr['rho_f'] = pd.Categorical(thr['rho_f'], categories=RHO_ORDER, ordered=True)
    pal = [mcolors.to_hex(c) for c in ggsci.pal_futurama("planetexpress")(12)]
    cdict = {'LAIV': pal[0], 'TIV': pal[3]}
    na = _n_per_arm(k)
    p = (ggplot(b, aes(x='x'))
         + geom_boxplot(aes(ymin='ymin', lower='lower', middle='middle', upper='upper',
                            ymax='ymax', fill='arm'), stat='identity', alpha=0.85, width=0.7)
         + geom_hline(thr, aes(yintercept='yint'), linetype='dashed', color='#555555')
         + facet_grid('rho_f ~ role_f', scales='free')
         + scale_fill_manual(values=cdict, name='vaccine arm')
         + theme_bw()
         + theme(axis_text_x=element_text(angle=45, vjust=1, hjust=1, size=7),
                 strip_text=element_text(size=8), figure_size=(12, 8))
         + labs(x='vaccine arm / strain (H3N2 = shared strain)',
                y='p(rho | x)  —  SPR (prob) / GMFR (fold)',
                title=f"SDY269 head-to-head (JOINT fit): LAIV (n={na.get('LAIV')}) vs TIV "
                      f"(n={na.get('TIV')}) HAI seroresponse at full accrual (dashed = SPR>0.70, GMFR>2.5)"))
    p.save(f"{OUT}/headhead_compare.pdf", verbose=False)
    p.save(f"{OUT}/headhead_compare.png", dpi=120, verbose=False)


def diffs_summary_plot(diffs):
    import ggsci, matplotlib.colors as mcolors
    from plotnine import (ggplot, aes, geom_boxplot, geom_hline, facet_wrap, theme_bw, theme,
                          labs, scale_fill_manual, element_text)
    ep = {3: 'SPR_diff', 4: 'GMFR_diff'}
    rows_s, rows_b = [], []
    for (rid, n), g in diffs.groupby(['rho_id', 'n_arm']):
        q = g['pps_rho_x'].quantile([0.025, 0.25, 0.5, 0.75, 0.975]).to_numpy()
        rows_b.append(dict(endpoint=ep[rid], n=n, ymin=q[0], lower=q[1], middle=q[2],
                           upper=q[3], ymax=q[4]))
        rows_s.append(dict(endpoint=ep[rid], n_arm=n, LAIV_median=float(g['group1'].median()),
                           TIV_median=float(g['group2'].median()), diff_median=q[2],
                           diff_lo=q[0], diff_hi=q[4], P_TIV_gt_LAIV=float((g['pps_rho_x'] > 0).mean())))
    pd.DataFrame(rows_s).sort_values(['endpoint', 'n_arm']).to_csv(
        f"{OUT}/headhead_shared_strain_diff.csv", index=False)
    b = pd.DataFrame(rows_b); ns = sorted(b['n'].unique())
    b['n_f'] = pd.Categorical(b['n'].astype(str), categories=[str(x) for x in ns], ordered=True)
    unit = {'SPR_diff': 'SPR rate difference', 'GMFR_diff': 'GMFR fold difference'}
    b['panel'] = b['endpoint'] + '  (' + b['endpoint'].map(unit) + ')'
    pal = [mcolors.to_hex(c) for c in ggsci.pal_futurama("planetexpress")(12)]
    cdict = {p: pal[i] for i, p in enumerate(sorted(b['panel'].unique()))}
    p = (ggplot(b, aes(x='n_f'))
         + geom_hline(yintercept=0, linetype='dashed', color='#555555')     # 0 = arms equal
         + geom_boxplot(aes(ymin='ymin', lower='lower', middle='middle', upper='upper',
                            ymax='ymax', fill='panel'), stat='identity', alpha=0.85, width=0.6)
         + facet_wrap('panel', ncol=2, scales='free_y')
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
    levels = load(['SPR', 'GMFR']); diffs = load(['SPR_diff', 'GMFR_diff'])
    s = summarise_final(levels); s.to_csv(f"{OUT}/headhead_final_summary.csv", index=False)
    plot_levels(levels); diffs_summary_plot(diffs)
    pd.set_option('display.width', 180, 'display.max_columns', 20)
    print("=== per-arm final seroresponse (median [95% CI], P(H1)) ===")
    print(s.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print("\n=== shared strain A/Uruguay H3N2: exact per-draw TIV-LAIV difference ===")
    print(pd.read_csv(f"{OUT}/headhead_shared_strain_diff.csv").to_string(
        index=False, float_format=lambda v: f"{v:.3f}"))
    print(f"\n-> {OUT}")
