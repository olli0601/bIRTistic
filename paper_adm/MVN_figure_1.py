#!/usr/bin/env python3
"""Composite results figure for the MVN benchmark (§13).

Assembles the standalone MVN diagnostic panels into ONE multi-panel figure:

    +-------------------------------+-------------------+
    | (A) per-component recovery    | (B) empirical cov |   <- top part
    |     effect size x component j |     matrix J=100  |
    +-------------------------------+-------------------+
    | (C) PPS heatmap  J=20 | J=100  (full width)       |   <- more below later
    +---------------------------------------------------+

The panels are built with plotnine from the saved simulation pickles (no re-simulation),
rendered to PNG at their FINAL composite size (so point font sizes are consistent across
panels -> homogeneous), then placed with matplotlib imshow at equal aspect (no warping).

Run:  pixi run python paper_adm/MVN_figure_1.py   (figure -> sandbox/paper-adm/mvn_figure_1.{pdf,png})
"""
import os
import sys
import glob
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'python'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from plotnine import (
    ggplot, aes, geom_tile, geom_errorbarh, geom_vline, geom_point,
    geom_boxplot, geom_hline, geom_col, geom_text, position_dodge, coord_flip,
    facet_grid, facet_wrap, coord_equal, scale_x_continuous, scale_x_discrete,
    scale_y_continuous, scale_y_discrete,
    scale_fill_gradient2, scale_fill_cmap, scale_fill_manual, scale_colour_manual, theme_bw, theme,
    element_text, element_blank, labs, guides, guide_legend,
)
from mizani.transforms import trans_new
from utils import _futurama_palette                             # compare_methods method palette
from amortiser_io import load_interim_data

# ---- shared method universe, short labels, palette + the signed-sqrt timing axis ------
_SHORT = {
    'analytic': 'analytic (closed-form)',
    'nested-MC using HMC for each (x,z)': 'nested-SVI\nof theta',
    'IS reweighting of theta|x': 'importance-sampling\nof theta',
    'Regression of endpt-x on w(z) - mquantile': 'rho-regression\n(known features)',
    'Amortiser-features-fixed-idcomp-qpsi-MLP-loss-multiquantilehead': 'ADM (known features,\nfixed N and J)',
    'Amortiser-features-MLP-xcomp-qpsi-MLP-loss-multiquantilehead': 'ADM (MLP features,\nfixed N and J)',
    'deepsetXcompAtt plain + affine (§14.4.6)': 'ADM (XAttention features,\nany N, any J)',
    'deepsetXcompAtt A power-law + BvM (§14.4.7)': 'ADM (XAttention features,\nBvM architecture, any N, any J)',
    # removed from panels (labels retained for completeness / other tooling):
    'Regression of endpt-x on w(z) using Gaussian approx': 'regression (Gauss)',
    'Regression of endpt-x on w(z) - quantile': 'regression (quantile)',
    'Amortiser-features-itemScompAtt-qpsi-MLP-loss-multiquantilehead': 'itemScompAtt',
    'Amortiser-features-itemXcompAtt-qpsi-MLP-loss-multiquantilehead': 'itemXcompAtt',
    'Amortiser-features-deepsetScompAtt-qpsi-MLP-loss-multiquantilehead': 'deepsetScompAtt',
    'Amortiser-features-deepsetXcompAtt-qpsi-MLP-loss-multiquantilehead': 'deepsetXcompAtt',
    'deepsetXcompAtt C floor + BvM (§14.4.7)': 'deepset C floor+BvM',
}
_ALL_METHODS = [                                               # matches compare_methods palette order
    'nested-MC using HMC for each (x,z)', 'IS reweighting of theta|x',
    'Regression of endpt-x on w(z) using Gaussian approx',
    'Regression of endpt-x on w(z) - quantile', 'Regression of endpt-x on w(z) - mquantile',
    'Amortiser-features-fixed-idcomp-qpsi-MLP-loss-multiquantilehead',
    'Amortiser-features-MLP-xcomp-qpsi-MLP-loss-multiquantilehead',
    'Amortiser-features-itemScompAtt-qpsi-MLP-loss-multiquantilehead',
    'Amortiser-features-itemXcompAtt-qpsi-MLP-loss-multiquantilehead',
    'Amortiser-features-deepsetScompAtt-qpsi-MLP-loss-multiquantilehead',
    'Amortiser-features-deepsetXcompAtt-qpsi-MLP-loss-multiquantilehead',
    'deepsetXcompAtt plain + affine (§14.4.6)', 'deepsetXcompAtt A power-law + BvM (§14.4.7)',
    'deepsetXcompAtt C floor + BvM (§14.4.7)',
]
_METHOD_COLOURS = dict(zip(_ALL_METHODS, _futurama_palette(len(_ALL_METHODS))))
_REMOVE = {                                                    # dropped from panels d/e/f
    'Regression of endpt-x on w(z) using Gaussian approx',
    'Regression of endpt-x on w(z) - quantile',
    'Amortiser-features-MLP-xcomp-qpsi-MLP-loss-multiquantilehead',
    'Amortiser-features-itemScompAtt-qpsi-MLP-loss-multiquantilehead',
    'Amortiser-features-itemXcompAtt-qpsi-MLP-loss-multiquantilehead',
    'Amortiser-features-deepsetScompAtt-qpsi-MLP-loss-multiquantilehead',
    'Amortiser-features-deepsetXcompAtt-qpsi-MLP-loss-multiquantilehead',
    'deepsetXcompAtt C floor + BvM (§14.4.7)',
}
_PANEL_METHODS = [m for m in _ALL_METHODS if m not in _REMOVE]  # kept methods for d/e/f
_signed_sqrt = trans_new(
    'signed_sqrt',
    lambda x: np.sign(x) * np.sqrt(np.abs(x)),
    lambda x: np.sign(x) * np.asarray(x, dtype=float) ** 2,
)

warnings.filterwarnings('ignore')

SB = os.environ.get('MVN_SB', "/Users/or105/sandbox/bIRTistic")
DIR_SIM = os.path.join(SB, "py-mvn-interim-simulations-260609")
DIR_CMP = os.path.join(SB, "py-mvn-interim-compare-methods-260609")
SCRATCH = os.environ.get('MVN_FIG_SCRATCH',
                         "/private/tmp/claude-501/-Users-or105-git-bIRTistic/"
                         "49372e22-d11b-4882-a442-d1c60bcbdfb0/scratchpad")
os.makedirs(SCRATCH, exist_ok=True)
DPI = 200

# --- shared style so every panel's fonts match (homogeneous) ---------------------------
BASE = (theme_bw()
        + theme(text=element_text(size=11),
                axis_title=element_text(size=11),
                axis_text=element_text(size=8),
                strip_text=element_text(size=9, face='bold'),
                strip_background=element_blank(),
                plot_title=element_text(size=13),
                legend_title=element_text(size=10),
                legend_text=element_text(size=8),
                panel_grid_major=element_blank(),
                panel_grid_minor=element_blank()))

# --- load the saved diagnostics (no re-simulation) -------------------------------------
sim = pd.read_pickle(os.path.join(DIR_SIM, 'mvn_sim_data.pkl'))
MU0 = float(sim['simu_params']['mu_0_baseline'])                 # 1.0 -> the "no effect" line
J_GRID = list(sim['simu_params']['J_grid'])
cf = pd.read_pickle(os.path.join(DIR_SIM, 'mvn_pps_closed_form.pkl'))
post, pps_cf = cf['p_h1_x'].copy(), cf['pps_cf'].copy()
diag_a = pd.read_pickle(os.path.join(DIR_SIM, 'mvn_diag_per_component.pkl'))
heat = pd.read_pickle(os.path.join(DIR_SIM, 'mvn_diag_cov_heatmap.pkl'))

interim_order = (post[['interim_date', 'interim_month_year']].drop_duplicates()
                 .sort_values('interim_date')['interim_month_year'].tolist())
_io_d = interim_order                                          # shared interim (x) order across panels
_N = post.drop_duplicates('interim_month_year').set_index('interim_month_year')['n_obs'].to_dict()


def _ylab_n(cats):                                              # interim y label with participant count
    return [f"{c}  (n={int(_N[c])})" if c in _N else str(c) for c in cats]


def _xlab_n(cats):                                             # interim x label: date + participant count
    return [f"{c}\nn={int(_N[c])}" if c in _N else str(c) for c in cats]


def _render(p, name, w, h):                                    # ggplot -> PNG at final size, return image
    path = os.path.join(SCRATCH, name)
    p.save(path, width=w, height=h, dpi=DPI, verbose=False, limitsize=False)
    return np.asarray(Image.open(path))


# =======================================================================================
# Panel A -- per-component recovery. Edits: red mu_true points removed; error bars +
# posterior mean in dark grey; empirical mean black; dashed "no effect" line at x = MU0.
# =======================================================================================
_PC_J = [J for J in sorted(J_GRID) if J in (20, 100)]
dA = diag_a[diag_a['J'].isin(_PC_J)].copy()
dA['J_label'] = pd.Categorical('J=' + dA['J'].astype(str),
                               categories=[f'J={J}' for J in _PC_J], ordered=True)
_pc_lab = {c: f"Interim\n{c}\n(n={int(_N[c])})" for c in interim_order}
dA['interim_lab'] = pd.Categorical(dA['interim_month_year'].map(_pc_lab),
                                   categories=[_pc_lab[c] for c in interim_order], ordered=True)
pA = (
    ggplot(dA, aes(y='j'))
    + geom_vline(xintercept=MU0, linetype='dashed', colour='#404040', size=0.4)   # no-effect line (x=1)
    + geom_errorbarh(aes(xmin='crI_lo', xmax='crI_hi'), colour='#4d4d4d', height=0.4, size=0.35)
    + geom_point(aes(x='mu_post'), colour='#4d4d4d', size=1.2)                     # posterior mean (dark grey)
    + geom_point(aes(x='y_bar'), colour='black', size=0.8)                         # empirical mean (black)
    + facet_grid('J_label ~ interim_lab', scales='free_y', space='free_y')
    + scale_y_continuous(expand=(0, 0.6))
    + BASE + theme(panel_spacing=0.012)
    + labs(x='effect size', y='component j')
)

# =======================================================================================
# Panel B -- empirical J=100 covariance matrix (square).
# =======================================================================================
hB = heat[(heat['J'] == 100) & (heat['method'] == 'empirical')]
pB = (
    ggplot(hB, aes(x='j', y='i', fill='value'))
    + geom_tile()
    + scale_fill_gradient2(low='#3b4cc0', mid='white', high='#b40426', midpoint=0.0)
    + scale_x_continuous(expand=(0, 0)) + scale_y_continuous(expand=(0, 0))
    + coord_equal()
    + BASE
    + labs(x='component j', y='component i', fill='cov(y_i, y_j)',
           title='empirical covariance (J=100)')
)

# =======================================================================================
# Panel C -- PPS heatmap, J=20 (left, narrow) | J=100 (right, wide); no whitespace.
# =======================================================================================
_PPS_J = [J for J in sorted(J_GRID) if J in (20, 100)]
dC = pps_cf[pps_cf['J'].isin(_PPS_J)].copy()
dC['J_label'] = pd.Categorical('J=' + dC['J'].astype(str),
                               categories=[f'J={J}' for J in _PPS_J], ordered=True)
_c_lab = {c: f"{c}\nn={int(_N[c])}" for c in interim_order}
# first interim at the TOP, last at the bottom -> reverse the factor levels (first level = bottom)
dC['interim_lab'] = pd.Categorical(dC['interim_month_year'].map(_c_lab),
                                   categories=[_c_lab[c] for c in reversed(interim_order)], ordered=True)
pC = (
    ggplot(dC, aes(x='j', y='interim_lab', fill='pps'))
    + geom_tile()
    + facet_grid('. ~ J_label', scales='free_x', space='free_x')
    + scale_x_continuous(expand=(0, 0)) + scale_y_discrete(expand=(0, 0))
    + scale_fill_gradient2(low='#009392', mid='#F1EAC8', high='#A5006A',
                           midpoint=0.5, limits=[0.0, 1.0],
                           breaks=[0.0, 0.25, 0.5, 0.75, 1.0],
                           labels=['0%', '25%', '50%', '75%', '100%'])
    + BASE + theme(panel_spacing=0.02)
    + labs(x='component j', y='interim analysis',
           fill='predictive\nprobability\nof success\n(analytic solution)')
)

# =======================================================================================
# Panel D -- p(H1 | x, z) boxplots for three components (26, 51, 75) at J=100, across
# interims, for the "all-amortised" method set (from mvn_J100_..._p_h1_xz_all_amortised.pdf).
# Facet columns = component; full width.
# =======================================================================================
_ETA = float(sim['simu_params'].get('pps_ProbH1_target_lwr_quantile',
                                    sim['simu_params'].get('pps_ProbH1_thresh', 0.89)))
# panel D shows the per-z p(H1|x,z) DISTRIBUTION; IS is excluded here because its importance
# weights degenerate (ESS~1) so every per-z value collapses to a hard 0/1 -> a boxplot of that
# Bernoulli spans [0,1] and misrepresents it. IS stays in panels E/F (its PPS/MSE are valid).
_D_DROP = {'IS reweighting of theta|x'}
_AMORT_D = ['analytic'] + [m for m in _PANEL_METHODS if m not in _D_DROP]  # + analytic ref
_D_SHORT = {m: _SHORT[m] for m in _AMORT_D}
_D_COLOURS = {_D_SHORT[m]: _METHOD_COLOURS.get(m, '#FFFFFF') for m in _AMORT_D}  # analytic -> white (transparent) ref
_D_OUTLINE = {_D_SHORT[m]: ('#000000' if m == 'analytic' else '#404040') for m in _AMORT_D}  # analytic = black contour
COMPS = [26, 51, 75]
box = pd.read_pickle(os.path.join(DIR_CMP, 'mvn_J100_compare_methods_p_h1_xz_all.pkl'))
box = box[box['response_label'].astype(str).isin([f'response mu_{c}' for c in COMPS])
          & box['method'].astype(str).isin(_AMORT_D)].copy()
box['component'] = pd.Categorical(
    box['response_label'].astype(str).map(lambda s: f"component {s.split('_')[-1]}"),
    categories=[f'component {c}' for c in COMPS], ordered=True)
box['method'] = pd.Categorical(box['method'].astype(str).map(_D_SHORT),
                               categories=[_D_SHORT[m] for m in _AMORT_D], ordered=True)
box['interim_month_year'] = pd.Categorical(box['interim_month_year'], categories=_io_d, ordered=True)
# order the dodge grouping by (interim, method) so bars follow the legend/colour order
_glev = [f"{it}||{_D_SHORT[m]}" for it in _io_d for m in _AMORT_D]
box['grp'] = box['interim_month_year'].astype(str) + '||' + box['method'].astype(str)
box['grp'] = pd.Categorical(box['grp'], categories=[g for g in _glev if g in set(box['grp'])],
                            ordered=True)
pD = (
    ggplot(box, aes(x='interim_month_year', fill='method', colour='method'))
    + geom_boxplot(aes(ymin='q025', lower='q25', middle='q50', upper='q75', ymax='q975', group='grp'),
                   stat='identity', position=position_dodge(width=0.8, preserve='single'),
                   width=0.72, size=0.25)
    + geom_hline(yintercept=_ETA, colour='black', size=0.8)
    + facet_wrap('~ component', ncol=3)
    + scale_fill_manual(values=_D_COLOURS, name='method')
    + scale_colour_manual(values=_D_OUTLINE, name='method')     # analytic drawn as transparent + black contour
    + scale_y_continuous(limits=[0, 1], expand=(0, 0, 0.02, 0),  # small top pad so degenerate boxes at 1.0 show
                         breaks=[0, 0.25, 0.5, 0.75, 1.0],
                         labels=['0%', '25%', '50%', '75%', '100%'])
    + scale_x_discrete(labels=_xlab_n)
    + BASE + theme(axis_text_x=element_text(size=7, angle=45, vjust=1, hjust=1),
                   legend_position='right', legend_text=element_text(ma='left'), panel_spacing=0.03)
    + labs(x='interim analysis', y=r'$p(H_1^{\,j} \mid x,\, z^{s})$', fill='method', colour='method')
)

# =======================================================================================
# Panel E -- MSE heatmap (method x interim), full width; viridis. A subset of methods.
# =======================================================================================
DIR_CMP = os.path.join(SB, "py-mvn-interim-compare-methods-260609")
REMOVE_METHODS = _REMOVE                                        # drop the 7 methods from panel e
mse = pd.read_pickle(os.path.join(DIR_CMP, 'mvn_compare_methods_mse.pkl'))
_meth_order = [str(m) for m in mse['method'].cat.categories if str(m) not in REMOVE_METHODS]
mse = mse[mse['method'].astype(str).isin(_meth_order)].copy()
# keep J=20 | J=100 as two facet columns (not collapsed); log10 fill for nuance.
dD = mse[mse['J'].isin([20, 100])].copy()
dD['J_label'] = pd.Categorical('J=' + dD['J'].astype(str),
                               categories=['J=20', 'J=100'], ordered=True)
_io_d = (dD[['interim_date', 'interim_month_year']].drop_duplicates()
         .sort_values('interim_date')['interim_month_year'].tolist())
dD['interim_month_year'] = pd.Categorical(dD['interim_month_year'], categories=_io_d, ordered=True)
dD['method'] = pd.Categorical(dD['method'].astype(str),
                              categories=list(reversed(_meth_order)), ordered=True)  # HMC at top (match panel f)
dD['mse_c'] = dD['mse'].clip(lower=1e-4)                        # floor for log10
pE = (
    ggplot(dD, aes(x='interim_month_year', y='method', fill='mse_c'))
    + geom_tile()
    + facet_grid('. ~ J_label')
    + scale_x_discrete(expand=(0, 0), labels=_xlab_n)
    + scale_y_discrete(expand=(0, 0), labels=lambda ms: [_SHORT.get(m, m) for m in ms])
    + scale_fill_cmap(cmap_name='viridis', trans='log10')
    + BASE + theme(axis_text_x=element_text(size=7, angle=45, vjust=1, hjust=1),
                   axis_text_y=element_text(ma='left', ha='right'), panel_spacing=0.02)
    + labs(x='interim analysis', y='method',
           fill='MSE (log10)\n(mean over component-\nspecific squared errors)')
)

# =======================================================================================
# Panel F -- timing pyramid: one-off TRAIN (blue, LEFT) vs total DEPLOY summed over interims
# (orange, RIGHT), one method per row, J=20 (left facet) | J=100 (right facet). Signed-sqrt
# x so large walltimes keep length; bars emanate from 0; positive minute tick labels.
# =======================================================================================
_ts = pd.read_pickle(os.path.join(DIR_CMP, 'mvn_compare_methods_timing_summary.pkl'))
_ts = _ts[_ts['J'].isin([20, 100]) & _ts['method'].isin(_PANEL_METHODS)].copy()
_trows = []
for _, _r in _ts.iterrows():
    _trows.append({'J': int(_r['J']), 'method': _r['method'], 'segment': 'train (one-off)',
                   'signed': -float(_r['train']), 'mins': float(_r['train'])})
    _trows.append({'J': int(_r['J']), 'method': _r['method'], 'segment': 'deploy (all interims)',
                   'signed': float(_r['deploy']), 'mins': float(_r['deploy'])})
tF = pd.DataFrame(_trows)
tF['mlabel'] = pd.Categorical(
    tF['method'].map(_SHORT),
    categories=list(reversed([_SHORT[m] for m in _PANEL_METHODS])), ordered=True)
tF['J_label'] = pd.Categorical('J=' + tF['J'].astype(str), categories=['J=20', 'J=100'], ordered=True)
tF['segment'] = pd.Categorical(tF['segment'],
                               ['train (one-off)', 'deploy (all interims)'], ordered=True)
tF['vlabel'] = tF['mins'].map(lambda v: f"{v:.3g}" if v > 0 else "")  # minutes printed outside bars
_fmx = float(tF['mins'].max()) if len(tF) else 1.0
_fcand = [1, 3, 10, 30, 100, 300, 1000]
_fbrk = sorted(set([-c for c in _fcand if c <= _fmx * 1.3] + [0.0]
                   + [c for c in _fcand if c <= _fmx * 1.3]))
_TD_COL = {'train (one-off)': '#1f77b4', 'deploy (all interims)': '#ff7f0e'}
_fL = tF[tF['signed'] < 0]                                       # train -> label to the LEFT
_fR = tF[tF['signed'] > 0]                                       # deploy -> label to the RIGHT
pF = (
    ggplot(tF, aes(x='mlabel', y='signed', fill='segment'))
    + geom_col(width=0.72)
    + geom_hline(yintercept=0, colour='#666666', size=0.4)
    + geom_text(_fL, aes(label='vlabel'), ha='right', va='center', size=6, colour='black')
    + geom_text(_fR, aes(label='vlabel'), ha='left', va='center', size=6, colour='black')
    + facet_grid('. ~ J_label')
    + scale_fill_manual(values=_TD_COL, breaks=list(_TD_COL), name='')
    + scale_y_continuous(trans=_signed_sqrt, breaks=_fbrk,
                         labels=lambda bs: [f"{abs(b):g}" for b in bs],
                         expand=(0.12, 0, 0.12, 0))
    + coord_flip()
    + BASE + theme(axis_text_x=element_text(size=7),
                   axis_text_y=element_text(ma='left', ha='right'), panel_spacing=0.03,
                   legend_position='right')
    + labs(x='method', y='time (mins, signed-sqrt)', fill='')
)

# =======================================================================================
# Compose. Each panel is rendered at its final cell size (points therefore identical
# across panels), then imshow'd at equal aspect (no stretch/warp).
# =======================================================================================
# cell sizes (inches): top row A|B, then C (PPS), D (p_h1_xz), E (MSE), F (timing) full width.
WA, WB, HTOP, HC, HD, HE, HF = 9.6, 6.0, 6.0, 3.6, 3.64, 3.64, 3.64  # d at 70%; e,f match d
imgA = _render(pA, 'fig_A_percomponent.png', WA, HTOP)
imgB = _render(pB, 'fig_B_cov.png', WB, HTOP)
imgC = _render(pC, 'fig_C_pps.png', WA + WB, HC)
imgD = _render(pD, 'fig_D_ph1.png', WA + WB, HD)
imgE = _render(pE, 'fig_E_mse.png', WA + WB, HE)
imgF = _render(pF, 'fig_F_timing.png', WA + WB, HF)

# full-width rows below the A|B top row, in order.
fw_rows = [(imgC, HC, 'c'), (imgD, HD, 'd'), (imgE, HE, 'e'), (imgF, HF, 'f')]

_heights = [HTOP] + [h for _, h, _ in fw_rows]
fig = plt.figure(figsize=(WA + WB, sum(_heights) + 0.7))
gs = fig.add_gridspec(len(_heights), 2, height_ratios=_heights, width_ratios=[WA, WB],
                      hspace=0.06, wspace=0.03,
                      left=0.005, right=0.995, top=0.997, bottom=0.004)
_placements = [(fig.add_subplot(gs[0, 0]), imgA, 'a'),
               (fig.add_subplot(gs[0, 1]), imgB, 'b')]
for _r, (_img, _h, _lett) in enumerate(fw_rows, start=1):
    _placements.append((fig.add_subplot(gs[_r, :]), _img, _lett))
for ax, img, lett in _placements:
    ax.imshow(img)                       # default aspect='equal' -> no warp
    ax.axis('off')
    ax.text(0.003, 0.997, lett, transform=ax.transAxes, fontsize=18, fontweight='bold',
            va='top', ha='left')

PAPER_OUT = os.environ.get('PAPER_ADM_OUT', "/Users/or105/sandbox/bIRTistic/paper-adm")
os.makedirs(PAPER_OUT, exist_ok=True)
out = os.path.join(PAPER_OUT, 'mvn_figure_1')
fig.savefig(out + '.pdf', dpi=DPI, bbox_inches='tight', pad_inches=0.02)
fig.savefig(out + '.png', dpi=DPI, bbox_inches='tight', pad_inches=0.02)
plt.close(fig)
print(f"Saved figure 1 -> {out}.pdf / .png")
