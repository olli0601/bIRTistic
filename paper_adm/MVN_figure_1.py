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
    geom_boxplot, geom_hline, position_dodge,
    facet_grid, facet_wrap, coord_equal, scale_x_continuous, scale_x_discrete,
    scale_y_continuous, scale_y_discrete,
    scale_fill_gradient2, scale_fill_cmap, scale_fill_manual, scale_colour_manual, theme_bw, theme,
    element_text, element_blank, labs,
)
from utils import _futurama_palette                             # compare_methods method palette
from amortiser_io import load_interim_data

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
_D_SHORT = {
    'analytic': 'analytic (closed-form)',                        # per-z^s ground-truth reference
    'nested-MC using HMC for each (x,z)': 'nested-MC (HMC)',
    'Regression of endpt-x on w(z) using Gaussian approx': 'regression (Gauss)',
    'Amortiser-features-fixed-idcomp-qpsi-MLP-loss-multiquantilehead': 'fixed-idcomp',
    'Amortiser-features-MLP-xcomp-qpsi-MLP-loss-multiquantilehead': 'MLP-xcomp',
    'Amortiser-features-itemScompAtt-qpsi-MLP-loss-multiquantilehead': 'itemScompAtt',
    'Amortiser-features-itemXcompAtt-qpsi-MLP-loss-multiquantilehead': 'itemXcompAtt',
    'Amortiser-features-deepsetScompAtt-qpsi-MLP-loss-multiquantilehead': 'deepsetScompAtt',
    'Amortiser-features-deepsetXcompAtt-qpsi-MLP-loss-multiquantilehead': 'deepsetXcompAtt',
}
_AMORT_D = list(_D_SHORT)                                       # order of the all-amortised methods
# reproduce compare_methods.method_colours (futurama over the full _all_methods order) so panel d
# colours match the compare-methods plots, then key by the short label.
_ALL_METHODS = [
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
pD = (
    ggplot(box, aes(x='interim_month_year', fill='method', colour='method'))
    + geom_boxplot(aes(ymin='q025', lower='q25', middle='q50', upper='q75', ymax='q975', group='grp'),
                   stat='identity', position=position_dodge(width=0.8, preserve='single'),
                   width=0.72, size=0.25)
    + geom_hline(yintercept=_ETA, colour='black', size=0.8)
    + facet_wrap('~ component', ncol=3)
    + scale_fill_manual(values=_D_COLOURS, name='method')
    + scale_colour_manual(values=_D_OUTLINE, name='method')     # analytic drawn as transparent + black contour
    + scale_y_continuous(limits=[0, 1], expand=(0, 0), breaks=[0, 0.25, 0.5, 0.75, 1.0],
                         labels=['0%', '25%', '50%', '75%', '100%'])
    + scale_x_discrete(labels=_xlab_n)
    + BASE + theme(axis_text_x=element_text(size=7, angle=45, vjust=1, hjust=1),
                   legend_position='right', panel_spacing=0.03)
    + labs(x='interim analysis', y=r'$p(H_1^{\,j} \mid x,\, z^{s})$', fill='method', colour='method')
)

# =======================================================================================
# Panel E -- MSE heatmap (method x interim), full width; viridis. A subset of methods.
# =======================================================================================
DIR_CMP = os.path.join(SB, "py-mvn-interim-compare-methods-260609")
REMOVE_METHODS = {
    'Regression of endpt-x on w(z) using Gaussian approx',
    'Regression of endpt-x on w(z) - quantile',
    'Amortiser-features-itemScompAtt-qpsi-MLP-loss-multiquantilehead',
    'Amortiser-features-deepsetScompAtt-qpsi-MLP-loss-multiquantilehead',
    'Amortiser-features-deepsetXcompAtt-qpsi-MLP-loss-multiquantilehead',
}
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
dD['method'] = pd.Categorical(dD['method'].astype(str), categories=_meth_order, ordered=True)
dD['mse_c'] = dD['mse'].clip(lower=1e-4)                        # floor for log10
pE = (
    ggplot(dD, aes(x='interim_month_year', y='method', fill='mse_c'))
    + geom_tile()
    + facet_grid('. ~ J_label')
    + scale_x_discrete(expand=(0, 0), labels=_xlab_n) + scale_y_discrete(expand=(0, 0))
    + scale_fill_cmap(cmap_name='viridis', trans='log10')
    + BASE + theme(axis_text_x=element_text(size=7, angle=45, vjust=1, hjust=1),
                   panel_spacing=0.02)
    + labs(x='interim analysis', y='method',
           fill='MSE (log10)\n(mean over component-\nspecific squared errors)')
)

# =======================================================================================
# Panel F -- PIT-KS heatmap (same grammar as E). PIT-KS is only defined for the amortiser
# methods: the RAW architectures from by-architecture (per-interim), and the deploy-calibrated
# deepset (RGXP/RGXA, J=20) from the BvM driver's per-interim detail. Non-amortiser rows blank.
# DEACTIVATED by default (PIT-KS is not available for every method, so the panel is mostly
# blank). Re-enable with MVN_FIG_PIT_KS=1; the code below is kept as an option for later.
# =======================================================================================
INCLUDE_PIT_KS = os.environ.get('MVN_FIG_PIT_KS', '0') == '1'
pF = None
if INCLUDE_PIT_KS:
    SUF2LABEL = {
        'RGEA': 'Amortiser-features-fixed-idcomp-qpsi-MLP-loss-multiquantilehead',
        'RGEC': 'Amortiser-features-MLP-xcomp-qpsi-MLP-loss-multiquantilehead',
        'RGEF': 'Amortiser-features-itemXcompAtt-qpsi-MLP-loss-multiquantilehead',
        'RGED': 'Amortiser-features-itemScompAtt-qpsi-MLP-loss-multiquantilehead',
        'RGDS': 'Amortiser-features-deepsetScompAtt-qpsi-MLP-loss-multiquantilehead',
        'RGDX': 'Amortiser-features-deepsetXcompAtt-qpsi-MLP-loss-multiquantilehead',
        'RGXP': 'deepsetXcompAtt plain + affine (§14.4.6)',
        'RGXA': 'deepsetXcompAtt A power-law + BvM (§14.4.7)',
        'RGXC': 'deepsetXcompAtt C floor + BvM (§14.4.7)',
    }
    DIR_ARCH = os.path.join(SB, "py-mvn-interim-diagnostics-260916")
    DIR_BVM = os.path.join(SB, "py-mvn-interim-amortise-deepsetXcompAtt-ragged-bvm-260914")
    # (1) raw architectures, per-interim PIT (all J)
    _ap = pd.read_csv(os.path.join(DIR_ARCH, 'mvn_arch_pit_by_interim.csv'))
    _ap['method'] = _ap['suf'].map(SUF2LABEL)
    _ap = _ap[['method', 'J', 'interim_month_year', 'interim_date', 'pit_ks']]
    # (2) deploy-calibrated deepset (RGXP/RGXA), per-interim PIT from the BvM detail (J-labelled files)
    _bd = pd.concat([pd.read_csv(f) for f in
                     glob.glob(os.path.join(DIR_BVM, 'mvn_J*_calibration_detail.csv'))], ignore_index=True)
    _bp = _bd.groupby(['suf', 'J', 'interim_id'], as_index=False)['pit_ks'].mean()
    _bp['method'] = _bp['suf'].map(SUF2LABEL)
    _jref = load_interim_data(os.path.join(DIR_SIM, 'mvn_J20_interim_data.pkl'))    # shared interim schedule
    _bp['interim_month_year'] = _bp['interim_id'].map({k: _jref[k]['interim_month_year'] for k in _jref})
    _bp['interim_date'] = _bp['interim_id'].map({k: _jref[k]['interim_date'] for k in _jref})
    _bp = _bp[['method', 'J', 'interim_month_year', 'interim_date', 'pit_ks']]

    dE = pd.concat([_ap, _bp], ignore_index=True)
    dE = dE[dE['method'].isin(_meth_order) & dE['J'].isin([20, 100])].copy()
    dE['J_label'] = pd.Categorical('J=' + dE['J'].astype(str),
                                   categories=['J=20', 'J=100'], ordered=True)
    dE['interim_month_year'] = pd.Categorical(dE['interim_month_year'], categories=_io_d, ordered=True)
    dE['method'] = pd.Categorical(dE['method'].astype(str), categories=_meth_order, ordered=True)
    pF = (
        ggplot(dE, aes(x='interim_month_year', y='method', fill='pit_ks'))
        + geom_tile()
        + facet_grid('. ~ J_label')
        + scale_x_discrete(expand=(0, 0), labels=_xlab_n)
        + scale_y_discrete(expand=(0, 0), drop=False)                 # keep the (blank) non-amortiser rows
        + scale_fill_cmap(cmap_name='viridis')
        + BASE + theme(axis_text_x=element_text(size=7, angle=45, vjust=1, hjust=1),
                       panel_spacing=0.02)
        + labs(x='interim analysis', y='method',
               fill='PIT-KS\n(conditional calibration\nvs the reference)')
    )

# =======================================================================================
# Compose. Each panel is rendered at its final cell size (points therefore identical
# across panels), then imshow'd at equal aspect (no stretch/warp).
# =======================================================================================
# cell sizes (inches): top row A|B, then C (PPS), D (p_h1_xz), E (MSE), F (PIT-KS) full width.
WA, WB, HTOP, HC, HD, HE, HF = 9.6, 6.0, 6.0, 3.6, 5.2, 4.8, 4.8
imgA = _render(pA, 'fig_A_percomponent.png', WA, HTOP)
imgB = _render(pB, 'fig_B_cov.png', WB, HTOP)
imgC = _render(pC, 'fig_C_pps.png', WA + WB, HC)
imgD = _render(pD, 'fig_D_ph1.png', WA + WB, HD)
imgE = _render(pE, 'fig_E_mse.png', WA + WB, HE)

# full-width rows below the A|B top row, in order; F (PIT-KS) optional (see INCLUDE_PIT_KS).
fw_rows = [(imgC, HC, 'c'), (imgD, HD, 'd'), (imgE, HE, 'e')]
if INCLUDE_PIT_KS:
    fw_rows.append((_render(pF, 'fig_F_pit.png', WA + WB, HF), HF, 'f'))

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
