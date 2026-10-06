#!/usr/bin/env python3
"""Composite results figure 3 for the HVTN 505 BAMA amortised-PPS application (§3.19).

    +---------------------------------------------------------------+
    | (a) participant-accrual timeline (vaccine red | placebo blue) |
    +---------------------------------------------------------------+
    | (b) PCM fit interim 1        (b) PCM fit interim 10            |   <- dashed arrows to (a)
    +---------------------------------------------------------------+
    | (c) rho vs eta0 thresholds   Con6 gp120 | A244 V1V2 | Gag p24  |
    +---------------------------------------------------------------+
    | (d) amortised PPS trajectory Con6 gp120 | A244 V1V2 | Gag p24  |
    +---------------------------------------------------------------+

Parts c/d are reproduced inline from the federated eta0 CSVs (items on COLUMNS,
unlike FederatedDiagnostics which puts them on rows); part b reuses the existing
prob_by_question_fit PNGs; part a is built from the per-interim subject tables.

Run:  pixi run python paper_adm/hvtn505_figure_3.py
"""
import os
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'python'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from PIL import Image
from plotnine import (
    ggplot, aes, geom_area, geom_boxplot, geom_hline, geom_line, geom_point,
    facet_wrap, scale_x_continuous, scale_y_continuous, scale_fill_manual,
    scale_fill_gradient2, scale_colour_manual, theme_bw, theme, element_text,
    element_blank, labs,
)
from mizani.bounds import squish
from amortiser_diag_plots import FederatedDiagnostics

warnings.filterwarnings('ignore')

SB = os.environ.get('CAVD_SB', "/Users/or105/sandbox/bIRTistic")
SVI = os.path.join(SB, "py-cavd-vtn505-bama_260918")
FED = os.path.join(SB, "py-cavd-vtn505-bama-amortise-deepsetXcompAtt-itemamortise-J64-"
                       "ftheadexpand-bvm-federated_260923", "vaccine_vs_placebo")
PFX = "pcm_1_interim"
SCRATCH = os.environ.get('HVTN_FIG_SCRATCH',
                         "/private/tmp/claude-501/-Users-or105-git-bIRTistic/"
                         "49372e22-d11b-4882-a442-d1c60bcbdfb0/scratchpad")
os.makedirs(SCRATCH, exist_ok=True)
DPI = 200
N_INTERIMS = 10
TARGET3 = ['Con6 gp120 (B)', 'A244 V1V2 (AE)', 'Gag p24']         # items shown in c/d (columns)
VAC_COL, PLA_COL = '#c1272d', '#1f5fa8'                            # vaccine red, placebo blue

BASE = (theme_bw()
        + theme(text=element_text(size=11), axis_title=element_text(size=11),
                axis_text=element_text(size=8), strip_text=element_text(size=9, face='bold'),
                strip_background=element_blank(),
                panel_grid_minor=element_blank()))


def _render(p, name, w, h):
    path = os.path.join(SCRATCH, name)
    p.save(path, width=w, height=h, dpi=DPI, verbose=False, limitsize=False)
    return np.asarray(Image.open(path))


# =======================================================================================
# Interim n (total, from the per-interim subject tables) + per-arm cumulative counts.
# =======================================================================================
_rows = []
for k in range(1, N_INTERIMS + 1):
    dk = pd.read_csv(os.path.join(SVI, f"{PFX}_{k}_data_dp1.csv")).drop_duplicates('pid')
    n_pla = int((dk['group'] == 0).sum())
    n_vac = int((dk['group'] == 1).sum())
    _rows.append({'interim': k, 'placebo': n_pla, 'vaccine': n_vac, 'n': n_pla + n_vac})
acc = pd.DataFrame(_rows)
NTOT = dict(zip(acc['interim'], acc['n']))
_ilab = {k: f"interim {k}\n(n={NTOT[k]})" for k in range(1, N_INTERIMS + 1)}
_iorder = [_ilab[k] for k in range(1, N_INTERIMS + 1)]


def _xcat(df):
    df = df.copy()
    df['xlab'] = pd.Categorical(df['interim'].astype(int).map(_ilab),
                                categories=_iorder, ordered=True)
    return df


# =======================================================================================
# Panel A -- participant-accrual timeline. Stacked area: placebo (blue) + vaccine (red)
# to the total-participant line; x = interim (labelled by total n).
# =======================================================================================
accL = acc.melt(id_vars=['interim', 'n'], value_vars=['placebo', 'vaccine'],
                var_name='arm', value_name='count')
accL['arm'] = pd.Categorical(accL['arm'], categories=['placebo', 'vaccine'], ordered=True)
pA = (
    ggplot(accL, aes('interim', 'count', fill='arm'))
    + geom_area(position='stack', colour='white', size=0.2)
    + scale_fill_manual(values={'vaccine': VAC_COL, 'placebo': PLA_COL}, name='arm')
    + scale_x_continuous(breaks=list(range(1, N_INTERIMS + 1)),
                         labels=[str(NTOT[k]) for k in range(1, N_INTERIMS + 1)],
                         expand=(0.01, 0))
    + scale_y_continuous(expand=(0, 0, 0.02, 0))
    + BASE + theme(legend_position='right', panel_grid_major_x=element_blank())
    + labs(x='cumulative participants at interim analysis', y='participants')
)


# =======================================================================================
# Panels C & D -- read the federated eta0 CSVs, filter to the 3 target items (COLUMNS).
# =======================================================================================
rho = pd.read_csv(os.path.join(FED, f"{PFX}_pps_RAGD_eta0_rho.csv"))       # interim,date,item,rho_pct
pps = pd.read_csv(os.path.join(FED, f"{PFX}_pps_RAGD_eta0_pps.csv"))       # interim,date,item,eta0_pct,pps
rho = rho[rho['item'].isin(TARGET3)].rename(columns={'item': 'item_label'})
pps = pps[pps['item'].isin(TARGET3)].rename(columns={'item': 'item_label'})
for d in (rho, pps):
    d['item_label'] = pd.Categorical(d['item_label'], categories=TARGET3, ordered=True)
rho = _xcat(rho)
pps = _xcat(pps)

# eta0 threshold colours (viridis, darker = harder), shared by c (lines) and d (trajectory)
_thr = sorted(pps['eta0_pct'].unique())
_vir = cm.get_cmap('viridis'); _nn = max(1, len(_thr) - 1)
_ETA_COL = {str(int(t)): mcolors.to_hex(_vir(1 - i / _nn)) for i, t in enumerate(_thr)}
pps['eta0'] = pd.Categorical(pps['eta0_pct'].astype(int).astype(str),
                             categories=[str(int(t)) for t in _thr], ordered=True)
# threshold levels drawn as horizontal rho lines in panel c
ldf = pps[['item_label', 'eta0_pct']].drop_duplicates().copy()
ldf['rho_pct'] = ldf['eta0_pct'].astype(float)
ldf['eta0'] = pd.Categorical(ldf['eta0_pct'].astype(int).astype(str),
                             categories=[str(int(t)) for t in _thr], ordered=True)
ldf = _xcat(ldf) if 'interim' in ldf.columns else ldf

# box = 25-75%, whiskers = 2.5-97.5% of the per-draw rho (no Tukey outlier dots).
qc = (rho.groupby(['item_label', 'xlab'], observed=True)['rho_pct']
      .agg(q025=lambda s: s.quantile(0.025), q25=lambda s: s.quantile(0.25),
           q50='median', q75=lambda s: s.quantile(0.75), q975=lambda s: s.quantile(0.975))
      .reset_index())
pC = (
    ggplot(qc, aes('xlab'))
    + geom_boxplot(aes(ymin='q025', lower='q25', middle='q50', upper='q75', ymax='q975'),
                   stat='identity', fill='#d9d9d9', size=0.3)
    + geom_hline(ldf, aes(yintercept='rho_pct', colour='eta0'), size=0.55)
    + scale_colour_manual(values=_ETA_COL, name='success threshold\neta0 (rel. change x100)',
                          labels=lambda l: [f'{v}%' for v in l])
    + facet_wrap('~ item_label', ncol=3, scales='free_y')
    + BASE + theme(axis_text_x=element_text(size=7, angle=45, vjust=1, hjust=1),
                   legend_position='right', legend_title=element_text(ma='left', ha='left'), panel_spacing=0.03)
    + labs(x='interim analysis',
           y='posterior effect size at interim\n$p(\\rho \\mid \\mathrm{current\\ data\\ at\\ interim})$')
)

pD = (
    ggplot(pps, aes('xlab', 'pps', colour='eta0', group='eta0'))
    + geom_line(size=0.5) + geom_point(size=1.3)
    + scale_colour_manual(values=_ETA_COL, name='success threshold\neta0 (rel. change x100)',
                          labels=lambda l: [f'{v}%' for v in l])
    + scale_y_continuous(limits=[0, 1], labels=lambda l: [f'{v:.0%}' for v in l])
    + facet_wrap('~ item_label', ncol=3)
    + BASE + theme(axis_text_x=element_text(size=7, angle=45, vjust=1, hjust=1),
                   legend_position='right', legend_title=element_text(ma='left', ha='left'), panel_spacing=0.03)
    + labs(x='interim analysis', y='PPS')
)


# =======================================================================================
# Panel E -- conditional-calibration PIT box per interim, 3 items on COLUMNS. Fill = PIT
# median - 0.5 (blue = amortiser median too low, red = too high). Data from the diagnostics
# in-memory store (no CSV) via FederatedDiagnostics.load_plotdata().
# =======================================================================================
_AMORT = dict(
    svi='py-cavd-vtn505-bama_260918',
    fed='py-cavd-vtn505-bama-amortise-deepsetXcompAtt-itemamortise-J64-ftheadexpand-bvm-federated_260923',
    title='CAVD HVTN 505 BAMA', item_kind='antigen', mixed_units=False,
    calib_prefix='cavd_bama_amortiser', eta0_units='relative change x100',
    ctag_prefix='cavd-vtn505-fed',
    endpoints=[dict(rho='vaccine_vs_placebo', instance='S3 rel-change', build='svi',
                    net='scale-feat', widetok=0, case_c=2, warp='none',
                    eta0=0.5, grid='0.25,0.5,0.75,1.0')])
pit = FederatedDiagnostics(SB, _AMORT, verbose=False).load_plotdata()['pitbox']['RAGD'].copy()
pit = pit[pit['item_label'].isin(TARGET3)].copy()
pit['item_label'] = pd.Categorical(pit['item_label'], categories=TARGET3, ordered=True)
pit['xlab'] = pd.Categorical(pit['interim_id'].astype(int).map(_ilab), categories=_iorder, ordered=True)
pit['dev'] = pit['middle'] - 0.5
pE = (
    ggplot(pit, aes('xlab', ymin='ymin', lower='lower', middle='middle',
                    upper='upper', ymax='ymax', fill='dev'))
    + geom_hline(yintercept=[0.10, 0.25, 0.50, 0.75, 0.90], linetype='dashed',
                 colour='#9e9e9e', size=0.3)
    + geom_boxplot(stat='identity', alpha=0.9, size=0.3, width=0.7)
    + scale_fill_gradient2(low='#2166ac', mid='#f7f7f7', high='#b2182b', midpoint=0.0,
                           limits=[-0.3, 0.3], oob=squish,
                           name='Median of\nPIT diagnostic\n(target = 0.0)')
    + facet_wrap('~ item_label', ncol=3)
    + scale_y_continuous(limits=[0, 1])
    + BASE + theme(axis_text_x=element_text(size=7, angle=45, vjust=1, hjust=1),
                   legend_position='right', legend_title=element_text(ma='left', ha='left'), panel_spacing=0.03)
    + labs(x='interim analysis', y='PIT  $u = F_{amortiser}(\\rho_{SVI} \\mid x,\\, z)$')
)


# =======================================================================================
# Compose. a (timeline) | b (two PCM fits, dashed arrows up to a) | c | d | e.
# =======================================================================================
WFULL, HA, HB, HC, HD, HE = 15.6, 2.0, 5.2, 4.4, 4.0, 4.0
imgA = _render(pA, 'hv_A_accrual.png', WFULL, HA)
imgC = _render(pC, 'hv_C_rho.png', WFULL, HC)
imgD = _render(pD, 'hv_D_traj.png', WFULL, HD)
imgE = _render(pE, 'hv_E_pit.png', WFULL, HE)
# part b: split the prob-fit PNGs into plot (above) + the shared antigen legend (below);
# keep ONE legend, placed TOP-LEFT above the two plots (native size), plots aligned.
_b1 = Image.open(os.path.join(SVI, f"{PFX}_1_prob_by_question_fit_part1.png"))
_b10 = Image.open(os.path.join(SVI, f"{PFX}_10_prob_by_question_fit_part1.png"))
_Wp, _Hp = _b1.size
_SPLIT = 1300                                                    # row: below the x-title, above the legend
imgB1 = np.asarray(_b1.crop((0, 0, _Wp, _SPLIT)))               # plot only
imgB10 = np.asarray(_b10.crop((0, 0, _Wp, _SPLIT)))
imgBleg = np.asarray(_b1.crop((0, _SPLIT, _Wp, _Hp)))          # single antigen legend (wide, short)

fig = plt.figure(figsize=(WFULL, HA + HB + HC + HD + HE + 0.9))
gs = fig.add_gridspec(5, 1, height_ratios=[HA, HB, HC, HD, HE],
                      hspace=0.10, left=0.005, right=0.995, top=0.995, bottom=0.004)

# --- vertical alignment: map every panel's plot-area (detected via the tall panel-border
# columns) to the same figure-x frame [TL, TR]; the y-titles fall left of TL, legends right of TR.
TL, TR = 0.055, 0.85


def _plot_bounds(img):
    g = np.asarray(Image.fromarray(img).convert('L')); H, W = g.shape
    band = (g[int(0.20 * H):int(0.80 * H)] < 130).mean(axis=0)
    cols = np.where(band > 0.6)[0]
    return (cols[0] / W, cols[-1] / W) if len(cols) > 1 else (0.05, 0.95)


def _place(ax, img, t0, t1):
    ax.imshow(img, aspect='auto'); ax.axis('off')
    L, R = _plot_bounds(img)
    pos = ax.get_position()
    w = (t1 - t0) / (R - L)
    ax.set_position([t0 - w * L, pos.y0, w, pos.height])


axA = fig.add_subplot(gs[0]); _place(axA, imgA, TL, TR)
# part-b block: single legend centred above; the two plots split [TL, TR] below.
gsB = gs[1].subgridspec(2, 2, height_ratios=[0.22, 0.78], width_ratios=[1, 1],
                        hspace=0.02, wspace=0.03)
axBleg = fig.add_subplot(gsB[0, :]); axBleg.imshow(imgBleg); axBleg.axis('off')
_mid = TL + 0.5 * (TR - TL)
axB1 = fig.add_subplot(gsB[1, 0]); _place(axB1, imgB1, TL, _mid - 0.008)
axB10 = fig.add_subplot(gsB[1, 1]); _place(axB10, imgB10, _mid + 0.008, TR)
axC = fig.add_subplot(gs[2]); _place(axC, imgC, TL, TR)
axD = fig.add_subplot(gs[3]); _place(axD, imgD, TL, TR)
axE = fig.add_subplot(gs[4]); _place(axE, imgE, TL, TR)

# dashed arrows: top-centre of each part-b panel -> the interim's x-position on the (a) strip.
# after _place, panel a's plot-area maps to [TL, TR]; its x is continuous 1..N with expand=(0.01,0),
# so interim k sits at fraction (k-1+e)/(N-1+2e) of [TL, TR], e = 0.01*(N-1).
pA_bb = axA.get_position()
_e = 0.01 * (N_INTERIMS - 1)


def _acc_x_fig(interim):
    frac = (interim - 1 + _e) / (N_INTERIMS - 1 + 2 * _e)
    return TL + (TR - TL) * frac


for ax, interim in ((axB1, 1), (axB10, 10)):
    pb = ax.get_position()
    x_from, y_from = 0.5 * (pb.x0 + pb.x1), pb.y1
    x_to, y_to = _acc_x_fig(interim), pA_bb.y0 + 0.10 * (pA_bb.y1 - pA_bb.y0)
    arr = matplotlib.patches.FancyArrowPatch(
        (x_from, y_from), (x_to, y_to), arrowstyle='-|>', linestyle='--',
        color='black', lw=1.2, mutation_scale=16, shrinkA=3, shrinkB=3)
    arr.set_transform(fig.transFigure)
    fig.add_artist(arr)

# panel letters, all at a common figure-x so 'b' lines up with a/c/d/e
_LETX = 0.010
for ax, lett in ((axA, 'a'), (axBleg, 'b'), (axC, 'c'), (axD, 'd'), (axE, 'e')):
    fig.text(_LETX, ax.get_position().y1 - 0.002, lett, fontsize=18, fontweight='bold',
             va='top', ha='left')

PAPER_OUT = os.environ.get('PAPER_ADM_OUT', os.path.dirname(os.path.abspath(__file__)))
os.makedirs(PAPER_OUT, exist_ok=True)
out = os.path.join(PAPER_OUT, 'hvtn505_figure_3')
fig.savefig(out + '.pdf', dpi=DPI, pad_inches=0.02)
fig.savefig(out + '.png', dpi=DPI, pad_inches=0.02)
plt.close(fig)
print(f"Saved figure 3 -> {out}.pdf / .png")
