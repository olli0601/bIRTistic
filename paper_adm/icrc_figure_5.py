#!/usr/bin/env python3
"""Composite results figure 5 for the ICRC DASS-21 (DRC) amortised-PPS application.

    +-----------------------------+-------------------------------------------+
    | (a) DASS-21 severity bands  | (b) PCM fit, interim 3 (Baseline|Endline)  |
    +-----------------------------+-------------------------------------------+
    | (c) marginal p(rho|x):      | (d) contraction law: SD of marginal        |
    |     SVI vs ADM (BvM / no-BvM)|     posterior vs sqrt(n)                    |
    +-----------------------------+-------------------------------------------+
    | (e) amortised PPS trajectory   DASS-21 Anxiety | Stress | Depression     |
    +-----------------------------------------------------------------------+

All subfigures are rebuilt from the federated CSVs (distress_reduction/) + the SVI
interim-3 prob-by-question CSV, with the display edits requested for the paper.

Run:  pixi run python paper_adm/icrc_figure_5.py
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
from PIL import Image
from plotnine import (
    ggplot, aes, geom_col, geom_errorbar, geom_boxplot, geom_point, geom_line,
    facet_wrap, facet_grid, position_dodge, scale_x_continuous, scale_y_continuous,
    scale_fill_manual, scale_colour_manual, theme_bw, theme, element_text,
    element_blank, labs, labeller, guides, guide_legend,
)

warnings.filterwarnings('ignore')

SB = os.environ.get('ICRC_SB', "/Users/or105/sandbox/bIRTistic")
SVI = os.path.join(SB, "py-icrc-dass-drc_260902")
FED = os.path.join(SB, "py-icrc-dass-drc-amortise-deepsetScompAtt-itemamortise-J64-"
                       "perpart-hyper-logit-bvm-federated_261018", "distress_reduction")
# parallel deploy: head fine-tune + affine shift but NO BvM reshape (bvm=0) -> its RAGD cdf
NEW_FED = os.path.join(SB, "py-icrc-dass-drc-amortise-deepsetScompAtt-itemamortise-J64-"
                           "perpart-hyper-logit-affine-noBvM-federated_261019", "distress_reduction")
PFX = "pcm_1_interim"
SCRATCH = os.environ.get('ICRC_FIG_SCRATCH',
                         "/private/tmp/claude-501/-Users-or105-git-bIRTistic/"
                         "49372e22-d11b-4882-a442-d1c60bcbdfb0/scratchpad")
os.makedirs(SCRATCH, exist_ok=True)
DPI = 200

# --- colours ---------------------------------------------------------------------------
CSVI = '#2166ac'       # SVI reference (blue)
CBVM = '#b2182b'       # ADM: fine-tune + affine + BvM (dark red)
CAFF = '#984ea3'       # ADM: fine-tune + affine, no BvM (purple)
CRAW = '#ef8a62'       # ADM: raw net (orange)
SVI_LAB = 'SVI fit'
BVM_LAB = 'ADM with\nfine tune + BvM calibration'          # RAGD (main)     = head-ft + affine + BvM
AFF_LAB = 'ADM with\nfine tune, no BvM calibration'        # RAGD (NEW_FED)  = head-ft + affine, no BvM
RAW_LAB = 'ADM without\nfine tune + BvM calibration'       # RAGDbase        = raw net
SRC_ORDER = [SVI_LAB, BVM_LAB, AFF_LAB, RAW_LAB]
SRC_COL = {SVI_LAB: CSVI, BVM_LAB: CBVM, AFF_LAB: CAFF, RAW_LAB: CRAW}

# DASS-21 subscale display labels + facet orders ----------------------------------------
SCALE_LAB = {'DASS_Anxiety': 'DASS-21 Anxiety', 'DASS_Stress': 'DASS-21 Stress',
             'DASS_Depression': 'DASS-21 Depression'}
ROW_ORDER = ['DASS-21 Anxiety', 'DASS-21 Stress', 'DASS-21 Depression']   # (c)/(d) rows, (e) columns

XBREAKS = [1, 5, 10, 15, 20]        # interim x-ticks shown only at these (every 5th)


def _ilab_map():
    n = (pd.read_csv(os.path.join(FED, f"{PFX}_pps_RAGD_contraction.csv"))
         .groupby('interim_id')['n'].first())
    return {int(k): f"interim {int(k)}\nn={int(v)}" for k, v in n.items()}


ILAB = _ilab_map()


def _xlab(breaks=XBREAKS):
    return [ILAB.get(b, str(b)) for b in breaks]

BASE = (theme_bw()
        + theme(text=element_text(size=10), axis_title=element_text(size=10),
                axis_text=element_text(size=8), strip_text=element_text(size=9, face='bold'),
                strip_background=element_blank(), panel_grid_minor=element_blank(),
                legend_title=element_text(size=9), legend_text=element_text(size=8)))


def _render(p, name, w, h):
    path = os.path.join(SCRATCH, name)
    p.save(path, width=w, height=h, dpi=DPI, verbose=False, limitsize=False)
    return np.asarray(Image.open(path))


def _scale_cat(s, order=ROW_ORDER):
    return pd.Categorical(s.map(SCALE_LAB).fillna(s), categories=order, ordered=True)


# =======================================================================================
# Panel A -- DASS-21 severity bands, 3 subscale columns (static reference table).
# =======================================================================================
DASS_BANDS = [
    ('Normal',            '0-9',  '0-7',  '0-14'),
    ('Mild',              '10-13', '8-9',  '15-18'),
    ('Moderate',          '14-20', '10-14', '19-25'),
    ('Severe',            '21-27', '15-19', '26-33'),
    ('Extremely severe',  '28+',  '20+',  '34+'),
]
DASS_COLS = ['Severity', 'DASS-21\nDepression', 'DASS-21\nAnxiety', 'DASS-21\nStress']


def panel_a(ax, tl=0.05):
    ax.axis('off')
    pos = ax.get_position()
    x0 = max(0.0, (tl - pos.x0) / pos.width)   # table/title left-aligned to figure-x `tl` (= c/e plot-left)
    # title sits just below the panel letter 'a' (added in the compose step), table below the title.
    ax.text(x0, 0.80, 'DASS-21 severity thresholds (sum score per subscale)',
            transform=ax.transAxes, va='top', ha='left', fontsize=9)
    tab = ax.table(cellText=[list(r) for r in DASS_BANDS], colLabels=DASS_COLS,
                   cellLoc='center', colLoc='center', bbox=[x0, 0.14, 1.0 - x0, 0.58],
                   colWidths=[0.40, 0.20, 0.20, 0.20])   # wider first column so "Extremely severe" fits
    tab.auto_set_font_size(False); tab.set_fontsize(9)
    for (r, c), cell in tab.get_celld().items():
        cell.set_edgecolor('#d0d0d0')
        if r == 0:
            cell.set_text_props(fontweight='bold'); cell.set_facecolor('#f2eef0')
        if c == 0 and r > 0:
            cell.set_text_props(ha='left'); cell.PAD = 0.04


# =======================================================================================
# Panel B -- PCM fit, interim 3: category distribution per subscale, Baseline vs Endline.
# =======================================================================================
TIME_LAB = {'Baseline': 'Baseline\nprior to group psychological\nsupport intervention',
            'Endline': 'Endline\nafter group psychological\nsupport intervention'}


def panel_b():
    df = pd.read_csv(os.path.join(SVI, f"{PFX}_3_prob_by_question_fit.csv"))
    df['scale'] = _scale_cat(df['item_label'])
    df['time'] = pd.Categorical(df['group_label'].map(TIME_LAB),   # migration: time_label -> group_label
                                categories=[TIME_LAB['Baseline'], TIME_LAB['Endline']], ordered=True)
    # strip the leading severity index ("0 Normal" -> "Normal", "4 Extremely severe" -> "Extremely severe")
    cats = sorted(df['y_label'].unique())
    strip = {c: c.split(' ', 1)[1] if c[:1].isdigit() else c for c in cats}
    df['cat'] = pd.Categorical(df['y_label'].map(strip),
                               categories=[strip[c] for c in cats], ordered=True)
    dodge = position_dodge(width=0.9, preserve='single')
    p = (ggplot(df, aes('cat', group='scale'))
         + geom_col(aes(y='p_emp', fill='scale'), position=dodge, width=0.8, alpha=0.8)   # empirical bars
         + geom_boxplot(aes(ymin='q_lower', lower='iqr_lower', middle='median',           # SVI fit (outline)
                            upper='iqr_upper', ymax='q_upper'),
                        position=dodge, stat='identity', alpha=0, width=0.3, size=0.3)
         + facet_wrap('~ time', ncol=2)
         + scale_fill_manual(values={'DASS-21 Anxiety': '#e08214', 'DASS-21 Stress': '#4d4d4d',
                                     'DASS-21 Depression': '#b3b3b3'},
                             name='discretised DASS-21 scale')
         + scale_y_continuous(labels=lambda l: [f'{v:.0%}' for v in l])
         + BASE + theme(axis_text_x=element_text(size=7, angle=30, vjust=1, hjust=1),
                        legend_position='right')
         + labs(x='', y='responses at 3rd interim'))
    return p


# =======================================================================================
# Panels C & D -- contraction (marginal p(rho|x) box) and contraction law (SD vs sqrt n).
# Colour = SVI / ADM-with-BvM / ADM-no-BvM. Rows = the 3 subscales.
# =======================================================================================
def _load_cdf():
    ragd = pd.read_csv(os.path.join(FED, f"{PFX}_pps_RAGD_contraction_cdf_by_item.csv"))
    base = pd.read_csv(os.path.join(FED, f"{PFX}_pps_RAGDbase_contraction_cdf_by_item.csv"))
    rows = []
    r = ragd.copy()
    r['grp'] = np.where(r['source'] == 'SVI', SVI_LAB, BVM_LAB)
    rows.append(r)
    b = base[base['source'] == 'amortiser'].copy()
    b['grp'] = RAW_LAB
    rows.append(b)
    aff_csv = os.path.join(NEW_FED, f"{PFX}_pps_RAGD_contraction_cdf_by_item.csv")
    if os.path.exists(aff_csv):                         # head-ft + affine, no BvM
        a = pd.read_csv(aff_csv)
        a = a[a['source'] == 'amortiser'].copy(); a['grp'] = AFF_LAB
        rows.append(a)
    df = pd.concat(rows, ignore_index=True)
    df['scale'] = _scale_cat(df['item_label'])
    df['grp'] = pd.Categorical(df['grp'], categories=SRC_ORDER, ordered=True)
    return df


def panel_c():
    df = _load_cdf()
    p = (ggplot(df, aes('interim_id', fill='grp'))
         + geom_boxplot(aes(ymin='ymin', lower='lower', middle='middle', upper='upper', ymax='ymax',
                            group='grp'),
                        stat='identity', position=position_dodge(width=0.8), size=0.25, width=0.7, alpha=0.85)
         + facet_wrap('~ scale', ncol=1, scales='free_y')
         + scale_fill_manual(values=SRC_COL, name='')
         + scale_x_continuous(breaks=XBREAKS, labels=_xlab())
         + BASE + theme(legend_position='bottom', legend_title=element_blank())
         + labs(x='',
                y='posterior effect size at interim\n$p(\\rho \\mid \\mathrm{current\\ data\\ at\\ interim})$'))
    return p


def _law_df():
    # posterior SD of the marginal p(rho|x) per interim = IQR(box)/1.349 (matches the deploy's contlaw sd);
    # interim -> n from the contraction table. Three sources: SVI, ADM+BvM (RAGD), ADM no-BvM (RAGDbase).
    nmap = (pd.read_csv(os.path.join(FED, f"{PFX}_pps_RAGD_contraction.csv"))
            .groupby('interim_id')['n'].first())
    ragd = pd.read_csv(os.path.join(FED, f"{PFX}_pps_RAGD_contraction_cdf_by_item.csv"))
    base = pd.read_csv(os.path.join(FED, f"{PFX}_pps_RAGDbase_contraction_cdf_by_item.csv"))
    srcs = [(ragd, 'main'), (base, 'raw')]
    aff_csv = os.path.join(NEW_FED, f"{PFX}_pps_RAGD_contraction_cdf_by_item.csv")
    if os.path.exists(aff_csv):
        srcs.append((pd.read_csv(aff_csv), 'aff'))
    recs = []
    for df_src, kind in srcs:
        for _, r in df_src.iterrows():
            if kind in ('raw', 'aff') and r['source'] != 'amortiser':
                continue
            grp = ({'raw': RAW_LAB, 'aff': AFF_LAB}.get(kind) if kind != 'main'
                   else (SVI_LAB if r['source'] == 'SVI' else BVM_LAB))
            n = nmap.get(r['interim_id'], np.nan)
            recs.append(dict(item_label=r['item_label'], n=n, sqrt_n=np.sqrt(n),
                             sd=(r['upper'] - r['lower']) / 1.349, grp=grp))
    df = pd.DataFrame(recs)
    df['scale'] = _scale_cat(df['item_label'])
    df['grp'] = pd.Categorical(df['grp'], categories=SRC_ORDER, ordered=True)
    return df


def panel_d():
    df = _law_df()
    p = (ggplot(df, aes('sqrt_n', 'sd', colour='grp'))
         + geom_point(size=1.2)
         + geom_line(size=0.5)
         + facet_wrap('~ scale', ncol=1, scales='free_y')
         + scale_colour_manual(values=SRC_COL, name='')
         + BASE + theme(legend_position='bottom', legend_title=element_blank())
         + labs(x='sqrt(#participants at interim)',
                y='std dev of marginal posterior\neffect size at interim'))
    return p


# =======================================================================================
# Panel E -- amortised PPS trajectory. Columns = the 3 subscales; colour = eta0 threshold.
# =======================================================================================
def panel_e():
    pps = pd.read_csv(os.path.join(FED, f"{PFX}_pps_RAGD_eta0_pps.csv"))   # interim,date,item,eta0_pct,pps
    pps['scale'] = _scale_cat(pps['item'])
    thr = sorted(pps['eta0_pct'].unique())
    import matplotlib.cm as cm, matplotlib.colors as mcolors
    vir = plt.get_cmap('viridis'); nn = max(1, len(thr) - 1)
    eta_col = {str(int(t)): mcolors.to_hex(vir(1 - i / nn)) for i, t in enumerate(thr)}
    pps['eta0'] = pd.Categorical(pps['eta0_pct'].astype(int).astype(str),
                                 categories=[str(int(t)) for t in thr], ordered=True)
    p = (ggplot(pps, aes('interim', 'pps', colour='eta0', group='eta0'))
         + geom_line(size=0.5) + geom_point(size=1.2)
         + facet_wrap('~ scale', ncol=3)
         + scale_colour_manual(values=eta_col, name='success threshold\neta0 (rel change x100)',
                               labels=lambda l: [f'{v}%' for v in l])
         + scale_x_continuous(breaks=XBREAKS, labels=_xlab())
         + scale_y_continuous(limits=[0, 1], labels=lambda l: [f'{v:.0%}' for v in l])
         + BASE + theme(legend_position='right')
         + labs(x='', y='PPS'))
    return p


# =======================================================================================
# Compose: row1 = a | b, row2 = c | d, row3 = e (full width).
# =======================================================================================
def _plot_bounds(img):
    """Figure-x fractions of the tall panel borders (plot area), ignoring y-axis labels / legends."""
    g = np.asarray(Image.fromarray(img).convert('L')); H, W = g.shape
    band = (g[int(0.20 * H):int(0.80 * H)] < 130).mean(axis=0)
    cols = np.where(band > 0.6)[0]
    return (cols[0] / W, cols[-1] / W) if len(cols) > 1 else (0.05, 0.95)


def _place(ax, img, t0, t1):
    """Reposition `ax` so the image's plot area (not its margins) maps to figure-x [t0, t1]."""
    ax.imshow(img, aspect='auto'); ax.axis('off')
    L, R = _plot_bounds(img)
    pos = ax.get_position()
    w = (t1 - t0) / (R - L)
    ax.set_position([t0 - w * L, pos.y0, w, pos.height])


def main():
    imgC = _render(panel_c(), 'icrc_c.png', 7.4, 7.0)
    imgD = _render(panel_d(), 'icrc_d.png', 7.4, 7.0)
    imgE = _render(panel_e(), 'icrc_e.png', 15.0, 3.6)

    TL = 0.05                                   # common left frame: (a) table, (c) & (e) plot axes
    FIGW, FIGH = 15.2, 14.4
    fig = plt.figure(figsize=(FIGW, FIGH))
    gs = fig.add_gridspec(3, 12, height_ratios=[3.3, 7.0, 3.8],
                          hspace=0.08, wspace=0.03, left=0.004, right=0.996, top=0.985, bottom=0.004)
    axA = fig.add_subplot(gs[0, 0:5]); panel_a(axA, TL)    # wider so the table fits
    axB = fig.add_subplot(gs[0, 6:12]); axB.axis('off')    # fill to right edge, no whitespace
    # Render panel (b) at its destination cell's EXACT inch dimensions so plotnine's fonts are
    # not horizontally squashed by imshow's stretch-to-fill (aspect='auto' matches render==cell aspect).
    _pB = axB.get_position()
    imgB = _render(panel_b(), 'icrc_b.png', _pB.width * FIGW, _pB.height * FIGH)
    axB.imshow(imgB, aspect='auto')
    axC = fig.add_subplot(gs[1, 0:6]); _place(axC, imgC, TL, 0.485)
    axD = fig.add_subplot(gs[1, 6:12]); _place(axD, imgD, 0.55, 0.985)
    axE = fig.add_subplot(gs[2, 0:12]); _place(axE, imgE, TL, 0.86)
    _LETX = {'a': 0.006, 'c': 0.006, 'e': 0.006}      # left-column letters share one x
    for ax, lett in ((axA, 'a'), (axB, 'b'), (axC, 'c'), (axD, 'd'), (axE, 'e')):
        pos = ax.get_position()
        fig.text(_LETX.get(lett, pos.x0 + 0.002), pos.y1 - 0.004, lett,
                 fontsize=17, fontweight='bold', va='top', ha='left')

    out = os.path.join(os.environ.get('PAPER_ADM_OUT', os.path.dirname(os.path.abspath(__file__))),
                       'icrc_figure_5')
    fig.savefig(out + '.pdf', dpi=DPI, pad_inches=0.02)
    fig.savefig(out + '.png', dpi=DPI, pad_inches=0.02)
    plt.close(fig)
    print(f"Saved figure 5 -> {out}.pdf / .png")


if __name__ == "__main__":
    main()
