#!/usr/bin/env python3
"""Supplementary figure for the ICRC DASS-21 (DRC) amortised-PPS application.

Conditional-calibration PIT box (like pcm_1_interim_pps_RAGD_pit_box_by_item.pdf), but for the
FOUR algorithm options of figure-5 c/d, laid out 2x2. Each panel = one algorithm, with the 3 DASS
subscales on ROW facets; x = interim, fill = PIT median - 0.5 (blue = amortiser too low, red = too
high; white = calibrated). Dashed references at 0.10/0.25/0.5/0.75/0.90.

    +-----------------------------+-----------------------------+
    | SVI fit (reference)         | ADM fine-tune + BvM         |
    +-----------------------------+-----------------------------+
    | ADM fine-tune, no BvM       | ADM without fine-tune + BvM |
    +-----------------------------+-----------------------------+

Run:  pixi run python paper_adm/icrc_supp_figure.py
"""
import os
import sys
import pickle
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'python'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
from PIL import Image
from plotnine import (ggplot, aes, geom_boxplot, geom_hline, facet_wrap, scale_x_continuous,
                      scale_y_continuous, scale_fill_gradient2, theme_bw, theme, element_text,
                      element_blank, labs, ggtitle)
from mizani.bounds import squish

warnings.filterwarnings('ignore')

SB = os.environ.get('ICRC_SB', "/Users/or105/sandbox/bIRTistic")
MAIN = os.path.join(SB, "py-icrc-dass-drc-amortise-deepsetXcompAtt-itemamortise-J64-"
                        "ftheadexpand-bvm-federated_260924", "distress_reduction")
AFFN = os.path.join(SB, "py-icrc-dass-drc-amortise-deepsetXcompAtt-itemamortise-J64-"
                        "scale-feat-ftheadexpand-affine-noBvM-federated_261012", "distress_reduction")
PFX = "pcm_1_interim"
SCRATCH = os.environ.get('ICRC_FIG_SCRATCH', "/private/tmp/claude-501/-Users-or105-git-bIRTistic/"
                         "49372e22-d11b-4882-a442-d1c60bcbdfb0/scratchpad")
os.makedirs(SCRATCH, exist_ok=True)
DPI = 200

SCALE_LAB = {'DASS_Anxiety': 'DASS-21 Anxiety', 'DASS_Stress': 'DASS-21 Stress',
             'DASS_Depression': 'DASS-21 Depression'}
ROW_ORDER = ['DASS-21 Anxiety', 'DASS-21 Stress', 'DASS-21 Depression']
XBREAKS = [1, 5, 10, 15, 20]
LOW, MID, HIGH = '#2166ac', '#f7f7f7', '#b2182b'

BASE = (theme_bw()
        + theme(text=element_text(size=10), axis_title=element_text(size=10),
                axis_text=element_text(size=8), strip_text=element_text(size=9, face='bold'),
                strip_background=element_blank(), panel_grid_minor=element_blank(),
                plot_title=element_text(size=11, face='bold', ha='center')))


def _ilab_map():
    n = (pd.read_csv(os.path.join(MAIN, f"{PFX}_pps_RAGD_contraction.csv"))
         .groupby('interim_id')['n'].first())
    return {int(k): f"interim {int(k)}\nn={int(v)}" for k, v in n.items()}


ILAB = _ilab_map()
NI = len(ILAB)


def _pitbox(pu):
    """Rebuild per (item, interim) box from a pit/marg u frame (item_label, u); rows interim-major.
    Also carry the per-cell KS of u vs Uniform(0,1) (= pit-KS / marg-KS for that interim x item)."""
    recs = []
    for it, g in pu.groupby('item_label'):
        u = g['u'].values; nd = len(u) // NI
        U = u[:NI * nd].reshape(NI, nd)
        for k in range(NI):
            a = np.sort(U[k]); ec = np.arange(1, len(a) + 1) / len(a)
            recs.append(dict(interim_id=k + 1, item_label=it,
                             ymin=np.quantile(a, 0.025), lower=np.quantile(a, 0.25),
                             middle=np.quantile(a, 0.5), upper=np.quantile(a, 0.75),
                             ymax=np.quantile(a, 0.975), ks=float(np.max(np.abs(ec - a)))))
    return pd.DataFrame(recs)


def _svi_reference(items):
    """SVI self-PIT is Uniform(0,1) by construction -> the calibrated target (median 0.5, KS 0)."""
    q = dict(ymin=0.025, lower=0.25, middle=0.5, upper=0.75, ymax=0.975, ks=0.0)
    return pd.DataFrame([dict(interim_id=k + 1, item_label=it, **q)
                         for it in items for k in range(NI)])


def _prep(df):
    df = df.copy()
    df['scale'] = pd.Categorical(df['item_label'].map(SCALE_LAB), categories=ROW_ORDER, ordered=True)
    df['dev'] = df['middle'] - 0.5
    return df


def _panel(df, title, ylab, kind):
    from plotnine import scale_fill_gradient
    df = _prep(df)
    # PIT (conditional): colour = signed median bias (affine-shift target). marg: colour = marg-KS
    # magnitude (spread mismatch; what BvM fixes), since the median is ~0.5 regardless.
    fill_aes = 'dev' if kind == 'pit' else 'ks'
    scale = (scale_fill_gradient2(low=LOW, mid=MID, high=HIGH, midpoint=0.0, limits=[-0.3, 0.3], oob=squish)
             if kind == 'pit' else
             scale_fill_gradient(low=MID, high=HIGH, limits=[0, 0.3], oob=squish))
    return (ggplot(df, aes('interim_id', ymin='ymin', lower='lower', middle='middle',
                           upper='upper', ymax='ymax', fill=fill_aes))
            + geom_hline(yintercept=[0.10, 0.25, 0.50, 0.75, 0.90], linetype='dashed',
                         colour='#9e9e9e', size=0.3)
            + geom_boxplot(aes(group='interim_id'), stat='identity', alpha=0.9, size=0.3, width=0.8)
            + scale
            + facet_wrap('~ scale', ncol=1)
            + scale_x_continuous(breaks=XBREAKS, labels=[ILAB.get(b, str(b)) for b in XBREAKS])
            + scale_y_continuous(limits=[0, 1])
            + BASE + theme(legend_position='none')
            + ggtitle(title)
            + labs(x='', y=ylab))


def _render(p, name, w, h):
    path = os.path.join(SCRATCH, name)
    p.save(path, width=w, height=h, dpi=DPI, verbose=False, limitsize=False)
    return np.asarray(Image.open(path))


YLAB = {'pit': 'PIT  $u = F_{amortiser}(\\rho_{SVI} \\mid x,\\, z)$',          # conditional calibration
        'marg': 'marginal PIT  $u = F^{\\mathrm{marg}}_{amortiser}(\\rho_{SVI})$'}  # marginal calibration


def _build(kind):
    main_pd = pickle.load(open(os.path.join(MAIN, f"{PFX}_pps_RAGD_plotdata.pkl"), 'rb'))
    affn_pd = pickle.load(open(os.path.join(AFFN, f"{PFX}_pps_RAGD_plotdata.pkl"), 'rb'))
    key = f'{kind}_u'
    items = list(main_pd[f'RAGD::{key}']['item_label'].unique())
    panels = [
        ('SVI fit (reference)', _svi_reference(items)),
        ('ADM with fine tune + BvM calibration', _pitbox(main_pd[f'RAGD::{key}'])),
        ('ADM with fine tune, no BvM calibration', _pitbox(affn_pd[f'RAGD::{key}'])),
        ('ADM without fine tune + BvM calibration', _pitbox(main_pd[f'RAGDbase::{key}'])),
    ]
    imgs = [_render(_panel(df, t, YLAB[kind], kind), f"supp_{kind}_{i}.png", 6.6, 6.4)
            for i, (t, df) in enumerate(panels)]

    fig = plt.figure(figsize=(14.2, 13.4))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.08], hspace=0.06, wspace=0.03,
                          left=0.004, right=0.9, top=0.99, bottom=0.01)
    for ax, img in zip([fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]),
                        fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])], imgs):
        ax.imshow(img, aspect='auto'); ax.axis('off')
    cax = fig.add_axes([0.92, 0.33, 0.015, 0.34])
    if kind == 'pit':
        cmap = mcolors.LinearSegmentedColormap.from_list('pit', [LOW, MID, HIGH])
        cb = fig.colorbar(ScalarMappable(norm=mcolors.Normalize(-0.3, 0.3), cmap=cmap), cax=cax, extend='both')
        cb.set_label('Median of\nPIT diagnostic\n(target = 0.0)', fontsize=10)
    else:
        cmap = mcolors.LinearSegmentedColormap.from_list('marg', [MID, HIGH])
        cb = fig.colorbar(ScalarMappable(norm=mcolors.Normalize(0, 0.3), cmap=cmap), cax=cax, extend='max')
        cb.set_label('marginal KS\n(SVI vs amortiser\nmarginal; target 0)', fontsize=10)
    cb.ax.tick_params(labelsize=8)

    name = 'icrc_supp_figure' if kind == 'pit' else 'icrc_supp_figure_marg'
    out = os.path.join(os.environ.get('PAPER_ADM_OUT', os.path.dirname(os.path.abspath(__file__))), name)
    fig.savefig(out + '.pdf', dpi=DPI, pad_inches=0.02); fig.savefig(out + '.png', dpi=DPI, pad_inches=0.02)
    plt.close(fig)
    print(f"Saved {name} -> {out}.pdf / .png")


def main():
    for kind in ('pit', 'marg'):
        _build(kind)


if __name__ == "__main__":
    main()
