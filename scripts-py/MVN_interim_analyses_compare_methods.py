#!/usr/bin/env python3
"""
MVN interim analyses: cross-method comparison (§3.3 of
``dev/amortised_decision_making.md``).

Standalone. Loads per-J per-method artifacts from:

  - HMC nested-MC : ``MVN_interim_analyses_with_nested_MC_hmc.py``
  - IS reweighting : ``MVN_interim_analyses_with_IS_from_x.py``
  - Regression endpt-x (Gaussian approx / quantile / mquantile) : the RGE
    family of scripts.
  - Amortiser (features-fixed / features-MLP / xcomp / xcompAtt): the four
    ``MVN_interim_analysis_amortise_...`` scripts.

Family split (matches the Binomial compare-methods layout):

  Regression-family plots (amortiser dropped to reduce clutter):
    1. ``mvn_J{J}_compare_methods_p_h1_xz_all.pdf``
    2. ``mvn_J{J}_compare_methods_pps_all.pdf``

  Amortiser-focus plots (best regression baseline + amortiser variants):
    3. ``mvn_J{J}_compare_methods_p_h1_xz_all_amortised.pdf``
    4. ``mvn_J{J}_compare_methods_pps_all_amortised.pdf``

  Shared timing plots (both amortiser + non-amortiser methods, one bar per
  method; leftmost x-axis slot ``training`` shows the one-off training
  time for amortiser methods; per-interim slots show deployment time only):
    5. ``mvn_J{J}_compare_methods_timing.pdf``
    6. ``mvn_J{J}_compare_methods_timing_amortised.pdf``  (train + deploy
       stacked bar per method with two fill colours, coord_flip)

  Cross-J summary:
    7. ``mvn_compare_methods_mse.pdf``
    8. ``mvn_compare_methods_rge_stats_box.pdf``

  IS diagnostic:
    9. ``mvn_J{J}_compare_methods_ess.pdf``

Legend layouts single-column bottom; method → colour mapping is consistent
across every plot in the file. Analytic PPS bar drawn WHITE with thin
BLACK outline (size = 1.1) and forced first (leftmost) inside each dodge
group so the ground-truth reference reads as distinct from Monte-Carlo
estimates.

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/MVN_interim_analyses_compare_methods.py
"""

# %%

import os
import pickle
import sys
import warnings
from pathlib import Path

try:
    project_root = Path(__file__).resolve().parent.parent
except NameError:
    project_root = Path.cwd()
    if project_root.name == 'scripts-py':
        project_root = project_root.parent

sys.path.insert(0, str(project_root / 'python'))

import numpy as np
import pandas as pd
from plotnine import (
    ggplot, aes, coord_flip, geom_boxplot, geom_col, geom_errorbar, geom_hline,
    geom_text, facet_grid, facet_wrap,
    scale_colour_manual, scale_fill_manual, scale_y_continuous, scale_y_sqrt,
    position_dodge, position_stack, theme_bw, theme,
    element_text, element_blank, labs, guides, guide_legend,
)

warnings.filterwarnings('ignore')

from utils import _futurama_palette

print("Imports successful")

# %%

# =============================================================================
# Configuration.
# =============================================================================

_sandbox = "/Users/or105/sandbox/bIRTistic"
DIR_SIM = os.path.join(_sandbox, "py-mvn-interim-simulations-260609")
DIR_HMC = os.path.join(_sandbox, "py-mvn-interim-with-nested-MC-hmc-260609")
DIR_IS  = os.path.join(_sandbox, "py-mvn-interim-with-IS-260611")
DIR_RGE  = os.path.join(_sandbox, "py-mvn-interim-with-regression-on-endptx-wz-260611")
DIR_RGEQ = os.path.join(_sandbox, "py-mvn-interim-with-regression-on-endptx-wz-quantile-regr-260702")
DIR_RGEM = os.path.join(_sandbox, "py-mvn-interim-with-regression-on-endptx-wz-mquantile-regr-260702")
DIR_RGEA = os.path.join(
    _sandbox,
    "py-mvn-interim-amortise-endptx-on-wz-with-features-fixed-"
    "idcomp-qpsi-MLP-loss-multiquantilehead-260714",
)
DIR_RGEB = os.path.join(
    _sandbox,
    "py-mvn-interim-amortise-endptx-on-wz-with-features-MLP-"
    "idcomp-qpsi-MLP-loss-multiquantilehead-260714",
)
DIR_RGEC = os.path.join(
    _sandbox,
    "py-mvn-interim-amortise-endptx-on-wz-with-features-MLP-"
    "xcomp-qpsi-MLP-loss-multiquantilehead-260715",
)
DIR_RGED = os.path.join(
    _sandbox,
    "py-mvn-interim-amortise-endptx-on-wz-with-features-"
    "itemScompAtt-qpsi-MLP-loss-multiquantilehead_15k_260716",
)
DIR_RGEF = os.path.join(
    _sandbox,
    "py-mvn-interim-amortise-endptx-on-wz-with-features-"
    "itemXcompAtt-qpsi-MLP-loss-multiquantilehead_15k_260716",
)
DIR_RGDS = os.path.join(
    _sandbox,
    "py-mvn-interim-amortise-endptx-on-wz-with-features-"
    "deepsetScompAtt-qpsi-MLP-loss-multiquantilehead_260803",
)
DIR_RGDX = os.path.join(
    _sandbox,
    "py-mvn-interim-amortise-endptx-on-wz-with-features-"
    "deepsetXcompAtt-qpsi-MLP-loss-multiquantilehead_260803",
)
DIR_OUT  = os.path.join(_sandbox, "py-mvn-interim-compare-methods-260609")
os.makedirs(DIR_OUT, exist_ok=True)

# Method labels. Amortiser labels use the filename convention so the legend
# is unambiguous when the family plot vs the amortiser-focus plot are read
# side-by-side.
LBL_HMC   = 'nested-MC using HMC for each (x,z)'
LBL_IS    = 'IS reweighting of theta|x'
LBL_RGE   = 'Regression of endpt-x on w(z) using Gaussian approx'
LBL_RGEQ  = 'Regression of endpt-x on w(z) - quantile'
LBL_RGEM  = 'Regression of endpt-x on w(z) - mquantile'
LBL_RGEA  = 'Amortiser-features-fixed-idcomp-qpsi-MLP-loss-multiquantilehead'
LBL_RGEB  = 'Amortiser-features-MLP-idcomp-qpsi-MLP-loss-multiquantilehead'
LBL_RGEC  = 'Amortiser-features-MLP-xcomp-qpsi-MLP-loss-multiquantilehead'
LBL_RGED  = 'Amortiser-features-itemScompAtt-qpsi-MLP-loss-multiquantilehead'
LBL_RGEF  = 'Amortiser-features-itemXcompAtt-qpsi-MLP-loss-multiquantilehead'
LBL_RGDS  = 'Amortiser-features-deepsetScompAtt-qpsi-MLP-loss-multiquantilehead'
LBL_RGDX  = 'Amortiser-features-deepsetXcompAtt-qpsi-MLP-loss-multiquantilehead'

non_amort_methods = [LBL_HMC, LBL_IS, LBL_RGE, LBL_RGEQ, LBL_RGEM]
# RGEB (features-MLP idcomp) dropped: undertrained on MVN scale (§13.2).
amort_focus_methods = [LBL_RGE, LBL_RGEA, LBL_RGEC, LBL_RGED, LBL_RGEF, LBL_RGDS, LBL_RGDX]
timing_methods = non_amort_methods + [LBL_RGEA, LBL_RGEC, LBL_RGED, LBL_RGEF, LBL_RGDS, LBL_RGDX]

# Ordered union used to derive the global palette so the same method reuses
# the same colour on every plot.
_all_methods = non_amort_methods + [LBL_RGEA, LBL_RGEC, LBL_RGED, LBL_RGEF, LBL_RGDS, LBL_RGDX]
method_colours = dict(zip(_all_methods, _futurama_palette(len(_all_methods))))
# Analytic drawn as white fill + black outline so it reads as ground-truth
# reference, not a competing estimate.
pps_colours = {'analytic': '#FFFFFF', **method_colours}
pps_outline_colours = {'analytic': '#000000',
                       **{m: c for m, c in method_colours.items()}}

# Fill colours for the stacked train / deploy timing bar.
_TRAIN_DEPLOY_COLOURS = {
    'train (one-off)':      '#1f77b4',
    'deploy (all interims)': '#ff7f0e',
}


def _load_training_mins(ckpt_path):
    if not os.path.exists(ckpt_path):
        return float('nan')
    with open(ckpt_path, 'rb') as f:
        payload = pickle.load(f)
    return float(payload.get('training_mins', float('nan')))


training_mins_by_method = {
    LBL_RGEA: _load_training_mins(
        os.path.join(DIR_RGEA, 'mvn_interim_amortised_pps_net.pkl')
    ),
    LBL_RGEC: _load_training_mins(
        os.path.join(DIR_RGEC, 'mvn_interim_amortised_pps_net.pkl')
    ),
    LBL_RGED: _load_training_mins(
        os.path.join(DIR_RGED, 'mvn_interim_amortised_pps_net.pkl')
    ),
    LBL_RGEF: _load_training_mins(
        os.path.join(DIR_RGEF, 'mvn_interim_amortised_pps_net.pkl')
    ),
    LBL_RGDS: _load_training_mins(
        os.path.join(DIR_RGDS, 'mvn_interim_amortised_pps_net.pkl')
    ),
    LBL_RGDX: _load_training_mins(
        os.path.join(DIR_RGDX, 'mvn_interim_amortised_pps_net.pkl')
    ),
}
for lbl in non_amort_methods:
    training_mins_by_method.setdefault(lbl, 0.0)
print(f"Training mins per method: "
      f"{ {k: round(v, 3) for k, v in training_mins_by_method.items()} }")

# %%

# =============================================================================
# Load upstream artifacts.
# =============================================================================

sim_data = pd.read_pickle(os.path.join(DIR_SIM, 'mvn_sim_data.pkl'))
simu_params = sim_data['simu_params']
di = sim_data['interim_grid']
J_GRID = list(simu_params['J_grid'])
K_levels = simu_params['K_levels']
pps_ProbH1_target_lwr_quantile = float(simu_params.get(
    'pps_ProbH1_target_lwr_quantile',
    simu_params.get('pps_ProbH1_thresh', 0.89),
))

pps_cf = pd.read_pickle(
    os.path.join(DIR_SIM, 'mvn_pps_closed_form.pkl'),
)['pps_cf']

hmc_artifacts = pd.read_pickle(
    os.path.join(DIR_HMC, 'mvn_pps_nested_mc.pkl'),
)
hmc_p_h1_xz_all = hmc_artifacts['p_h1_xz']        # has J column
hmc_pps_all     = hmc_artifacts['pps']
hmc_timing_all  = hmc_artifacts['timing']

print(f"Loaded analytic ({len(pps_cf)}) + HMC nested-MC artifacts."
      f" J_GRID = {J_GRID}.")


# =============================================================================
# Plot helpers (parallel to the Binomial compare-methods script).
# =============================================================================


def _boxplot_p_h1_xz(box_stats, methods, response_label_cats,
                     n_per_level, pdf_path):
    """Faceted dodged boxplot of ``p_h1_xz`` per (response, interim, method)."""
    keep = [m for m in methods if m != LBL_RGEQ]
    d = box_stats[box_stats['method'].isin(keep)].copy()
    d['method'] = pd.Categorical(
        d['method'].astype(str), categories=keep, ordered=True,
    )
    p = (
        ggplot(d, aes(x='interim_month_year', fill='method'))
        + geom_boxplot(
            aes(ymin='q025', lower='q25', middle='q50', upper='q75',
                ymax='q975', group='grp'),
            stat='identity',
            position=position_dodge(width=0.8, preserve='single'),
            width=0.7, colour='#404040', size=0.3,
        )
        + geom_hline(yintercept=pps_ProbH1_target_lwr_quantile,
                     colour='black', size=1.0)
        + scale_fill_manual(values=method_colours,
                            breaks=keep, limits=keep)
        + scale_y_continuous(
            limits=[0, 1],
            breaks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            labels=['0%', '20%', '40%', '60%', '80%', '100%'],
        )
        + facet_wrap('~ response_label', ncol=n_per_level, dir='h')
        + guides(fill=guide_legend(ncol=1))
        + theme_bw()
        + theme(
            axis_text_x=element_text(angle=45, vjust=1, hjust=1),
            legend_position='bottom',
            legend_direction='vertical',
            figure_size=(3.0 * n_per_level, 2.5 * K_levels + 2),
            strip_background=element_blank(),
            strip_text=element_text(face='bold'),
        )
        + labs(x='Interim', y='p(H_1j | x, z) for predicted z samples',
               fill='method')
    )
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"Saved p(H_1 | x, z) boxplot to {pdf_path}")


def _pps_bars(pps_g_ci, methods, n_per_level, pdf_path, width_scale=1.0):
    """Faceted PPS bars per (response, interim, method) with bootstrap CI.

    Analytic bar white + thin black outline (size = 1.1), forced first in
    the dodge order so it lands leftmost inside each group. ``width_scale``
    multiplies the figure width (default 1.0; caller passes >1 when a
    single facet needs to fit more dodged bars).
    """
    keep_with_analytic = ['analytic', *methods]
    d = pps_g_ci[pps_g_ci['method'].isin(keep_with_analytic)].copy()
    # Re-apply Categorical AFTER any pandas ops so the level ordering is
    # preserved end-to-end.
    d['method'] = pd.Categorical(
        d['method'].astype(str),
        categories=keep_with_analytic, ordered=True,
    )
    d = d.sort_values(
        ['response_label', 'interim_month_year', 'method']
    ).reset_index(drop=True)
    p = (
        ggplot(d, aes(x='interim_month_year', y='pps', fill='method',
                      colour='method', group='method'))
        + geom_col(position=position_dodge(width=0.8, preserve='single'),
                   width=0.7, size=1.1)
        + geom_errorbar(
            aes(x='interim_month_year', ymin='q025_bs', ymax='q975_bs',
                group='method'),
            position=position_dodge(width=0.8, preserve='single'),
            width=0.3, colour='black', size=0.4, inherit_aes=False,
        )
        + geom_hline(yintercept=pps_ProbH1_target_lwr_quantile,
                     colour='black', size=1.0, linetype='dashed')
        # Fill + colour scales share the same ``name`` so plotnine merges
        # them into ONE legend showing both aesthetics per entry -- the
        # legend swatch for 'analytic' picks up its black outline.
        + scale_fill_manual(values=pps_colours,
                            breaks=keep_with_analytic,
                            limits=keep_with_analytic,
                            name='method')
        + scale_colour_manual(values=pps_outline_colours,
                              breaks=keep_with_analytic,
                              limits=keep_with_analytic,
                              name='method')
        + scale_y_continuous(
            limits=[0, 1.05],
            breaks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            labels=['0%', '20%', '40%', '60%', '80%', '100%'],
        )
        + facet_wrap('~ response_label', ncol=n_per_level, dir='h')
        + guides(
            fill=guide_legend(ncol=1),
            colour=guide_legend(ncol=1),
        )
        + theme_bw()
        + theme(
            axis_text_x=element_text(angle=45, vjust=1, hjust=1),
            legend_position='bottom',
            legend_direction='vertical',
            figure_size=(3.0 * n_per_level * width_scale,
                         2.5 * K_levels + 2),
            strip_background=element_blank(),
            strip_text=element_text(face='bold'),
        )
        + labs(x='Interim',
               y='PPS = int P(p(H_1j | x, z) > eta) dz', fill='method')
    )
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"Saved PPS bars plot to {pdf_path}")


def _timing_bars(timing_df, methods, interim_order_local, pdf_path,
                 figure_size=(15, 8)):
    """Two-panel-like layout on ONE x-axis:

      - Leftmost slot ``training``: one-off amortiser training time.
        Non-amortised methods have zero here (bar suppressed).
      - Interim slots: per-interim DEPLOY time only.
    """
    d = timing_df[timing_df['method'].isin(methods)][
        ['interim_month_year', 'method', 'mins']
    ].copy()
    d['slot'] = d['interim_month_year'].astype(str)
    d = d[['slot', 'method', 'mins']]
    train_rows = pd.DataFrame({
        'slot':   ['training'] * len(methods),
        'method': methods,
        'mins':   [training_mins_by_method.get(m, 0.0) for m in methods],
    })
    slot_order = ['training'] + interim_order_local
    long = pd.concat([train_rows, d], ignore_index=True)
    full = pd.MultiIndex.from_product(
        [methods, slot_order], names=['method', 'slot'],
    ).to_frame(index=False)
    long = full.merge(long, on=['method', 'slot'], how='left').fillna({'mins': 0.0})
    long['slot'] = pd.Categorical(
        long['slot'], categories=slot_order, ordered=True,
    )
    long['method'] = pd.Categorical(
        long['method'], categories=methods, ordered=True,
    )
    long['value_label'] = long['mins'].map(
        lambda v: f"{v:.2f}" if v > 0 else "",
    )
    p = (
        ggplot(long, aes(x='slot', y='mins', fill='method'))
        + geom_col(position=position_dodge(width=0.8, preserve='single'),
                   width=0.7)
        + geom_text(
            aes(label='value_label'),
            position=position_dodge(width=0.8, preserve='single'),
            size=6, angle=30, ha='left', va='bottom',
        )
        + scale_fill_manual(values=method_colours,
                            breaks=methods, limits=methods)
        + scale_y_sqrt(expand=(0, 0, 0.15, 0))
        + guides(fill=guide_legend(ncol=1))
        + theme_bw()
        + theme(
            axis_text_x=element_text(angle=45, vjust=1, hjust=1),
            legend_position='bottom',
            legend_direction='vertical',
            figure_size=figure_size,
        )
        + labs(x='training (one-off) | deployment (per interim)',
               y='time (mins)',
               fill='method')
    )
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"Saved timing plot to {pdf_path}")


def _train_deploy_stacked_bars(timing_df, methods, pdf_path,
                               figure_size=(11, 5)):
    """Horizontal stacked bar per method: train segment + total deploy
    across all interims. Two fill colours (train vs deploy). coord_flip
    so long labels stay readable."""
    d = timing_df[timing_df['method'].isin(methods)].copy()
    deploy_total = (
        d.groupby('method', observed=True)['mins'].sum().reset_index()
        .rename(columns={'mins': 'value'})
    )
    deploy_total['segment'] = 'deploy (all interims)'
    train_total = pd.DataFrame({
        'method':  methods,
        'value':   [training_mins_by_method.get(m, 0.0) for m in methods],
        'segment': ['train (one-off)'] * len(methods),
    })
    stacked = pd.concat([train_total, deploy_total], ignore_index=True)
    stacked['method'] = pd.Categorical(
        stacked['method'], categories=list(reversed(methods)), ordered=True,
    )
    stacked['segment'] = pd.Categorical(
        stacked['segment'],
        categories=['train (one-off)', 'deploy (all interims)'],
        ordered=True,
    )
    stacked['value_label'] = stacked['value'].map(
        lambda v: f"{v:.2f}" if v > 0 else "",
    )
    p = (
        ggplot(stacked, aes(x='method', y='value', fill='segment'))
        + geom_col(position=position_stack(reverse=True), width=0.6)
        + geom_text(
            aes(label='value_label'),
            position=position_stack(vjust=0.5, reverse=True),
            size=8, colour='black',
        )
        + scale_fill_manual(values=_TRAIN_DEPLOY_COLOURS)
        + guides(fill=guide_legend(ncol=1))
        + coord_flip()
        + theme_bw()
        + theme(
            legend_position='bottom',
            legend_direction='vertical',
            figure_size=figure_size,
        )
        + labs(x='method',
               y='wall-clock time (mins) across the full interim schedule',
               fill='segment')
    )
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"Saved train + deploy stacked timing plot to {pdf_path}")


# %%

# =============================================================================
# Per-J cross-method comparison.
# =============================================================================

mse_rows = []                  # per (J, method, interim) MSE rows
for J in J_GRID:
    block_size_J = J // K_levels
    n_per_level = min(5, block_size_J)
    j_indices_all = np.array([
        k * block_size_J + c for k in range(K_levels)
        for c in range(n_per_level)
    ])
    response_label_cats = [f'response mu_{j}' for j in j_indices_all]
    print(f"\n{'#' * 70}\n# J = {J} (responses shown: {len(j_indices_all)})"
          f"\n{'#' * 70}")

    # ---- Analytic per-z ----
    an_p = pd.read_pickle(
        os.path.join(DIR_HMC, f'mvn_J{J}_pps_p_h1_xz_analytic.pkl'),
    )[['interim_id', 'interim_date', 'interim_month_year', 'j', 's',
       'p_h1_xz']].copy()

    # ---- HMC per-J slices ----
    hmc_p = hmc_p_h1_xz_all[hmc_p_h1_xz_all['J'] == J][[
        'interim_id', 'interim_date', 'interim_month_year', 'j', 's',
        'p_h1_xz',
    ]].copy()
    hmc_pps = hmc_pps_all[hmc_pps_all['J'] == J][[
        'interim_id', 'interim_date', 'interim_month_year', 'j', 'pps',
    ]].copy()
    hmc_timing = hmc_timing_all[hmc_timing_all['J'] == J][[
        'interim_id', 'interim_date', 'interim_month_year', 'mins',
    ]].drop_duplicates().copy()

    # ---- IS per-J slices ----
    is_p = pd.read_pickle(
        os.path.join(DIR_IS, f'mvn_J{J}_pps_IS_p_h1_xz.pkl'),
    )[['interim_id', 'interim_date', 'interim_month_year', 'j', 's',
       'p_h1_xz']].copy()
    is_pps = pd.read_csv(
        os.path.join(DIR_IS, f'mvn_J{J}_pps_IS.csv'),
    )[['interim_id', 'interim_date', 'interim_month_year', 'j',
       'pps']].copy()
    is_timing = pd.read_csv(
        os.path.join(DIR_IS, f'mvn_J{J}_pps_IS_timing.csv'),
    )[['interim_id', 'interim_date', 'interim_month_year',
       'mins_interim_id']].rename(columns={'mins_interim_id': 'mins'}).copy()
    is_perf = pd.read_csv(
        os.path.join(DIR_IS, f'mvn_J{J}_pps_IS_perf.csv'),
    )[['interim_id', 'interim_date', 'interim_month_year', 's',
       'ess', 'ess_over_n']].copy()

    def _read_slice(d, suffix, is_pkl=True):
        p_ = pd.read_pickle(
            os.path.join(d, f'mvn_J{J}_pps_{suffix}_p_h1_xz.pkl'),
        )[['interim_id', 'interim_date', 'interim_month_year', 'j', 's',
           'p_h1_xz']].copy()
        pps_ = pd.read_csv(
            os.path.join(d, f'mvn_J{J}_pps_{suffix}.csv'),
        )[['interim_id', 'interim_date', 'interim_month_year', 'j',
           'pps']].copy()
        t_ = pd.read_csv(
            os.path.join(d, f'mvn_J{J}_pps_{suffix}_timing.csv'),
        )[['interim_id', 'interim_date', 'interim_month_year',
           'mins_interim_id']].rename(
               columns={'mins_interim_id': 'mins'},
        ).copy()
        return p_, pps_, t_

    rge_p,  rge_pps,  rge_timing  = _read_slice(DIR_RGE,  'RGE')
    rgeq_p, rgeq_pps, rgeq_timing = _read_slice(DIR_RGEQ, 'RGEQ')
    rgem_p, rgem_pps, rgem_timing = _read_slice(DIR_RGEM, 'RGEM')
    rgea_p, rgea_pps, rgea_timing = _read_slice(DIR_RGEA, 'RGEA')
    rgeb_p, rgeb_pps, rgeb_timing = _read_slice(DIR_RGEB, 'RGEB')
    rgec_p, rgec_pps, rgec_timing = _read_slice(DIR_RGEC, 'RGEC')
    rged_p, rged_pps, rged_timing = _read_slice(DIR_RGED, 'RGED')
    rgef_p, rgef_pps, rgef_timing = _read_slice(DIR_RGEF, 'RGEF')

    def _read_slice_optional(d, suffix):
        try:
            return _read_slice(d, suffix)
        except (FileNotFoundError, OSError):
            return None, None, None

    rgds_p, rgds_pps, rgds_timing = _read_slice_optional(DIR_RGDS, 'RGDS')
    rgdx_p, rgdx_pps, rgdx_timing = _read_slice_optional(DIR_RGDX, 'RGDX')

    p_h1_xz = pd.concat(
        [
            an_p  .assign(method='analytic'),
            hmc_p .assign(method=LBL_HMC),
            is_p  .assign(method=LBL_IS),
            rge_p .assign(method=LBL_RGE),
            rgeq_p.assign(method=LBL_RGEQ),
            rgem_p.assign(method=LBL_RGEM),
            rgea_p.assign(method=LBL_RGEA),
            rgeb_p.assign(method=LBL_RGEB),
            rgec_p.assign(method=LBL_RGEC),
            rged_p.assign(method=LBL_RGED),
            rgef_p.assign(method=LBL_RGEF),
            *([rgds_p.assign(method=LBL_RGDS)] if rgds_p is not None else []),
            *([rgdx_p.assign(method=LBL_RGDX)] if rgdx_p is not None else []),
        ],
        ignore_index=True,
    )
    pps = pd.concat(
        [
            pps_cf[pps_cf['J'] == J][[
                'interim_id', 'interim_date', 'interim_month_year', 'j', 'pps',
            ]].assign(method='analytic'),
            hmc_pps .assign(method=LBL_HMC),
            is_pps  .assign(method=LBL_IS),
            rge_pps .assign(method=LBL_RGE),
            rgeq_pps.assign(method=LBL_RGEQ),
            rgem_pps.assign(method=LBL_RGEM),
            rgea_pps.assign(method=LBL_RGEA),
            rgeb_pps.assign(method=LBL_RGEB),
            rgec_pps.assign(method=LBL_RGEC),
            rged_pps.assign(method=LBL_RGED),
            rgef_pps.assign(method=LBL_RGEF),
            *([rgds_pps.assign(method=LBL_RGDS)] if rgds_pps is not None else []),
            *([rgdx_pps.assign(method=LBL_RGDX)] if rgdx_pps is not None else []),
        ],
        ignore_index=True,
    )
    timing = pd.concat(
        [
            hmc_timing .assign(method=LBL_HMC),
            is_timing  .assign(method=LBL_IS),
            rge_timing .assign(method=LBL_RGE),
            rgeq_timing.assign(method=LBL_RGEQ),
            rgem_timing.assign(method=LBL_RGEM),
            rgea_timing.assign(method=LBL_RGEA),
            rgeb_timing.assign(method=LBL_RGEB),
            rgec_timing.assign(method=LBL_RGEC),
            rged_timing.assign(method=LBL_RGED),
            rgef_timing.assign(method=LBL_RGEF),
            *([rgds_timing.assign(method=LBL_RGDS)] if rgds_timing is not None else []),
            *([rgdx_timing.assign(method=LBL_RGDX)] if rgdx_timing is not None else []),
        ],
        ignore_index=True,
    )

    # Normalise mixed Timestamp/str dtypes.
    for df in (p_h1_xz, pps, timing, is_perf):
        df['interim_date'] = pd.to_datetime(df['interim_date'])

    interim_order = (
        p_h1_xz[['interim_date', 'interim_month_year']]
        .drop_duplicates()
        .sort_values('interim_date')['interim_month_year']
        .tolist()
    )

    for df in (p_h1_xz, pps, timing, is_perf):
        df['interim_month_year'] = pd.Categorical(
            df['interim_month_year'], categories=interim_order, ordered=True,
        )
    p_h1_xz['method'] = pd.Categorical(
        p_h1_xz['method'],
        categories=['analytic'] + _all_methods, ordered=True,
    )
    pps['method'] = pd.Categorical(
        pps['method'], categories=['analytic'] + _all_methods, ordered=True,
    )
    timing['method'] = pd.Categorical(
        timing['method'], categories=_all_methods, ordered=True,
    )

    p_h1_xz_g = p_h1_xz[p_h1_xz['j'].isin(j_indices_all)].copy()
    pps_g     = pps[pps['j'].isin(j_indices_all)].copy()
    for df in (p_h1_xz_g, pps_g):
        df['response_label'] = pd.Categorical(
            'response mu_' + df['j'].astype(str),
            categories=response_label_cats, ordered=True,
        )

    # ---- Box-stats + PPS bootstrap CI computed once, used across the
    # ---- regression-family + amortiser-focus variants of the plot.
    p_h1_xz_g_box = p_h1_xz_g[p_h1_xz_g['method'] != LBL_RGEQ].copy()
    p_h1_xz_g_box['method'] = (
        p_h1_xz_g_box['method'].cat.remove_unused_categories()
    )
    box_stats = (
        p_h1_xz_g_box
        .groupby(['response_label', 'interim_month_year', 'method'],
                 observed=True)['p_h1_xz']
        .agg(
            q025=lambda s: s.quantile(0.025),
            q25 =lambda s: s.quantile(0.25),
            q50 =lambda s: s.quantile(0.50),
            q75 =lambda s: s.quantile(0.75),
            q975=lambda s: s.quantile(0.975),
        )
        .reset_index()
    )
    box_stats['grp'] = (
        box_stats['interim_month_year'].astype(str) + '|'
        + box_stats['method'].astype(str)
    )
    box_stats.to_pickle(os.path.join(
        DIR_OUT, f'mvn_J{J}_compare_methods_p_h1_xz_all.pkl',
    ))

    bs_B = 2000
    rng_bs = np.random.default_rng(simu_params['seed'])
    bs_parts = []
    p_h1_xz_g_mc = p_h1_xz_g[p_h1_xz_g['method'] != 'analytic']
    for m_label, g in p_h1_xz_g_mc.groupby('method', observed=True):
        wide = (
            g.pivot_table(index=['j', 'interim_id', 'interim_date'],
                          columns='s', values='p_h1_xz')
            .sort_index()
        )
        Parr = wide.to_numpy()
        n_rows, S = Parr.shape
        idx = rng_bs.integers(0, S, size=(n_rows, bs_B, S))
        Pbs = Parr[np.arange(n_rows)[:, None, None], idx]
        pps_bs = (Pbs > pps_ProbH1_target_lwr_quantile).mean(axis=2)
        q = np.quantile(pps_bs, [0.025, 0.975], axis=1).T
        part = pd.DataFrame(q, columns=['q025_bs', 'q975_bs'])
        part['j']          = wide.index.get_level_values('j').to_numpy()
        part['interim_id'] = wide.index.get_level_values('interim_id').to_numpy()
        part['method'] = m_label
        bs_parts.append(part)
    bs = pd.concat(bs_parts, ignore_index=True)
    pps_g_ci = pps_g.merge(
        bs, on=['j', 'interim_id', 'method'], how='left',
    )
    pps_g_ci.to_pickle(os.path.join(
        DIR_OUT, f'mvn_J{J}_compare_methods_pps_all.pkl',
    ))

    # ---- Regression-family plots (amortiser excluded) ----
    _boxplot_p_h1_xz(
        box_stats, non_amort_methods, response_label_cats, n_per_level,
        os.path.join(DIR_OUT, f'mvn_J{J}_compare_methods_p_h1_xz_all.pdf'),
    )
    _pps_bars(
        pps_g_ci, non_amort_methods, n_per_level,
        os.path.join(DIR_OUT, f'mvn_J{J}_compare_methods_pps_all.pdf'),
        width_scale=1.4,
    )

    # ---- Amortiser-focus plots (best regression + amortiser variants) ----
    _boxplot_p_h1_xz(
        box_stats, amort_focus_methods, response_label_cats, n_per_level,
        os.path.join(
            DIR_OUT, f'mvn_J{J}_compare_methods_p_h1_xz_all_amortised.pdf',
        ),
    )
    _pps_bars(
        pps_g_ci, amort_focus_methods, n_per_level,
        os.path.join(
            DIR_OUT, f'mvn_J{J}_compare_methods_pps_all_amortised.pdf',
        ),
        width_scale=1.4,
    )

    # ---- Timing: leftmost training slot + per-interim deploy only ----
    _timing_bars(
        timing, timing_methods, interim_order,
        os.path.join(DIR_OUT, f'mvn_J{J}_compare_methods_timing.pdf'),
    )
    _train_deploy_stacked_bars(
        timing, amort_focus_methods,
        os.path.join(
            DIR_OUT, f'mvn_J{J}_compare_methods_timing_amortised.pdf',
        ),
    )

    # ---- IS ESS / particle per interim ----
    is_ess = (
        is_perf.groupby('interim_month_year', observed=True)
        .agg(value=('ess_over_n', 'median'))
        .reset_index()
        .assign(method=LBL_IS)
    )
    is_ess['method'] = pd.Categorical(
        is_ess['method'], categories=_all_methods, ordered=True,
    )
    p = (
        ggplot(is_ess, aes(x='interim_month_year', y='value', fill='method'))
        + geom_col(width=0.7)
        + scale_fill_manual(values=method_colours)
        + guides(fill=guide_legend(ncol=1))
        + theme_bw()
        + theme(
            axis_text_x=element_text(angle=45, vjust=1, hjust=1),
            legend_position='bottom',
            legend_direction='vertical',
            figure_size=(12, 6),
            strip_background=element_blank(),
            strip_text=element_text(face='bold'),
        )
        + labs(x='Interim', y='ESS / particle', fill='method')
    )
    pdf_path = os.path.join(DIR_OUT, f'mvn_J{J}_compare_methods_ess.pdf')
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved IS ESS / particle plot to {pdf_path}")

    # ---- MSE: per (interim, method) average over all J components ----
    pps_pivot = (
        pps.pivot_table(
            index=['interim_id', 'interim_date', 'interim_month_year', 'j'],
            columns='method', values='pps',
        ).reset_index()
    )
    for method in _all_methods:
        if method not in pps_pivot.columns:
            continue
        se = (pps_pivot[method] - pps_pivot['analytic']) ** 2
        mse = (
            pd.DataFrame({
                'interim_id':         pps_pivot['interim_id'],
                'interim_date':       pps_pivot['interim_date'],
                'interim_month_year': pps_pivot['interim_month_year'],
                'se':                 se,
            })
            .groupby(['interim_id', 'interim_date', 'interim_month_year'],
                     observed=True)['se']
            .mean().reset_index(name='mse')
        )
        mse['J'] = J
        mse['method'] = method
        mse_rows.append(mse)

# %%

# =============================================================================
# Cross-J MSE plot: x = interim date, y = MSE (over responses j),
# fill = method, facet = J.
# =============================================================================

mse_all = pd.concat(mse_rows, ignore_index=True)
mse_all['method'] = pd.Categorical(
    mse_all['method'], categories=_all_methods, ordered=True,
)
mse_all['J_label'] = pd.Categorical(
    'J=' + mse_all['J'].astype(str),
    categories=[f'J={J}' for J in J_GRID], ordered=True,
)
_mse_interim_order = (
    mse_all[['interim_date', 'interim_month_year']]
    .drop_duplicates()
    .sort_values('interim_date')['interim_month_year'].tolist()
)
mse_all['interim_month_year'] = pd.Categorical(
    mse_all['interim_month_year'],
    categories=_mse_interim_order, ordered=True,
)
mse_all.to_pickle(os.path.join(DIR_OUT, 'mvn_compare_methods_mse.pkl'))

p = (
    ggplot(mse_all,
           aes(x='interim_month_year', y='mse', fill='method'))
    + geom_col(position=position_dodge(width=0.8, preserve='single'),
               width=0.7)
    + facet_wrap('~ J_label', ncol=1, scales='free_y')
    + scale_fill_manual(values=method_colours,
                        breaks=_all_methods, limits=_all_methods)
    + scale_y_sqrt(expand=(0, 0, 0.05, 0))
    + guides(fill=guide_legend(ncol=1))
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        legend_position='bottom',
        legend_direction='vertical',
        figure_size=(14, 4 * len(J_GRID) + 2),
        strip_background=element_blank(),
        strip_text=element_text(face='bold'),
    )
    + labs(x='Interim date',
           y='MSE of PPS vs analytic (mean over responses j; sqrt-scale)',
           fill='method')
)
pdf_path = os.path.join(DIR_OUT, 'mvn_compare_methods_mse.pdf')
p.save(pdf_path, verbose=False, limitsize=False)
print(f"Saved cross-J MSE plot to {pdf_path}")

# %%

# =============================================================================
# Regression-only: violin plot of per-response R^2 (top row) and empirical
# rho (bottom row) across all J responses j, faceted by J in columns.
# =============================================================================

rge_stats_parts = []
for J in J_GRID:
    pkl = os.path.join(DIR_RGE, f'mvn_J{J}_pps_RGE_stats.pkl')
    if not os.path.exists(pkl):
        print(f"[rge_stats] missing {pkl}; skipping J={J}.")
        continue
    rge_stats_parts.append(pd.read_pickle(pkl))

if rge_stats_parts:
    rge_stats_all = pd.concat(rge_stats_parts, ignore_index=True)
    rge_stats_all['interim_date'] = pd.to_datetime(
        rge_stats_all['interim_date'],
    )
    _stats_interim_order = (
        rge_stats_all[['interim_date', 'interim_month_year']]
        .drop_duplicates()
        .sort_values('interim_date')['interim_month_year'].tolist()
    )
    rge_stats_all['interim_month_year'] = pd.Categorical(
        rge_stats_all['interim_month_year'],
        categories=_stats_interim_order, ordered=True,
    )
    rge_stats_all['J_label'] = pd.Categorical(
        'J=' + rge_stats_all['J'].astype(str),
        categories=[f'J={J}' for J in J_GRID], ordered=True,
    )

    long = rge_stats_all.melt(
        id_vars=['J', 'J_label', 'interim_id', 'interim_date',
                 'interim_month_year', 'j'],
        value_vars=['r2', 'rho'],
        var_name='metric', value_name='value',
    )
    long['metric'] = pd.Categorical(
        long['metric'].map({'r2': 'R^2', 'rho': 'rho'}),
        categories=['R^2', 'rho'], ordered=True,
    )
    rge_box_stats = (
        long.groupby(['J_label', 'metric', 'interim_month_year'],
                     observed=True)['value']
        .agg(
            q025=lambda s: s.quantile(0.025),
            q25 =lambda s: s.quantile(0.25),
            q50 =lambda s: s.quantile(0.50),
            q75 =lambda s: s.quantile(0.75),
            q975=lambda s: s.quantile(0.975),
        )
        .reset_index()
    )
    rge_box_stats['grp'] = (
        rge_box_stats['J_label'].astype(str) + '|'
        + rge_box_stats['metric'].astype(str) + '|'
        + rge_box_stats['interim_month_year'].astype(str)
    )
    rge_box_stats.to_pickle(os.path.join(
        DIR_OUT, 'mvn_compare_methods_rge_stats_box.pkl',
    ))

    p = (
        ggplot(rge_box_stats,
               aes(x='interim_month_year',
                   ymin='q025', lower='q25', middle='q50',
                   upper='q75', ymax='q975', group='grp'))
        + geom_boxplot(stat='identity',
                       fill='#197EC0', alpha=0.6,
                       colour='#404040', size=0.3)
        + facet_grid('metric ~ J_label', scales='free_y')
        + theme_bw()
        + theme(
            axis_text_x=element_text(angle=45, vjust=1, hjust=1),
            legend_position='none',
            figure_size=(5 * len(J_GRID), 8),
            strip_background=element_blank(),
            strip_text=element_text(face='bold'),
        )
        + labs(
            x='Interim date',
            y='regression statistic (across responses j)',
        )
    )
    pdf_path = os.path.join(
        DIR_OUT, 'mvn_compare_methods_rge_stats_box.pdf',
    )
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"Saved RGE rho / R^2 boxplot to {pdf_path}")

print(f"\nCross-method comparison complete. Outputs in: {DIR_OUT}")
