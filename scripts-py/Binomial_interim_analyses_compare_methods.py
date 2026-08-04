#!/usr/bin/env python3
"""
Binomial interim analyses: cross-method comparison.

Loads the per-method artifacts produced by

  - HMC nested-MC (refit)        : ``Binomial_interim_analyses_with_nested_MC_hmc.py``
  - SVI nested-MC (refit)        : ``Binomial_interim_analyses_with_nested_MC_svi_autodiagnormal.py``
  - IS (reweight)                : ``Binomial_interim_analyses_with_IS_from_x.py``
  - Regression H1x (Strong-Oakley): ``Binomial_interim_analysis_regression_H1x_on_wz.py``
  - Regression endptx (Strong-Oakley): ``Binomial_interim_analysis_regression_endptx_on_wz_with_gauss_approx.py``
  - Regression endptx quantile (Strong-Oakley): ``Binomial_interim_analysis_regression_endptx_on_wz_with_quantile_regr.py``
  - Regression endptx mquantile (Strong-Oakley): ``Binomial_interim_analysis_regression_endptx_on_wz_with_mquantile_regr.py``
  - Regression endptx amortised (features-fixed): ``Binomial_interim_analysis_amortise_endptx_on_wz_with_features_fixed_qpsi_MLP_loss_multiquantilehead.py``
  - Regression endptx amortised (features-MLP)  : ``Binomial_interim_analysis_amortise_endptx_on_wz_with_features_MLP_qpsi_MLP_loss_multiquantilehead.py``

and emits two families of figures compared against the closed-form analytic
PPS (common across every run).

Regression-family plots (amortised methods removed to reduce clutter):
  1. ``binomial_compare_methods_p_h1_xz_boxplot.pdf``
  2. ``binomial_compare_methods_pps.pdf``
  3. ``binomial_compare_methods_timing.pdf``  (train + test times; for the
     non-amortised family train time is 0, so this is per-interim inference)

Amortised-focus plots (best regression + amortised variants):
  4. ``binomial_compare_methods_p_h1_xz_boxplot_amortised.pdf``
  5. ``binomial_compare_methods_pps_amortised.pdf``
  6. ``binomial_compare_methods_timing_amortised.pdf``  (train + deploy times,
     stacked bar per method with one fill colour per train / deploy segment,
     one bar per method summing per-interim deploy times across the schedule)

Plus the IS diagnostic:
  7. ``binomial_compare_methods_ess.pdf``

Legend layouts single-column; method → colour mapping is consistent across
every plot in the file.

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/Binomial_interim_analyses_compare_methods.py
"""

# %%

import os
import sys
import pickle
from pathlib import Path

try:
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
except NameError:
    project_root = Path.cwd()
    if project_root.name == 'scripts-py':
        project_root = project_root.parent

python_path = str(project_root / 'python')
if python_path not in sys.path:
    sys.path.insert(0, python_path)

import numpy as np
import pandas as pd
from plotnine import (
    ggplot, aes, coord_flip, geom_boxplot, geom_col, geom_errorbar, geom_hline,
    geom_text, scale_colour_manual, scale_fill_manual, scale_y_continuous,
    scale_y_sqrt, position_dodge, position_stack,
    theme_bw, theme, element_text, labs,
    guides, guide_legend,
)
import warnings
warnings.filterwarnings('ignore')

from utils import _futurama_palette

print("✓ Imports successful")

# %%

# =============================================================================
# Configuration
# =============================================================================

_sandbox = "/Users/or105/sandbox/bIRTistic"
dir_hmc = os.path.join(_sandbox, "py-binomial-interim-with-nested-MC-hmc-260605")
dir_svi = os.path.join(_sandbox, "py-binomial-interim-with-nested-MC-svi-autodiagnormal-260604")
dir_is = os.path.join(_sandbox, "py-binomial-interim-with-IS-260605")
dir_rg = os.path.join(_sandbox, "py-binomial-interim-with-regression-on-H1x-wz-260606")
dir_rge = os.path.join(_sandbox, "py-binomial-interim-with-regression-on-endptx-wz-260606")
dir_rgeq = os.path.join(_sandbox, "py-binomial-interim-with-regression-on-endptx-wz-quantile-regr-260702")
dir_rgem = os.path.join(_sandbox, "py-binomial-interim-with-regression-on-endptx-wz-mquantile-regr-260702")
dir_rgea = os.path.join(_sandbox, "py-binomial-interim-amortise-endptx-on-wz-with-features-fixed-qpsi-MLP-loss-multiquantilehead-260714")
dir_rgeb = os.path.join(_sandbox, "py-binomial-interim-amortise-endptx-on-wz-with-features-MLP-qpsi-MLP-loss-multiquantilehead-260711")
dir_out = os.path.join(_sandbox, "py-binomial-interim-compare-methods-260606")
os.makedirs(dir_out, exist_ok=True)

# Method labels. Amortised labels include the filename variant marker
# (``features-fixed`` / ``features-MLP``) so the legend is unambiguous when
# the family plot vs the amortised-focus plot are read side-by-side.
LBL_HMC     = 'nested-MC using HMC for each (x,z)'
LBL_SVI     = 'nested-MC using SVI for each (x,z)'
LBL_IS      = 'IS reweighting of theta|x'
LBL_RG      = 'Regression of H1-x on w(z)'
LBL_RGE     = 'Regression of endpt-x on w(z) using Gaussian approx'
LBL_RGEQ    = 'Regression of endpt-x on w(z) - quantile'
LBL_RGEM    = 'Regression of endpt-x on w(z) - mquantile'
LBL_RGEA    = 'Amortiser-features-fixed-qpsi-MLP-loss-multiquantilehead'
LBL_RGEB    = 'Amortiser-features-MLP-qpsi-MLP-loss-multiquantilehead'

# Regression-family method order (amortised excluded from PPS + boxplot;
# both amortised variants ARE added back to the shared timing plot per
# user request, since train + test times are comparable there).
non_amort_methods = [
    LBL_HMC, LBL_SVI, LBL_IS, LBL_RG,
    LBL_RGE, LBL_RGEQ, LBL_RGEM,
]
# Amortised-focus method order: best regression baseline (Gaussian approx)
# + the two amortised variants.
amort_focus_methods = [LBL_RGE, LBL_RGEA, LBL_RGEB]
# Timing plot includes non-amortised methods PLUS the two amortised variants
# (train + test comparable across all methods).
timing_methods = non_amort_methods + [LBL_RGEA, LBL_RGEB]

# Global palette: assign one colour per method up-front so the same method
# reuses the same colour on every plot in the script.
_all_methods = non_amort_methods + [LBL_RGEA, LBL_RGEB]
method_colours = dict(zip(_all_methods, _futurama_palette(len(_all_methods))))
# Analytic drawn as white fill + black outline so it reads as
# ground-truth reference, not a competing estimate.
pps_colours = {'analytic': '#FFFFFF', **method_colours}
# Outline colour per method: only analytic gets a visible black outline;
# every other method uses matching-fill colour so the border blends in.
pps_outline_colours = {'analytic': '#000000',
                       **{m: c for m, c in method_colours.items()}}

# Fill colours for the stacked train / deploy timing bar.
_TRAIN_DEPLOY_COLOURS = {
    'train (one-off)':      '#1f77b4',
    'deploy (all interims)': '#ff7f0e',
}

pps_ProbH1_target_lwr_quantile = 0.89

# %%

# =============================================================================
# Load per-method artifacts.
# =============================================================================


def _load_p_h1_xz(d, fname):
    df = pd.read_pickle(os.path.join(d, fname))
    return df[['interim_id', 'interim_date', 'interim_month_year', 's', 'p_h1_xz']]


def _load_pps(d, fname):
    df = pd.read_pickle(os.path.join(d, fname))
    return df[['interim_id', 'interim_date', 'interim_month_year', 'pps']]


p_h1_xz_hmc = _load_p_h1_xz(dir_hmc, 'binomial_interim_pps_p_h1_xz.pkl')
p_h1_xz_svi = _load_p_h1_xz(dir_svi, 'binomial_interim_pps_p_h1_xz.pkl')
p_h1_xz_is = _load_p_h1_xz(dir_is, 'pps_p_h1_xz_IS.pkl')
p_h1_xz_rg = _load_p_h1_xz(dir_rg, 'binomial_interim_pps_RG_p_h1_xz.pkl')
p_h1_xz_rge = _load_p_h1_xz(dir_rge, 'binomial_interim_pps_RGE_p_h1_xz.pkl')
p_h1_xz_rgeq = _load_p_h1_xz(dir_rgeq, 'binomial_interim_pps_RGEQ_p_h1_xz.pkl')
p_h1_xz_rgem = _load_p_h1_xz(dir_rgem, 'binomial_interim_pps_RGEM_p_h1_xz.pkl')
p_h1_xz_rgea = _load_p_h1_xz(dir_rgea, 'binomial_interim_pps_RGEA_p_h1_xz.pkl')
p_h1_xz_rgeb = _load_p_h1_xz(dir_rgeb, 'binomial_interim_pps_RGEB_p_h1_xz.pkl')

p_h1_xz = pd.concat(
    [
        p_h1_xz_hmc.assign(method=LBL_HMC),
        p_h1_xz_svi.assign(method=LBL_SVI),
        p_h1_xz_is.assign(method=LBL_IS),
        p_h1_xz_rg.assign(method=LBL_RG),
        p_h1_xz_rge.assign(method=LBL_RGE),
        p_h1_xz_rgeq.assign(method=LBL_RGEQ),
        p_h1_xz_rgem.assign(method=LBL_RGEM),
        p_h1_xz_rgea.assign(method=LBL_RGEA),
        p_h1_xz_rgeb.assign(method=LBL_RGEB),
    ],
    ignore_index=True,
)

pps_hmc = _load_pps(dir_hmc, 'binomial_interim_pps_nested_mc.pkl')
pps_svi = _load_pps(dir_svi, 'binomial_interim_pps_nested_mc.pkl')
pps_is = _load_pps(dir_is, 'pps_IS.pkl')
# Regression PPS tables are saved as CSV; load + project.
pps_rg = pd.read_csv(os.path.join(dir_rg, 'binomial_interim_pps_RG.csv'))[
    ['interim_id', 'interim_date', 'interim_month_year', 'pps']
]
pps_rge = pd.read_csv(os.path.join(dir_rge, 'binomial_interim_pps_RGE.csv'))[
    ['interim_id', 'interim_date', 'interim_month_year', 'pps']
]
pps_rgeq = pd.read_csv(os.path.join(dir_rgeq, 'binomial_interim_pps_RGEQ.csv'))[
    ['interim_id', 'interim_date', 'interim_month_year', 'pps']
]
pps_rgem = pd.read_csv(os.path.join(dir_rgem, 'binomial_interim_pps_RGEM.csv'))[
    ['interim_id', 'interim_date', 'interim_month_year', 'pps']
]
pps_rgea = pd.read_csv(os.path.join(dir_rgea, 'binomial_interim_pps_RGEA.csv'))[
    ['interim_id', 'interim_date', 'interim_month_year', 'pps']
]
pps_rgeb = pd.read_csv(os.path.join(dir_rgeb, 'binomial_interim_pps_RGEB.csv'))[
    ['interim_id', 'interim_date', 'interim_month_year', 'pps']
]
ppsa = pd.read_pickle(
    os.path.join(dir_hmc, 'binomial_interim_pps_closed_form.pkl')
)[['interim_id', 'interim_date', 'interim_month_year', 'pps']]

pps = pd.concat(
    [
        ppsa.assign(method='analytic'),
        pps_hmc.assign(method=LBL_HMC),
        pps_svi.assign(method=LBL_SVI),
        pps_is.assign(method=LBL_IS),
        pps_rg.assign(method=LBL_RG),
        pps_rge.assign(method=LBL_RGE),
        pps_rgeq.assign(method=LBL_RGEQ),
        pps_rgem.assign(method=LBL_RGEM),
        pps_rgea.assign(method=LBL_RGEA),
        pps_rgeb.assign(method=LBL_RGEB),
    ],
    ignore_index=True,
)


# Timing: HMC / SVI carry mins as `interim_mins` on every row of p_h1_xz;
# IS / regression / amortised have dedicated timing csvs.
def _hmc_svi_timing(d, fname, label):
    df = pd.read_pickle(os.path.join(d, fname))
    return (
        df[['interim_id', 'interim_date', 'interim_month_year', 'interim_mins']]
        .drop_duplicates()
        .rename(columns={'interim_mins': 'mins'})
        .assign(method=label)
    )


timing_hmc = _hmc_svi_timing(dir_hmc, 'binomial_interim_pps_p_h1_xz.pkl', LBL_HMC)
timing_svi = _hmc_svi_timing(dir_svi, 'binomial_interim_pps_p_h1_xz.pkl', LBL_SVI)


def _csv_timing(d, fname, p_h1_xz_method, label):
    meta = (
        p_h1_xz_method[['interim_id', 'interim_date', 'interim_month_year']]
        .drop_duplicates()
    )
    csv = pd.read_csv(os.path.join(d, fname))
    return (
        meta
        .merge(csv[['interim_id', 'mins_interim_id']], on='interim_id')
        .rename(columns={'mins_interim_id': 'mins'})
        .assign(method=label)
    )


timing_is = _csv_timing(dir_is, 'pps_timing.csv', p_h1_xz_is, LBL_IS)
timing_rg = _csv_timing(dir_rg, 'binomial_interim_pps_RG_timing.csv', p_h1_xz_rg, LBL_RG)
timing_rge = _csv_timing(dir_rge, 'binomial_interim_pps_RGE_timing.csv', p_h1_xz_rge, LBL_RGE)
timing_rgeq = _csv_timing(dir_rgeq, 'binomial_interim_pps_RGEQ_timing.csv', p_h1_xz_rgeq, LBL_RGEQ)
timing_rgem = _csv_timing(dir_rgem, 'binomial_interim_pps_RGEM_timing.csv', p_h1_xz_rgem, LBL_RGEM)
timing_rgea = _csv_timing(dir_rgea, 'binomial_interim_pps_RGEA_timing.csv', p_h1_xz_rgea, LBL_RGEA)
timing_rgeb = _csv_timing(dir_rgeb, 'binomial_interim_pps_RGEB_timing.csv', p_h1_xz_rgeb, LBL_RGEB)

timing = pd.concat(
    [timing_hmc, timing_svi, timing_is, timing_rg, timing_rge, timing_rgeq,
     timing_rgem, timing_rgea, timing_rgeb],
    ignore_index=True,
)


def _load_training_mins(ckpt_path):
    if not os.path.exists(ckpt_path):
        return float('nan')
    with open(ckpt_path, 'rb') as f:
        payload = pickle.load(f)
    return float(payload.get('training_mins', float('nan')))


training_mins_by_method = {
    LBL_RGEA: _load_training_mins(
        os.path.join(dir_rgea, 'binomial_interim_amortised_pps_net.pkl')
    ),
    LBL_RGEB: _load_training_mins(
        os.path.join(dir_rgeb, 'binomial_interim_amortised_pps_net.pkl')
    ),
}
# Non-amortised methods have no training step (train = 0).
for lbl in non_amort_methods:
    training_mins_by_method.setdefault(lbl, 0.0)
print(f"Training mins per method: "
      f"{ {k: round(v, 3) for k, v in training_mins_by_method.items()} }")

# IS perf: ESS + E(w^2) per (interim, s); attach interim metadata for facets.
_meta_is = (
    p_h1_xz_is[['interim_id', 'interim_date', 'interim_month_year']]
    .drop_duplicates()
)
is_perf = pd.read_csv(os.path.join(dir_is, 'pps_is_perf.csv'))
is_perf = is_perf.merge(_meta_is, on='interim_id')

# %%

# =============================================================================
# Establish a single interim_month_year ordering shared across all panels.
# =============================================================================

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


# %%

# =============================================================================
# Helper: dodged boxplot of p(H_1 | x, z) quantiles per (method, interim).
# =============================================================================


def _boxplot_p_h1_xz(df, methods, pdf_path, figure_size=(14, 8)):
    """Emit a dodged boxplot of ``p_h1_xz`` per (method, interim).

    Quantile method is dropped: its p_h1_xz is a hard 0/1 label, not a
    probability, so the boxplot would mix semantically-different quantities.
    """
    keep = [m for m in methods if m != LBL_RGEQ]
    d = df[df['method'].isin(keep)].copy()
    d['method'] = pd.Categorical(d['method'], categories=keep, ordered=True)
    stats = (
        d.groupby(['interim_month_year', 'method'], observed=True)['p_h1_xz']
        .agg(
            q025=lambda s: s.quantile(0.025),
            q25=lambda s: s.quantile(0.25),
            q50=lambda s: s.quantile(0.50),
            q75=lambda s: s.quantile(0.75),
            q975=lambda s: s.quantile(0.975),
        )
        .reset_index()
    )
    stats['grp'] = (
        stats['interim_month_year'].astype(str) + '|'
        + stats['method'].astype(str)
    )
    p = (
        ggplot(stats, aes(x='interim_month_year', fill='method'))
        + geom_boxplot(
            aes(ymin='q025', lower='q25', middle='q50', upper='q75', ymax='q975',
                group='grp'),
            stat='identity', position=position_dodge(width=0.8), width=0.7,
            colour='#808080', size=0.3,
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
        + guides(fill=guide_legend(ncol=1))
        + theme_bw()
        + theme(
            axis_text_x=element_text(angle=45, vjust=1, hjust=1),
            legend_position='bottom',
            legend_direction='vertical',
            figure_size=figure_size,
        )
        + labs(x='Interim', y='p(H_1 | x, z) for predicted z samples',
               fill='method')
    )
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"Saved p(H_1 | x, z) boxplot to: {pdf_path}")


# =============================================================================
# Helper: dodged PPS bars with 95% bootstrap CI per (method, interim).
# =============================================================================


def _bootstrap_pps_ci(p_h1_xz_df, seed=123, B=2000):
    """Return per-(interim_id, method) 95% bootstrap CI on the MC PPS
    estimate (analytic method has no MC noise -> skipped)."""
    parts = []
    rng = np.random.default_rng(seed)
    for m_label, g in p_h1_xz_df.groupby('method', observed=True):
        wide = (
            g.pivot_table(index=['interim_id', 'interim_date'],
                          columns='s', values='p_h1_xz')
            .sort_index()
        )
        P = wide.to_numpy()
        n_interim, S = P.shape
        idx = rng.integers(0, S, size=(n_interim, B, S))
        P_bs = P[np.arange(n_interim)[:, None, None], idx]
        pps_bs = (P_bs > pps_ProbH1_target_lwr_quantile).mean(axis=2)
        q = np.quantile(pps_bs, [0.025, 0.975], axis=1).T
        part = pd.DataFrame(q, columns=['q025_bs', 'q975_bs'])
        part['interim_id'] = wide.index.get_level_values('interim_id').to_numpy()
        part['method'] = m_label
        parts.append(part)
    return pd.concat(parts, ignore_index=True)


def _pps_bars(pps_df, methods, pdf_path, figure_size=(14, 8)):
    """Emit a per-(method, interim) PPS bar chart with bootstrap CI.

    Analytic bar is drawn WHITE with a thin BLACK outline (size = 1.1)
    so the ground-truth reference is visually distinct from the
    Monte-Carlo estimates (which use fill = method colour,
    matching-fill outline). Analytic is forced first in the dodge
    order so it always sits at the LEFT of each interim's group. The
    legend swatch for `analytic` picks up the same black outline via
    ``override_aes`` on the fill guide.
    """
    keep_with_analytic = ['analytic', *methods]
    d = pps_df[pps_df['method'].isin(keep_with_analytic)].copy()
    ci = _bootstrap_pps_ci(
        p_h1_xz[p_h1_xz['method'].isin(methods)].copy(),
    )
    # Merge BEFORE forcing Categorical: merge on object dtype avoids the
    # pandas quirk that drops the level ordering on the join key.
    d = d.merge(ci, on=['interim_id', 'method'], how='left')
    d['method'] = pd.Categorical(
        d['method'].astype(str),
        categories=keep_with_analytic, ordered=True,
    )
    # Sort so ``analytic`` is drawn first (and lands leftmost inside each
    # interim's dodge group).
    d = d.sort_values(['interim_month_year', 'method']).reset_index(drop=True)
    # Per-legend override: only 'analytic' gets a black outline; every
    # other row gets its own fill colour again so the outline blends in.
    _legend_outline = [pps_outline_colours[m] for m in keep_with_analytic]
    p = (
        ggplot(d, aes(x='interim_month_year', y='pps', fill='method',
                      colour='method', group='method'))
        + geom_col(position=position_dodge(width=0.8, preserve='single'),
                   width=0.7, size=1.1)
        + geom_errorbar(
            data=d,
            mapping=aes(x='interim_month_year', ymin='q025_bs',
                        ymax='q975_bs', group='method'),
            position=position_dodge(width=0.8, preserve='single'), width=0.3,
            colour='black', size=0.4, inherit_aes=False,
        )
        + geom_hline(yintercept=pps_ProbH1_target_lwr_quantile,
                     colour='black', size=1.0, linetype='dashed')
        + scale_fill_manual(values=pps_colours,
                            breaks=keep_with_analytic,
                            limits=keep_with_analytic)
        + scale_colour_manual(values=pps_outline_colours,
                              breaks=keep_with_analytic,
                              limits=keep_with_analytic,
                              guide=None)
        + guides(fill=guide_legend(
            ncol=1,
            override_aes={'colour': _legend_outline, 'size': 1.1},
        ))
        + theme_bw()
        + theme(
            axis_text_x=element_text(angle=45, vjust=1, hjust=1),
            legend_position='bottom',
            legend_direction='vertical',
            figure_size=figure_size,
        )
        + labs(x='Interim',
               y='PPS = int P(p(H_1 | x, z) > eta) dz', fill='method')
    )
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"Saved PPS comparison plot to: {pdf_path}")


# =============================================================================
# Helper: per-interim inference time bar chart with sqrt y-axis. When
# ``add_training_time`` is True the training minutes are added to each
# per-interim bar (non-amortised methods carry train = 0 so the visualisation
# is unchanged for them).
# =============================================================================


def _timing_bars(timing_df, methods, pdf_path, figure_size=(15, 8)):
    """Two-panel-like layout on ONE x-axis:

      - Leftmost slot ``training``: one-off amortiser training time.
        Non-amortised methods have zero here (bar suppressed).
      - Interim slots (2024-Jan .. 2024-Nov): per-interim DEPLOY time
        only (no training added).

    Bar widths stay consistent across the whole x-axis because every
    (method, slot) pair carries an explicit row (zeros where absent).
    """
    d = timing_df[timing_df['method'].isin(methods)][
        ['interim_month_year', 'method', 'mins']
    ].copy()
    d['slot'] = d['interim_month_year'].astype(str)
    d = d[['slot', 'method', 'mins']]

    # Prepend a ``training`` slot per method with training_mins_by_method.
    train_rows = pd.DataFrame({
        'slot':   ['training'] * len(methods),
        'method': methods,
        'mins':   [training_mins_by_method.get(m, 0.0) for m in methods],
    })

    slot_order = ['training'] + interim_order
    long = pd.concat([train_rows, d], ignore_index=True)

    # Fill zeros for any missing (method, slot) so dodge width is uniform.
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
    print(f"Saved timing comparison plot to: {pdf_path}")


# =============================================================================
# Helper: total-wall-clock stacked bar per method (train segment + total
# deploy segment) for the amortised-focus plot.
# =============================================================================


def _train_deploy_stacked_bars(timing_df, methods, pdf_path,
                               figure_size=(11, 4)):
    """One horizontal bar per method summing per-interim deploy times
    across the full schedule, stacked on the one-off training time. Two
    fill colours per bar (train vs deploy). coord_flip so long method
    labels stay readable."""
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
    # coord_flip reverses factor order along the flipped x-axis; feed
    # reversed levels so the top bar is the first method listed.
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
    print(f"Saved train + deploy stacked timing plot to: {pdf_path}")


# %%

# =============================================================================
# Regression-family (amortised excluded) plots.
# =============================================================================

_boxplot_p_h1_xz(
    p_h1_xz, non_amort_methods,
    os.path.join(dir_out, 'binomial_compare_methods_p_h1_xz_boxplot.pdf'),
)
_pps_bars(
    pps, non_amort_methods,
    os.path.join(dir_out, 'binomial_compare_methods_pps.pdf'),
)
_timing_bars(
    timing, timing_methods,
    os.path.join(dir_out, 'binomial_compare_methods_timing.pdf'),
)

# %%

# =============================================================================
# Amortised-focus plots: best regression baseline + amortised variants.
# =============================================================================

_boxplot_p_h1_xz(
    p_h1_xz, amort_focus_methods,
    os.path.join(
        dir_out, 'binomial_compare_methods_p_h1_xz_boxplot_amortised.pdf',
    ),
)
_pps_bars(
    pps, amort_focus_methods,
    os.path.join(dir_out, 'binomial_compare_methods_pps_amortised.pdf'),
)
_train_deploy_stacked_bars(
    timing, amort_focus_methods,
    os.path.join(dir_out, 'binomial_compare_methods_timing_amortised.pdf'),
)

# %%

# =============================================================================
# IS ESS / particle plot (unchanged).
# =============================================================================

is_ess = (
    is_perf
    .groupby(['interim_month_year'], observed=True)
    .agg(value=('ess_over_n', 'median'))
    .reset_index()
    .assign(method=LBL_IS)
)
is_ess['method'] = pd.Categorical(
    is_ess['method'], categories=non_amort_methods, ordered=True,
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
    )
    + labs(x='Interim', y='ESS / particle', fill='method')
)
pdf_path = os.path.join(dir_out, 'binomial_compare_methods_ess.pdf')
p.save(pdf_path, verbose=False, limitsize=False)
print(f"Saved IS ESS / particle plot to: {pdf_path}")

# %%

print("\n✓ Cross-method comparison complete.")
print(f"Outputs in: {dir_out}")
