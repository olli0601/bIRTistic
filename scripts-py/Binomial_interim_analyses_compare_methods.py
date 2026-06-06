#!/usr/bin/env python3
"""
Binomial interim analyses: cross-method comparison.

Loads the per-method artifacts produced by

  - HMC nested-MC (refit)        : ``Binomial_interim_analyses_with_nested_MC_hmc.py``
  - SVI nested-MC (refit)        : ``Binomial_interim_analyses_with_nested_MC_svi_autodiagnormal.py``
  - IS (reweight)                : ``Binomial_interim_analyses_with_IS_from_x.py``
  - Regression H1x (Strong-Oakley): ``Binomial_interim_analysis_regression_H1x_on_wz.py``
  - Regression endptx (Strong-Oakley): ``Binomial_interim_analysis_regression_endptx_on_wz.py``

and emits four figures comparing them against the closed-form analytic PPS
(common across all five runs):

  1. ``binomial_compare_methods_p_h1_xz_boxplot.pdf``: per-interim
     p(H_1 | x, z) quantile boxplot dodged by method.
  2. ``binomial_compare_methods_pps.pdf``: PPS bar chart per interim
     (closed-form analytic in black + each method in a Futurama colour).
  3. ``binomial_compare_methods_timing.pdf``: per-interim inference time
     (mins) per method, Futurama colours, sqrt y-axis.
  4. ``binomial_compare_methods_ess.pdf``: ESS / particle + E(w^2) for the
     IS-from-x method (the nested-MC methods refit so have no IS weights).

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/Binomial_interim_analyses_compare_methods.py
"""

# %%

import os
import sys
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
    ggplot, aes, geom_boxplot, geom_col, geom_errorbar, geom_hline, geom_text,
    scale_fill_manual, scale_y_continuous, scale_y_sqrt,
    position_dodge, theme_bw, theme, element_text, labs,
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
dir_out = os.path.join(_sandbox, "py-binomial-interim-compare-methods-260606")
os.makedirs(dir_out, exist_ok=True)

method_order = [
    'nested-MC using HMC for each (x,z)',
    'nested-MC using SVI for each (x,z)',
    'IS reweighting of theta|x',
    'Regression of H1-x on w(z)',
    'Regression of endpt-x on w(z)',
]
method_colours = dict(zip(method_order, _futurama_palette(len(method_order))))

pps_ProbH1_thresh = 0.89

# %%

# =============================================================================
# Load per-method artifacts.
#   - p_h1_xz: long form, one row per (interim, s) with column ``p_h1_xz``.
#   - pps:     scalar PPS per interim per method.
#   - timing:  mins per interim per method.
#   - is_perf: ESS / particle + E(w^2) per (interim, s) for the IS method.
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

p_h1_xz = pd.concat(
    [
        p_h1_xz_hmc.assign(method='nested-MC using HMC for each (x,z)'),
        p_h1_xz_svi.assign(method='nested-MC using SVI for each (x,z)'),
        p_h1_xz_is.assign(method='IS reweighting of theta|x'),
        p_h1_xz_rg.assign(method='Regression of H1-x on w(z)'),
        p_h1_xz_rge.assign(method='Regression of endpt-x on w(z)'),
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
ppsa = pd.read_pickle(
    os.path.join(dir_hmc, 'binomial_interim_pps_closed_form.pkl')
)[['interim_id', 'interim_date', 'interim_month_year', 'pps']]

pps = pd.concat(
    [
        ppsa.assign(method='analytic'),
        pps_hmc.assign(method='nested-MC using HMC for each (x,z)'),
        pps_svi.assign(method='nested-MC using SVI for each (x,z)'),
        pps_is.assign(method='IS reweighting of theta|x'),
        pps_rg.assign(method='Regression of H1-x on w(z)'),
        pps_rge.assign(method='Regression of endpt-x on w(z)'),
    ],
    ignore_index=True,
)

# Timing: HMC / SVI carry mins as `interim_mins` on every row of p_h1_xz;
# IS has a dedicated timing csv.
def _hmc_svi_timing(d, fname, label):
    df = pd.read_pickle(os.path.join(d, fname))
    return (
        df[['interim_id', 'interim_date', 'interim_month_year', 'interim_mins']]
        .drop_duplicates()
        .rename(columns={'interim_mins': 'mins'})
        .assign(method=label)
    )


timing_hmc = _hmc_svi_timing(dir_hmc, 'binomial_interim_pps_p_h1_xz.pkl',
                             'nested-MC using HMC for each (x,z)')
timing_svi = _hmc_svi_timing(dir_svi, 'binomial_interim_pps_p_h1_xz.pkl',
                             'nested-MC using SVI for each (x,z)')

def _csv_timing(d, fname, p_h1_xz_method, label):
    """Regression / IS scripts save mins per interim in a small csv keyed by
    interim_id; merge in the (interim_date, interim_month_year) metadata from
    that method's p_h1_xz frame."""
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


timing_is = _csv_timing(dir_is, 'pps_timing.csv',
                        p_h1_xz_is, 'IS reweighting of theta|x')
timing_rg = _csv_timing(dir_rg, 'binomial_interim_pps_RG_timing.csv',
                        p_h1_xz_rg, 'Regression of H1-x on w(z)')
timing_rge = _csv_timing(dir_rge, 'binomial_interim_pps_RGE_timing.csv',
                         p_h1_xz_rge, 'Regression of endpt-x on w(z)')

timing = pd.concat(
    [timing_hmc, timing_svi, timing_is, timing_rg, timing_rge],
    ignore_index=True,
)

# Reuse _meta_is for IS-perf merge below.
_meta_is = (
    p_h1_xz_is[['interim_id', 'interim_date', 'interim_month_year']]
    .drop_duplicates()
)

# IS perf: ESS + E(w^2) per (interim, s); attach interim metadata for facets.
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

p_h1_xz['method'] = pd.Categorical(
    p_h1_xz['method'], categories=method_order, ordered=True,
)
pps['method'] = pd.Categorical(
    pps['method'], categories=['analytic', *method_order], ordered=True,
)
timing['method'] = pd.Categorical(
    timing['method'], categories=method_order, ordered=True,
)

# %%

# =============================================================================
# Plot 1: per-interim p(H_1 | x, z) distribution, methods dodged with
# Futurama fill colours. Quantiles q025/q25/q50/q75/q975 per (interim, method).
# =============================================================================

box_stats = (
    p_h1_xz
    .groupby(['interim_month_year', 'method'], observed=True)['p_h1_xz']
    .agg(
        q025=lambda s: s.quantile(0.025),
        q25=lambda s: s.quantile(0.25),
        q50=lambda s: s.quantile(0.50),
        q75=lambda s: s.quantile(0.75),
        q975=lambda s: s.quantile(0.975),
    )
    .reset_index()
)
box_stats['grp'] = (
    box_stats['interim_month_year'].astype(str)
    + '|' + box_stats['method'].astype(str)
)

p = (
    ggplot(box_stats, aes(x='interim_month_year', fill='method'))
    + geom_boxplot(
        aes(ymin='q025', lower='q25', middle='q50', upper='q75', ymax='q975',
            group='grp'),
        stat='identity', position=position_dodge(width=0.8), width=0.7,
        colour='#808080', size=0.3,
    )
    + geom_hline(yintercept=pps_ProbH1_thresh, colour='black', size=1.0)
    + scale_fill_manual(values=method_colours)
    + scale_y_continuous(
        limits=[0, 1],
        breaks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        labels=['0%', '20%', '40%', '60%', '80%', '100%'],
    )
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        legend_position='top',
        figure_size=(14, 6),
    )
    + labs(x='Interim', y='p(H_1 | x, z) for predicted z samples', fill='method')
)
pdf_path = os.path.join(
    dir_out, 'binomial_compare_methods_p_h1_xz_boxplot.pdf',
)
p.save(pdf_path, verbose=False, limitsize=False)
print(f"Saved p(H_1 | x, z) boxplot to: {pdf_path}")

# %%

# =============================================================================
# Plot 2: PPS per interim. Black = closed-form analytic; Futurama colours per
# method. Dodged side-by-side. Vertical bars = 95% bootstrap interval on the
# Monte-Carlo PPS estimate (B resamples with replacement of the S per-interim
# p(H_1 | x, z) samples; analytic has no MC noise so no errorbar).
# =============================================================================

_pps_colours = {'analytic': '#000000', **method_colours}

bs_B = 2000
bs_parts = []
rng_bs = np.random.default_rng(123)
for m_label, g in p_h1_xz.groupby('method', observed=True):
    wide = (
        g.pivot_table(index=['interim_id', 'interim_date'],
                      columns='s', values='p_h1_xz')
        .sort_index()
    )
    P = wide.to_numpy()
    n_interim, S = P.shape
    idx = rng_bs.integers(0, S, size=(n_interim, bs_B, S))
    P_bs = P[np.arange(n_interim)[:, None, None], idx]
    pps_bs = (P_bs > pps_ProbH1_thresh).mean(axis=2)
    q = np.quantile(pps_bs, [0.025, 0.975], axis=1).T
    part = pd.DataFrame(q, columns=['q025_bs', 'q975_bs'])
    part['interim_id'] = wide.index.get_level_values('interim_id').to_numpy()
    part['method'] = m_label
    bs_parts.append(part)
bs = pd.concat(bs_parts, ignore_index=True)
pps = pps.merge(bs, on=['interim_id', 'method'], how='left')

p = (
    ggplot(pps, aes(x='interim_month_year', y='pps', fill='method'))
    + geom_col(position=position_dodge(width=0.8), width=0.7)
    + geom_errorbar(
        data=pps,
        mapping=aes(x='interim_month_year', ymin='q025_bs', ymax='q975_bs',
                    group='method'),
        position=position_dodge(width=0.8), width=0.3,
        colour='black', size=0.4, inherit_aes=False,
    )
    + geom_hline(yintercept=pps_ProbH1_thresh, colour='black', size=1.0,
                 linetype='dashed')
    + scale_fill_manual(values=_pps_colours)
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        legend_position='top',
        figure_size=(14, 6),
    )
    + labs(
        x='Interim', y='PPS = int P(p(H_1 | x, z) > eta) dz', fill='method',
    )
)
pdf_path = os.path.join(dir_out, 'binomial_compare_methods_pps.pdf')
p.save(pdf_path, verbose=False, limitsize=False)
print(f"Saved PPS comparison plot to: {pdf_path}")

# %%

# =============================================================================
# Plot 3: per-interim inference time (mins), Futurama colours per method.
# sqrt y-axis since the IS path is ~1000x faster than HMC nested-MC.
# =============================================================================

timing['value_label'] = timing['mins'].map(lambda v: f"{v:.2f} mins")
p = (
    ggplot(timing, aes(x='interim_month_year', y='mins', fill='method'))
    + geom_col(position=position_dodge(width=0.8), width=0.7)
    + geom_text(
        aes(label='value_label'),
        position=position_dodge(width=0.8),
        size=6, angle=30, ha='left', va='bottom',
    )
    + scale_fill_manual(values=method_colours)
    + scale_y_sqrt(expand=(0, 0, 0.15, 0))
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        legend_position='top',
        figure_size=(14, 6),
    )
    + labs(x='Interim', y='time (mins)', fill='method')
)
pdf_path = os.path.join(dir_out, 'binomial_compare_methods_timing.pdf')
p.save(pdf_path, verbose=False, limitsize=False)
print(f"Saved timing comparison plot to: {pdf_path}")

# %%

# =============================================================================
# Plot 4: ESS / particle per interim for the IS-reweight method
# (nested-MC methods refit -> no IS weights, no ESS).
# =============================================================================

is_ess = (
    is_perf
    .groupby(['interim_month_year'], observed=True)
    .agg(value=('ess_over_n', 'median'))
    .reset_index()
    .assign(method='IS reweighting of theta|x')
)
is_ess['method'] = pd.Categorical(
    is_ess['method'], categories=method_order, ordered=True,
)

p = (
    ggplot(is_ess, aes(x='interim_month_year', y='value', fill='method'))
    + geom_col(width=0.7)
    + scale_fill_manual(values=method_colours)
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        legend_position='top',
        figure_size=(12, 5),
    )
    + labs(
        x='Interim', y='ESS / particle', fill='method',        
    )
)
pdf_path = os.path.join(dir_out, 'binomial_compare_methods_ess.pdf')
p.save(pdf_path, verbose=False, limitsize=False)
print(f"Saved IS ESS / E(w^2) plot to: {pdf_path}")

# %%

print("\n✓ Cross-method comparison complete.")
print(f"Outputs in: {dir_out}")
