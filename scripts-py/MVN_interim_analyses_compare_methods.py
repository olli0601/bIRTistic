#!/usr/bin/env python3
"""
MVN interim analyses: closed-form vs nested-MC HMC, per-J cross-comparison.

Loads the per-J artifacts produced by
``MVN_interim_analyses_with_nested_MC_hmc.py`` (one
``py-mvn-interim-with-nested-MC-hmc-J{J}-260609`` directory per
``J in {5, 20, 50, 100}``) and emits:

  1. ``mvn_compare_methods_pps.pdf``: per-component PPS per interim,
     closed-form analytic (black) vs nested-MC HMC (Futurama blue),
     faceted by (J, item j). Vertical bars on the HMC bar are 95%
     bootstrap intervals over the S Monte-Carlo p(H_1 | x, z) draws.
  2. ``mvn_compare_methods_p_h1_xz_boxplot.pdf``: per-interim per-J
     boxplot of the nested-MC p(H_1 | x, z) draws, showing a small set of
     components from each J (j = 0, J/4, J/2, 3J/4, J-1).
  3. ``mvn_compare_methods_timing.pdf``: per-interim nested-MC HMC wall
     time per J, sqrt y-axis.
  4. ``mvn_compare_methods_pps_scatter.pdf``: per-interim scatter of
     closed-form PPS vs nested-MC HMC PPS, faceted by J. Identity line
     dashed; close to identity = the HMC nested-MC estimator is unbiased.

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/MVN_interim_analyses_compare_methods.py
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
    ggplot, aes, geom_abline, geom_boxplot, geom_col, geom_errorbar,
    geom_hline, geom_point, facet_grid, facet_wrap,
    scale_fill_manual, scale_colour_manual, scale_y_sqrt,
    position_dodge, theme_bw, theme, element_text, labs,
)
import warnings
warnings.filterwarnings('ignore')

from utils import _futurama_palette

print("Imports successful")

# %%

# =============================================================================
# Configuration
# =============================================================================

SANDBOX = "/Users/or105/sandbox/bIRTistic"
J_GRID = [5, 20, 50, 100]
METHOD_HMC = 'nested-MC HMC'
METHOD_ANA = 'analytic'
PPS_ProbH1_THRESH = 0.89

dir_out = os.path.join(SANDBOX, "py-mvn-interim-compare-methods-260609")
os.makedirs(dir_out, exist_ok=True)


def _dir_for(J: int) -> str:
    return os.path.join(SANDBOX, f"py-mvn-interim-with-nested-MC-hmc-J{J}-260609")


# %%

# =============================================================================
# Load per-J artifacts.
# =============================================================================


def _load_J(J: int):
    d = _dir_for(J)
    if not os.path.isdir(d):
        raise FileNotFoundError(
            f"Missing per-J artifact dir for J={J}: {d}\n"
            f"Run MVN_interim_analyses_with_nested_MC_hmc.py first."
        )
    ppsa = pd.read_pickle(os.path.join(d, 'mvn_interim_pps_closed_form.pkl'))
    pps_nm = pd.read_pickle(os.path.join(d, 'mvn_interim_pps_nested_mc.pkl'))
    dp_h1_xz = pd.read_pickle(os.path.join(d, 'mvn_interim_pps_p_h1_xz.pkl'))
    # j is on ppsa already; recover j on pps_nm + dp_h1_xz from item_label = 'mu_<j>'.
    pps_nm = pps_nm.assign(
        j=pps_nm['item_label'].str.replace('mu_', '').astype(int),
    )
    dp_h1_xz = dp_h1_xz.assign(
        j=dp_h1_xz['item_label'].str.replace('mu_', '').astype(int),
    )
    for df in (ppsa, pps_nm, dp_h1_xz):
        df['J'] = J
    return ppsa, pps_nm, dp_h1_xz


ppsa_parts = []
pps_nm_parts = []
dp_h1_xz_parts = []
for J in J_GRID:
    a, b, c = _load_J(J)
    ppsa_parts.append(a)
    pps_nm_parts.append(b)
    dp_h1_xz_parts.append(c)

ppsa = pd.concat(ppsa_parts, ignore_index=True)
pps_nm = pd.concat(pps_nm_parts, ignore_index=True)
dp_h1_xz = pd.concat(dp_h1_xz_parts, ignore_index=True)

print(f"Loaded: ppsa rows = {len(ppsa)}; pps_nm rows = {len(pps_nm)};"
      f" dp_h1_xz rows = {len(dp_h1_xz)}")

# Common interim ordering across all J (interims share the same grid).
interim_order = (
    pps_nm[['interim_date', 'interim_month_year']]
    .drop_duplicates()
    .sort_values('interim_date')['interim_month_year']
    .tolist()
)

# %%

# =============================================================================
# Plot 1: per-(J, j) PPS bars, closed-form vs nested-MC HMC.
# Show one component j per facet column; J on the row.
# =============================================================================


def _js_show(J: int) -> list[int]:
    return sorted(set([0, J // 4, J // 2, (3 * J) // 4, J - 1]))


def _facet_subset(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for J in J_GRID:
        js = _js_show(J)
        sub = df[(df['J'] == J) & (df['j'].isin(js))].copy()
        sub['j_rank'] = sub['j'].map(lambda j: js.index(j))
        sub['j_label'] = sub['j'].map(lambda j: f'mu_{j} / J={J}')
        rows.append(sub)
    return pd.concat(rows, ignore_index=True)


pps_long = pd.concat(
    [
        _facet_subset(ppsa).assign(method=METHOD_ANA)
                          .rename(columns={'pps': 'pps_value'}),
        _facet_subset(pps_nm).assign(method=METHOD_HMC)
                             .rename(columns={'pps': 'pps_value'}),
    ],
    ignore_index=True,
)

# Bootstrap CI on the nested-MC PPS estimate per (J, j, interim).
bs_B = 2000
rng_bs = np.random.default_rng(123)
bs_rows = []
for (J, j), g in dp_h1_xz.groupby(['J', 'j']):
    wide = (
        g.pivot_table(index=['interim_id', 'interim_date'],
                      columns='s', values='p_h1_xz')
         .sort_index()
    )
    P = wide.to_numpy()
    n_interim, S = P.shape
    idx = rng_bs.integers(0, S, size=(n_interim, bs_B, S))
    P_bs = P[np.arange(n_interim)[:, None, None], idx]
    pps_bs = (P_bs > PPS_ProbH1_THRESH).mean(axis=2)
    q = np.quantile(pps_bs, [0.025, 0.975], axis=1).T
    part = pd.DataFrame(q, columns=['q025_bs', 'q975_bs'])
    part['interim_id'] = wide.index.get_level_values('interim_id').to_numpy()
    part['J'] = J
    part['j'] = j
    bs_rows.append(part)
bs = pd.concat(bs_rows, ignore_index=True)
pps_long = pps_long.merge(
    bs[['J', 'j', 'interim_id', 'q025_bs', 'q975_bs']],
    on=['J', 'j', 'interim_id'], how='left',
)
# Analytic has no MC noise; suppress its errorbars.
mask_ana = pps_long['method'] == METHOD_ANA
pps_long.loc[mask_ana, ['q025_bs', 'q975_bs']] = np.nan

pps_long['interim_month_year'] = pd.Categorical(
    pps_long['interim_month_year'], categories=interim_order, ordered=True,
)
pps_long['method'] = pd.Categorical(
    pps_long['method'], categories=[METHOD_ANA, METHOD_HMC], ordered=True,
)
fut2 = _futurama_palette(2)
pps_colours = {METHOD_ANA: '#000000', METHOD_HMC: fut2[0]}

p = (
    ggplot(pps_long, aes(x='interim_month_year', y='pps_value', fill='method'))
    + geom_col(position=position_dodge(width=0.8), width=0.7)
    + geom_errorbar(
        mapping=aes(x='interim_month_year', ymin='q025_bs', ymax='q975_bs',
                    group='method'),
        position=position_dodge(width=0.8), width=0.3,
        colour='black', size=0.4, inherit_aes=False,
    )
    + geom_hline(yintercept=PPS_ProbH1_THRESH, colour='black',
                 linetype='dashed', size=0.5)
    + facet_grid('J ~ j_rank', labeller='label_both', scales='free_y')
    + scale_fill_manual(values=pps_colours)
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        legend_position='top',
        figure_size=(16, 12),
    )
    + labs(
        x='interim', y='PPS = P( P(H_1j | x,z) > eta )', fill='method',
        title=('MVN per-component PPS: closed-form Phi-tail vs nested-MC '
               'HMC, faceted by J (rows) and component rank (cols)'),
    )
)
pdf_path = os.path.join(dir_out, 'mvn_compare_methods_pps.pdf')
p.save(pdf_path, verbose=False, limitsize=False)
print(f"Saved per-component PPS comparison to: {pdf_path}")

# %%

# =============================================================================
# Plot 2: per-(J, j) boxplot of nested-MC p(H_1 | x, z).
# =============================================================================

box_sub = _facet_subset(dp_h1_xz)
box_stats = (
    box_sub.groupby(['J', 'j', 'j_rank', 'interim_month_year'],
                    observed=True)['p_h1_xz']
    .agg(
        q025=lambda s: s.quantile(0.025),
        q25=lambda s: s.quantile(0.25),
        q50=lambda s: s.quantile(0.50),
        q75=lambda s: s.quantile(0.75),
        q975=lambda s: s.quantile(0.975),
    )
    .reset_index()
)
box_stats['interim_month_year'] = pd.Categorical(
    box_stats['interim_month_year'], categories=interim_order, ordered=True,
)

p = (
    ggplot(
        box_stats,
        aes(x='interim_month_year', ymin='q025', lower='q25',
            middle='q50', upper='q75', ymax='q975'),
    )
    + geom_boxplot(stat='identity', fill=fut2[0], alpha=0.45,
                   colour='#404040', size=0.3)
    + geom_hline(yintercept=PPS_ProbH1_THRESH, colour='black', size=1.0)
    + facet_grid('J ~ j_rank', labeller='label_both')
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        figure_size=(16, 12),
    )
    + labs(
        x='interim', y='p(H_1j | x, z) for predicted z samples',
        title=('MVN nested-MC HMC p(H_1j | x, z) distribution across the S '
               'predicted z samples, faceted by J (rows) and component '
               'rank (cols)'),
    )
)
pdf_path = os.path.join(dir_out, 'mvn_compare_methods_p_h1_xz_boxplot.pdf')
p.save(pdf_path, verbose=False, limitsize=False)
print(f"Saved nested-MC p(H_1j | x, z) boxplot to: {pdf_path}")

# %%

# =============================================================================
# Plot 3: per-J nested-MC HMC wall time per interim.
# =============================================================================

timing = (
    dp_h1_xz[['J', 'interim_id', 'interim_date', 'interim_month_year',
              'interim_mins']]
    .drop_duplicates()
    .rename(columns={'interim_mins': 'mins'})
)
timing['interim_month_year'] = pd.Categorical(
    timing['interim_month_year'], categories=interim_order, ordered=True,
)
timing['J_label'] = timing['J'].map(lambda J: f'J={J}')
J_colours = dict(zip(
    [f'J={J}' for J in J_GRID],
    _futurama_palette(len(J_GRID)),
))

p = (
    ggplot(timing, aes(x='interim_month_year', y='mins', fill='J_label'))
    + geom_col(position=position_dodge(width=0.8), width=0.7)
    + scale_fill_manual(values=J_colours)
    + scale_y_sqrt(expand=(0, 0, 0.15, 0))
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        legend_position='top',
        figure_size=(14, 5),
    )
    + labs(
        x='interim', y='nested-MC HMC time (mins)', fill='dimension',
        title=('MVN nested-MC HMC wall time per interim by problem '
               'dimension J (sqrt y-axis)'),
    )
)
pdf_path = os.path.join(dir_out, 'mvn_compare_methods_timing.pdf')
p.save(pdf_path, verbose=False, limitsize=False)
print(f"Saved nested-MC HMC timing plot to: {pdf_path}")

# %%

# =============================================================================
# Plot 4: closed-form vs nested-MC PPS scatter, faceted by J.
# =============================================================================

scatter = ppsa[['J', 'interim_id', 'j', 'pps']].rename(
    columns={'pps': 'pps_analytic'},
).merge(
    pps_nm[['J', 'interim_id', 'j', 'pps']].rename(
        columns={'pps': 'pps_hmc'},
    ),
    on=['J', 'interim_id', 'j'], how='inner',
)
scatter['J_label'] = scatter['J'].map(lambda J: f'J={J}')

p = (
    ggplot(scatter, aes(x='pps_analytic', y='pps_hmc', colour='J_label'))
    + geom_abline(slope=1, intercept=0, linetype='dashed', colour='black')
    + geom_point(alpha=0.7, size=2)
    + facet_wrap('~ J_label', ncol=2)
    + scale_colour_manual(values=J_colours)
    + theme_bw()
    + theme(
        legend_position='none',
        figure_size=(10, 9),
    )
    + labs(
        x='PPS closed-form (analytic Phi-tail)',
        y='PPS nested-MC HMC',
        title=('MVN PPS: analytic Phi-tail vs nested-MC HMC, per '
               '(interim, component j). Dashed line = identity.'),
    )
)
pdf_path = os.path.join(dir_out, 'mvn_compare_methods_pps_scatter.pdf')
p.save(pdf_path, verbose=False, limitsize=False)
print(f"Saved closed-form vs HMC scatter to: {pdf_path}")

# %%

print(f"\nCross-J comparison complete.")
print(f"Outputs in: {dir_out}")
