#!/usr/bin/env python3
"""
MVN interim analyses: cross-method comparison (§3.3 of
``dev/amortised_decision_making.md``).

Standalone. Loads per-J per-method artifacts from:

  - HMC nested-MC : ``MVN_interim_analyses_with_nested_MC_hmc.py``
  - IS reweighting : ``MVN_interim_analyses_with_IS_from_x.py``
  - Regression endpt-x (Gaussian approx)  : ``MVN_interim_analysis_regression_endptx_on_wz_with_gauss_approx.py``
  - Regression endpt-x (quantile regr)    : ``MVN_interim_analysis_regression_endptx_on_wz_with_quantile_regr.py``
  - Regression endpt-x (multi-quantile)   : ``MVN_interim_analysis_regression_endptx_on_wz_with_mquantile_regr.py``

and the closed-form analytic PPS from the simulations directory. Emits
four PDFs per J:

  1. ``mvn_J{J}_compare_methods_p_h1_xz_all.pdf`` -- per-(response j,
     interim) p(H_1j | x, z) quantile boxplots faceted in a
     ``K_levels x n_per_level`` grid (same layout as
     ``mvn_J{J}_pps_RGE_p_h1_xz_all.pdf``), dodged by method, Futurama
     fill per method.
  2. ``mvn_J{J}_compare_methods_pps_all.pdf`` -- PPS bar chart in the
     same facet grid, analytic (black) + each method in a Futurama
     colour; 95% bootstrap CIs on the Monte-Carlo bars.
  3. ``mvn_J{J}_compare_methods_timing.pdf`` -- per-interim inference
     time (mins) per method, Futurama colours, sqrt y-axis.
  4. ``mvn_J{J}_compare_methods_ess.pdf`` -- median ESS / particle per
     interim for the IS method (other methods don't have IS weights).

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/MVN_interim_analyses_compare_methods.py
"""

# %%

import os
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
    ggplot, aes, geom_boxplot, geom_col, geom_errorbar, geom_hline,
    geom_text,
    facet_grid, facet_wrap,
    scale_fill_manual, scale_y_continuous, scale_y_sqrt,
    position_dodge, theme_bw, theme, element_text, element_blank, labs,
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
DIR_OUT  = os.path.join(_sandbox, "py-mvn-interim-compare-methods-260609")
os.makedirs(DIR_OUT, exist_ok=True)

METHOD_ORDER = [
    'nested-MC using HMC for each (x,z)',
    'IS reweighting of theta|x',
    'Regression of endpt-x on w(z) using Gaussian approx',
    'Regression of endpt-x on w(z) - quantile',
    'Regression of endpt-x on w(z) - mquantile',
]
method_colours = dict(zip(METHOD_ORDER, _futurama_palette(len(METHOD_ORDER))))
pps_colours = {'analytic': '#000000', **method_colours}

# %%

# =============================================================================
# Load upstream artifacts.
# =============================================================================

sim_data = pd.read_pickle(os.path.join(DIR_SIM, 'mvn_sim_data.pkl'))
simu_params = sim_data['simu_params']
di = sim_data['interim_grid']
J_GRID = list(simu_params['J_grid'])
K_levels = simu_params['K_levels']
pps_ProbH1_thresh = simu_params['pps_ProbH1_thresh']

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

    # ---- Analytic per-z (loaded from nested-MC dir, computed by
    # ---- model.fit_closed_form_pH1 in nested_MC_hmc.py Phase 4.5) ----
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

    # ---- RGE per-J slices (Gaussian approx) ----
    rge_p = pd.read_pickle(
        os.path.join(DIR_RGE, f'mvn_J{J}_pps_RGE_p_h1_xz.pkl'),
    )[['interim_id', 'interim_date', 'interim_month_year', 'j', 's',
       'p_h1_xz']].copy()
    rge_pps = pd.read_csv(
        os.path.join(DIR_RGE, f'mvn_J{J}_pps_RGE.csv'),
    )[['interim_id', 'interim_date', 'interim_month_year', 'j',
       'pps']].copy()
    rge_timing = pd.read_csv(
        os.path.join(DIR_RGE, f'mvn_J{J}_pps_RGE_timing.csv'),
    )[['interim_id', 'interim_date', 'interim_month_year',
       'mins_interim_id']].rename(columns={'mins_interim_id': 'mins'}).copy()

    # ---- RGEQ per-J slices (quantile regression) ----
    rgeq_p = pd.read_pickle(
        os.path.join(DIR_RGEQ, f'mvn_J{J}_pps_RGEQ_p_h1_xz.pkl'),
    )[['interim_id', 'interim_date', 'interim_month_year', 'j', 's',
       'p_h1_xz']].copy()
    rgeq_pps = pd.read_csv(
        os.path.join(DIR_RGEQ, f'mvn_J{J}_pps_RGEQ.csv'),
    )[['interim_id', 'interim_date', 'interim_month_year', 'j',
       'pps']].copy()
    rgeq_timing = pd.read_csv(
        os.path.join(DIR_RGEQ, f'mvn_J{J}_pps_RGEQ_timing.csv'),
    )[['interim_id', 'interim_date', 'interim_month_year',
       'mins_interim_id']].rename(columns={'mins_interim_id': 'mins'}).copy()

    # ---- RGEM per-J slices (multi-quantile regression) ----
    rgem_p = pd.read_pickle(
        os.path.join(DIR_RGEM, f'mvn_J{J}_pps_RGEM_p_h1_xz.pkl'),
    )[['interim_id', 'interim_date', 'interim_month_year', 'j', 's',
       'p_h1_xz']].copy()
    rgem_pps = pd.read_csv(
        os.path.join(DIR_RGEM, f'mvn_J{J}_pps_RGEM.csv'),
    )[['interim_id', 'interim_date', 'interim_month_year', 'j',
       'pps']].copy()
    rgem_timing = pd.read_csv(
        os.path.join(DIR_RGEM, f'mvn_J{J}_pps_RGEM_timing.csv'),
    )[['interim_id', 'interim_date', 'interim_month_year',
       'mins_interim_id']].rename(columns={'mins_interim_id': 'mins'}).copy()

    # ---- Stack methods ----
    p_h1_xz = pd.concat(
        [
            an_p  .assign(method='analytic'),
            hmc_p .assign(method='nested-MC using HMC for each (x,z)'),
            is_p  .assign(method='IS reweighting of theta|x'),
            rge_p .assign(method='Regression of endpt-x on w(z) using Gaussian approx'),
            rgeq_p.assign(method='Regression of endpt-x on w(z) - quantile'),
            rgem_p.assign(method='Regression of endpt-x on w(z) - mquantile'),
        ],
        ignore_index=True,
    )
    pps = pd.concat(
        [
            pps_cf[pps_cf['J'] == J][[
                'interim_id', 'interim_date', 'interim_month_year', 'j', 'pps',
            ]].assign(method='analytic'),
            hmc_pps .assign(method='nested-MC using HMC for each (x,z)'),
            is_pps  .assign(method='IS reweighting of theta|x'),
            rge_pps .assign(method='Regression of endpt-x on w(z) using Gaussian approx'),
            rgeq_pps.assign(method='Regression of endpt-x on w(z) - quantile'),
            rgem_pps.assign(method='Regression of endpt-x on w(z) - mquantile'),
        ],
        ignore_index=True,
    )
    timing = pd.concat(
        [
            hmc_timing .assign(method='nested-MC using HMC for each (x,z)'),
            is_timing  .assign(method='IS reweighting of theta|x'),
            rge_timing .assign(method='Regression of endpt-x on w(z) using Gaussian approx'),
            rgeq_timing.assign(method='Regression of endpt-x on w(z) - quantile'),
            rgem_timing.assign(method='Regression of endpt-x on w(z) - mquantile'),
        ],
        ignore_index=True,
    )

    # Normalise mixed Timestamp/str dtypes (CSV loads as str).
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
        categories=['analytic'] + METHOD_ORDER, ordered=True,
    )
    pps['method'] = pd.Categorical(
        pps['method'], categories=['analytic'] + METHOD_ORDER, ordered=True,
    )
    timing['method'] = pd.Categorical(
        timing['method'], categories=METHOD_ORDER, ordered=True,
    )

    p_h1_xz_g = p_h1_xz[p_h1_xz['j'].isin(j_indices_all)].copy()
    pps_g     = pps[pps['j'].isin(j_indices_all)].copy()
    for df in (p_h1_xz_g, pps_g):
        df['response_label'] = pd.Categorical(
            'response mu_' + df['j'].astype(str),
            categories=response_label_cats, ordered=True,
        )

    # ---- Plot 1: p(H_1j | x, z) boxplot, faceted by response, dodged by method.
    # ---- Drop quantile method: its p_h1_xz is a hard 0/1 label, not a
    # ---- probability, so the boxplot would mix semantically-different quantities.
    # ---- Quantile method still appears in PPS bars, timing, and MSE below.
    p_h1_xz_g_box = p_h1_xz_g[
        p_h1_xz_g['method'] != 'Regression of endpt-x on w(z) - quantile'
    ].copy()
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

    p = (
        ggplot(box_stats, aes(x='interim_month_year', fill='method'))
        + geom_boxplot(
            aes(ymin='q025', lower='q25', middle='q50', upper='q75',
                ymax='q975', group='grp'),
            stat='identity', position=position_dodge(width=0.8), width=0.7,
            colour='#404040', size=0.3,
        )
        + geom_hline(yintercept=pps_ProbH1_thresh, colour='black', size=1.0)
        + scale_fill_manual(values=pps_colours)
        + scale_y_continuous(
            limits=[0, 1],
            breaks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            labels=['0%', '20%', '40%', '60%', '80%', '100%'],
        )
        + facet_wrap('~ response_label', ncol=n_per_level, dir='h')
        + theme_bw()
        + theme(
            axis_text_x=element_text(angle=45, vjust=1, hjust=1),
            legend_position='top',
            figure_size=(3.0 * n_per_level, 2.5 * K_levels),
            strip_background=element_blank(),
            strip_text=element_text(face='bold'),
        )
        + labs(x='Interim', y='p(H_1j | x, z) for predicted z samples',
               fill='method')
    )
    pdf_path = os.path.join(
        DIR_OUT, f'mvn_J{J}_compare_methods_p_h1_xz_all.pdf',
    )
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved p(H_1 | x, z) boxplot to {pdf_path}")

    # ---- Plot 2: PPS bars per (response, interim), dodged by method,
    # ---- bootstrap CI on Monte-Carlo bars ----
    bs_B = 2000
    rng_bs = np.random.default_rng(simu_params['seed'])
    bs_parts = []
    # Skip 'analytic' for bootstrap: it has no Monte-Carlo noise.
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
        pps_bs = (Pbs > pps_ProbH1_thresh).mean(axis=2)
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

    p = (
        ggplot(pps_g_ci, aes(x='interim_month_year', y='pps', fill='method'))
        + geom_col(position=position_dodge(width=0.8), width=0.7)
        + geom_errorbar(
            aes(x='interim_month_year',
                ymin='q025_bs', ymax='q975_bs', group='method'),
            position=position_dodge(width=0.8), width=0.3,
            colour='black', size=0.4, inherit_aes=False,
        )
        + geom_hline(yintercept=pps_ProbH1_thresh, colour='black', size=1.0,
                     linetype='dashed')
        + scale_fill_manual(values=pps_colours)
        + scale_y_continuous(
            limits=[0, 1.05],
            breaks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            labels=['0%', '20%', '40%', '60%', '80%', '100%'],
        )
        + facet_wrap('~ response_label', ncol=n_per_level, dir='h')
        + theme_bw()
        + theme(
            axis_text_x=element_text(angle=45, vjust=1, hjust=1),
            legend_position='top',
            figure_size=(3.0 * n_per_level, 2.5 * K_levels),
            strip_background=element_blank(),
            strip_text=element_text(face='bold'),
        )
        + labs(x='Interim',
               y='PPS = int P(p(H_1j | x, z) > eta) dz',
               fill='method')
    )
    pdf_path = os.path.join(DIR_OUT, f'mvn_J{J}_compare_methods_pps_all.pdf')
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved PPS comparison plot to {pdf_path}")

    # ---- Plot 3: Timing, methods dodged, sqrt y-axis, futurama colours ----
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
            strip_background=element_blank(),
            strip_text=element_text(face='bold'),
        )
        + labs(x='Interim', y='time (mins)', fill='method')
    )
    pdf_path = os.path.join(DIR_OUT, f'mvn_J{J}_compare_methods_timing.pdf')
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved timing comparison plot to {pdf_path}")

    # ---- Plot 4: IS ESS / particle per interim ----
    is_ess = (
        is_perf.groupby('interim_month_year', observed=True)
        .agg(value=('ess_over_n', 'median'))
        .reset_index()
        .assign(method='IS reweighting of theta|x')
    )
    is_ess['method'] = pd.Categorical(
        is_ess['method'], categories=METHOD_ORDER, ordered=True,
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
            strip_background=element_blank(),
            strip_text=element_text(face='bold'),
        )
        + labs(x='Interim', y='ESS / particle', fill='method')
    )
    pdf_path = os.path.join(DIR_OUT, f'mvn_J{J}_compare_methods_ess.pdf')
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved IS ESS / particle plot to {pdf_path}")

    # ---- MSE: per (interim, method) average over all J components of
    # ---- (pps_method - pps_analytic)^2. Uses the full per-J pps table
    # ---- (not the j_indices_all subset) so MSE is over all responses.
    pps_pivot = (
        pps.pivot_table(
            index=['interim_id', 'interim_date', 'interim_month_year', 'j'],
            columns='method', values='pps',
        ).reset_index()
    )
    for method in METHOD_ORDER:
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
# fill = algorithm, facet = J. Standard Futurama palette.
# =============================================================================

mse_all = pd.concat(mse_rows, ignore_index=True)
mse_all['method'] = pd.Categorical(
    mse_all['method'], categories=METHOD_ORDER, ordered=True,
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
    + geom_col(position=position_dodge(width=0.8), width=0.7)
    + facet_wrap('~ J_label', ncol=1, scales='free_y')
    + scale_fill_manual(values=method_colours)
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        legend_position='top',
        figure_size=(14, 4 * len(J_GRID)),
        strip_background=element_blank(),
        strip_text=element_text(face='bold'),
    )
    + labs(x='Interim date',
           y='MSE of PPS vs analytic (mean over responses j)',
           fill='method')
)
pdf_path = os.path.join(DIR_OUT, 'mvn_compare_methods_mse.pdf')
p.save(pdf_path, verbose=False, limitsize=False)
print(f"Saved cross-J MSE plot to {pdf_path}")

# %%

# =============================================================================
# Regression-only: violin plot of per-response R^2 (top row) and
# empirical rho (bottom row) across all J responses j, faceted by J in
# columns ordered J=20, 60, 100. Source: rge_stats pkls written by the
# regression-endptx script.
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

    # Pre-compute (q025, q25, q50, q75, q975) per (J, metric, interim) for
    # the stat='identity' boxplot pattern used elsewhere in this codebase.
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
