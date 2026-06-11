#!/usr/bin/env python3
"""
MVN nested-MC HMC PPS (§3.3 of dev/amortised_decision_making.md).

Standalone: loads the simulated x cohort from
``mvn_sim_data.pkl``, the per-(J, interim) future-data blocks from
``mvn_J{J}_interim_data.pkl``, and the closed-form analytic PPS from
``mvn_pps_closed_form.pkl`` -- all produced upstream in the simulations
directory. Runs nested-MC HMC PPS in a flat loop over (J, interim) and
emits per-J plots adapted from the Binomial nested-MC layout.

If ``mvn_pps_nested_mc.pkl`` already exists in the nested-MC output
directory the Phase 5 loop is skipped and per-J plots are regenerated
from the cached artifacts.

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/MVN_interim_analyses_with_nested_MC_hmc.py
"""

# %%

import os
import sys
import time
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
    ggplot, aes, geom_abline, geom_boxplot, geom_col, geom_errorbar,
    geom_hline, geom_point,
    facet_wrap, scale_color_manual, scale_fill_manual,
    scale_x_continuous, scale_x_datetime, scale_y_continuous,
    position_dodge, theme_bw, theme,
    element_text, element_blank, labs,
)

warnings.filterwarnings('ignore')

from model_mvn import MVNModel
from fit_interim import fit_interim_posterior_xz_with_nested_monte_carlo
from utils import _material_palette

print("Imports successful")

# %%

# =============================================================================
# Configuration. The simulations directory holds the upstream artifacts
# (``mvn_sim_data.pkl``, ``mvn_J{J}_interim_data.pkl``,
# ``mvn_pps_closed_form.pkl``); ``DIR_OUT`` is the dedicated nested-MC
# HMC output directory.
# =============================================================================

DIR_SIM = os.path.join(
    "/Users/or105/sandbox/bIRTistic", "py-mvn-interim-simulations-260609",
)
DIR_OUT = os.path.join(
    "/Users/or105/sandbox/bIRTistic", "py-mvn-interim-with-nested-MC-hmc-260609",
)
os.makedirs(DIR_OUT, exist_ok=True)

PPS_Z_TOTAL = 200
HMC_CHAINS = 2
HMC_ITER_WARMUP = 200
HMC_ITER_SAMPLING = 300

# %%

# =============================================================================
# Load upstream artifacts.
# =============================================================================

sim_data = pd.read_pickle(os.path.join(DIR_SIM, 'mvn_sim_data.pkl'))
simu_params = sim_data['simu_params']
di = sim_data['interim_grid']
J_GRID = list(simu_params['J_grid'])
pps_cf = pd.read_pickle(os.path.join(DIR_SIM, 'mvn_pps_closed_form.pkl'))['pps_cf']

interim_data_by_J = {}
for J in J_GRID:
    pkl_path = os.path.join(DIR_SIM, f'mvn_J{J}_interim_data.pkl')
    if not os.path.exists(pkl_path):
        raise FileNotFoundError(
            f"{pkl_path} missing; "
            f"run MVN_interim_analyses_make_interim_data.py first."
        )
    interim_data_by_J[J] = pd.read_pickle(pkl_path)

print(f"Loaded sim_data + {len(J_GRID)} interim_data J cells + pps_cf "
      f"({len(pps_cf)} rows).")

nm_pkl = os.path.join(DIR_OUT, 'mvn_pps_nested_mc.pkl')

# %%

# =============================================================================
# Phase 4.5. Analytic per-component P(H_1j | x, z_s) using the cached
# zi for fair comparison with HMC. One pkl per J.
# =============================================================================

an_pkls_by_J = {
    J: os.path.join(DIR_OUT, f'mvn_J{J}_pps_p_h1_xz_analytic.pkl')
    for J in J_GRID
}
an_by_J = {}
for J in J_GRID:
    if os.path.exists(an_pkls_by_J[J]):
        print(f"[analytic] J={J} loading cached {an_pkls_by_J[J]}")
        an_by_J[J] = pd.read_pickle(an_pkls_by_J[J])
        continue
    cell = sim_data['cells'][J]
    dit, R_chol = cell['dit'], cell['R_chol']
    interim_data = interim_data_by_J[J]
    an_parts = []
    for interim_id, blk in interim_data.items():
        model_x = MVNModel(
            dit=dit, dcati=blk['dpi'], seed=simu_params['seed'],
            J=J, K_chol=R_chol, sigma=simu_params['sigma'],
            prior_tau=simu_params['prior_tau'],
            mu_0_baseline=simu_params['mu_0_baseline'],
        )
        an = model_x.fit_closed_form_pH1(blk['zi'])
        an['J'] = J
        an['interim_id'] = interim_id
        an['interim_date'] = blk['interim_date']
        an['interim_month_year'] = blk['interim_month_year']
        an_parts.append(an)
        print(f"[analytic] J={J} interim {interim_id}: computed"
              f" {len(an)} rows ({an['s'].nunique()} samples"
              f" x {an['j'].nunique()} components).")
    an_full = pd.concat(an_parts, ignore_index=True)
    pd.to_pickle(an_full, an_pkls_by_J[J])
    print(f"[analytic] J={J} saved {len(an_full)} rows to"
          f" {an_pkls_by_J[J]}")
    an_by_J[J] = an_full

# %%

# =============================================================================
# Phase 5. Nested-MC HMC PPS, flat loop over (J, interim). For each cell
# uses the cached zi block; the inner per-sample refits construct fresh
# s_models on (xi + z_s) inside the driver.
# =============================================================================

if os.path.exists(nm_pkl):
    print(f"\n[nested-MC] loading cached artifacts from {nm_pkl}")
    _cached = pd.read_pickle(nm_pkl)
    nm        = _cached['p_h1_xz']
    nm_timing = _cached['timing']
    nm_pps    = _cached['pps']
else:
    nm_parts = []
    nm_timing_parts = []
    t_all0 = time.time()
    for J in J_GRID:
        cell = sim_data['cells'][J]
        dit, R_chol = cell['dit'], cell['R_chol']
        model_init_kwargs = dict(
            J=J, K_chol=R_chol, sigma=simu_params['sigma'],
            prior_tau=simu_params['prior_tau'],
            mu_0_baseline=simu_params['mu_0_baseline'],
        )
        interim_data = interim_data_by_J[J]
        for interim_id, blk in interim_data.items():
            interim_date = blk['interim_date']
            interim_month_year = blk['interim_month_year']
            dpi = blk['dpi']
            zi = blk['zi']
            n_obs = blk['n_obs']
            interim_m = blk['interim_m']

            print(f"\n{'=' * 70}\nJ={J} PPS interim {interim_id}"
                  f" ({interim_date.date()}) n={n_obs} m={interim_m}"
                  f"\n{'=' * 70}")
            t0 = time.time()

            model_x = MVNModel(
                dit=dit, dcati=dpi, seed=simu_params['seed'],
                **model_init_kwargs,
            )
            tmp = fit_interim_posterior_xz_with_nested_monte_carlo(
                model_x, zi,
                interim_method_args={
                    'pps_z_total':        PPS_Z_TOTAL,
                    'pps_H1_def':         simu_params['pps_H1_def'],
                    'pps_ProbH1_thresh':  simu_params['pps_ProbH1_thresh'],
                    'fit_method':         'fit_pyro_hmc',
                    'seed':               simu_params['seed'],
                    'save_to_file':       False,
                    'verbose':            False,
                    'cpu_n':              1,
                    'output_file_prefix': None,
                    'model_init_kwargs':  model_init_kwargs,
                },
                fit_method_args={
                    'x_formula':     '~ 1',
                    'chains':        HMC_CHAINS,
                    'iter_warmup':   HMC_ITER_WARMUP,
                    'iter_sampling': HMC_ITER_SAMPLING,
                },
            )
            tmp['J'] = J
            tmp['interim_id'] = interim_id
            tmp['interim_date'] = interim_date
            tmp['interim_month_year'] = interim_month_year
            interim_mins = (time.time() - t0) / 60.0
            tmp['interim_mins'] = round(interim_mins, 3)
            nm_parts.append(tmp)
            nm_timing_parts.append({
                'J':                  J,
                'interim_id':         interim_id,
                'interim_date':       interim_date,
                'interim_month_year': interim_month_year,
                'n_obs':              n_obs,
                'm_future':           interim_m,
                'mins':               interim_mins,
            })
            print(f"  J={J} interim {interim_id} PPS done in"
                  f" {interim_mins:.2f} min")

    if not nm_parts:
        raise RuntimeError("[nested-MC] no interims with m > 0.")
    nm = pd.concat(nm_parts, ignore_index=True)
    nm = nm.assign(
        j=nm['item_label'].str.replace('mu_', '').astype(int),
    )
    nm_timing = pd.DataFrame(nm_timing_parts)
    nm_pps = (
        nm.groupby(['J', 'interim_id', 'interim_date', 'interim_month_year',
                    'j', 'item_label'], observed=True)['p_h1_xz']
          .apply(lambda p: float(
              (p > simu_params['pps_ProbH1_thresh']).mean()
          ))
          .reset_index(name='pps')
    )
    nm_pps['eta'] = simu_params['pps_ProbH1_thresh']
    nm_pps['S'] = PPS_Z_TOTAL
    pd.to_pickle({'p_h1_xz': nm, 'timing': nm_timing, 'pps': nm_pps}, nm_pkl)
    print(f"\nSaved nested-MC HMC artifacts to {nm_pkl};"
          f" total time {(time.time() - t_all0) / 60.0:.2f} min")

# %%

# =============================================================================
# Phase 6. Per-J plots (mirror of Binomial nested_MC_hmc.py line 475ff).
# Boxplot of p(H_1j | x, z) per (interim, j) on K_levels representative
# components; PPS bar chart per (interim, j) with closed-form analytic
# + HMC nested-MC dodged using material-teal palette for HMC bars
# (legend collapses to {analytic, HMC}); 4xn_per_level grid 'all'
# variants; scatter (analytic per-sample vs HMC per-sample).
# =============================================================================

bs_B = 2000
K_levels = simu_params['K_levels']

for J in sorted(set(nm['J'].unique()) & set(J_GRID)):
    block_size_J = J // K_levels
    j_repr = [block_size_J // 2 + block_size_J * k for k in range(K_levels)]
    nm_J = nm[(nm['J'] == J) & (nm['j'].isin(j_repr))].copy()
    cf_J = pps_cf[(pps_cf['J'] == J) & (pps_cf['j'].isin(j_repr))].copy()
    nm_pps_J = nm_pps[(nm_pps['J'] == J) & (nm_pps['j'].isin(j_repr))].copy()
    if nm_J.empty:
        continue
    nm_J['j_label']     = 'response mu_' + nm_J['j'].astype(str)
    cf_J['j_label']     = 'response mu_' + cf_J['j'].astype(str)
    nm_pps_J['j_label'] = 'response mu_' + nm_pps_J['j'].astype(str)
    j_label_cats = [f'response mu_{j}' for j in j_repr]

    order = (
        di[di['interim_id'].isin(nm_J['interim_id'].unique())]
        .sort_values('interim_date')['interim_month_year'].tolist()
    )
    for df in (nm_J, cf_J, nm_pps_J):
        df['interim_month_year'] = pd.Categorical(
            df['interim_month_year'], categories=order, ordered=True,
        )
        df['j_label'] = pd.Categorical(
            df['j_label'], categories=j_label_cats, ordered=True,
        )

    hmc_palette = _material_palette('teal', len(j_repr))
    fill_values = {'analytic': '#000000'}
    fill_values.update({f'mu_{j}': c for j, c in zip(j_repr, hmc_palette)})
    hmc_repr_key = f'mu_{j_repr[-1]}'

    # ---- Boxplot of p(H_1j | x, z): analytic vs HMC dodged ----
    an_J_box = an_by_J.get(J, pd.DataFrame())
    if not an_J_box.empty:
        an_J_box = an_J_box[an_J_box['j'].isin(j_repr)].copy()
        an_J_box['j_label'] = 'response mu_' + an_J_box['j'].astype(str)
        an_J_box['interim_month_year'] = pd.Categorical(
            an_J_box['interim_month_year'], categories=order, ordered=True,
        )
        an_J_box['j_label'] = pd.Categorical(
            an_J_box['j_label'], categories=j_label_cats, ordered=True,
        )

    def _quantile_table(df, method):
        q = (
            df.groupby(['j_label', 'interim_month_year', 'j'],
                       observed=True)['p_h1_xz']
              .agg(
                  q025=lambda s: s.quantile(0.025),
                  q25 =lambda s: s.quantile(0.25),
                  q50 =lambda s: s.quantile(0.50),
                  q75 =lambda s: s.quantile(0.75),
                  q975=lambda s: s.quantile(0.975),
              ).reset_index()
        )
        q['method'] = method
        return q

    box = pd.concat(
        [_quantile_table(an_J_box, 'analytic'),
         _quantile_table(nm_J, 'HMC')],
        ignore_index=True,
    )
    box['method'] = pd.Categorical(
        box['method'], categories=['analytic', 'HMC'], ordered=True,
    )
    box['fill_key'] = np.where(
        box['method'].astype(str) == 'analytic',
        'analytic',
        'mu_' + box['j'].astype(str),
    )
    box['fill_key'] = pd.Categorical(
        box['fill_key'],
        categories=['analytic'] + [f'mu_{j}' for j in j_repr],
        ordered=True,
    )
    box['grp'] = (
        box['interim_month_year'].astype(str) + '|'
        + box['fill_key'].astype(str)
    )
    box.to_pickle(os.path.join(DIR_OUT, f'mvn_J{J}_pps_p_h1_xz.pkl'))

    p = (
        ggplot(box, aes(x='interim_month_year',
                        ymin='q025', lower='q25', middle='q50',
                        upper='q75', ymax='q975',
                        fill='fill_key', group='grp'))
        + geom_boxplot(stat='identity',
                       position=position_dodge(width=0.8), width=0.7,
                       colour='#404040', size=0.3)
        + geom_hline(yintercept=simu_params['pps_ProbH1_thresh'],
                     colour='black', size=1.0)
        + scale_fill_manual(values=fill_values,
                            breaks=['analytic', hmc_repr_key],
                            labels=['analytic', 'HMC'],
                            name='method')
        + scale_y_continuous(
            breaks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            labels=['0%', '20%', '40%', '60%', '80%', '100%'],
        )
        + facet_wrap('~ j_label', ncol=K_levels)
        + theme_bw()
        + theme(axis_text_x=element_text(angle=45, vjust=1, hjust=1),
                figure_size=(14, 5),
                legend_position='top',
                strip_background=element_blank(),
                strip_text=element_text(face='bold'))
        + labs(x='Interim',
               y='p(H_1j | x, z) for predicted z samples')
    )
    pdf_path = os.path.join(DIR_OUT, f'mvn_J{J}_pps_p_h1_xz.pdf')
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved p(H_1 | x, z) boxplot to {pdf_path}")

    # ---- PPS bars: closed-form + HMC dodged + bootstrap CI ----
    wide = (
        nm_J.pivot_table(index=['j', 'j_label', 'interim_id', 'interim_date'],
                         columns='s', values='p_h1_xz')
        .sort_index()
    )
    bs = wide.to_numpy()
    n_rows, S = bs.shape
    rng = np.random.default_rng(simu_params['seed'])
    idx = rng.integers(0, S, size=(n_rows, bs_B, S))
    bs = bs[np.arange(n_rows)[:, None, None], idx]
    bs = (bs > simu_params['pps_ProbH1_thresh']).mean(axis=2)
    ci = np.quantile(bs, [0.025, 0.975], axis=1).T
    ci_df = pd.DataFrame(ci, columns=['q025_bs', 'q975_bs'])
    ci_df['j']          = wide.index.get_level_values('j').to_numpy()
    ci_df['interim_id'] = wide.index.get_level_values('interim_id').to_numpy()

    hmc_bar = nm_pps_J[['interim_id', 'interim_date',
                        'interim_month_year', 'j', 'j_label',
                        'pps']].assign(method='HMC').merge(
        ci_df, on=['j', 'interim_id'], how='left',
    )
    cf_bar = cf_J[['interim_id', 'interim_date',
                   'interim_month_year', 'j', 'j_label',
                   'pps']].assign(method='analytic',
                                  q025_bs=np.nan, q975_bs=np.nan)
    bar = pd.concat([cf_bar, hmc_bar], ignore_index=True)
    bar['interim_date'] = pd.to_datetime(bar['interim_date'])
    bar['method'] = pd.Categorical(
        bar['method'], categories=['analytic', 'HMC'], ordered=True,
    )
    bar['j_label'] = pd.Categorical(
        bar['j_label'], categories=j_label_cats, ordered=True,
    )
    bar['fill_key'] = np.where(
        bar['method'].astype(str) == 'analytic',
        'analytic',
        'mu_' + bar['j'].astype(str),
    )
    bar['fill_key'] = pd.Categorical(
        bar['fill_key'],
        categories=['analytic'] + [f'mu_{j}' for j in j_repr],
        ordered=True,
    )
    bar.to_pickle(os.path.join(DIR_OUT, f'mvn_J{J}_pps_bars.pkl'))

    p = (
        ggplot(bar, aes(x='interim_date', y='pps', fill='fill_key'))
        + geom_col(position=position_dodge(width=20), width=18)
        + geom_errorbar(
            aes(x='interim_date',
                ymin='q025_bs', ymax='q975_bs', group='method'),
            position=position_dodge(width=20), width=8,
            colour='black', size=0.5, inherit_aes=False,
        )
        + scale_fill_manual(values=fill_values,
                            breaks=['analytic', hmc_repr_key],
                            labels=['analytic', 'HMC'],
                            name='method')
        + scale_x_datetime(date_breaks='1 month', date_labels='%Y-%b')
        + scale_y_continuous(
            breaks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            labels=['0%', '20%', '40%', '60%', '80%', '100%'],
        )
        + facet_wrap('~ j_label', ncol=K_levels)
        + theme_bw()
        + theme(axis_text_x=element_text(angle=45, vjust=1, hjust=1),
                figure_size=(14, 5),
                legend_position='top',
                strip_background=element_blank(),
                strip_text=element_text(face='bold'))
        + labs(x='interim date',
               y='PPS = int P(p(H_1j | x, z) > eta) dz')
    )
    pdf_path = os.path.join(DIR_OUT, f'mvn_J{J}_pps_bars.pdf')
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved PPS bar chart to {pdf_path}")

    # ---- "all" PPS bars (K_levels x n_per_level grid) ----
    n_per_level = min(5, block_size_J)
    j_indices_all = np.array([
        k * block_size_J + c for k in range(K_levels)
        for c in range(n_per_level)
    ])
    nm_all_J = nm[(nm['J'] == J) & (nm['j'].isin(j_indices_all))].copy()
    cf_all_J = pps_cf[(pps_cf['J'] == J) & (pps_cf['j'].isin(j_indices_all))].copy()
    nm_pps_all_J = nm_pps[(nm_pps['J'] == J) & (nm_pps['j'].isin(j_indices_all))].copy()

    for df in (nm_all_J, cf_all_J, nm_pps_all_J):
        df['k_level'] = ('block ' + (df['j'] // block_size_J).astype(str))
        df['within_block_idx'] = (df['j'] % block_size_J).astype(int)
        df['response_label'] = 'response mu_' + df['j'].astype(str)

    k_level_cats = [f'block {k}' for k in range(K_levels)]
    response_label_cats = [f'response mu_{j}' for j in j_indices_all]

    wide_all = (
        nm_all_J.pivot_table(index=['j', 'interim_id', 'interim_date'],
                             columns='s', values='p_h1_xz')
        .sort_index()
    )
    bs_all = wide_all.to_numpy()
    n_rows_a, S_a = bs_all.shape
    rng = np.random.default_rng(simu_params['seed'])
    idx = rng.integers(0, S_a, size=(n_rows_a, bs_B, S_a))
    bs_all = bs_all[np.arange(n_rows_a)[:, None, None], idx]
    bs_all = (bs_all > simu_params['pps_ProbH1_thresh']).mean(axis=2)
    ci_all = np.quantile(bs_all, [0.025, 0.975], axis=1).T
    ci_df_all = pd.DataFrame(ci_all, columns=['q025_bs', 'q975_bs'])
    ci_df_all['j']          = wide_all.index.get_level_values('j').to_numpy()
    ci_df_all['interim_id'] = wide_all.index.get_level_values('interim_id').to_numpy()

    hmc_bar_all = nm_pps_all_J[[
        'interim_id', 'interim_date', 'interim_month_year',
        'j', 'k_level', 'within_block_idx', 'response_label', 'pps',
    ]].assign(method='HMC').merge(
        ci_df_all, on=['j', 'interim_id'], how='left',
    )
    cf_bar_all = cf_all_J[[
        'interim_id', 'interim_date', 'interim_month_year',
        'j', 'k_level', 'within_block_idx', 'response_label', 'pps',
    ]].assign(method='analytic', q025_bs=np.nan, q975_bs=np.nan)
    bar_all = pd.concat([cf_bar_all, hmc_bar_all], ignore_index=True)
    bar_all['interim_date'] = pd.to_datetime(bar_all['interim_date'])
    bar_all['method'] = pd.Categorical(
        bar_all['method'], categories=['analytic', 'HMC'], ordered=True,
    )
    bar_all['k_level'] = pd.Categorical(
        bar_all['k_level'], categories=k_level_cats, ordered=True,
    )
    bar_all['response_label'] = pd.Categorical(
        bar_all['response_label'], categories=response_label_cats, ordered=True,
    )

    hmc_palette_all = _material_palette('teal', len(j_indices_all))
    fill_values_all = {'analytic': '#000000'}
    fill_values_all.update({
        f'mu_{j}': c for j, c in zip(j_indices_all, hmc_palette_all)
    })
    bar_all['fill_key'] = np.where(
        bar_all['method'].astype(str) == 'analytic',
        'analytic',
        'mu_' + bar_all['j'].astype(str),
    )
    bar_all['fill_key'] = pd.Categorical(
        bar_all['fill_key'],
        categories=['analytic'] + [f'mu_{j}' for j in j_indices_all],
        ordered=True,
    )
    hmc_repr_key_all = f'mu_{j_indices_all[-1]}'

    p = (
        ggplot(bar_all, aes(x='interim_date', y='pps', fill='fill_key'))
        + geom_col(position=position_dodge(width=20), width=18)
        + geom_errorbar(
            aes(x='interim_date',
                ymin='q025_bs', ymax='q975_bs', group='method'),
            position=position_dodge(width=20), width=8,
            colour='black', size=0.5, inherit_aes=False,
        )
        + scale_fill_manual(values=fill_values_all,
                            breaks=['analytic', hmc_repr_key_all],
                            labels=['analytic', 'HMC'],
                            name='method')
        + scale_x_datetime(date_breaks='2 months', date_labels='%Y-%b')
        + scale_y_continuous(
            breaks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            labels=['0%', '20%', '40%', '60%', '80%', '100%'],
        )
        + facet_wrap('~ response_label',
                     ncol=n_per_level, dir='h')
        + theme_bw()
        + theme(axis_text_x=element_text(angle=45, vjust=1, hjust=1),
                figure_size=(3.0 * n_per_level,
                             2.5 * K_levels),
                legend_position='top',
                strip_background=element_blank(),
                strip_text=element_text(face='bold'))
        + labs(x='interim date',
               y='PPS = int P(p(H_1j | x, z) > eta) dz')
    )
    pdf_path = os.path.join(DIR_OUT, f'mvn_J{J}_pps_bars_all.pdf')
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved 'all' PPS bar chart to {pdf_path}")

    # ---- "all" boxplot ----
    an_J_box_all = an_by_J.get(J, pd.DataFrame())
    if an_J_box_all.empty:
        continue
    an_J_box_all = an_J_box_all[an_J_box_all['j'].isin(j_indices_all)].copy()
    an_J_box_all['response_label'] = 'response mu_' + an_J_box_all['j'].astype(str)
    an_J_box_all['interim_month_year'] = pd.Categorical(
        an_J_box_all['interim_month_year'], categories=order, ordered=True,
    )
    an_J_box_all['response_label'] = pd.Categorical(
        an_J_box_all['response_label'],
        categories=response_label_cats, ordered=True,
    )

    nm_all_box = nm_all_J.copy()
    nm_all_box['interim_month_year'] = pd.Categorical(
        nm_all_box['interim_month_year'], categories=order, ordered=True,
    )

    def _q_table_all(df, method):
        q = (
            df.groupby(['response_label', 'interim_month_year', 'j'],
                       observed=True)['p_h1_xz']
              .agg(
                  q025=lambda s: s.quantile(0.025),
                  q25 =lambda s: s.quantile(0.25),
                  q50 =lambda s: s.quantile(0.50),
                  q75 =lambda s: s.quantile(0.75),
                  q975=lambda s: s.quantile(0.975),
              ).reset_index()
        )
        q['method'] = method
        return q

    box_all = pd.concat(
        [_q_table_all(an_J_box_all, 'analytic'),
         _q_table_all(nm_all_box, 'HMC')],
        ignore_index=True,
    )
    box_all['method'] = pd.Categorical(
        box_all['method'], categories=['analytic', 'HMC'], ordered=True,
    )
    box_all['fill_key'] = np.where(
        box_all['method'].astype(str) == 'analytic',
        'analytic',
        'mu_' + box_all['j'].astype(str),
    )
    box_all['fill_key'] = pd.Categorical(
        box_all['fill_key'],
        categories=['analytic'] + [f'mu_{j}' for j in j_indices_all],
        ordered=True,
    )
    box_all['response_label'] = pd.Categorical(
        box_all['response_label'],
        categories=response_label_cats, ordered=True,
    )
    box_all['grp'] = (
        box_all['interim_month_year'].astype(str) + '|'
        + box_all['fill_key'].astype(str)
    )
    box_all.to_pickle(os.path.join(DIR_OUT, f'mvn_J{J}_pps_p_h1_xz_all.pkl'))

    p = (
        ggplot(box_all, aes(x='interim_month_year',
                            ymin='q025', lower='q25', middle='q50',
                            upper='q75', ymax='q975',
                            fill='fill_key', group='grp'))
        + geom_boxplot(stat='identity',
                       position=position_dodge(width=0.8), width=0.7,
                       colour='#404040', size=0.3)
        + geom_hline(yintercept=simu_params['pps_ProbH1_thresh'],
                     colour='black', size=1.0)
        + scale_fill_manual(values=fill_values_all,
                            breaks=['analytic', hmc_repr_key_all],
                            labels=['analytic', 'HMC'],
                            name='method')
        + scale_y_continuous(
            breaks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            labels=['0%', '20%', '40%', '60%', '80%', '100%'],
        )
        + facet_wrap('~ response_label',
                     ncol=n_per_level, dir='h')
        + theme_bw()
        + theme(axis_text_x=element_text(angle=45, vjust=1, hjust=1),
                figure_size=(3.0 * n_per_level,
                             2.5 * K_levels),
                legend_position='top',
                strip_background=element_blank(),
                strip_text=element_text(face='bold'))
        + labs(x='Interim',
               y='p(H_1j | x, z) for predicted z samples')
    )
    pdf_path = os.path.join(DIR_OUT, f'mvn_J{J}_pps_p_h1_xz_all.pdf')
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved 'all' p(H_1 | x, z) boxplot to {pdf_path}")

    # ---- Scatter analytic vs HMC per-sample ----
    from ggsci import pal_simpsons
    n_interims_J = nm_all_J['interim_id'].nunique()
    simpsons_colours = pal_simpsons()(n_interims_J)
    scatter = an_J_box_all[[
        'interim_id', 'interim_month_year', 'j', 's', 'p_h1_xz',
    ]].rename(columns={'p_h1_xz': 'p_h1_xz_analytic'}).merge(
        nm_all_box[['interim_id', 'j', 's', 'p_h1_xz']].rename(
            columns={'p_h1_xz': 'p_h1_xz_hmc'},
        ),
        on=['interim_id', 'j', 's'], how='inner',
    )
    scatter['response_label'] = 'response mu_' + scatter['j'].astype(str)
    scatter['response_label'] = pd.Categorical(
        scatter['response_label'],
        categories=response_label_cats, ordered=True,
    )
    scatter['interim_month_year'] = pd.Categorical(
        scatter['interim_month_year'], categories=order, ordered=True,
    )
    scatter.to_pickle(
        os.path.join(DIR_OUT, f'mvn_J{J}_pps_p_h1_xz_scatter.pkl'),
    )

    p = (
        ggplot(scatter, aes(x='p_h1_xz_analytic', y='p_h1_xz_hmc',
                            color='interim_month_year'))
        + geom_abline(slope=1, intercept=0,
                      colour='#404040', linetype='dashed', size=0.4)
        + geom_point(alpha=0.5, size=0.6)
        + scale_color_manual(values=simpsons_colours, name='interim')
        + scale_x_continuous(
            breaks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            labels=['0%', '20%', '40%', '60%', '80%', '100%'],
        )
        + scale_y_continuous(
            breaks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            labels=['0%', '20%', '40%', '60%', '80%', '100%'],
        )
        + facet_wrap('~ response_label',
                     ncol=n_per_level, dir='h')
        + theme_bw()
        + theme(axis_text_x=element_text(angle=45, vjust=1, hjust=1),
                figure_size=(3.0 * n_per_level, 2.5 * K_levels),
                legend_position='top',
                strip_background=element_blank(),
                strip_text=element_text(face='bold'))
        + labs(x='p(H_1j | x, z) analytic',
               y='p(H_1j | x, z) HMC')
    )
    pdf_path = os.path.join(DIR_OUT, f'mvn_J{J}_pps_p_h1_xz_scatter.pdf')
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved analytic vs HMC scatter to {pdf_path}")
