#!/usr/bin/env python3
"""
Build the MVN interim simulation cohorts + every closed-form artifact and
diagnostic plot used by ``MVN_interim_analyses_with_nested_MC_hmc.py``
(§3.3 of ``dev/amortised_decision_making.md``).

Flat-cell layout. Running the file end-to-end (or stepping through it)
persists to ``DIR_OUT``:

  - ``mvn_sim_data.pkl``                (Phase 1)
  - ``mvn_diag_per_component.pdf``      (Phase 2a)
  - ``mvn_diag_cov_heatmap.pdf``        (Phase 2b)
  - ``mvn_diag_cov_bars.pdf``           (Phase 2c)
  - ``mvn_pps_closed_form.pkl``         (Phase 3)
  - ``mvn_p_h1_x.pdf``                  (Phase 4a)
  - ``mvn_pps_closed_form.pdf``         (Phase 4b)

The companion HMC script imports ``simu_params``, ``sim_data``, ``di`` and
``interim_order`` from this module; top-level cells run as a side-effect
of import.

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/MVN_interim_analyses_make_sim_data.py
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
from scipy.stats import norm as _norm
from plotnine import (
    ggplot, aes, geom_col, geom_errorbar, geom_hline, geom_point, geom_tile,
    facet_grid, facet_wrap,
    scale_fill_gradient2, scale_fill_manual,
    scale_x_continuous, scale_y_continuous, scale_y_discrete,
    position_dodge2, theme_bw, theme,
    element_text, element_blank, labs,
)

warnings.filterwarnings('ignore')

from model_mvn import MVNModel, make_block_equicorr_chol

print("Imports successful")

# %%

# =============================================================================
# Configuration. All simulation + closed-form parameters live in the
# ``simu_params`` dict; downstream code accesses them via
# ``simu_params['key']``.
# =============================================================================

simu_params = {
    'seed':              123,
    'J_grid':            [20, 60, 100],
    'K_levels':          4,
    'rho_w':             0.8,
    'rho_b':             0.1,
    'sigma':             1.0,
    'delta':             0.075,
    'mu_0_baseline':     1.0,
    'prior_tau':         10.0,
    'pps_H1_min_effect_size_thresh':        0.0,
    'pps_ProbH1_target_lwr_quantile': 0.89,
    'N_full':            1050,
    'N_interims':        7,
}

DIR_OUT = os.path.join("/Users/or105/sandbox/bIRTistic", "py-mvn-interim-simulations-260609")
os.makedirs(DIR_OUT, exist_ok=True)

J_LABEL_CATS = [f'J={J}' for J in sorted(simu_params['J_grid'])]

# %%

# =============================================================================
# Phase 1. Simulate one cohort per J. Persist everything in one pkl.
# =============================================================================

dates = pd.date_range(
    '2025-01-01', '2025-12-31', periods=simu_params['N_full'],
)
interim_dates = pd.date_range(
    '2025-01-01', '2025-12-31', periods=simu_params['N_interims'] + 2,
)[1:-1]
di = pd.DataFrame({
    'interim_id':         np.arange(1, simu_params['N_interims'] + 1, dtype=int),
    'interim_date':       interim_dates,
    'interim_month_year': pd.to_datetime(interim_dates).strftime('%Y-%b'),
})
interim_order = di['interim_month_year'].tolist()

# K_levels block-mean levels (above-baseline first, descending toward 1;
# below-baseline mirror), one block per level for every J.
n_up = simu_params['K_levels'] // 2
level_means = np.concatenate([
    simu_params['mu_0_baseline'] + np.arange(n_up, 0, -1) * simu_params['delta'],
    simu_params['mu_0_baseline']
    - np.arange(1, simu_params['K_levels'] - n_up + 1) * simu_params['delta'],
])

sim_data = {
    'simu_params':   simu_params,
    'level_means':   level_means,
    'interim_grid':  di,
    'interim_order': interim_order,
    'cells':         {},
}

for J in simu_params['J_grid']:
    if J % simu_params['K_levels'] != 0:
        raise ValueError(
            f"J={J} not divisible by K_levels={simu_params['K_levels']}."
        )
    block_size_J = J // simu_params['K_levels']
    R_chol = make_block_equicorr_chol(
        J, block_size_J, simu_params['rho_w'], simu_params['rho_b'],
    )
    mu_true = np.repeat(level_means, block_size_J)
    Sigma_true = simu_params['sigma'] ** 2 * (R_chol @ R_chol.T)

    rng_J = np.random.default_rng(simu_params['seed'] + J)
    Y = mu_true[None, :] + simu_params['sigma'] * (
        rng_J.standard_normal(size=(simu_params['N_full'], J)) @ R_chol.T
    )                                                    # (N, J)
    pid = np.repeat(np.arange(1, simu_params['N_full'] + 1), J)
    j = np.tile(np.arange(J), simu_params['N_full'])
    dp = pd.DataFrame({
        'pid':              pid,
        'j':                j,
        'submission_date':  np.repeat(dates, J),
        'y':                Y.reshape(-1),
        'oid':              np.arange(1, simu_params['N_full'] * J + 1,
                                      dtype=int),
    })
    dit = pd.DataFrame({
        'item_label':       [f'mu_{k}' for k in range(J)],
        'item_type':        ['mvn'] * J,
        'item_high_label':  ['higher_is_better'] * J,
    })
    sim_data['cells'][J] = {
        'J':           J,
        'R_chol':      R_chol,
        'mu_true':     mu_true,
        'Sigma_true':  Sigma_true,
        'Y':           Y,
        'dit':         dit,
        'dp':          dp,
    }
    print(f"  J={J}: simulated {simu_params['N_full']} obs in R^{J};"
          f" block-equicorr (K_levels={simu_params['K_levels']} blocks,"
          f" block_size={block_size_J}, rho_w={simu_params['rho_w']},"
          f" rho_b={simu_params['rho_b']});"
          f" true mu levels = {level_means.round(3).tolist()}")

sim_pkl = os.path.join(DIR_OUT, 'mvn_sim_data.pkl')
pd.to_pickle(sim_data, sim_pkl)
print(f"\nSaved sim data ({len(simu_params['J_grid'])} cells) to {sim_pkl}")

# %%

# =============================================================================
# Phase 2a. Per-component diagnostic per (J, interim): true mu_j (red x),
# empirical interim cohort mean (black o), closed-form posterior + 95% CrI
# (blue point + errorbar). Cols by J, rows by interim.
# =============================================================================

diag_a_parts = []
for J in simu_params['J_grid']:
    cell = sim_data['cells'][J]
    for i in range(len(di)):
        interim_date = di.iloc[i]['interim_date']
        dpi = MVNModel.get_interim_data_x(
            cell['dp'][cell['dp']['submission_date'] <= interim_date],
        )
        n_obs = int(dpi['pid'].nunique())
        if n_obs == 0:
            continue
        model = MVNModel(
            dit=cell['dit'], dcati=dpi, seed=simu_params['seed'],
            J=J, K_chol=cell['R_chol'], sigma=simu_params['sigma'],
            prior_tau=simu_params['prior_tau'],
            mu_0_baseline=simu_params['mu_0_baseline'],
        )
        mom = model._posterior_moments()
        sd = np.sqrt(np.diag(mom['cov_n']))
        y_bar = cell['Y'][:n_obs].mean(axis=0)
        diag_a_parts.append(pd.DataFrame({
            'J':                  J,
            'interim_id':         int(di.iloc[i]['interim_id']),
            'interim_date':       interim_date,
            'interim_month_year': di.iloc[i]['interim_month_year'],
            'n_obs':              n_obs,
            'j':                  np.arange(J),
            'mu_true':            cell['mu_true'],
            'y_bar':              y_bar,
            'mu_post':            mom['mu_n'],
            'crI_lo':             mom['mu_n'] - 1.96 * sd,
            'crI_hi':             mom['mu_n'] + 1.96 * sd,
        }))
diag_a = pd.concat(diag_a_parts, ignore_index=True)
diag_a['J_label'] = pd.Categorical(
    'J=' + diag_a['J'].astype(str),
    categories=J_LABEL_CATS, ordered=True,
)
diag_a['interim_month_year'] = pd.Categorical(
    diag_a['interim_month_year'],
    categories=interim_order, ordered=True,
)
diag_a.to_pickle(os.path.join(DIR_OUT, 'mvn_diag_per_component.pkl'))

p_a = (
    ggplot(diag_a, aes(x='j'))
    + geom_errorbar(aes(ymin='crI_lo', ymax='crI_hi'),
                    colour='#1f77b4', width=0.4, size=0.4)
    + geom_point(aes(y='mu_post'), colour='#1f77b4', size=1.4)
    + geom_point(aes(y='y_bar'), colour='black', fill='black', size=1.0)
    + geom_point(aes(y='mu_true'), colour='red', shape='x', size=2.0)
    + geom_hline(yintercept=simu_params['mu_0_baseline'], linetype='dashed',
                 colour='#808080', size=0.3)
    + facet_grid('interim_month_year ~ J_label', scales='free_x')
    + theme_bw()
    + theme(figure_size=(5 * len(simu_params['J_grid']),
                         2.4 * simu_params['N_interims']),
            axis_text_x=element_text(angle=0),
            strip_background=element_blank(),
            strip_text=element_text(face='bold'))
    + labs(x='component j', y='effect size')
)
p_a.save(os.path.join(DIR_OUT, 'mvn_diag_per_component.pdf'),
         verbose=False, limitsize=False)
print(f"Saved per-component diagnostic to "
      f"{os.path.join(DIR_OUT, 'mvn_diag_per_component.pdf')}")

# %%

# =============================================================================
# Phase 2b. True Sigma vs empirical covariance heatmap, faceted by J + method.
# =============================================================================

emp_cov_by_J = {J: np.cov(sim_data['cells'][J]['Y'].T, ddof=1)
                for J in simu_params['J_grid']}

heat_parts = []
for J in simu_params['J_grid']:
    for label, M in (('true',      sim_data['cells'][J]['Sigma_true']),
                     ('empirical', emp_cov_by_J[J])):
        ii, jj = np.meshgrid(np.arange(J), np.arange(J), indexing='ij')
        heat_parts.append(pd.DataFrame({
            'J':      J,
            'i':      ii.reshape(-1),
            'j':      jj.reshape(-1),
            'value':  M.reshape(-1),
            'method': label,
        }))
heat = pd.concat(heat_parts, ignore_index=True)
heat['J_label'] = pd.Categorical(
    'J=' + heat['J'].astype(str),
    categories=J_LABEL_CATS, ordered=True,
)
heat['method'] = pd.Categorical(
    heat['method'], categories=['true', 'empirical'], ordered=True,
)
heat.to_pickle(os.path.join(DIR_OUT, 'mvn_diag_cov_heatmap.pkl'))

p_b = (
    ggplot(heat, aes(x='j', y='i', fill='value'))
    + geom_tile()
    + facet_grid('J_label ~ method', scales='free', space='free')
    + scale_fill_gradient2(low='#3b4cc0', mid='white', high='#b40426',
                           midpoint=0.0)
    + scale_x_continuous(expand=(0, 0))
    + scale_y_continuous(expand=(0, 0))
    + theme_bw()
    + theme(figure_size=(10, 18),
            axis_text_x=element_text(angle=0),
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            strip_background=element_blank(),
            strip_text=element_text(face='bold'))
    + labs(x='component j', y='component i', fill='cov(y_i, y_j)')
)
p_b.save(os.path.join(DIR_OUT, 'mvn_diag_cov_heatmap.pdf'),
         verbose=False, limitsize=False)
print(f"Saved covariance heatmap to "
      f"{os.path.join(DIR_OUT, 'mvn_diag_cov_heatmap.pdf')}")

# %%

# =============================================================================
# Phase 2c. Element-wise covariance barplot (true vs empirical), dodge2,
# upper triangle including diagonal so each unique entry appears once.
# =============================================================================

bar_parts = []
for J in simu_params['J_grid']:
    iu, ju = np.triu_indices(J)
    pair_idx = np.arange(len(iu))
    for label, M in (('true',      sim_data['cells'][J]['Sigma_true']),
                     ('empirical', emp_cov_by_J[J])):
        bar_parts.append(pd.DataFrame({
            'J':         J,
            'pair_idx':  pair_idx,
            'value':     M[iu, ju],
            'method':    label,
        }))
bar = pd.concat(bar_parts, ignore_index=True)
bar['J_label'] = pd.Categorical(
    'J=' + bar['J'].astype(str),
    categories=J_LABEL_CATS, ordered=True,
)
bar['method'] = pd.Categorical(
    bar['method'], categories=['true', 'empirical'], ordered=True,
)
bar.to_pickle(os.path.join(DIR_OUT, 'mvn_diag_cov_bars.pkl'))

p_c = (
    ggplot(bar, aes(x='pair_idx', y='value', fill='method'))
    + geom_col(position=position_dodge2(width=0.9, preserve='single'),
               width=0.8)
    + facet_wrap('~ J_label', ncol=1, scales='free')
    + scale_fill_manual(values={'true': '#000000', 'empirical': '#d62728'})
    + theme_bw()
    + theme(figure_size=(14, 16),
            axis_text_x=element_blank(),
            legend_position='top',
            strip_background=element_blank(),
            strip_text=element_text(face='bold'))
    + labs(x='upper-triangle covariance entry (i, j)',
           y='cov(y_i, y_j)', fill='method')
)
p_c.save(os.path.join(DIR_OUT, 'mvn_diag_cov_bars.pdf'),
         verbose=False, limitsize=False)
print(f"Saved per-element covariance bar plot to "
      f"{os.path.join(DIR_OUT, 'mvn_diag_cov_bars.pdf')}")

# %%

# =============================================================================
# Phase 3. Closed-form posterior P(H_1j | x) + closed-form PPS per
# (J, interim). One pkl across all J in J_grid.
# =============================================================================

post_parts = []
pps_cf_parts = []
for J in simu_params['J_grid']:
    cell = sim_data['cells'][J]
    for i in range(len(di)):
        interim_id = int(di.iloc[i]['interim_id'])
        interim_date = di.iloc[i]['interim_date']
        interim_month_year = di.iloc[i]['interim_month_year']
        dpi = MVNModel.get_interim_data_x(
            cell['dp'][cell['dp']['submission_date'] <= interim_date],
        )
        n_obs = int(dpi['pid'].nunique())
        if n_obs == 0:
            continue
        m_future = simu_params['N_full'] - n_obs
        model = MVNModel(
            dit=cell['dit'], dcati=dpi, seed=simu_params['seed'],
            J=J, K_chol=cell['R_chol'], sigma=simu_params['sigma'],
            prior_tau=simu_params['prior_tau'],
            mu_0_baseline=simu_params['mu_0_baseline'],
        )
        mom = model._posterior_moments()
        sd = np.sqrt(np.diag(mom['cov_n']))
        post_parts.append(pd.DataFrame({
            'J':                  J,
            'interim_id':         interim_id,
            'interim_date':       interim_date,
            'interim_month_year': interim_month_year,
            'n_obs':              n_obs,
            'm_future':           m_future,
            'j':                  np.arange(J),
            'mu_n':               mom['mu_n'],
            'sd_n':               sd,
            'p_h1_x':             1.0 - _norm.cdf(
                                     (simu_params['mu_0_baseline']
                                      - mom['mu_n']) / sd),
            'mu_true':            cell['mu_true'],
        }))
        if m_future > 0:
            cf = model.fit_closed_form_pps(
                m=m_future,
                pps_H1_min_effect_size_thresh=simu_params['pps_H1_min_effect_size_thresh'],
                pps_ProbH1_target_lwr_quantile=simu_params['pps_ProbH1_target_lwr_quantile'],
            ).assign(
                J=J,
                interim_id=interim_id,
                interim_date=interim_date,
                interim_month_year=interim_month_year,
                n_obs=n_obs,
                mu_true=cell['mu_true'],
            )
            pps_cf_parts.append(cf)

post = pd.concat(post_parts, ignore_index=True)
pps_cf = pd.concat(pps_cf_parts, ignore_index=True)

for df in (post, pps_cf):
    df['interim_month_year'] = pd.Categorical(
        df['interim_month_year'], categories=interim_order, ordered=True,
    )
    df['J_label'] = pd.Categorical(
        'J=' + df['J'].astype(str),
        categories=J_LABEL_CATS, ordered=True,
    )

cf_pkl = os.path.join(DIR_OUT, 'mvn_pps_closed_form.pkl')
pd.to_pickle({'p_h1_x': post, 'pps_cf': pps_cf}, cf_pkl)
print(f"Saved closed-form posterior + PPS to {cf_pkl}: "
      f"{len(post)} P(H_1|x) rows, {len(pps_cf)} PPS rows.")

# %%

# =============================================================================
# Phase 4a. P(H_1j | x) heatmap: x=j (shared scale, identical cell width),
# y=interim, fill=P(H_1j | x). TealRose divergingx midpoint 0.5.
# =============================================================================

p_4a = (
    ggplot(post, aes(x='j', y='interim_month_year', fill='p_h1_x'))
    + geom_tile()
    + facet_wrap('~ J_label', ncol=1)
    + scale_x_continuous(expand=(0, 0))
    + scale_y_discrete(expand=(0, 0))
    + scale_fill_gradient2(low='#009392', mid='#F1EAC8', high='#A5006A',
                           midpoint=0.5, limits=[0.0, 1.0])
    + theme_bw()
    + theme(figure_size=(14, 3 * len(simu_params['J_grid'])),
            axis_text_y=element_text(size=8),
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            strip_background=element_blank(),
            strip_text=element_text(face='bold'))
    + labs(x='component j', y='interim', fill='P(H_1j | x)')
)
p_4a.save(os.path.join(DIR_OUT, 'mvn_p_h1_x.pdf'),
          verbose=False, limitsize=False)
print(f"Saved P(H_1j | x) heatmap to "
      f"{os.path.join(DIR_OUT, 'mvn_p_h1_x.pdf')}")

# %%

# =============================================================================
# Phase 4b. PPS heatmap, same layout as 4a.
# =============================================================================

p_4b = (
    ggplot(pps_cf, aes(x='j', y='interim_month_year', fill='pps'))
    + geom_tile()
    + facet_wrap('~ J_label', ncol=1)
    + scale_x_continuous(expand=(0, 0))
    + scale_y_discrete(expand=(0, 0))
    + scale_fill_gradient2(low='#009392', mid='#F1EAC8', high='#A5006A',
                           midpoint=0.5, limits=[0.0, 1.0])
    + theme_bw()
    + theme(figure_size=(14, 3 * len(simu_params['J_grid'])),
            axis_text_y=element_text(size=8),
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            strip_background=element_blank(),
            strip_text=element_text(face='bold'))
    + labs(x='component j', y='interim', fill='PPS')
)
p_4b.save(os.path.join(DIR_OUT, 'mvn_pps_closed_form.pdf'),
          verbose=False, limitsize=False)
print(f"Saved closed-form PPS heatmap to "
      f"{os.path.join(DIR_OUT, 'mvn_pps_closed_form.pdf')}")
