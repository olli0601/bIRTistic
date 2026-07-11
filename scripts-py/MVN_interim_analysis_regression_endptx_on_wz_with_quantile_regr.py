#!/usr/bin/env python3
"""
MVN interim analysis: regression-based EVSI labels (Strong & Oakley,
2014), continuous-endpoint variant with **conditional-quantile regression**
instead of the Gaussian mean + variance approximation used by
``MVN_interim_analysis_regression_endptx_on_wz_with_gauss_approx.py``.

Per interim, per response j, fit

    Q_{1 - eta_H}(mu_j - mu_0 | x, w_j(z))  ~  pinball loss on S rows

with a linear ``statsmodels.QuantReg`` at tau = 1 - eta_H, so
``P(H_1j | x, z) > eta_H`` is equivalent to ``Q_hat > pps_H1_min_effect_size_thresh`` — no
Gaussian step, no plug-in variance.

Outputs per J (``_RGEQ_`` suffix so compare-methods can load them next to
the Gaussian ``_RGE_`` variant) plus the same set of plots as the Gaussian
variant, with the diagnostic scatter drawing a QuantReg line at
tau = 1 - eta_H rather than a Gaussian GLM smoother.

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/MVN_interim_analysis_regression_endptx_on_wz_with_quantile_regr.py
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
import xarray as xr
import arviz as az
from plotnine import (
    ggplot, aes, geom_abline, geom_boxplot, geom_col, geom_errorbar,
    geom_hline, geom_point, geom_smooth, geom_text,
    facet_grid, facet_wrap,
    scale_color_manual, scale_fill_manual,
    scale_x_continuous, scale_x_datetime, scale_y_continuous,
    position_dodge, theme_bw, theme,
    element_text, element_blank, labs,
)

warnings.filterwarnings('ignore')

import statsmodels.api as sm
from model_mvn import MVNModel
from fit_interim import (
    fit_interim_regress_endptx_on_wz_with_quantile_regr,
    fit_interim_regress_H1x_on_wz_per_item_summary,
)
from utils import _material_palette

print("Imports successful")

# %%

# =============================================================================
# Configuration. ``DIR_SIM`` provides the upstream artifacts;
# ``DIR_OUT`` is the dedicated regression-endptx output directory.
# =============================================================================

DIR_SIM = os.path.join(
    "/Users/or105/sandbox/bIRTistic", "py-mvn-interim-simulations-260609",
)
DIR_OUT = os.path.join(
    "/Users/or105/sandbox/bIRTistic",
    "py-mvn-interim-with-regression-on-endptx-wz-quantile-regr-260702",
)
os.makedirs(DIR_OUT, exist_ok=True)

PPS_Z_TOTAL = 4000

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

K_levels  = simu_params['K_levels']
seed      = simu_params['seed']
pps_H1_min_effect_size_thresh        = simu_params['pps_H1_min_effect_size_thresh']
pps_ProbH1_target_lwr_quantile = simu_params['pps_ProbH1_target_lwr_quantile']

# %%

# =============================================================================
# Phase 5. Per (J, interim) regression compute + per-(j, s) analytic
# label for downstream comparison plots.
# =============================================================================

rge_artifacts = {}
for J in J_GRID:
    cell = sim_data['cells'][J]
    dit, R_chol = cell['dit'], cell['R_chol']

    print(f"\n{'#' * 70}\n# J = {J} (n_items = {J})\n{'#' * 70}")

    p_h1_xz_rows = []
    perf_rows = []
    timing_rows = []
    wa_rows = []
    an_rows = []
    t_J0 = time.time()
    for interim_id, blk in interim_data_by_J[J].items():
        t0 = time.time()
        dpi = blk['dpi']
        zi_full = blk['zi']
        mu_draws = blk['mu_draws']                  # (S_cache, J)
        interim_date = blk['interim_date']
        interim_month_year = blk['interim_month_year']
        interim_m = blk['interim_m']
        n_obs = blk['n_obs']

        ypred_cols = sorted(
            [c for c in zi_full.columns if c.startswith('ypred_')],
            key=lambda c: int(c.split('_')[1]),
        )
        if PPS_Z_TOTAL > len(ypred_cols):
            raise ValueError(
                f"PPS_Z_TOTAL={PPS_Z_TOTAL} > cached S={len(ypred_cols)}."
            )
        keep_cols = ypred_cols[:PPS_Z_TOTAL]
        meta_cols = [c for c in zi_full.columns if not c.startswith('ypred_')]
        zi = zi_full[meta_cols + keep_cols].copy()
        mu_S = mu_draws[:PPS_Z_TOTAL, :]            # (S, J)

        print(f"\n{'=' * 70}\nJ={J} RGE interim {interim_id}"
              f" ({interim_date.date()}) n={n_obs} m={interim_m}"
              f" S={PPS_Z_TOTAL}\n{'=' * 70}")

        model_x = MVNModel(
            dit=dit, dcati=dpi, seed=seed,
            J=J, K_chol=R_chol, sigma=simu_params['sigma'],
            prior_tau=simu_params['prior_tau'],
            mu_0_baseline=simu_params['mu_0_baseline'],
        )

        # Build an arviz idata from the cached mu draws so
        # ``get_endpoints_per_draw`` can be reused unchanged.
        fake_idata = az.InferenceData(posterior=xr.Dataset({
            'mu': (('chain', 'draw', 'j'), mu_S[None, :, :]),
        }))
        wa = model_x.get_w(zi)
        x_ratio = model_x.get_endpoints_per_draw(draws=fake_idata)
        x_ratio = x_ratio.rename(columns={'ratio': 'pps_ratio_x'})
        x_ratio['pps_H1_x'] = (
            x_ratio['pps_ratio_x'] > pps_H1_min_effect_size_thresh
        ).astype(int)
        wa = wa.merge(
            x_ratio[['draw', 'item_label', 'item_type',
                     'pps_ratio_x', 'pps_H1_x']],
            on=['draw', 'item_label', 'item_type'], how='inner',
        )
        wa['J'] = J
        wa['interim_id'] = interim_id
        wa['interim_month_year'] = interim_month_year
        wa['interim_date'] = interim_date
        wa_rows.append(wa)

        p_h1_xz_interim, perf_interim = fit_interim_regress_endptx_on_wz_with_quantile_regr(
            wa.drop(columns=['J', 'interim_id', 'interim_month_year',
                             'interim_date']),
            pps_H1_min_effect_size_thresh=pps_H1_min_effect_size_thresh,
            pps_ProbH1_target_lwr_quantile=pps_ProbH1_target_lwr_quantile,
        )
        p_h1_xz_interim['J'] = J
        p_h1_xz_interim['interim_id'] = interim_id
        p_h1_xz_interim['interim_date'] = interim_date
        p_h1_xz_interim['interim_month_year'] = interim_month_year
        p_h1_xz_rows.append(p_h1_xz_interim)

        perf_interim['J'] = J
        perf_interim['interim_id'] = interim_id
        perf_interim['interim_date'] = interim_date
        perf_interim['interim_month_year'] = interim_month_year
        perf_rows.append(perf_interim)

        # Analytic per-(s, j) for comparison plots: reuse model_x.
        an = model_x.fit_closed_form_pH1(zi)
        an['J'] = J
        an['interim_id'] = interim_id
        an['interim_date'] = interim_date
        an['interim_month_year'] = interim_month_year
        an_rows.append(an)

        mins = (time.time() - t0) / 60.0
        timing_rows.append({
            'J':                  J,
            'interim_id':         interim_id,
            'interim_date':       interim_date,
            'interim_month_year': interim_month_year,
            'n_obs':              n_obs,
            'm_future':           interim_m,
            'mins_interim_id':    round(mins, 3),
        })
        print(f"  J={J} interim {interim_id} regression-endptx done in"
              f" {mins:.2f} min")

    if not p_h1_xz_rows:
        print(f"J={J}: no interims with m > 0; skipping save.")
        continue

    dp_h1_xz = pd.concat(p_h1_xz_rows, ignore_index=True)
    dp_h1_xz['pps_H1_min_effect_size_thresh'] = pps_H1_min_effect_size_thresh
    dp_h1_xz['pps_ProbH1_target_lwr_quantile'] = pps_ProbH1_target_lwr_quantile
    dp_h1_xz['S'] = PPS_Z_TOTAL
    dp_h1_xz['j'] = dp_h1_xz['item_label'].str.replace('mu_', '').astype(int)
    perf_all = pd.concat(perf_rows, ignore_index=True)
    perf_all['j'] = perf_all['item_label'].str.replace('mu_', '').astype(int)
    timing_all = pd.DataFrame(timing_rows)
    timing_all['mins_J_total'] = round((time.time() - t_J0) / 60.0, 3)
    wa_all = pd.concat(wa_rows, ignore_index=True)
    wa_all['j'] = wa_all['item_label'].str.replace('mu_', '').astype(int)
    an_all = pd.concat(an_rows, ignore_index=True)

    pps_df = (
        dp_h1_xz.groupby(['J', 'interim_id', 'interim_date',
                          'interim_month_year', 'item_label', 'item_type',
                          'item_high_label', 'j'], observed=True)['p_h1_xz']
        .apply(lambda p: float(
            (p > pps_ProbH1_target_lwr_quantile).mean()
        ))
        .reset_index(name='pps')
    )
    pps_df['eta'] = pps_ProbH1_target_lwr_quantile
    pps_df['S'] = PPS_Z_TOTAL

    dp_h1_xz.to_pickle(os.path.join(DIR_OUT, f'mvn_J{J}_pps_RGEQ_p_h1_xz.pkl'))
    perf_all.to_csv(os.path.join(DIR_OUT, f'mvn_J{J}_pps_RGEQ_perf.csv'),
                    index=False)
    timing_all.to_csv(os.path.join(DIR_OUT, f'mvn_J{J}_pps_RGEQ_timing.csv'),
                      index=False)
    pps_df.to_csv(os.path.join(DIR_OUT, f'mvn_J{J}_pps_RGEQ.csv'),
                  index=False)
    wa_all.to_pickle(os.path.join(DIR_OUT, f'mvn_J{J}_pps_RGEQ_wa.pkl'))
    an_all.to_pickle(os.path.join(DIR_OUT,
                                  f'mvn_J{J}_pps_RGEQ_p_h1_xz_analytic.pkl'))
    rge_artifacts[J] = {
        'dp_h1_xz': dp_h1_xz, 'perf': perf_all,
        'timing':   timing_all, 'pps': pps_df, 'wa': wa_all,
        'an':       an_all,
    }
    print(f"J={J}: saved RGE artifacts."
          f" Total time {(time.time() - t_J0) / 60.0:.2f} min."
          f" {len(dp_h1_xz)} p_h1_xz rows; PPS table has {len(pps_df)} rows.")

# %%

# =============================================================================
# Phase 6. Per-J plots: PPS bars / boxplot / scatter / "all" variants /
# Binomial-style diagnostic.
# =============================================================================

bs_B = 2000
METHOD_LABEL = 'Regression of endpt-x on w(z) - quantile'

for J in J_GRID:
    if J not in rge_artifacts:
        continue
    art = rge_artifacts[J]
    dp_h1_xz = art['dp_h1_xz']
    pps_df   = art['pps']
    perf_all = art['perf']
    wa_all   = art['wa']
    an_all   = art['an']

    block_size_J = J // K_levels
    j_repr = [block_size_J // 2 + block_size_J * k for k in range(K_levels)]
    n_per_level = min(5, block_size_J)
    j_indices_all = np.array([
        k * block_size_J + c for k in range(K_levels)
        for c in range(n_per_level)
    ])

    order = (
        di[di['interim_id'].isin(dp_h1_xz['interim_id'].unique())]
        .sort_values('interim_date')['interim_month_year'].tolist()
    )

    # Per-J "small" slices.
    rge_small      = dp_h1_xz[dp_h1_xz['j'].isin(j_repr)].copy()
    cf_small       = pps_cf[(pps_cf['J'] == J) & (pps_cf['j'].isin(j_repr))].copy()
    rge_pps_small  = pps_df[pps_df['j'].isin(j_repr)].copy()
    an_small       = an_all[an_all['j'].isin(j_repr)].copy()

    j_label_cats = [f'response mu_{j}' for j in j_repr]
    for df in (rge_small, cf_small, rge_pps_small, an_small):
        df['j_label'] = 'response mu_' + df['j'].astype(str)
        df['interim_month_year'] = pd.Categorical(
            df['interim_month_year'], categories=order, ordered=True,
        )
        df['j_label'] = pd.Categorical(
            df['j_label'], categories=j_label_cats, ordered=True,
        )

    hmc_palette = _material_palette('teal', len(j_repr))
    fill_values = {'analytic': '#000000'}
    fill_values.update({f'mu_{j}': c for j, c in zip(j_repr, hmc_palette)})
    method_repr_key = f'mu_{j_repr[-1]}'

    # ---- 6a. Boxplot of p(H_1j | x, z): analytic vs RGE dodged ----
    def _qtab(df, method):
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
        [_qtab(an_small, 'analytic'), _qtab(rge_small, METHOD_LABEL)],
        ignore_index=True,
    )
    box['method'] = pd.Categorical(
        box['method'], categories=['analytic', METHOD_LABEL], ordered=True,
    )
    box['fill_key'] = np.where(
        box['method'].astype(str) == 'analytic',
        'analytic', 'mu_' + box['j'].astype(str),
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
    box.to_pickle(os.path.join(DIR_OUT, f'mvn_J{J}_pps_RGEQ_p_h1_xz_box.pkl'))

    p = (
        ggplot(box, aes(x='interim_month_year',
                        ymin='q025', lower='q25', middle='q50',
                        upper='q75', ymax='q975',
                        fill='fill_key', group='grp'))
        + geom_boxplot(stat='identity',
                       position=position_dodge(width=0.8), width=0.7,
                       colour='#404040', size=0.3)
        + geom_hline(yintercept=pps_ProbH1_target_lwr_quantile,
                     colour='black', size=1.0)
        + scale_fill_manual(values=fill_values,
                            breaks=['analytic', method_repr_key],
                            labels=['analytic', METHOD_LABEL],
                            name='method')
        + scale_y_continuous(
            breaks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            labels=['0%', '20%', '40%', '60%', '80%', '100%'],
        )
        + facet_wrap('~ j_label', ncol=K_levels)
        + theme_bw()
        + theme(axis_text_x=element_text(angle=45, vjust=1, hjust=1),
                figure_size=(14, 5), legend_position='top',
                strip_background=element_blank(),
                strip_text=element_text(face='bold'))
        + labs(x='Interim',
               y='p(H_1j | x, z) for predicted z samples')
    )
    pdf_path = os.path.join(DIR_OUT, f'mvn_J{J}_pps_RGEQ_p_h1_xz.pdf')
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved RGE p(H_1 | x, z) boxplot to {pdf_path}")

    # ---- 6b. PPS bars: analytic + RGE dodged + bootstrap CI ----
    wide = (
        rge_small.pivot_table(index=['j', 'j_label', 'interim_id',
                                     'interim_date'],
                              columns='s', values='p_h1_xz')
        .sort_index()
    )
    bs = wide.to_numpy()
    n_rows, S = bs.shape
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, S, size=(n_rows, bs_B, S))
    bs = bs[np.arange(n_rows)[:, None, None], idx]
    bs = (bs > pps_ProbH1_target_lwr_quantile).mean(axis=2)
    ci = np.quantile(bs, [0.025, 0.975], axis=1).T
    ci_df = pd.DataFrame(ci, columns=['q025_bs', 'q975_bs'])
    ci_df['j']          = wide.index.get_level_values('j').to_numpy()
    ci_df['interim_id'] = wide.index.get_level_values('interim_id').to_numpy()

    rge_bar = rge_pps_small[['interim_id', 'interim_date',
                             'interim_month_year', 'j', 'j_label',
                             'pps']].assign(method=METHOD_LABEL).merge(
        ci_df, on=['j', 'interim_id'], how='left',
    )
    cf_bar = cf_small[['interim_id', 'interim_date',
                       'interim_month_year', 'j', 'pps']].copy()
    cf_bar['j_label'] = pd.Categorical(
        'response mu_' + cf_bar['j'].astype(str),
        categories=j_label_cats, ordered=True,
    )
    cf_bar = cf_bar.assign(method='analytic',
                           q025_bs=np.nan, q975_bs=np.nan)
    bar = pd.concat([cf_bar, rge_bar], ignore_index=True)
    bar['interim_date'] = pd.to_datetime(bar['interim_date'])
    bar['method'] = pd.Categorical(
        bar['method'], categories=['analytic', METHOD_LABEL], ordered=True,
    )
    bar['fill_key'] = np.where(
        bar['method'].astype(str) == 'analytic',
        'analytic', 'mu_' + bar['j'].astype(str),
    )
    bar['fill_key'] = pd.Categorical(
        bar['fill_key'],
        categories=['analytic'] + [f'mu_{j}' for j in j_repr],
        ordered=True,
    )
    bar.to_pickle(os.path.join(DIR_OUT, f'mvn_J{J}_pps_RGEQ_bars.pkl'))

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
                            breaks=['analytic', method_repr_key],
                            labels=['analytic', METHOD_LABEL],
                            name='method')
        + scale_x_datetime(date_breaks='1 month', date_labels='%Y-%b')
        + scale_y_continuous(
            breaks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            labels=['0%', '20%', '40%', '60%', '80%', '100%'],
        )
        + facet_wrap('~ j_label', ncol=K_levels)
        + theme_bw()
        + theme(axis_text_x=element_text(angle=45, vjust=1, hjust=1),
                figure_size=(14, 5), legend_position='top',
                strip_background=element_blank(),
                strip_text=element_text(face='bold'))
        + labs(x='interim date',
               y='PPS = int P(p(H_1j | x, z) > eta) dz')
    )
    pdf_path = os.path.join(DIR_OUT, f'mvn_J{J}_pps_RGEQ_bars.pdf')
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved RGE PPS bar chart to {pdf_path}")

    # ---- 6c, 6d. "all" variants ----
    rge_all_pps = pps_df[pps_df['j'].isin(j_indices_all)].copy()
    cf_all  = pps_cf[(pps_cf['J'] == J) & (pps_cf['j'].isin(j_indices_all))].copy()
    rge_all = dp_h1_xz[dp_h1_xz['j'].isin(j_indices_all)].copy()
    an_all_J = an_all[an_all['j'].isin(j_indices_all)].copy()
    response_label_cats = [f'response mu_{j}' for j in j_indices_all]
    for df in (rge_all_pps, cf_all, rge_all, an_all_J):
        df['response_label'] = 'response mu_' + df['j'].astype(str)
        df['interim_month_year'] = pd.Categorical(
            df['interim_month_year'], categories=order, ordered=True,
        )
        df['response_label'] = pd.Categorical(
            df['response_label'], categories=response_label_cats, ordered=True,
        )

    hmc_palette_all = _material_palette('teal', len(j_indices_all))
    fill_values_all = {'analytic': '#000000'}
    fill_values_all.update({
        f'mu_{j}': c for j, c in zip(j_indices_all, hmc_palette_all)
    })
    method_repr_key_all = f'mu_{j_indices_all[-1]}'

    # bars all
    wide_all = (
        rge_all.pivot_table(index=['j', 'interim_id', 'interim_date'],
                            columns='s', values='p_h1_xz')
        .sort_index()
    )
    bs_a = wide_all.to_numpy()
    n_rows_a, S_a = bs_a.shape
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, S_a, size=(n_rows_a, bs_B, S_a))
    bs_a = bs_a[np.arange(n_rows_a)[:, None, None], idx]
    bs_a = (bs_a > pps_ProbH1_target_lwr_quantile).mean(axis=2)
    ci_a = np.quantile(bs_a, [0.025, 0.975], axis=1).T
    ci_df_a = pd.DataFrame(ci_a, columns=['q025_bs', 'q975_bs'])
    ci_df_a['j']          = wide_all.index.get_level_values('j').to_numpy()
    ci_df_a['interim_id'] = wide_all.index.get_level_values('interim_id').to_numpy()

    rge_bar_all = rge_all_pps[['interim_id', 'interim_date',
                               'interim_month_year', 'j', 'response_label',
                               'pps']].assign(method=METHOD_LABEL).merge(
        ci_df_a, on=['j', 'interim_id'], how='left',
    )
    cf_bar_all = cf_all[['interim_id', 'interim_date',
                         'interim_month_year', 'j', 'response_label',
                         'pps']].assign(method='analytic',
                                        q025_bs=np.nan, q975_bs=np.nan)
    bar_all = pd.concat([cf_bar_all, rge_bar_all], ignore_index=True)
    bar_all['interim_date'] = pd.to_datetime(bar_all['interim_date'])
    bar_all['method'] = pd.Categorical(
        bar_all['method'], categories=['analytic', METHOD_LABEL], ordered=True,
    )
    bar_all['fill_key'] = np.where(
        bar_all['method'].astype(str) == 'analytic',
        'analytic', 'mu_' + bar_all['j'].astype(str),
    )
    bar_all['fill_key'] = pd.Categorical(
        bar_all['fill_key'],
        categories=['analytic'] + [f'mu_{j}' for j in j_indices_all],
        ordered=True,
    )
    bar_all['response_label'] = pd.Categorical(
        bar_all['response_label'], categories=response_label_cats, ordered=True,
    )

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
                            breaks=['analytic', method_repr_key_all],
                            labels=['analytic', METHOD_LABEL],
                            name='method')
        + scale_x_datetime(date_breaks='2 months', date_labels='%Y-%b')
        + scale_y_continuous(
            breaks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            labels=['0%', '20%', '40%', '60%', '80%', '100%'],
        )
        + facet_wrap('~ response_label', ncol=n_per_level, dir='h')
        + theme_bw()
        + theme(axis_text_x=element_text(angle=45, vjust=1, hjust=1),
                figure_size=(3.0 * n_per_level, 2.5 * K_levels),
                legend_position='top',
                strip_background=element_blank(),
                strip_text=element_text(face='bold'))
        + labs(x='interim date',
               y='PPS = int P(p(H_1j | x, z) > eta) dz')
    )
    pdf_path = os.path.join(DIR_OUT, f'mvn_J{J}_pps_RGEQ_bars_all.pdf')
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved RGE 'all' PPS bar chart to {pdf_path}")

    # boxplot all
    box_all = pd.concat(
        [_qtab(an_all_J.rename(columns={'response_label': 'j_label'}),
               'analytic'),
         _qtab(rge_all.rename(columns={'response_label': 'j_label'}),
               METHOD_LABEL)],
        ignore_index=True,
    ).rename(columns={'j_label': 'response_label'})
    box_all['method'] = pd.Categorical(
        box_all['method'], categories=['analytic', METHOD_LABEL], ordered=True,
    )
    box_all['fill_key'] = np.where(
        box_all['method'].astype(str) == 'analytic',
        'analytic', 'mu_' + box_all['j'].astype(str),
    )
    box_all['fill_key'] = pd.Categorical(
        box_all['fill_key'],
        categories=['analytic'] + [f'mu_{j}' for j in j_indices_all],
        ordered=True,
    )
    box_all['response_label'] = pd.Categorical(
        box_all['response_label'], categories=response_label_cats, ordered=True,
    )
    box_all['grp'] = (
        box_all['interim_month_year'].astype(str) + '|'
        + box_all['fill_key'].astype(str)
    )

    p = (
        ggplot(box_all, aes(x='interim_month_year',
                            ymin='q025', lower='q25', middle='q50',
                            upper='q75', ymax='q975',
                            fill='fill_key', group='grp'))
        + geom_boxplot(stat='identity',
                       position=position_dodge(width=0.8), width=0.7,
                       colour='#404040', size=0.3)
        + geom_hline(yintercept=pps_ProbH1_target_lwr_quantile,
                     colour='black', size=1.0)
        + scale_fill_manual(values=fill_values_all,
                            breaks=['analytic', method_repr_key_all],
                            labels=['analytic', METHOD_LABEL],
                            name='method')
        + scale_y_continuous(
            breaks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            labels=['0%', '20%', '40%', '60%', '80%', '100%'],
        )
        + facet_wrap('~ response_label', ncol=n_per_level, dir='h')
        + theme_bw()
        + theme(axis_text_x=element_text(angle=45, vjust=1, hjust=1),
                figure_size=(3.0 * n_per_level, 2.5 * K_levels),
                legend_position='top',
                strip_background=element_blank(),
                strip_text=element_text(face='bold'))
        + labs(x='Interim',
               y='p(H_1j | x, z) for predicted z samples')
    )
    pdf_path = os.path.join(DIR_OUT, f'mvn_J{J}_pps_RGEQ_p_h1_xz_all.pdf')
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved RGE 'all' p(H_1 | x, z) boxplot to {pdf_path}")

    # ---- 6e. Scatter: analytic per-sample vs RGE per-sample ----
    scatter = an_all_J[['interim_id', 'interim_month_year', 'j', 's',
                        'p_h1_xz']].rename(
        columns={'p_h1_xz': 'p_h1_xz_analytic'},
    ).merge(
        rge_all[['interim_id', 'j', 's', 'p_h1_xz']].rename(
            columns={'p_h1_xz': 'p_h1_xz_rge'},
        ),
        on=['interim_id', 'j', 's'], how='inner',
    )
    scatter['response_label'] = pd.Categorical(
        'response mu_' + scatter['j'].astype(str),
        categories=response_label_cats, ordered=True,
    )
    scatter['interim_month_year'] = pd.Categorical(
        scatter['interim_month_year'], categories=order, ordered=True,
    )
    from ggsci import pal_simpsons
    simpsons_colours = pal_simpsons()(len(order))
    scatter.to_pickle(
        os.path.join(DIR_OUT, f'mvn_J{J}_pps_RGEQ_scatter.pkl'),
    )

    p = (
        ggplot(scatter, aes(x='p_h1_xz_analytic', y='p_h1_xz_rge',
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
        + facet_wrap('~ response_label', ncol=n_per_level, dir='h')
        + theme_bw()
        + theme(axis_text_x=element_text(angle=45, vjust=1, hjust=1),
                figure_size=(3.0 * n_per_level, 2.5 * K_levels),
                legend_position='top',
                strip_background=element_blank(),
                strip_text=element_text(face='bold'))
        + labs(x='p(H_1j | x, z) analytic',
               y='p(H_1j | x, z) RGE')
    )
    pdf_path = os.path.join(DIR_OUT, f'mvn_J{J}_pps_RGEQ_scatter.pdf')
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved RGE vs analytic scatter to {pdf_path}")

    # ---- 6f. Diagnostic: w_ratio vs pps_ratio_x per (response, interim).
    # ---- Build rge_stats first: per-(J, interim, j) rho + r2 across ALL
    # ---- responses (for the cross-method violin plot); annotate the
    # ---- diagnostic only on the j_repr subset. ----
    rge_stats_rows = []
    for (im_id, im, im_my, jj), g in wa_all.groupby(
            ['interim_id', 'interim_date', 'interim_month_year', 'j'],
            observed=True):
        s = fit_interim_regress_H1x_on_wz_per_item_summary(g)
        rge_stats_rows.append({
            'J':                  J,
            'interim_id':         int(im_id),
            'interim_date':       im,
            'interim_month_year': im_my,
            'j':                  int(jj),
            'item_label':         g['item_label'].iloc[0],
            'rho':                float(s['rho']),
            'r2':                 float(s['r2']),
            'x_left':             float(s['x_left']),
            'x_right':            float(s['x_right']),
            'y_top':              float(s['y_top']),
        })
    rge_stats = pd.DataFrame(rge_stats_rows)
    rge_stats.to_pickle(
        os.path.join(DIR_OUT, f'mvn_J{J}_pps_RGEQ_stats.pkl'),
    )

    wa_small = wa_all[wa_all['j'].isin(j_repr)].copy()
    wa_small['response_label'] = 'response mu_' + wa_small['j'].astype(str)
    wa_small['interim_month_year'] = pd.Categorical(
        wa_small['interim_month_year'], categories=order, ordered=True,
    )
    wa_small['response_label'] = pd.Categorical(
        wa_small['response_label'], categories=j_label_cats, ordered=True,
    )

    rge_stats_small = rge_stats[rge_stats['j'].isin(j_repr)].copy()
    rge_stats_small['response_label'] = pd.Categorical(
        'response mu_' + rge_stats_small['j'].astype(str),
        categories=j_label_cats, ordered=True,
    )
    rge_stats_small['interim_month_year'] = pd.Categorical(
        rge_stats_small['interim_month_year'],
        categories=order, ordered=True,
    )
    rge_stats_small['rho_label'] = rge_stats_small['rho'].map(
        lambda r: f"rho={r:.2f}",
    )
    rge_stats_small['r2_label'] = rge_stats_small['r2'].map(
        lambda r: f"R2={100 * r:.1f}%",
    )

    # Precompute per-facet QuantReg line at tau = 1 - eta_H so we can
    # overlay via geom_line (plotnine's geom_smooth has no 'quantile' method).
    tau_used = 1.0 - pps_ProbH1_target_lwr_quantile
    from plotnine import geom_line
    _line_parts = []
    for (imy, resp_lab), gg in wa_small.groupby(
            ['interim_month_year', 'response_label'], observed=True):
        mask = np.isfinite(gg['w_ratio']) & np.isfinite(gg['pps_ratio_x'])
        ggf = gg[mask]
        if len(ggf) < 3:
            continue
        w = ggf['w_ratio'].to_numpy()
        y = ggf['pps_ratio_x'].to_numpy()
        X_fit = sm.add_constant(w)
        try:
            res = sm.QuantReg(y, X_fit).fit(q=tau_used, max_iter=5000)
        except Exception:
            continue
        w_grid = np.linspace(float(w.min()), float(w.max()), 50)
        X_pred = sm.add_constant(w_grid)
        y_pred = np.asarray(res.predict(X_pred), dtype=float)
        _line_parts.append(pd.DataFrame({
            'interim_month_year': imy,
            'response_label':     resp_lab,
            'w_ratio':            w_grid,
            'pps_ratio_x':        y_pred,
        }))
    _qline = pd.concat(_line_parts, ignore_index=True) if _line_parts else pd.DataFrame()
    if not _qline.empty:
        _qline['interim_month_year'] = pd.Categorical(
            _qline['interim_month_year'], categories=order, ordered=True,
        )
        _qline['response_label'] = pd.Categorical(
            _qline['response_label'], categories=j_label_cats, ordered=True,
        )

    p = (
        ggplot(wa_small, aes(x='w_ratio', y='pps_ratio_x'))
        + geom_point(alpha=0.3, size=0.6, colour='#009392')
        + geom_line(_qline, aes(x='w_ratio', y='pps_ratio_x'),
                    colour='black', size=0.7, inherit_aes=False)
        + geom_hline(yintercept=pps_H1_min_effect_size_thresh, colour='#009392',
                     linetype='dashed', size=0.5)
        + geom_text(rge_stats_small,
                    aes(x='x_left', y='y_top', label='rho_label'),
                    ha='left', va='top', size=8, inherit_aes=False)
        + geom_text(rge_stats_small,
                    aes(x='x_right', y='y_top', label='r2_label'),
                    ha='right', va='top', size=8, inherit_aes=False)
        + facet_grid('response_label ~ interim_month_year', scales='free')
        + theme_bw()
        + theme(figure_size=(3.0 * len(order), 3.0 * K_levels),
                strip_background=element_blank(),
                strip_text=element_text(face='bold'),
                legend_position='none')
        + labs(x='empirical endpoint summary w(z) = mean(z_j) - mu_0',
               y=f'endpoint mu_j - mu_0 (quantile line at tau={tau_used:.3f})')
    )
    pdf_path = os.path.join(DIR_OUT, f'mvn_J{J}_endpointratio_vs_wratio.pdf')
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved RGEQ diagnostic scatter to {pdf_path}")
