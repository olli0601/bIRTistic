#!/usr/bin/env python3
"""
MVN interim PPS via self-normalised importance sampling from the x-fit
posterior (Case A, fixed x).

Standalone -- loads the simulated x cohort + per-(J, interim) cached
future-data block (``zi``) + first ``PPS_Z_TOTAL`` outer-HMC ``mu``
draws (``mu_draws``) from the simulations directory. No refits on
(x, z_s). Per (J, interim, s):

  1. ``theta_k = mu_draws[k]`` for ``k = 1..K`` (the K = PPS_Z_TOTAL
     cached posterior draws).
  2. ``log p(z_s | theta_k) = -1/(2 sigma^2) sum_i (z_{s,i} - theta_k)^T
     K^{-1} (z_{s,i} - theta_k)`` (constants drop in softmax). Evaluated
     in closed-form via sufficient-statistic expansion
     ``quad_k = T_2(z_s) - 2 theta_k^T K^{-1} T_1(z_s)
                + m * theta_k^T K^{-1} theta_k``.
  3. ``w_k = softmax_k(-quad_k / (2 sigma^2))``.
  4. Self-normalised IS estimate per component j:
     ``p(H_{1j} | x, z_s) = sum_k w_k * 1{theta_{k,j} - mu_0 > pps_H1_def}``.
  5. PPS = mean over s of ``1{p(H_{1j} | x, z_s) > pps_ProbH1_thresh}``.

Outputs per J: ``mvn_J{J}_pps_IS_p_h1_xz.pkl``, ``mvn_J{J}_pps_IS.csv``,
``mvn_J{J}_pps_IS_perf.csv`` (ESS / particle, E(w^2)),
``mvn_J{J}_pps_IS_timing.csv``, plus the standard quintet of plots
adapted from the nested-MC HMC layout.

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/MVN_interim_analyses_with_IS_from_x.py
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
from scipy.special import softmax as _softmax
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
from utils import _material_palette

print("Imports successful")

# %%

# =============================================================================
# Configuration.
# =============================================================================

DIR_SIM = os.path.join(
    "/Users/or105/sandbox/bIRTistic", "py-mvn-interim-simulations-260609",
)
DIR_OUT = os.path.join(
    "/Users/or105/sandbox/bIRTistic", "py-mvn-interim-with-IS-260611",
)
os.makedirs(DIR_OUT, exist_ok=True)

PPS_Z_TOTAL = 200
svi_algorithm_automultivariatenormal = 'AutoMultivariateNormal'
SVI_LR = 0.05
SVI_NUM_STEPS = 2000
SVI_OUTPUT_SAMPLES = 4000

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
pps_H1_def        = simu_params['pps_H1_def']
pps_ProbH1_thresh = simu_params['pps_ProbH1_thresh']
mu_0_baseline     = simu_params['mu_0_baseline']
sigma             = simu_params['sigma']
sigma2            = sigma ** 2

# %%

# =============================================================================
# Phase 5. Per (J, interim) self-normalised IS estimator + analytic
# per-(j, s) labels (in-memory) for downstream comparison plots.
# =============================================================================

is_artifacts = {}
for J in J_GRID:
    dit = sim_data['cells'][J]['dit']
    dp = sim_data['cells'][J]['dp']
    R_chol = sim_data['cells'][J]['R_chol']
    interim_data = interim_data_by_J[J]
    # Re-use cached HMC zarrs co-located with the interim_data pkl so
    # resume=True picks them up instantly.
    dir_J = os.path.join(DIR_SIM, f"nested_mc_J{J}")
    os.makedirs(dir_J, exist_ok=True)

    print(f"\n{'#' * 70}\n# J = {J} (n_items = {J})\n{'#' * 70}")

    p_h1_xz_rows = []
    is_perf_rows = []
    timing_rows = []
    an_rows = []
    t_J0 = time.time()
    # Index loop so the cell body can be stepped through interactively.
    for i in range(len(di)):
        interim_id = int(di.iloc[i]['interim_id'])
        interim_date = di.iloc[i]['interim_date']
        interim_month_year = di.iloc[i]['interim_month_year']
        if interim_id not in interim_data:
            continue
        blk = interim_data[interim_id]
        interim_m = int(blk['interim_m'])
        n_obs = int(blk['n_obs'])
        if interim_m <= 0:
            continue

        print(f"\n{'=' * 70}\nJ={J} IS interim {interim_id}"
              f" ({interim_date.date()}) n={n_obs} m={interim_m}"
              f" S={PPS_Z_TOTAL}\n{'=' * 70}")
        t0 = time.time()

        # a) fit model on today's data x via pyro_SVI with
        #    AutoMultivariateNormal autoguide. resume=False so timing
        #    matches a fresh run (comparable to fresh-run nested-MC HMC).
        xi = MVNModel.get_interim_data_x(
            dp[dp['submission_date'] <= interim_date],
        )
        interim_prefix = os.path.join(dir_J, f"mvn_pps_is_i{interim_id}")
        model_x = MVNModel(
            dit=dit, dcati=xi, seed=seed,
            J=J, K_chol=R_chol, sigma=sigma,
            prior_tau=simu_params['prior_tau'],
            mu_0_baseline=mu_0_baseline,
        )
        fit_x = model_x.fit_pyro_svi(
            output_file_prefix=f"{interim_prefix}_x",
            algorithm=svi_algorithm_automultivariatenormal,
            lr=SVI_LR, num_steps=SVI_NUM_STEPS,
            output_samples=SVI_OUTPUT_SAMPLES,
            save_to_file=True, resume=False, verbose=False,
        )

        # b) expand to future data keeping draw index linked to posterior
        #    draws (zi[ypred_s] aligned with theta draw s).
        zi = model_x.get_interim_z_from_ypredi(
            f"{interim_prefix}_x_draws.zarr", interim_m,
            pps_z_total=PPS_Z_TOTAL, seed=seed, keep_order=True,
        )

        # c) calculate endpoint on theta | x: ind_{k,j} = 1{ratio_{k,j}
        #    = mu_{k,j} - mu_0 > pps_H1_def} via the per-(item, draw) frame
        #    (matches Binomial_interim_analyses_with_IS_from_x.py line 280).
        x_ratio = model_x.get_endpoints_per_draw(draws=fit_x['draws'])
        ind_K_J = (
            x_ratio.pivot_table(index='draw', columns='j', values='ratio')
                   .sort_index().to_numpy() > pps_H1_def
        ).astype(np.float64)                                       # (K, J)
        theta = model_x.get_stacked_posterior(fit_x['draws'])
        mu_K = np.asarray(theta['mu'], dtype=np.float64)            # (K, J)
        K = mu_K.shape[0]
        K_inv = model_x.K_inv                                      # (J, J)
        mu_kinv  = mu_K @ K_inv                                    # (K, J)
        mu_kinv_mu = (mu_kinv * mu_K).sum(axis=1)                  # (K,)

        # d) IS loop. Reshape zi's ypred_s columns to (S, m, J) row-major
        #    over (s, pid, j), then evaluate the self-normalised IS weights
        #    via sufficient statistics:
        #      log p(z_s | theta_k) propto - quad_{k,s} / (2 sigma^2),
        #      quad_{k,s} = T_2(z_s) - 2 theta_k^T K^{-1} T_1(z_s)
        #                   + m theta_k^T K^{-1} theta_k.
        #    Vectorised across (k, s); scipy.special.softmax handles the
        #    log-sum-exp shift internally.
        zi_cols = sorted(
            [c for c in zi.columns if c.startswith('ypred_')],
            key=lambda c: int(c.split('_')[1]),
        )[:PPS_Z_TOTAL]
        zi_sorted = zi.sort_values(['pid', 'j']).reset_index(drop=True)
        zi_arr = zi_sorted[zi_cols].to_numpy().T.reshape(
            PPS_Z_TOTAL, interim_m, J,
        )                                                          # (S, m, J)
        T1_S = zi_arr.sum(axis=1)                               # (S, J)
        T2_S = np.einsum('sij,jk,sik->s', zi_arr, K_inv, zi_arr)
        quad = (
            T2_S[None, :]                                          # (1, S)
            - 2.0 * mu_kinv @ T1_S.T                               # (K, S)
            + interim_m * mu_kinv_mu[:, None]                      # (K, 1)
        )                                                          # (K, S)
        w = _softmax(-quad / (2.0 * sigma2), axis=0)               # (K, S)

        p_h1_xz = w.T @ ind_K_J                                    # (S, J)
        sum_w2 = (w ** 2).sum(axis=0)                              # (S,)
        ess        = 1.0 / sum_w2
        ess_over_n = 1.0 / (K * sum_w2)
        ew2        = sum_w2 / K

        # Long form (s, j) for p_h1_xz.
        s_idx_grid = np.repeat(np.arange(1, PPS_Z_TOTAL + 1), J)
        j_grid = np.tile(np.arange(J), PPS_Z_TOTAL)
        item_labels = dit['item_label'].to_numpy()
        tmp = pd.DataFrame({
            'J':                  J,
            'interim_id':         interim_id,
            'interim_date':       interim_date,
            'interim_month_year': interim_month_year,
            'item_label':         np.tile(item_labels, PPS_Z_TOTAL),
            'item_type':          np.tile(dit['item_type'].to_numpy(),
                                          PPS_Z_TOTAL),
            'item_high_label':    np.tile(dit['item_high_label'].to_numpy(),
                                          PPS_Z_TOTAL),
            'j':                  j_grid,
            's':                  s_idx_grid,
            'p_h1_xz':            p_h1_xz.reshape(-1),
        })
        p_h1_xz_rows.append(tmp)

        is_perf_rows.append(pd.DataFrame({
            'J':                  J,
            'interim_id':         interim_id,
            'interim_date':       interim_date,
            'interim_month_year': interim_month_year,
            's':                  np.arange(1, PPS_Z_TOTAL + 1),
            'N':                  K,
            'ess':                ess,
            'ess_over_n':         ess_over_n,
            'ew2':                ew2,
        }))

        # Record IS interim wall time BEFORE the analytic-only comparison
        # so the comparison-method timing reflects SVI fit + zi expand +
        # IS reweight only.
        mins = (time.time() - t0) / 60.0

        # Analytic comparison (NOT part of the IS timing). Used only by
        # downstream comparison plots.
        an = model_x.fit_closed_form_pH1(zi_sorted[
            [c for c in zi_sorted.columns if not c.startswith('ypred_')]
            + ypred_cols
        ])
        an['J'] = J
        an['interim_id'] = interim_id
        an['interim_date'] = interim_date
        an['interim_month_year'] = interim_month_year
        an_rows.append(an)

        timing_rows.append({
            'J':                  J,
            'interim_id':         interim_id,
            'interim_date':       interim_date,
            'interim_month_year': interim_month_year,
            'n_obs':              n_obs,
            'm_future':           interim_m,
            'mins_interim_id':    round(mins, 3),
        })
        print(f"  J={J} interim {interim_id} IS done in"
              f" {mins:.2f} min; K={K}; median ESS/K"
              f" = {np.median(ess_over_n):.3g}")

    if not p_h1_xz_rows:
        print(f"J={J}: no interims with m > 0; skipping save.")
        continue

    dp_h1_xz = pd.concat(p_h1_xz_rows, ignore_index=True)
    dp_h1_xz['pps_H1_def'] = pps_H1_def
    dp_h1_xz['pps_ProbH1_thresh'] = pps_ProbH1_thresh
    dp_h1_xz['S'] = PPS_Z_TOTAL
    is_perf = pd.concat(is_perf_rows, ignore_index=True)
    timing_all = pd.DataFrame(timing_rows)
    timing_all['mins_J_total'] = round((time.time() - t_J0) / 60.0, 3)
    an_all = pd.concat(an_rows, ignore_index=True)

    pps_df = (
        dp_h1_xz.groupby(['J', 'interim_id', 'interim_date',
                          'interim_month_year', 'item_label', 'item_type',
                          'item_high_label', 'j'], observed=True)['p_h1_xz']
        .apply(lambda p: float(
            (p > pps_ProbH1_thresh).mean()
        ))
        .reset_index(name='pps')
    )
    pps_df['eta'] = pps_ProbH1_thresh
    pps_df['S'] = PPS_Z_TOTAL

    dp_h1_xz.to_pickle(os.path.join(DIR_OUT, f'mvn_J{J}_pps_IS_p_h1_xz.pkl'))
    is_perf.to_csv(os.path.join(DIR_OUT, f'mvn_J{J}_pps_IS_perf.csv'),
                   index=False)
    timing_all.to_csv(os.path.join(DIR_OUT, f'mvn_J{J}_pps_IS_timing.csv'),
                      index=False)
    pps_df.to_csv(os.path.join(DIR_OUT, f'mvn_J{J}_pps_IS.csv'), index=False)
    an_all.to_pickle(os.path.join(DIR_OUT,
                                  f'mvn_J{J}_pps_IS_p_h1_xz_analytic.pkl'))
    is_artifacts[J] = {
        'dp_h1_xz': dp_h1_xz, 'perf': is_perf,
        'timing':   timing_all, 'pps': pps_df, 'an': an_all,
    }
    print(f"J={J}: saved IS artifacts."
          f" Total time {(time.time() - t_J0) / 60.0:.2f} min.")

# %%

# =============================================================================
# Phase 6. Per-J plots: PPS bars / boxplot / scatter / "all" / ESS bar.
# =============================================================================

bs_B = 2000
METHOD_LABEL = 'IS reweighting of theta|x'

for J in J_GRID:
    if J not in is_artifacts:
        continue
    art = is_artifacts[J]
    dp_h1_xz = art['dp_h1_xz']
    pps_df   = art['pps']
    is_perf  = art['perf']
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

    is_small      = dp_h1_xz[dp_h1_xz['j'].isin(j_repr)].copy()
    cf_small      = pps_cf[(pps_cf['J'] == J) & (pps_cf['j'].isin(j_repr))].copy()
    is_pps_small  = pps_df[pps_df['j'].isin(j_repr)].copy()
    an_small      = an_all[an_all['j'].isin(j_repr)].copy()

    j_label_cats = [f'response mu_{j}' for j in j_repr]
    for df in (is_small, cf_small, is_pps_small, an_small):
        df['j_label'] = 'response mu_' + df['j'].astype(str)
        df['interim_month_year'] = pd.Categorical(
            df['interim_month_year'], categories=order, ordered=True,
        )
        df['j_label'] = pd.Categorical(
            df['j_label'], categories=j_label_cats, ordered=True,
        )

    palette = _material_palette('teal', len(j_repr))
    fill_values = {'analytic': '#000000'}
    fill_values.update({f'mu_{j}': c for j, c in zip(j_repr, palette)})
    method_repr_key = f'mu_{j_repr[-1]}'

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

    # ---- 6a. Boxplot of p(H_1j | x, z) ----
    box = pd.concat(
        [_qtab(an_small, 'analytic'), _qtab(is_small, METHOD_LABEL)],
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
    box.to_pickle(os.path.join(DIR_OUT, f'mvn_J{J}_pps_IS_p_h1_xz_box.pkl'))

    p = (
        ggplot(box, aes(x='interim_month_year',
                        ymin='q025', lower='q25', middle='q50',
                        upper='q75', ymax='q975',
                        fill='fill_key', group='grp'))
        + geom_boxplot(stat='identity',
                       position=position_dodge(width=0.8), width=0.7,
                       colour='#404040', size=0.3)
        + geom_hline(yintercept=pps_ProbH1_thresh,
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
    pdf_path = os.path.join(DIR_OUT, f'mvn_J{J}_pps_IS_p_h1_xz.pdf')
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved IS p(H_1 | x, z) boxplot to {pdf_path}")

    # ---- 6b. PPS bars: analytic + IS + bootstrap CI ----
    wide = (
        is_small.pivot_table(index=['j', 'j_label', 'interim_id',
                                    'interim_date'],
                             columns='s', values='p_h1_xz')
        .sort_index()
    )
    bs = wide.to_numpy()
    n_rows, S = bs.shape
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, S, size=(n_rows, bs_B, S))
    bs = bs[np.arange(n_rows)[:, None, None], idx]
    bs = (bs > pps_ProbH1_thresh).mean(axis=2)
    ci = np.quantile(bs, [0.025, 0.975], axis=1).T
    ci_df = pd.DataFrame(ci, columns=['q025_bs', 'q975_bs'])
    ci_df['j']          = wide.index.get_level_values('j').to_numpy()
    ci_df['interim_id'] = wide.index.get_level_values('interim_id').to_numpy()

    is_bar = is_pps_small[['interim_id', 'interim_date',
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
    bar = pd.concat([cf_bar, is_bar], ignore_index=True)
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
    bar.to_pickle(os.path.join(DIR_OUT, f'mvn_J{J}_pps_IS_bars.pkl'))

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
    pdf_path = os.path.join(DIR_OUT, f'mvn_J{J}_pps_IS_bars.pdf')
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved IS PPS bar chart to {pdf_path}")

    # ---- 6c, 6d. "all" variants ----
    is_all_pps = pps_df[pps_df['j'].isin(j_indices_all)].copy()
    cf_all  = pps_cf[(pps_cf['J'] == J) & (pps_cf['j'].isin(j_indices_all))].copy()
    is_all = dp_h1_xz[dp_h1_xz['j'].isin(j_indices_all)].copy()
    an_all_J = an_all[an_all['j'].isin(j_indices_all)].copy()
    response_label_cats = [f'response mu_{j}' for j in j_indices_all]
    for df in (is_all_pps, cf_all, is_all, an_all_J):
        df['response_label'] = 'response mu_' + df['j'].astype(str)
        df['interim_month_year'] = pd.Categorical(
            df['interim_month_year'], categories=order, ordered=True,
        )
        df['response_label'] = pd.Categorical(
            df['response_label'], categories=response_label_cats, ordered=True,
        )

    palette_all = _material_palette('teal', len(j_indices_all))
    fill_values_all = {'analytic': '#000000'}
    fill_values_all.update({
        f'mu_{j}': c for j, c in zip(j_indices_all, palette_all)
    })
    method_repr_key_all = f'mu_{j_indices_all[-1]}'

    wide_all = (
        is_all.pivot_table(index=['j', 'interim_id', 'interim_date'],
                           columns='s', values='p_h1_xz')
        .sort_index()
    )
    bs_a = wide_all.to_numpy()
    n_rows_a, S_a = bs_a.shape
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, S_a, size=(n_rows_a, bs_B, S_a))
    bs_a = bs_a[np.arange(n_rows_a)[:, None, None], idx]
    bs_a = (bs_a > pps_ProbH1_thresh).mean(axis=2)
    ci_a = np.quantile(bs_a, [0.025, 0.975], axis=1).T
    ci_df_a = pd.DataFrame(ci_a, columns=['q025_bs', 'q975_bs'])
    ci_df_a['j']          = wide_all.index.get_level_values('j').to_numpy()
    ci_df_a['interim_id'] = wide_all.index.get_level_values('interim_id').to_numpy()

    is_bar_all = is_all_pps[['interim_id', 'interim_date',
                             'interim_month_year', 'j', 'response_label',
                             'pps']].assign(method=METHOD_LABEL).merge(
        ci_df_a, on=['j', 'interim_id'], how='left',
    )
    cf_bar_all = cf_all[['interim_id', 'interim_date',
                         'interim_month_year', 'j', 'response_label',
                         'pps']].assign(method='analytic',
                                        q025_bs=np.nan, q975_bs=np.nan)
    bar_all = pd.concat([cf_bar_all, is_bar_all], ignore_index=True)
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
    pdf_path = os.path.join(DIR_OUT, f'mvn_J{J}_pps_IS_bars_all.pdf')
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved IS 'all' PPS bar chart to {pdf_path}")

    # boxplot all
    box_all = pd.concat(
        [_qtab(an_all_J.rename(columns={'response_label': 'j_label'}),
               'analytic'),
         _qtab(is_all.rename(columns={'response_label': 'j_label'}),
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
        + geom_hline(yintercept=pps_ProbH1_thresh,
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
    pdf_path = os.path.join(DIR_OUT, f'mvn_J{J}_pps_IS_p_h1_xz_all.pdf')
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved IS 'all' p(H_1 | x, z) boxplot to {pdf_path}")

    # ---- 6e. Scatter analytic vs IS per-sample ----
    scatter = an_all_J[['interim_id', 'interim_month_year', 'j', 's',
                        'p_h1_xz']].rename(
        columns={'p_h1_xz': 'p_h1_xz_analytic'},
    ).merge(
        is_all[['interim_id', 'j', 's', 'p_h1_xz']].rename(
            columns={'p_h1_xz': 'p_h1_xz_IS'},
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
    scatter.to_pickle(os.path.join(DIR_OUT, f'mvn_J{J}_pps_IS_scatter.pkl'))

    p = (
        ggplot(scatter, aes(x='p_h1_xz_analytic', y='p_h1_xz_IS',
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
               y='p(H_1j | x, z) IS')
    )
    pdf_path = os.path.join(DIR_OUT, f'mvn_J{J}_pps_IS_scatter.pdf')
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved IS vs analytic scatter to {pdf_path}")

    # ---- 6f. IS ESS / K diagnostic ----
    is_perf['interim_month_year'] = pd.Categorical(
        is_perf['interim_month_year'], categories=order, ordered=True,
    )
    ess_bar = (
        is_perf.groupby('interim_month_year',
                        observed=True)['ess_over_n']
        .median().reset_index(name='median_ess_over_n')
    )
    p = (
        ggplot(ess_bar, aes(x='interim_month_year',
                            y='median_ess_over_n'))
        + geom_col(fill='#009392', width=0.7)
        + scale_y_continuous(
            breaks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            labels=['0%', '20%', '40%', '60%', '80%', '100%'],
        )
        + theme_bw()
        + theme(axis_text_x=element_text(angle=45, vjust=1, hjust=1),
                figure_size=(12, 4),
                strip_background=element_blank(),
                strip_text=element_text(face='bold'))
        + labs(x='Interim', y='median ESS / K (over S samples)')
    )
    pdf_path = os.path.join(DIR_OUT, f'mvn_J{J}_pps_IS_ess.pdf')
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"J={J}: saved IS ESS/K diagnostic to {pdf_path}")
