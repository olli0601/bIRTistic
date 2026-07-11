#!/usr/bin/env python3
"""
Binomial interim analyses, PPS via importance sampling (Case A, fixed x).

Beta-Bernoulli §3.1 simulation: ``N = 500`` Bernoulli observations with
true ``p = 0.4`` spread over one year, monthly interim cutoffs. At each
interim:

- fit the x-posterior on the accrued cohort xi with NumPyro NUTS
  (``BinomialModel.fit_pyro_hmc``);
- analytic Beta-Binomial PPS via ``model.fit_closed_form_pps`` (reference);
- IS-from-x PPS: draw ``pps_z_total`` posterior-predictive samples zi_s
  from xi's fit, then estimate ``p(H_1 | x, z_s)`` by importance-
  reweighting the x-posterior draws with ``w_k = softmax_k log p(zi_s |
  p_k)``. No refits on (x, z_s).

H_1 convention: ``1 - p / p_0 > pps_H1_min_effect_size_thresh`` (right tail of the
directional endpoint ``ratio = 1 - p/p_0`` defined in
``BinomialModel.get_endpoints_per_draw``).

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/Binomial_interim_analyses_with_IS_from_x.py
"""

# %%

import os
import sys
import time
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

import jax
import jax.numpy as jnp
import numpy as np
import pandas as pd
from plotnine import (
    ggplot, aes, geom_boxplot, geom_col, geom_errorbar, geom_hline,
    scale_fill_manual, scale_x_datetime, position_dodge, theme_bw, theme,
    element_text, labs,
)
import warnings
warnings.filterwarnings('ignore')

from model_binomial import BinomialModel

print("✓ Imports successful")

# %%

# =============================================================================
# Configuration
# =============================================================================

np.random.seed(42)
seed = 123

N = 500
TRUE_P = 0.4
PRIOR_A = 1.0           # Beta(1,1) uniform
PRIOR_B = 1.0
P_0 = 0.5

pps_z_total = 200
pps_H1_min_effect_size_thresh = 0.25
pps_ProbH1_target_lwr_quantile = 0.89

dir_out = "/Users/or105/sandbox/bIRTistic/py-binomial-interim-with-IS-260605"
os.makedirs(dir_out, exist_ok=True)

# %%

# =============================================================================
# Simulate Bernoulli data of size N = 500 spaced over one year
# =============================================================================

start = pd.Timestamp('2025-01-01')
end = pd.Timestamp('2025-12-31')
dates = pd.date_range(start, end, periods=N)
rng = np.random.default_rng(seed)
y = rng.binomial(1, TRUE_P, size=N)

dp = pd.DataFrame({
    'pid': np.arange(N),
    'submission_date': dates,
    'y': y.astype(int),
})

dit = pd.DataFrame({
    'item_label': ['intervention'],
    'item_type': ['binomial'],
    'item_high_label': ['higher_is_better'],
})

print(f"  Simulated {N} Bernoulli outcomes, true p={TRUE_P:.3f},"
      f" observed k={int(y.sum())}.")

# %%

# =============================================================================
# Interim grid: one cohort per calendar month
# =============================================================================

month_starts = pd.date_range(start, end, freq='MS')
di = pd.DataFrame({
    'interim_id': range(1, len(month_starts) + 1),
    'month_start': month_starts,
    'interim_date': month_starts + pd.offsets.MonthEnd(0),
})
di['interim_month_year'] = di['month_start'].dt.strftime('%Y-%b')

print(f"  {len(di)} monthly interims from {di['interim_date'].min().date()} "
      f"to {di['interim_date'].max().date()}.")

# %%

# =============================================================================
# Fit per-interim x-posterior with NumPyro NUTS (pyro_HMC). resume=True so
# reruns pick up cached posteriors from {dir_out}/binomial_i{interim_id}_x_*.
# =============================================================================

po = []   # per-interim fitted endpoint draws
for i in range(len(di)):
    interim_id = int(di.iloc[i]['interim_id'])
    interim_date = di.iloc[i]['interim_date']

    dpi = dp[dp["submission_date"] <= interim_date].copy()
    dpi["oid"] = np.arange(1, len(dpi) + 1, dtype=int)
    n_obs = len(dpi)
    if n_obs == 0:
        continue

    print(f"\n--- Interim {interim_id} ({interim_date.date()}): n={n_obs} ---")

    model = BinomialModel(
        dit=dit, dcati=dpi,
        prior_a=PRIOR_A, prior_b=PRIOR_B, p_0=P_0, seed=seed,
    )

    prefix = os.path.join(dir_out, f"binomial_i{interim_id}_pyro_hmc")
    fit = model.fit_closed_form_posterior(
        output_file_prefix=prefix,
        save_to_file=True, resume=True, verbose=False,
    )
    de = model.get_endpoints_per_draw(draws=fit['draws'])
    de['interim_id'] = interim_id
    de['interim_month_year'] = di.iloc[i]['interim_month_year']
    de['interim_date'] = interim_date
    de['method'] = 'pyro_HMC'
    po.append(de)

po = pd.concat(po, ignore_index=True)
print(f"\n  Collected {len(po):,} pyro_HMC posterior draws across "
      f"{po['interim_id'].nunique()} interims.")

# %%

# =============================================================================
# Persist tidy outputs
# =============================================================================

po.to_pickle(os.path.join(dir_out, 'posterior_draws.pkl'))

# %%

# =============================================================================
# PPS via closed-form Beta-Binomial (analytic integration over zi).
# For each interim with m = N - n_obs > 0 future trials, model.fit_closed_form_pps
# returns the analytic Beta-Binomial tail sum: P(P(H_1 | x, z) > eta).
# =============================================================================

n_full = len(dp)
ppsa = []
for i in range(len(di)):
    interim_id = int(di.iloc[i]['interim_id'])
    interim_date = di.iloc[i]['interim_date']
    interim_month_year = di.iloc[i]['interim_month_year']
    dpi = dp[dp["submission_date"] <= interim_date].copy()
    dpi["oid"] = np.arange(1, len(dpi) + 1, dtype=int)
    n_obs = len(dpi)
    m = n_full - n_obs
    if m <= 0:
        continue
    model = BinomialModel(
        dit=dit, dcati=dpi,
        prior_a=PRIOR_A, prior_b=PRIOR_B, p_0=P_0, seed=seed,
    )
    pps = model.fit_closed_form_pps(
        m,
        pps_H1_min_effect_size_thresh=pps_H1_min_effect_size_thresh,
        pps_ProbH1_target_lwr_quantile=pps_ProbH1_target_lwr_quantile,
    )
    ppsa.append({
        'interim_id': interim_id,
        'interim_month_year': interim_month_year,
        'interim_date': interim_date,
        'item_label':      dit.iloc[0]['item_label'],
        'item_type':       dit.iloc[0]['item_type'],
        'item_high_label': dit.iloc[0]['item_high_label'],
        'm': m,
        'pps': float(pps),
        'eta': pps_ProbH1_target_lwr_quantile,
    })

ppsa = pd.DataFrame(ppsa)
ppsa.to_pickle(os.path.join(dir_out, 'pps_closed_form.pkl'))
print(f"\nClosed-form PPS table saved to "
      f"{os.path.join(dir_out, 'pps_closed_form.pkl')}")
print(ppsa.to_string(index=False))

# %%

# =============================================================================
# PPS via importance sampling at every interim.
#
# For each interim with interim_m = n_full - n_obs > 0:
#   1. Fit x-posterior on xi via pyro_HMC (resume from cache).
#   2. Draw S=pps_z_total posterior-predictive samples zi_s for the
#      interim_m missing trials from xi's fit (drawn in posterior order via
#      get_interim_z_from_ypredi(keep_order=True), so column ypred_s
#      aligns with posterior draw s).
#   3. For each s, compute IS weights w_k = softmax_k log p(zi_s | p_k)
#      over the x-posterior draws p_k. p(zi_s | p_k) = product of
#      Bernoulli(p_k) at each future obs. Self-normalised IS estimate:
#      p(H_1 | x, z_s) = sum_k w_k * 1{ratio_k > pps_H1_min_effect_size_thresh}.
#   4. PPS = mean over s of 1{p(H_1 | x, z_s) > pps_ProbH1_target_lwr_quantile}.
# =============================================================================

dp_h1_xz_rows = []
is_perf_rows = []
pps_timing_rows = []
t_all0 = time.time()
for i in range(len(di)):
    interim_id = int(di.iloc[i]['interim_id'])
    interim_date = di.iloc[i]['interim_date']
    interim_month_year = di.iloc[i]['interim_month_year']
    dpi = dp[dp["submission_date"] <= interim_date].copy()
    dpi["oid"] = np.arange(1, len(dpi) + 1, dtype=int)
    n_obs = len(dpi)
    interim_m = n_full - n_obs
    if interim_m <= 0:
        continue
    print(f"\n{'='*70}\nIS-PPS for interim {interim_id} ({interim_date.date()})\n{'='*70}")
    print(f"  n_obs={n_obs:,} | m={interim_m}")

    t0 = time.time()

    # fit model on today's data x
    interim_prefix = os.path.join(dir_out, f"binomial_i{interim_id}_pyro_hmc")
    model = BinomialModel(
        dit=dit, dcati=dpi,
        prior_a=PRIOR_A, prior_b=PRIOR_B, p_0=P_0, seed=seed,
    )
    fit = model.fit_closed_form_posterior(
        output_file_prefix=interim_prefix,
        save_to_file=True, resume=True, verbose=False,
    )

    # expand to future data keeping draw index linked to posterior draws
    zi = model.get_interim_z_from_ypredi(
        f"{interim_prefix}_draws.zarr", interim_m,
        pps_z_total=pps_z_total, seed=seed, keep_order=True,
    )

    # calculate endpoint on theta | x
    x_ratio = model.get_endpoints_per_draw(draws=fit['draws'])
    ind = (x_ratio['ratio'].to_numpy() > pps_H1_min_effect_size_thresh).astype(float)  # (K,)

    # stacked x-posterior parameters for IS likelihood evaluation
    theta = model.get_stacked_posterior(fit['draws'])
    n_draw = int(theta['p'].shape[0])

    item_row = dit.iloc[0]
    # IS loop: one (K,)-vector of weights per s
    for s_idx in range(pps_z_total):
        s_label = s_idx + 1
        zcol = f'ypred_{s_idx}'
        z_y = zi[zcol].astype(int).to_numpy()                # (interim_m,)
        z_stan = {
            'y': jnp.asarray(z_y, dtype=jnp.int32),
            'N_total': int(z_y.shape[0]),
            'k': int(z_y.sum()),
            'a': model.prior_a,
            'b': model.prior_b,
        }
        # log p(z_s | theta_k) = sum_i log Bernoulli(z_s,i | p_k); vmap over draws.
        loglik = jax.vmap(lambda p: model.eval_loglik(z_stan, p))(theta)  # (K, interim_m)
        w = np.asarray(jax.nn.softmax(loglik.sum(axis=1)))               # (K,)

        sum_w2 = float(np.sum(w ** 2))
        is_perf_rows.append({
            'interim_id': interim_id,
            's': s_label,
            'N': n_draw,
            'ess': float(1.0 / sum_w2),
            'ess_over_n': float(1.0 / (n_draw * sum_w2)),
            'ew2': float(np.mean(w ** 2)),
        })

        # Self-normalised IS estimate of P(H_1 | x, z_s).
        dp_h1_xz_s = float((w * ind).sum())
        dp_h1_xz_rows.append({
            'interim_id':      interim_id,
            'interim_date':    interim_date,
            'interim_month_year': interim_month_year,
            'item_label':      item_row['item_label'],
            'item_type':       item_row['item_type'],
            'item_high_label': item_row['item_high_label'],
            's':               s_label,
            'p_h1_xz':         dp_h1_xz_s,
        })

    mins_interim_id = (time.time() - t0) / 60.0
    pps_timing_rows.append({
        'interim_id': interim_id,
        'mins_interim_id': round(mins_interim_id, 3),
    })
    print(f"  interim {interim_id} IS-PPS done in {mins_interim_id:.2f} min")

if not dp_h1_xz_rows:
    raise RuntimeError("[IS-PPS] no interims with missing participants to evaluate.")

dp_h1_xz = pd.DataFrame(dp_h1_xz_rows)
dp_h1_xz['pps_H1_min_effect_size_thresh'] = pps_H1_min_effect_size_thresh
dp_h1_xz['pps_ProbH1_target_lwr_quantile'] = pps_ProbH1_target_lwr_quantile
dp_h1_xz['S'] = pps_z_total
is_perf = pd.DataFrame(is_perf_rows)
pps_timing = pd.DataFrame(pps_timing_rows)
pps_timing['mins_total'] = round((time.time() - t_all0) / 60.0, 3)

dp_h1_xz.to_pickle(os.path.join(dir_out, 'pps_p_h1_xz_IS.pkl'))
is_perf.to_csv(os.path.join(dir_out, 'pps_is_perf.csv'), index=False)
pps_timing.to_csv(os.path.join(dir_out, 'pps_timing.csv'), index=False)
print(f"\nTotal IS-PPS time: {pps_timing['mins_total'].iloc[0]:.2f} min")

# %%

# =============================================================================
# Aggregate: PPS = mean over s of 1{p(H_1 | x, z_s) > eta}
# =============================================================================

pps_df = (
    dp_h1_xz
    .groupby(['interim_id', 'interim_date', 'interim_month_year',
              'item_label', 'item_type', 'item_high_label'])['p_h1_xz']
    .apply(lambda p: float((p > pps_ProbH1_target_lwr_quantile).mean()))
    .reset_index(name='pps')
)
pps_df['eta'] = pps_ProbH1_target_lwr_quantile
pps_df['S'] = pps_z_total
pps_df.to_pickle(os.path.join(dir_out, 'pps_IS.pkl'))
pps_df.to_csv(os.path.join(dir_out, 'pps_IS.csv'), index=False)
print(f"\nIS PPS:\n{pps_df.to_string(index=False)}")

# %%

# =============================================================================
# Boxplot: distribution of p(H_1 | x, z) across the S IS samples per interim
# =============================================================================

order = (
    di[di['interim_id'].isin(dp_h1_xz['interim_id'].unique())]
    .sort_values('interim_date')['interim_month_year'].tolist()
)
dp_h1_xz['interim_month_year'] = pd.Categorical(
    dp_h1_xz['interim_month_year'], categories=order, ordered=True,
)

tmp = (
    dp_h1_xz
    .groupby(['interim_month_year'], observed=True)['p_h1_xz']
    .agg(
        q025=lambda s: s.quantile(0.025),
        q25=lambda s: s.quantile(0.25),
        q50=lambda s: s.quantile(0.50),
        q75=lambda s: s.quantile(0.75),
        q975=lambda s: s.quantile(0.975),
    )
    .reset_index()
)
tmp.to_pickle(os.path.join(dir_out, 'binomial_interim_pps_p_h1_xz_boxplot.pkl'))

p = (
    ggplot(
        tmp,
        aes(x='interim_month_year',
            ymin='q025', lower='q25', middle='q50',
            upper='q75', ymax='q975',
            group='interim_month_year'),
    )
    + geom_boxplot(stat='identity', fill='#9467bd', alpha=0.4, colour='#808080')
    + geom_hline(yintercept=pps_ProbH1_target_lwr_quantile, colour='black', size=1.0)
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        figure_size=(14, 5),
    )
    + labs(
        x='Interim', y='p(H_1 | x, z)',
        title=(f'Binomial IS p(H_1 | x, z) distribution per interim '
               f'(S = {pps_z_total}, eta = {pps_ProbH1_target_lwr_quantile}).'),
    )
)
tmp = os.path.join(dir_out, 'binomial_interim_pps_p_h1_xz_boxplot.pdf')
p.save(tmp, verbose=False, limitsize=False)
print(f"Saved IS p(H_1 | x, z) boxplot to: {tmp}")

# %%

# =============================================================================
# Bar chart: PPS per interim_date. Black bar = closed-form analytic PPS;
# #9467bd bar = IS PPS. Dodged side-by-side. Vertical bars on the IS bar =
# 95% bootstrap interval (B resamples with replacement of the S MC samples
# per interim, q025/q975 over the resulting bootstrap PPS dist).
# =============================================================================

bs_B = 2000
quantiles = [0.025, 0.25, 0.5, 0.75, 0.975]
quantile_names = ['q025', 'q25', 'q50', 'q75', 'q975']

wide = (
    dp_h1_xz
    .pivot_table(index=['interim_id', 'interim_date'],
                 columns='s', values='p_h1_xz')
    .sort_index()
)
bs = wide.to_numpy()                                      # (n_interim, S)
n_interim, S = bs.shape
rng = np.random.default_rng(seed)
tmp = rng.integers(0, S, size=(n_interim, bs_B, S))         # (n_interim, B, S)
bs = bs[np.arange(n_interim)[:, None, None], tmp]       # (n_interim, B, S)
bs = (bs > pps_ProbH1_target_lwr_quantile).mean(axis=2)                # (n_interim, B)
bs = pd.DataFrame(np.quantile(bs, quantiles, axis=1).T,
                  columns=quantile_names
                  )
bs['interim_id'] = wide.index.get_level_values('interim_id').to_numpy()
bs['interim_date'] = pd.to_datetime(
    wide.index.get_level_values('interim_date'),
)
bs['method'] = 'IS'

tmp = (
    dp_h1_xz.groupby(
        ['interim_id', 'interim_date', 'interim_month_year'],
        observed=True,
    )['p_h1_xz']
    .apply(lambda s: float((s > pps_ProbH1_target_lwr_quantile).mean()))
    .reset_index(name='pps')
)
tmp = pd.concat(
    [
        ppsa[['interim_id', 'interim_date', 'pps']].assign(method='analytic'),
        tmp[['interim_id', 'interim_date', 'pps']].assign(method='IS'),
    ],
    ignore_index=True,
)
tmp['interim_date'] = pd.to_datetime(tmp['interim_date'])
tmp['method'] = pd.Categorical(
    tmp['method'], categories=['analytic', 'IS'], ordered=True,
)
tmp = tmp.merge(
    bs[['interim_id', 'method', 'q025', 'q975']],
    on=['method', 'interim_id'], how='left',
)
tmp.to_pickle(os.path.join(dir_out, 'binomial_interim_pps_bootstrap.pkl'))

_bar_colours = {'analytic': '#000000', 'IS': '#9467bd'}
p = (
    ggplot(
        tmp,
        aes(x='interim_date', y='pps', fill='method'),
    )
    + geom_col(position=position_dodge(width=20), width=18)
    + geom_errorbar(
        data=tmp,
        mapping=aes(x='interim_date', ymin='q025', ymax='q975',
                    group='method'),
        position=position_dodge(width=20), width=8,
        colour='black', size=0.5, inherit_aes=False,
    )
    + scale_fill_manual(values=_bar_colours)
    + scale_x_datetime(date_breaks='1 month', date_labels='%Y-%b')
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        figure_size=(14, 5),
        legend_position='top',
    )
    + labs(
        x='interim date', y='PPS = P(p(H_1 | x, z) > eta)',
        fill='method',
    )
)
tmp = os.path.join(dir_out, 'binomial_interim_pps_bars.pdf')
p.save(tmp, verbose=False, limitsize=False)
print(f"Saved PPS bar chart to: {tmp}")
