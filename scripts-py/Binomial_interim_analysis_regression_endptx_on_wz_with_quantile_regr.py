#!/usr/bin/env python3
"""
Binomial interim analysis: regression-based EVSI labels (Strong & Oakley,
2014), continuous-target variant with **conditional-quantile regression**
instead of the Gaussian mean + variance approximation used by
``Binomial_interim_analysis_regression_endptx_on_wz_with_gauss_approx.py``.

Fits, per interim,

    Q_{1 - eta_H}(rho | x, w(z))  ~  pinball loss on {(w^(s), rho^(s))}

with a linear ``statsmodels.QuantReg`` at tau = 1 - eta_H, so
``P(H_1 | x, z) > eta_H`` is equivalent to ``Q_hat > eta_0`` — no Gaussian
step, no plug-in variance.

Outputs (``_RGEQ_`` suffix so compare-methods can load them next to the
Gaussian ``_RGE_`` variant):
p_h1_xz pkl, PPS csv, box-stats pkl, p_h1_xz boxplot, PPS bars with 95%
bootstrap CI vs analytic, perf long-form pkl (rho / R^2 / time(min)),
timing csv, plus a per-interim diagnostic scatter (w_ratio vs pps_ratio_x
with quantile-regression line + rho / R^2 annotations).

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/Binomial_interim_analysis_regression_endptx_on_wz_with_quantile_regr.py
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

import numpy as np
import pandas as pd
from plotnine import (
    ggplot, aes, geom_boxplot, geom_col, geom_errorbar, geom_hline, geom_point,
    geom_smooth, geom_text, facet_wrap, scale_fill_manual,
    scale_x_datetime, scale_y_continuous, position_dodge,
    theme_bw, theme, element_text, element_rect, labs,
)
import warnings
warnings.filterwarnings('ignore')

from model_binomial import BinomialModel
from fit_interim import (
    fit_interim_regress_endptx_on_wz_with_quantile_regr,
    fit_interim_regress_H1x_on_wz_per_item_summary,
)

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

pps_H1_min_effect_size_thresh = 0.25       # endpoint = 1 - p / p_0 > 0.25 (i.e. p < 0.375)
pps_ProbH1_target_lwr_quantile = 0.89
pps_z_total = 200

hmc_chains = 4
hmc_iter_warmup = 500
hmc_iter_sampling = 1000

dir_out = "/Users/or105/sandbox/bIRTistic/py-binomial-interim-with-regression-on-endptx-wz-quantile-regr-260702"
os.makedirs(dir_out, exist_ok=True)

file_prefix = "binomial_interim"

print(f"Output dir: {dir_out}")

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

n_full = len(dp)
print(f"  n_full={n_full} | looping over {len(di)} monthly interim rounds.")

# %%

# =============================================================================
# PPS via closed-form Beta-Binomial (analytic integration over zi).
# =============================================================================

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
        m, pps_H1_min_effect_size_thresh=pps_H1_min_effect_size_thresh, pps_ProbH1_target_lwr_quantile=pps_ProbH1_target_lwr_quantile,
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
ppsa.to_pickle(os.path.join(dir_out, f"{file_prefix}_pps_closed_form.pkl"))
print(f"\nClosed-form PPS table:\n{ppsa.to_string(index=False)}")

# %%

# =============================================================================
# Per-interim loop:
#   1. fit x-posterior on xi (pyro_HMC, resume=True)
#   2. build zi keeping posterior-draw order
#   3. wa = model.get_w(zi) + model.get_endpoints_per_draw -> pps_ratio_x
#   4. fit_interim_regress_endptx_on_wz_with_quantile_regr(wa)
# =============================================================================

p_h1_xz_rows = []
perf_rows = []
pps_timing_rows = []
t_all0 = time.time()
for i in range(len(di)):
    interim_id = int(di.iloc[i]['interim_id'])
    interim_date = di.iloc[i]['interim_date']
    interim_month_year = di.iloc[i]['interim_month_year']
    dpi = dp[dp["submission_date"] <= interim_date].copy()
    dpi["oid"] = np.arange(1, len(dpi) + 1, dtype=int)
    n_obs = len(dpi)
    if n_obs == 0:
        continue
    interim_m = n_full - n_obs
    if interim_m <= 0:
        print(f"\n[interim {interim_id}] skipped: no missing participants"); continue
    print(f"\n{'='*70}\nRegression-endptx (quantile) for interim {interim_id} ({interim_date.date()})\n{'='*70}")
    print(f"  n_obs={n_obs:,} | m={interim_m}")

    t0 = time.time()

    # fit model on today's data x
    interim_prefix = os.path.join(dir_out, f"{file_prefix}_i{interim_id}_pyro_hmc")
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
    wa = model.get_w(zi)

    # calculate endpoint on theta | x
    x_ratio = model.get_endpoints_per_draw(draws=fit['draws'])

    # calculate Prob(H1|x) ie mean over samples theta^s
    tmp = (
        model.get_p_h1(x_ratio, pps_H1_min_effect_size_thresh=pps_H1_min_effect_size_thresh)
        .rename(columns={'p_h1': 'pps_ProbH1_x'})
    )
    tmp['pps_H1_yes'] = (tmp['pps_ProbH1_x'] > pps_ProbH1_target_lwr_quantile).astype(int)

    # calculate 1{theta^s in H1} per draw
    x_ratio = x_ratio.rename(columns={'ratio': 'pps_ratio_x'})
    x_ratio['pps_H1_x'] = (x_ratio['pps_ratio_x'] > pps_H1_min_effect_size_thresh).astype(int)
    x_ratio = x_ratio.merge(
        tmp[['item_label', 'pps_ProbH1_x', 'pps_H1_yes']],
        on='item_label', how='left',
    )
    wa = wa.merge(
        x_ratio[['draw', 'item_label', 'item_type', 'pps_ratio_x', 'pps_H1_x']],
        on=['draw', 'item_label', 'item_type'], how='inner',
    )

    # save data for regression
    pkl_path = os.path.join(dir_out, f"{file_prefix}_i{interim_id}_regression_training.pkl")
    wa.to_pickle(pkl_path)
    print(f"  saved training pkl: {pkl_path}")

    # quantile regression of endpoint on w(x) at tau = 1 - eta_H
    p_h1_xz_interim, perf_interim = fit_interim_regress_endptx_on_wz_with_quantile_regr(
        wa, pps_H1_min_effect_size_thresh=pps_H1_min_effect_size_thresh, pps_ProbH1_target_lwr_quantile=pps_ProbH1_target_lwr_quantile,
    )
    p_h1_xz_interim['interim_id'] = interim_id
    p_h1_xz_interim['interim_date'] = interim_date
    p_h1_xz_interim['interim_month_year'] = interim_month_year
    p_h1_xz_rows.append(p_h1_xz_interim)

    perf_interim['interim_id'] = interim_id
    perf_interim['interim_date'] = interim_date
    perf_rows.append(perf_interim)

    mins_interim_id = (time.time() - t0) / 60.0
    pps_timing_rows.append({'interim_id': interim_id,
                            'mins_interim_id': round(mins_interim_id, 3)})
    print(f"  interim {interim_id} regression-endptx (quantile) done in {mins_interim_id:.2f} min")

if not p_h1_xz_rows:
    raise RuntimeError("[RGEQ-PPS] no interims with missing participants to evaluate.")

# %%

# =============================================================================
# Aggregate
# =============================================================================

dp_h1_xz = pd.concat(p_h1_xz_rows, ignore_index=True)
dp_h1_xz['pps_H1_min_effect_size_thresh'] = pps_H1_min_effect_size_thresh
dp_h1_xz['pps_ProbH1_target_lwr_quantile'] = pps_ProbH1_target_lwr_quantile
dp_h1_xz['S'] = pps_z_total
perf_all = pd.concat(perf_rows, ignore_index=True)
mins_total = (time.time() - t_all0) / 60.0
pps_timing = pd.DataFrame(pps_timing_rows)
pps_timing['mins_total'] = round(mins_total, 3)
print(f"\nAll interims RGEQ-PPS done in {mins_total:.2f} min")

perf_all.to_csv(os.path.join(dir_out, f"{file_prefix}_pps_RGEQ_perf.csv"), index=False)
pkl_path = os.path.join(dir_out, f"{file_prefix}_pps_RGEQ_p_h1_xz.pkl")
dp_h1_xz.to_pickle(pkl_path)
print(f"Saved RGEQ P(H_1 | x, z) samples to: {pkl_path}")
pps_timing.to_csv(os.path.join(dir_out, f"{file_prefix}_pps_RGEQ_timing.csv"), index=False)

# PPS = mean(p_h1_xz > eta_H) = mean(p_h1_xz == 1) since p_h1_xz in {0, 1}.
pps_df = (
    dp_h1_xz.groupby(['interim_id', 'interim_date', 'interim_month_year',
                      'item_label', 'item_type', 'item_high_label'])['p_h1_xz']
    .apply(lambda p: float((p > pps_ProbH1_target_lwr_quantile).mean()))
    .reset_index(name='pps')
)
pps_df['eta'] = pps_ProbH1_target_lwr_quantile
pps_df['S'] = pps_z_total
pps_df.to_csv(os.path.join(dir_out, f"{file_prefix}_pps_RGEQ.csv"), index=False)
print(pps_df.to_string(index=False))

perf = perf_all.merge(di[['interim_id', 'interim_month_year']], on='interim_id', how='left')
perf = perf.merge(pps_timing[['interim_id', 'mins_interim_id']], on='interim_id', how='left')
perf['time_min'] = perf['mins_interim_id']
metric_labels = {'rho': 'rho', 'r2': 'R^2', 'time_min': 'time (min)'}
perf_long = perf.melt(
    id_vars=['interim_id', 'interim_month_year', 'item_label'],
    value_vars=list(metric_labels.keys()),
    var_name='metric', value_name='value',
)
perf_long['metric'] = pd.Categorical(
    perf_long['metric'].map(metric_labels),
    categories=list(metric_labels.values()), ordered=True,
)
perf_long['method'] = 'Regression endpt quantile (Strong-Oakley)'
perf_long.to_pickle(os.path.join(dir_out, f"{file_prefix}_pps_RGEQ_perf_long.pkl"))
print(f"Saved RGEQ performance long-form pkl.")

# %%

# =============================================================================
# Plot: per-interim p(H_1 | x, z) is a 0/1 label under quantile regression,
# so the "distribution" degenerates to a stacked bar of the success fraction.
# We keep the same box-shape output for consistency with the Gaussian pipeline
# but the interquartile range will collapse to a single value.
# =============================================================================

order = (
    di[di['interim_id'].isin(dp_h1_xz['interim_id'].unique())]
    .sort_values('interim_date')['interim_month_year'].tolist()
)
dp_h1_xz['interim_month_year'] = pd.Categorical(
    dp_h1_xz['interim_month_year'], categories=order, ordered=True,
)

box_stats = (
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
box_stats.to_pickle(os.path.join(dir_out, f"{file_prefix}_pps_RGEQ_p_h1_xz_boxplot.pkl"))

p = (
    ggplot(
        box_stats,
        aes(x='interim_month_year',
            ymin='q025', lower='q25', middle='q50',
            upper='q75', ymax='q975',
            group='interim_month_year'),
    )
    + geom_boxplot(stat='identity', fill='#8c564b', alpha=0.4, colour='#808080')
    + geom_hline(yintercept=pps_ProbH1_target_lwr_quantile, colour='black', size=1.0)
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        figure_size=(14, 5),
    )
    + labs(
        x='Interim', y='p(H_1 | x, z) per z draw (0/1 under quantile regr)',
    )
)
tmp = os.path.join(dir_out, f"{file_prefix}_pps_RGEQ_p_h1_xz_boxplot.pdf")
p.save(tmp, verbose=False, limitsize=False)
print(f"Saved RGEQ p(H_1 | x, z) boxplot to: {tmp}")

# %%

# =============================================================================
# Bar chart: PPS per interim_date. Black bar = closed-form analytic PPS;
# #8c564b bar = quantile-regression PPS. Errorbars = 95% bootstrap CI on the
# regression PPS (B resamples with replacement of the S MC samples per
# interim).
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
bs = bs[np.arange(n_interim)[:, None, None], tmp]
bs = (bs > pps_ProbH1_target_lwr_quantile).mean(axis=2)                # (n_interim, B)
bs = pd.DataFrame(np.quantile(bs, quantiles, axis=1).T, columns=quantile_names)
bs['interim_id'] = wide.index.get_level_values('interim_id').to_numpy()
bs['interim_date'] = pd.to_datetime(
    wide.index.get_level_values('interim_date'),
)
bs['method'] = 'Regression endptx quantile'

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
        tmp[['interim_id', 'interim_date', 'pps']].assign(method='Regression endptx quantile'),
    ],
    ignore_index=True,
)
tmp['interim_date'] = pd.to_datetime(tmp['interim_date'])
tmp['method'] = pd.Categorical(
    tmp['method'], categories=['analytic', 'Regression endptx quantile'], ordered=True,
)
tmp = tmp.merge(
    bs[['interim_id', 'method', 'q025', 'q975']],
    on=['method', 'interim_id'], how='left',
)
tmp.to_pickle(os.path.join(dir_out, f"{file_prefix}_pps_RGEQ_bootstrap.pkl"))

_bar_colours = {'analytic': '#000000', 'Regression endptx quantile': '#8c564b'}
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
        x='interim date', y='PPS = int P(p(H_1 | x, z) > eta) dz',
        fill='method',
    )
)
tmp = os.path.join(dir_out, f"{file_prefix}_pps_RGEQ_bars.pdf")
p.save(tmp, verbose=False, limitsize=False)
print(f"Saved RGEQ PPS bar chart to: {tmp}")

# %%

# =============================================================================
# Diagnostic-plot loop: per-interim w_ratio vs pps_ratio_x scatter with
# quantile-regression line at tau = 1 - eta_H + rho / R^2 annotations.
# =============================================================================

print("\nBuilding diagnostic plot (quantile fit on pps_ratio_x), one facet per interim...")
wa_parts = []
for interim_id in di['interim_id']:
    interim_id = int(interim_id)
    wa_pkl = os.path.join(dir_out, f"{file_prefix}_i{interim_id}_regression_training.pkl")
    if not os.path.exists(wa_pkl):
        continue
    wa_i = pd.read_pickle(wa_pkl).copy()
    wa_i['interim_id'] = interim_id
    wa_i['interim_month_year'] = (
        di.loc[di['interim_id'] == interim_id, 'interim_month_year'].iloc[0]
    )
    wa_parts.append(wa_i)
wa_all = pd.concat(wa_parts, ignore_index=True)

order = (
    di[di['interim_id'].isin(wa_all['interim_id'].unique())]
    .sort_values('interim_date')['interim_month_year'].tolist()
)
wa_all['interim_month_year'] = pd.Categorical(
    wa_all['interim_month_year'], categories=order, ordered=True,
)

_ann = (
    wa_all.groupby(['interim_month_year', 'item_label'], observed=True)
    .apply(fit_interim_regress_H1x_on_wz_per_item_summary)
    .reset_index()
)
_ann['rho_label'] = _ann['rho'].apply(lambda r: f"ρ={r:.2f}")
_ann['r2_label'] = _ann['r2'].apply(lambda r: f"R2={100 * r:.1f}%")

tau_used = 1.0 - pps_ProbH1_target_lwr_quantile

# Precompute the quantile-regression line per facet: fit QuantReg on
# (w_ratio, pps_ratio_x) per (interim_month_year, item_label), then evaluate on
# a 50-point grid spanning the observed w_ratio range within that facet.
import statsmodels.api as sm
_line_parts = []
for (imy, item), gg in wa_all.groupby(['interim_month_year', 'item_label'],
                                      observed=True):
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
        'item_label': item,
        'w_ratio': w_grid,
        'pps_ratio_x': y_pred,
    }))
_qline = pd.concat(_line_parts, ignore_index=True) if _line_parts else pd.DataFrame()
if not _qline.empty:
    _qline['interim_month_year'] = pd.Categorical(
        _qline['interim_month_year'], categories=order, ordered=True,
    )

from plotnine import geom_line
p = (
    ggplot(wa_all, aes(x='w_ratio', y='pps_ratio_x'))
    + geom_point(alpha=0.3, size=0.5, colour='#8c564b')
    + geom_line(_qline, aes(x='w_ratio', y='pps_ratio_x'),
                colour='black', size=0.8, inherit_aes=False)
    + geom_hline(yintercept=pps_H1_min_effect_size_thresh, colour='#8c564b', linetype='dashed', size=0.6)
    + geom_text(_ann, aes(x='x_left', y='y_top', label='rho_label'),
                ha='left', va='top', size=8, inherit_aes=False)
    + geom_text(_ann, aes(x='x_right', y='y_top', label='r2_label'),
                ha='right', va='top', size=8, inherit_aes=False)
    + facet_wrap('~ interim_month_year', ncol=4, scales='free')
    + theme_bw()
    + theme(
        figure_size=(14, 10),
        strip_background=element_rect(fill='none', colour='none'),
        panel_background=element_rect(fill='none', colour='none'),
        legend_position='none',
    )
    + labs(x='empirical endpoint summary w(z)',
           y=f'actual endpoint from posterior p(p | x) draws (quantile line at tau={tau_used:.3f})')
)
pdf_path = os.path.join(dir_out, f"{file_prefix}_endpointratio_vs_wratio.pdf")
p.save(pdf_path, verbose=False, limitsize=False)
print(f"  saved combined diagnostic plot: {pdf_path}")

print("\n✓ Binomial regression-endptx-quantile PPS + diagnostic plot complete.")
