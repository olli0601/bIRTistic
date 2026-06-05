#!/usr/bin/env python3
"""
Binomial interim analyses, Python (first-development-phase).

Simulates ``N = 500`` Bernoulli observations with true ``p = 0.4`` equally
spaced over one year, defines monthly interim analyses, and at each
interim fits the Beta-Bernoulli model two ways:
- Stan HMC (cmdstanpy + ``src/stan/binomial_v260603.stan``)
- NumPyro HMC (NUTS on ``src/numpyro/binomial_v260603.pyro``)

Per-interim density plot facets compare the two HMC posteriors against
the analytic Beta(a + k, b + n - k) posterior (grey40 background fill).

PPS comes in a follow-up phase.

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/Binomial_interim_analyses_with_nested_MC.py
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
    ggplot, aes, geom_density, geom_vline, geom_hline, geom_boxplot,
    facet_wrap, scale_color_manual, scale_fill_manual, theme_bw, theme,
    element_text, labs, position_dodge,
)
import warnings
warnings.filterwarnings('ignore')

from model_binomial import BinomialModel
from fit_interim import fit_interim_nested_monte_carlo_of_posterior_xz

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

dir_out = "/Users/or105/sandbox/bIRTistic/py-binomial-interim-260604"
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

# =============================================================================
# Fit Stan HMC + NumPyro HMC + NumPyro SVI (AutoDiagonalNormal,
# AutoLowRankMultivariateNormal) at each interim. resume=True so reruns
# pick up cached posteriors from {dir_out}/binomial_i{interim_id}_*.
# =============================================================================

# (method_label, fit_fn_name, fit_kwargs, file_suffix)
methods = [
    ('stan_HMC', 'fit_stan_hmc',
     dict(chains=4, iter_warmup=500, iter_sampling=1000),
     'stan_hmc'),
    ('pyro_HMC', 'fit_pyro_hmc',
     dict(chains=4, iter_warmup=500, iter_sampling=1000),
     'pyro_hmc'),
    ('pyro_SVI_AutoDiagonalNormal', 'fit_pyro_svi',
     dict(algorithm='AutoDiagonalNormal', lr=0.05, num_steps=4000,
          output_samples=4000),
     'pyro_svi_autodiagnormal'),
    ('pyro_SVI_AutoLowRankMultivariateNormal', 'fit_pyro_svi',
     dict(algorithm='AutoLowRankMultivariateNormal', lr=0.05,
          num_steps=4000, output_samples=4000),
     'pyro_svi_autolowrankmvn'),
]

# %%

po = []   # per-method, per-interim fitted endpoint draws
pa = []    # per-interim Closed_form (analytic) endpoint draws
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

    # Analytic samples from posterior given x 
    fit = model.fit_closed_form_posterior(
        output_samples=8000, verbose=False,
    )
    de = model.get_endpoints_per_draw(draws=fit['draws'])
    de['interim_id'] = interim_id
    de['interim_month_year'] = di.iloc[i]['interim_month_year']
    de['interim_date'] = interim_date
    de['method'] = 'Closed_form'    
    pa.append(de)
    
    # Numeric estimates from posterior given x 
    for m in range(len(methods)):
        method_label, fit_name, fit_kwargs, file_suffix = methods[m]
        prefix = os.path.join(
            dir_out, f"binomial_i{interim_id}_{file_suffix}",
        )
        fit = getattr(model, fit_name)(
            output_file_prefix=prefix,
            save_to_file=True, resume=True, verbose=False,
            **fit_kwargs,
        )
        de = model.get_endpoints_per_draw(draws=fit['draws'])
        de['interim_id'] = interim_id
        de['interim_month_year'] = di.iloc[i]['interim_month_year']
        de['interim_date'] = interim_date
        de['method'] = method_label
        po.append(de)

po = pd.concat(po, ignore_index=True)
pa = pd.concat(pa, ignore_index=True)
print(f"\n  Collected {len(po):,} fitted posterior draws across "
      f"{po['interim_id'].nunique()} interims and "
      f"{po['method'].nunique()} methods + "
      f"{len(pa):,} Closed_form draws.")

# %%

# =============================================================================
# Persist tidy outputs
# =============================================================================

po.to_pickle(os.path.join(dir_out, 'posterior_draws.pkl'))
pa.to_pickle(os.path.join(dir_out, 'analytic_posterior_draws.pkl'))

# %%

# =============================================================================
# Density plot per interim facet:
# grey fill = Closed_form (analytic) posterior;
# colour line for each fitted method.
# =============================================================================

interim_order = di['interim_month_year'].tolist()
po['interim_month_year'] = pd.Categorical(
    po['interim_month_year'],
    categories=interim_order, ordered=True,
)
pa['interim_month_year'] = pd.Categorical(
    pa['interim_month_year'],
    categories=interim_order, ordered=True,
)

method_colours = {
    'stan_HMC': '#1f77b4',
    'pyro_HMC': '#d62728',
    'pyro_SVI_AutoDiagonalNormal': '#2ca02c',
    'pyro_SVI_AutoLowRankMultivariateNormal': '#ff7f0e',
}
_fit_method_labels = list(method_colours.keys())
_all_method_labels = ['Closed_form', *_fit_method_labels]
po['method'] = pd.Categorical(
    po['method'], categories=_fit_method_labels, ordered=True,
)

p = (
    ggplot()
    + geom_density(
        pa, aes(x='p'),
        fill='#666666', colour='#666666', alpha=0.5, size=0,
    )
    + geom_density(
        po, aes(x='p', colour='method'),
        size=0.9,
    )
    + geom_vline(xintercept=TRUE_P, linetype='dashed',
                 colour='black', size=0.4)
    + facet_wrap('~ interim_month_year', ncol=4, scales='free_y')
    + scale_color_manual(values=method_colours)
    + theme_bw()
    + theme(
        figure_size=(14, 10),
        legend_position='top',
        axis_text_x=element_text(angle=0),
    )
    + labs(
        x='p', y='posterior density',
        colour='method',
        title=(f'Binomial posterior across {len(di)} monthly interims '
               f'(true p={TRUE_P}, prior Beta({PRIOR_A:.0f},{PRIOR_B:.0f}))'),
    )
)
pdf_path = os.path.join(dir_out, 'binomial_interim_hmc_vs_analytic.pdf')
p.save(pdf_path, verbose=False, limitsize=False)
print(f"\nSaved per-interim density plot to: {pdf_path}")

#%%
# =============================================================================
# Endpoint boxplot per interim:
# y = endpoint ratio (p / p_0), x = interim;
# boxes = posterior quantiles (2.5/25/50/75/97.5%); dashed line at true ratio.
# =============================================================================

quantiles = [0.025, 0.25, 0.5, 0.75, 0.975]
quantiles_names = ['q025', 'q25', 'q50', 'q75', 'q975']
tmp = pd.concat([pa, po], ignore_index=True)
tmp['method'] = tmp['method'].astype(str)
tmp = (
    tmp
    .groupby(['interim_id', 'interim_month_year', 'interim_date', 'method'], observed=True)
    .apply(lambda g: pd.Series(dict(zip(quantiles_names, g['ratio'].quantile(quantiles).to_list()))))
    .reset_index()
)
tmp['interim_month_year'] = pd.Categorical(
    tmp['interim_month_year'],
    categories=interim_order, ordered=True,
)
tmp['method'] = pd.Categorical(tmp['method'], categories=_all_method_labels, ordered=True)

_box_fill_colours = {'Closed_form': '#666666', **method_colours}
true_ratio = 1 - TRUE_P / P_0

p = (
    ggplot(
        tmp,
        aes(x='interim_month_year', ymin='q025', lower='q25',
            middle='q50', upper='q75', ymax='q975', fill='method'),
    )
    + geom_boxplot(stat='identity', position=position_dodge(width=0.8),
                   width=0.7, colour='black', size=0.3)
    + geom_hline(yintercept=true_ratio, linetype='dashed',
                 colour='black', size=0.4)
    + scale_fill_manual(values=_box_fill_colours)
    + theme_bw()
    + theme(
        figure_size=(14, 6),
        legend_position='top',
        axis_text_x=element_text(angle=45, hjust=1),
    )
    + labs(
        x='interim date', y='endpoint 1- p / p_0', fill='method',
        title=(f'Binomial endpoint posterior across {len(di)} monthly '
               f'interims (true endpoint={true_ratio:.2f})'),
    )
)
tmp = os.path.join(dir_out, 'binomial_interim_endpoint_boxplot.pdf')
p.save(tmp, verbose=False, limitsize=False)
print(f"Saved per-interim endpoint boxplot to: {tmp}")

# %%

# =============================================================================
# PPS via closed-form Beta-Binomial (analytic integration over zi).
# For each interim with m = N - n_obs > 0 future trials, model.fit_closed_form_pps
# returns the analytic Beta-Binomial tail sum: P(P(H_1 | x, z) > eta).
# =============================================================================

pps_H1_def = 0.25
pps_ProbH1_thresh = 0.89
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
        pps_H1_def=pps_H1_def,
        pps_ProbH1_thresh=pps_ProbH1_thresh,
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
        'eta': pps_ProbH1_thresh,
    })

ppsa = pd.DataFrame(ppsa)
ppsa.to_pickle(os.path.join(dir_out, 'pps_closed_form.pkl'))
print(f"\nClosed-form PPS table saved to "
      f"{os.path.join(dir_out, 'pps_closed_form.pkl')}")
print(ppsa.to_string(index=False))

# %%

# =============================================================================
# PPS via nested Monte-Carlo over zi (mirror of
# Ukraine_interim_analyses_with_nested_MC.py lines 537-690).
#
# For each interim with interim_m = n_full - n_obs > 0:
#   1. Fit interim posterior on xi (closed-form Beta(a+k, b+n-k)).
#   2. Draw S=pps_z_total posterior-predictive samples for the interim_m
#      missing trials: p_s ~ posterior(p | xi); k_m_s ~ Binomial(interim_m, p_s).
#   3. For each s, refit the posterior on (xi, zi_s) via NumPyro SVI
#      (AutoDiagonalNormal) and record p(H_1 | x, z_s) = fraction of posterior
#      draws with endpoint 1 - p/p_0 > pps_H1_def (equivalently ratio
#      = p/p_0 < 1 - pps_H1_def).
#   4. PPS = mean over s of 1{p(H_1 | x, z_s) > pps_ProbH1_thresh}.
# =============================================================================

pps_z_total = 200
fitting_output_samples = 4000

p_h1_xz_rows = []
t_all0 = time.time()
for i in range(len(di)):
    interim_id = int(di.iloc[i]['interim_id'])
    interim_date = di.iloc[i]['interim_date']
    interim_month_year = di.iloc[i]['interim_month_year']
    dpi = dp[dp["submission_date"] <= interim_date].copy()
    dpi["oid"] = np.arange(1, len(dpi) + 1, dtype=int)
    n_obs = len(dpi)
    interim_m = n_full - n_obs  # participants missing to reach full data
    if interim_m <= 0:
        continue
    print(f"\n{'=' * 70}\nPPS for interim {interim_id} ({interim_date.date()})\n{'=' * 70}")

    t0 = time.time()
    interim_prefix = os.path.join(dir_out, f"binomial_pps_i{interim_id}")

    # Outer interim fit (SVI) -> zarr with posterior['ypred'] feeding zi.
    model_x = BinomialModel(
        dit=dit, dcati=dpi,
        prior_a=PRIOR_A, prior_b=PRIOR_B, p_0=P_0, seed=seed,
    )
    fit_x = model_x.fit_pyro_svi(
        output_file_prefix=f"{interim_prefix}_x",
        algorithm='AutoDiagonalNormal',
        lr=0.05, num_steps=4000, output_samples=fitting_output_samples,
        save_to_file=True, resume=True, verbose=False,
    )
    zi = model_x.get_interim_z_from_ypredi(
        f"{interim_prefix}_x_draws.zarr", interim_m,
        pps_z_total=pps_z_total, seed=seed,
    )
    sample_df = fit_interim_nested_monte_carlo_of_posterior_xz(
        xi=model_x.dcati,
        zi=zi,
        dit=dit,
        output_file_prefix=interim_prefix,
        fitting_method_args={
            'algorithm': 'AutoDiagonalNormal',
            'x_formula': "~ 1",
            'lr': 0.05,
            'num_steps': 4000,
            'output_samples': fitting_output_samples,
        },
        pps_z_total=pps_z_total,
        pps_H1_def=pps_H1_def,
        pps_ProbH1_thresh=pps_ProbH1_thresh,
        categorical_threshold=2,
        model_cls=BinomialModel,
        seed=seed,
        save_to_file=False,
        verbose=False,
        cpu_n=1,
    )
    sample_df['interim_id'] = interim_id
    sample_df['interim_date'] = interim_date
    sample_df['interim_month_year'] = interim_month_year
    mins_interim_id = (time.time() - t0) / 60.0
    sample_df['interim_mins'] = round(mins_interim_id, 3)
    p_h1_xz_rows.append(sample_df)
    print(f"  interim {interim_id} PPS done in {mins_interim_id:.2f} min")

if not p_h1_xz_rows:
    raise RuntimeError("[PPS] no interims with missing participants to evaluate.")

p_h1_xz_df = pd.concat(p_h1_xz_rows, ignore_index=True)
p_h1_xz_df.to_pickle(os.path.join(dir_out, 'pps_p_h1_xz.pkl'))

# 4. Aggregate: PPS = mean over s of 1{p(H_1 | x, z_s) > eta}.
pps_df = (
    p_h1_xz_df
    .groupby(['interim_id', 'interim_date', 'interim_month_year',
              'item_label', 'item_type', 'item_high_label'])['p_h1_xz']
    .apply(lambda p: float((p > pps_ProbH1_thresh).mean()))
    .reset_index(name='pps')
)
pps_df['eta'] = pps_ProbH1_thresh
pps_df['S'] = pps_z_total
pps_df.to_pickle(os.path.join(dir_out, 'pps_nested_mc.pkl'))
print(f"\nNested-MC PPS table saved to {os.path.join(dir_out, 'pps_nested_mc.pkl')}")
print(pps_df.to_string(index=False))
print(f"\nTotal nested-MC PPS time: {(time.time() - t_all0) / 60.0:.2f} min")

# %%

# =============================================================================
# Boxplot: distribution of p(H_1 | x, z) across the S nested-MC samples per interim
# =============================================================================

order = (
    di[di['interim_id'].isin(p_h1_xz_df['interim_id'].unique())]
    .sort_values('interim_date')['interim_month_year'].tolist()
)
p_h1_xz_df['interim_month_year'] = pd.Categorical(
    p_h1_xz_df['interim_month_year'], categories=order, ordered=True,
)

box_stats = (
    p_h1_xz_df
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
box_stats.to_pickle(os.path.join(dir_out, 'pps_p_h1_xz_boxplot.pkl'))

p = (
    ggplot(
        box_stats,
        aes(x='interim_month_year',
            ymin='q025', lower='q25', middle='q50',
            upper='q75', ymax='q975',
            group='interim_month_year'),
    )
    + geom_boxplot(stat='identity', fill='#1f77b4')
    + geom_hline(yintercept=pps_ProbH1_thresh, colour='black', size=1.0)
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        figure_size=(14, 5),
    )
    + labs(
        x='Interim', y='p(H_1 | x, z)',
        title=(f'Binomial nested-MC p(H_1 | x, z) distribution per interim '
               f'(S = {pps_z_total}, eta = {pps_ProbH1_thresh})'),
    )
)
tmp = os.path.join(dir_out, 'binomial_interim_pps_p_h1_xz_boxplot.pdf')
p.save(tmp, verbose=False, limitsize=False)
print(f"Saved nested-MC p(H_1 | x, z) boxplot to: {tmp}")


