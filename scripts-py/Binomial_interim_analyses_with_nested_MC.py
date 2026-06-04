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
    ggplot, aes, geom_density, geom_vline, facet_wrap, scale_color_manual,
    scale_fill_manual, theme_bw, theme, element_text, labs,
)
from scipy.stats import beta as _scipy_beta
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

# %%

# =============================================================================
# Fit Stan HMC + NumPyro HMC at each interim; collect posterior draws
# =============================================================================

posterior_rows = []   # (interim_id, method, p_draw)
analytic_rows = []    # (interim_id, a_post, b_post, k, n)

for _, row in di.iterrows():
    interim_id = int(row['interim_id'])
    interim_date = row['interim_date']

    cohort = dp[dp['submission_date'] <= interim_date].copy()
    n_obs = len(cohort)
    if n_obs == 0:
        continue

    print(f"\n--- Interim {interim_id} ({interim_date.date()}): n={n_obs} ---")

    model = BinomialModel(
        dit=dit, dcati=cohort,
        prior_a=PRIOR_A, prior_b=PRIOR_B, p_0=P_0, seed=seed,
    )

    k = int(cohort['y'].sum())
    a_post = PRIOR_A + k
    b_post = PRIOR_B + (n_obs - k)
    analytic_rows.append({
        'interim_id': interim_id, 'interim_month_year': row['interim_month_year'],
        'a_post': a_post, 'b_post': b_post, 'k': k, 'n': n_obs,
    })

    # Stan HMC
    stan_fit = model.fit_stan_hmc(
        chains=4, iter_warmup=500, iter_sampling=1000, verbose=False,
    )
    for p_val in stan_fit['posterior_samples']['p']:
        posterior_rows.append({
            'interim_id': interim_id,
            'interim_month_year': row['interim_month_year'],
            'method': 'stan_HMC', 'p': float(p_val),
        })

    # Pyro HMC (NumPyro NUTS)
    pyro_fit = model.fit_pyro_hmc(
        chains=4, iter_warmup=500, iter_sampling=1000, verbose=False,
    )
    for p_val in pyro_fit['posterior_samples']['p']:
        posterior_rows.append({
            'interim_id': interim_id,
            'interim_month_year': row['interim_month_year'],
            'method': 'pyro_HMC', 'p': float(p_val),
        })

posterior_df = pd.DataFrame(posterior_rows)
analytic_df = pd.DataFrame(analytic_rows)
print(f"\n  Collected {len(posterior_df):,} HMC posterior draws across "
      f"{posterior_df['interim_id'].nunique()} interims.")

# %%

# =============================================================================
# Analytic Beta posterior samples (grey background) per interim
# =============================================================================

_n_analytic_draws = 8000
analytic_samples_rows = []
analytic_rng = np.random.default_rng(seed + 1)
for _, ar in analytic_df.iterrows():
    p_samples = analytic_rng.beta(
        ar['a_post'], ar['b_post'], size=_n_analytic_draws,
    )
    for v in p_samples:
        analytic_samples_rows.append({
            'interim_id': ar['interim_id'],
            'interim_month_year': ar['interim_month_year'],
            'p': float(v),
        })
analytic_samples_df = pd.DataFrame(analytic_samples_rows)

# Persist tidy outputs.
posterior_df.to_csv(os.path.join(dir_out, 'posterior_hmc_draws.csv'), index=False)
analytic_df.to_csv(os.path.join(dir_out, 'analytic_posterior_params.csv'), index=False)
analytic_samples_df.to_csv(os.path.join(dir_out, 'analytic_posterior_draws.csv'), index=False)

# %%

# =============================================================================
# Density plot per interim facet:
# grey40 fill = analytic Beta posterior;
# colour line for stan_HMC and pyro_HMC.
# =============================================================================

interim_order = di['interim_month_year'].tolist()
posterior_df['interim_month_year'] = pd.Categorical(
    posterior_df['interim_month_year'],
    categories=interim_order, ordered=True,
)
analytic_samples_df['interim_month_year'] = pd.Categorical(
    analytic_samples_df['interim_month_year'],
    categories=interim_order, ordered=True,
)

method_colours = {'stan_HMC': '#1f77b4', 'pyro_HMC': '#d62728'}

p = (
    ggplot()
    + geom_density(
        analytic_samples_df, aes(x='p'),
        fill='#666666', colour='#666666', alpha=0.5, size=0,
    )
    + geom_density(
        posterior_df, aes(x='p', colour='method'),
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
        colour='HMC method',
        title=(f'Binomial posterior across {len(di)} monthly interims '
               f'(true p={TRUE_P}, prior Beta({PRIOR_A:.0f},{PRIOR_B:.0f}))'),
    )
)
pdf_path = os.path.join(dir_out, 'binomial_interim_hmc_vs_analytic.pdf')
p.save(pdf_path, verbose=False, limitsize=False)
print(f"\nSaved per-interim density plot to: {pdf_path}")

# Sanity print: max abs diff of posterior mean (HMC vs analytic) per interim.
print("\nPer-interim posterior-mean comparison (|HMC_mean - analytic_mean|):")
for _, ar in analytic_df.iterrows():
    analytic_mean = ar['a_post'] / (ar['a_post'] + ar['b_post'])
    for method in ('stan_HMC', 'pyro_HMC'):
        sub = posterior_df[
            (posterior_df['interim_id'] == ar['interim_id'])
            & (posterior_df['method'] == method)
        ]
        hmc_mean = sub['p'].mean()
        diff = abs(hmc_mean - analytic_mean)
        print(f"  interim {int(ar['interim_id']):2d} {method}: "
              f"HMC={hmc_mean:.4f} analytic={analytic_mean:.4f} diff={diff:.4f}")

print("\n✓ Binomial interim analysis complete.")
