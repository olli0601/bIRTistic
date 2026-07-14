#!/usr/bin/env python3
"""
Binomial interim analysis: **amortised** PPS via a DeepSets multi-quantile
network (§11.6 of ``dev/amortised_decision_making.md``).

Prototype outlined in §11:

- Model = Binomial (§3.1), hardcoded sufficient statistic
  ``T(x_i) = x_i`` so the DeepSets sum is just the total number of
  successes.
- Prior-predictive training sampler: ``p ~ Beta(a, b);
  kn ~ Binomial(n, p); km ~ Binomial(m, p); rho = 1 - p / p_0``.
- Multi-quantile head at 11 tau levels (0.05 ... 0.95). Continuous
  ``P(H_1 | x, z)`` obtained by monotone-corrected CDF interpolation at
  ``pps_H1_min_effect_size_thresh``.
- Deployment: at each interim, fit the x-posterior once (pyro HMC,
  ``resume=True``), draw ``pps_z_total`` posterior-predictive z samples,
  and forward-pass features ``(kn, km, n, m)`` through the trained net.

Outputs (``_RGEA_`` suffix; matches the compare-methods loader):
p_h1_xz pkl, PPS csv, box-stats pkl, p_h1_xz boxplot, PPS bars with 95%
bootstrap CI vs analytic, perf long-form pkl (rho / R^2 / time(min)),
timing csv, plus a per-interim diagnostic scatter.

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/Binomial_interim_analysis_amortise_endptx_on_wz_with_features_fixed_qpsi_MLP_loss_multiquantilehead.py
"""

# %%

import os
import sys
import time
from functools import partial
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
from fit_interim import fit_interim_regress_H1x_on_wz_per_item_summary
from amortiser_common import (
    train,
    save_trained_model,
    load_fitted_model,
)
from amortiser_pps_features_fixed_qpsi_MLP_loss_multiquantilehead import (
    Amortiser_PPS_features_fixed_qpsi_MLP_loss_multiquantilehead,
    predict_amortised_p_h1_for_many_xz,
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

# Amortised net training config. Override any of the knobs below by
# setting the environment variable ``AMORTISER_VARIANT`` to one of the
# entries in ``_VARIANTS``; leaves the default deployment untouched.
NET_HIDDEN = (64, 64)               # §12.4 ablation: 64x64 halves MSE
                                    # at 3x-faster training vs
                                    # (256, 256, 128).
NET_TAUS = (0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
NET_STEPS = 40000
NET_BATCH = 8192
NET_LR = 1e-3
NET_SEED = seed
NET_N_MAX = N
N_DIST = 'uniform'

_ROOT_DIR = (
    "/Users/or105/sandbox/bIRTistic/"
    "py-binomial-interim-amortise-endptx-on-wz-with-features-fixed-"
    "qpsi-MLP-loss-multiquantilehead"
)
dir_out = f"{_ROOT_DIR}-260714"          # §12.4 default swap to 64x64.

_VARIANTS = {
    '64x64':                 dict(NET_HIDDEN=(64, 64),
                                  dir_out=f"{_ROOT_DIR}-64x64_260713"),
    'num_quantile_levels_5': dict(
        NET_TAUS=(0.05, 0.25, 0.5, 0.75, 0.95),
        dir_out=f"{_ROOT_DIR}-num_quantile_levels_5_260713",
    ),
    'num_quantile_levels_21': dict(
        NET_TAUS=tuple(round(0.025 + 0.05 * i, 4) for i in range(20)) + (0.975,),
        dir_out=f"{_ROOT_DIR}-num_quantile_levels_21_260713",
    ),
    'S_2000':                dict(pps_z_total=2000,
                                  dir_out=f"{_ROOT_DIR}-S_2000_260713"),
    'log_uniform_n':         dict(N_DIST='log_uniform',
                                  dir_out=f"{_ROOT_DIR}-log_uniform_n_260713"),
    'combo_64x64_qlv5_S2000': dict(
        NET_HIDDEN=(64, 64),
        NET_TAUS=(0.05, 0.25, 0.5, 0.75, 0.95),
        pps_z_total=2000,
        dir_out=f"{_ROOT_DIR}-combo_64x64_qlv5_S2000_260714",
    ),
}

_variant = os.environ.get('AMORTISER_VARIANT', '')
if _variant:
    if _variant not in _VARIANTS:
        raise ValueError(
            f"unknown AMORTISER_VARIANT={_variant!r}; "
            f"available: {list(_VARIANTS)}"
        )
    for k, v in _VARIANTS[_variant].items():
        # Rebind the module-level knob for this run.
        globals()[k] = v
    print(f"[amortiser] variant = {_variant}")

os.makedirs(dir_out, exist_ok=True)

file_prefix = "binomial_interim"
NET_CKPT = os.path.join(dir_out, f"{file_prefix}_amortised_pps_net.pkl")

# Cached cohort + amortiser training data written by
# ``Binomial_interim_analyses_make_sim_data.py``.
SIM_DATA_DIR = (
    "/Users/or105/sandbox/bIRTistic/py-binomial-interim-make-sim-data-260702"
)
COHORT_PKL = os.path.join(SIM_DATA_DIR, 'binomial_sim_cohort.pkl')
TRAINING_PKL = os.path.join(SIM_DATA_DIR, 'binomial_amortiser_training_data.pkl')

print(f"Output dir: {dir_out}")

# %%

# =============================================================================
# Load the cached simulated Bernoulli cohort (dp, dit, di).
# =============================================================================

if not os.path.exists(COHORT_PKL):
    raise FileNotFoundError(
        f"{COHORT_PKL} missing; run "
        "scripts-py/Binomial_interim_analyses_make_sim_data.py first."
    )
cohort = pd.read_pickle(COHORT_PKL)
dp = cohort['dp']
dit = cohort['dit']
di = cohort['di']
assert cohort['N'] == N and cohort['TRUE_P'] == TRUE_P, (
    "cohort config drift; regenerate the sim-data cache."
)
n_full = len(dp)
print(f"  Loaded cohort: N={n_full}, "
      f"observed k={int(dp['y'].sum())}, "
      f"{len(di)} monthly interim rounds.")

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
# Train amortised net once (or load).
# =============================================================================

if os.path.exists(NET_CKPT):
    print(f"\nLoading cached amortised net from {NET_CKPT}")
    fit = load_fitted_model(NET_CKPT)
else:
    print(f"\nTraining amortised net "
          f"({NET_STEPS} steps, batch {NET_BATCH}) ...")

    # Feed batches from the model's training-data sampler. The BinomialModel
    # instance carries prior_a / prior_b / p_0 so the sampler only needs
    # ``n_max``. This keeps the amortiser data-generation model-agnostic.
    _amortiser_model = BinomialModel(
        dit=dit,
        dcati=dp.assign(oid=np.arange(1, N + 1, dtype=int)),
        prior_a=PRIOR_A, prior_b=PRIOR_B, p_0=P_0, seed=seed,
    )
    sample_fn = partial(
        _amortiser_model.make_training_data_with_features,
        n_max=NET_N_MAX, n_dist=N_DIST,
    )
    fit = train(
        sample_fn,
        Amortiser_PPS_features_fixed_qpsi_MLP_loss_multiquantilehead,
        net_kwargs={'hidden_dims': NET_HIDDEN},
        pps_ProbH1_lwr_quantiles_mesh=NET_TAUS,
        num_steps=NET_STEPS,
        batch_size=NET_BATCH,
        lr=NET_LR,
        seed=NET_SEED,
        verbose=True,
    )
    save_trained_model(fit, NET_CKPT)
    print(f"Saved amortised net to {NET_CKPT}")
mins_train = float(fit.get('training_mins', float('nan')))
print(f"Amortised net trained in {mins_train:.2f} min "
      f"(measured inside train loop; excludes JIT compilation, "
      f"data-loading, and load-cache time).")

# %%

# =============================================================================
# Per-interim loop:
#   1. fit x-posterior on xi (pyro_HMC, resume=True)
#   2. build zi keeping posterior-draw order
#   3. compute per-draw features (kn, km, n, m) and forward-pass net
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
    print(f"\n{'='*70}\nAmortised PPS for interim {interim_id} "
          f"({interim_date.date()})\n{'='*70}")
    print(f"  n_obs={n_obs:,} | m={interim_m}")

    t0 = time.time()
    
    # in order to generate future data z
    # fit model on today's data x (for the posterior-predictive z draws +
    # per-draw endpoint frame; the amortised net itself does not use HMC).
    interim_prefix = os.path.join(
        dir_out, f"{file_prefix}_i{interim_id}_pyro_hmc",
    )
    model = BinomialModel(
        dit=dit, dcati=dpi,
        prior_a=PRIOR_A, prior_b=PRIOR_B, p_0=P_0, seed=seed,
    )
    hmc_fit = model.fit_closed_form_posterior(
        output_file_prefix=interim_prefix,
        save_to_file=True, resume=True, verbose=False,
    )

    zi = model.get_interim_z_from_ypredi(
        f"{interim_prefix}_draws.zarr", interim_m,
        pps_z_total=pps_z_total, seed=seed, keep_order=True,
    )
    wa = model.get_w(zi)

    # calculate endpoint on theta | x (needed for the diagnostic perf table).
    x_ratio = model.get_endpoints_per_draw(draws=hmc_fit['draws'])
    x_ratio = x_ratio.rename(columns={'ratio': 'pps_ratio_x'})
    x_ratio['pps_H1_x'] = (
        x_ratio['pps_ratio_x'] > pps_H1_min_effect_size_thresh
    ).astype(int)
    wa = wa.merge(
        x_ratio[['draw', 'item_label', 'item_type',
                 'pps_ratio_x', 'pps_H1_x']],
        on=['draw', 'item_label', 'item_type'], how='inner',
    )

    # ---- Amortised feature construction (hardcoded T(x_i) = x_i) ----
    # ``zi`` was built above from ``model.fit_closed_form_posterior``, so
    # the ``ypred`` columns already contain analytic (joint) posterior-
    # predictive samples ``km_s = sum_i ypred_{s, i}``. ``wa['w']`` is the
    # mean over rows, so ``km_s = w_s * m``.
    kn = int(dpi['y'].sum())
    n_val = int(n_obs)
    m_val = int(interim_m)
    wa['km'] = np.rint(wa['w'].to_numpy() * m_val).astype(np.int32)
    # DeepSets pooling for hardcoded T(x_i) = x_i: sum the sufficient
    # statistics of x and z before feeding the head. Feature scaling
    # divides by NET_N_MAX so the training distribution is matched.
    scale = float(NET_N_MAX)
    wa['k_total'] = (kn + wa['km']) / scale
    wa['n_total'] = (n_val + m_val) / scale

    p_h1_xz_interim, perf_interim = predict_amortised_p_h1_for_many_xz(
        wa, fit,
        pps_H1_min_effect_size_thresh=pps_H1_min_effect_size_thresh,
        pps_ProbH1_target_lwr_quantile=pps_ProbH1_target_lwr_quantile,
        feature_cols=('k_total', 'n_total'),
    )
    p_h1_xz_interim['interim_id'] = interim_id
    p_h1_xz_interim['interim_date'] = interim_date
    p_h1_xz_interim['interim_month_year'] = interim_month_year
    p_h1_xz_rows.append(p_h1_xz_interim)

    perf_interim['interim_id'] = interim_id
    perf_interim['interim_date'] = interim_date
    perf_rows.append(perf_interim)

    # save data for downstream diagnostics.
    pkl_path = os.path.join(
        dir_out, f"{file_prefix}_i{interim_id}_regression_training.pkl",
    )
    wa.to_pickle(pkl_path)

    mins_interim_id = (time.time() - t0) / 60.0
    pps_timing_rows.append({
        'interim_id': interim_id,
        'mins_interim_id': round(mins_interim_id, 3),
    })
    print(f"  interim {interim_id} amortised PPS done in "
          f"{mins_interim_id:.2f} min")

if not p_h1_xz_rows:
    raise RuntimeError("[RGEA-PPS] no interims with missing participants to evaluate.")

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
print(f"\nAll interims RGEA-PPS done in {mins_total:.2f} min")

perf_all.to_csv(
    os.path.join(dir_out, f"{file_prefix}_pps_RGEA_perf.csv"), index=False,
)
pkl_path = os.path.join(dir_out, f"{file_prefix}_pps_RGEA_p_h1_xz.pkl")
dp_h1_xz.to_pickle(pkl_path)
print(f"Saved RGEA P(H_1 | x, z) samples to: {pkl_path}")
pps_timing.to_csv(
    os.path.join(dir_out, f"{file_prefix}_pps_RGEA_timing.csv"), index=False,
)

# PPS = mean(p_h1_xz > eta_H). p_h1_xz is continuous in [0, 1].
pps_df = (
    dp_h1_xz.groupby(['interim_id', 'interim_date', 'interim_month_year',
                      'item_label', 'item_type', 'item_high_label'])['p_h1_xz']
    .apply(lambda p: float((p > pps_ProbH1_target_lwr_quantile).mean()))
    .reset_index(name='pps')
)
pps_df['eta'] = pps_ProbH1_target_lwr_quantile
pps_df['S'] = pps_z_total
pps_df.to_csv(os.path.join(dir_out, f"{file_prefix}_pps_RGEA.csv"), index=False)
print(pps_df.to_string(index=False))

# Merge closed-form for comparison print.
tab = pps_df[['interim_id', 'interim_month_year', 'pps']].rename(
    columns={'pps': 'pps_amortised'},
).merge(
    ppsa[['interim_id', 'pps']].rename(columns={'pps': 'pps_analytic'}),
    on='interim_id', how='left',
)
tab['abs_err'] = (tab['pps_amortised'] - tab['pps_analytic']).abs()
print("\nAmortised vs analytic PPS (per interim):")
print(tab.to_string(index=False))
print(f"\nMax abs error: {tab['abs_err'].max():.4f}, mean: {tab['abs_err'].mean():.4f}")

# %%

# =============================================================================
# Simple perf long-form pkl for compare-methods integration.
# =============================================================================

perf = perf_all.merge(
    di[['interim_id', 'interim_month_year']], on='interim_id', how='left',
)
perf = perf.merge(
    pps_timing[['interim_id', 'mins_interim_id']], on='interim_id', how='left',
)
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
perf_long['method'] = 'Regression endpt amortised (Strong-Oakley)'
perf_long.to_pickle(
    os.path.join(dir_out, f"{file_prefix}_pps_RGEA_perf_long.pkl"),
)
print("Saved RGEA performance long-form pkl.")

# %%

# =============================================================================
# Plot: per-interim p(H_1 | x, z) boxplot (q025/q25/q50/q75/q975 over the
# S Monte-Carlo samples).
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
box_stats.to_pickle(
    os.path.join(dir_out, f"{file_prefix}_pps_RGEA_p_h1_xz_boxplot.pkl"),
)

p = (
    ggplot(
        box_stats,
        aes(x='interim_month_year',
            ymin='q025', lower='q25', middle='q50',
            upper='q75', ymax='q975',
            group='interim_month_year'),
    )
    + geom_boxplot(stat='identity', fill='#17becf', alpha=0.5,
                   colour='#808080')
    + geom_hline(yintercept=pps_ProbH1_target_lwr_quantile, colour='black', size=1.0)
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        figure_size=(14, 5),
    )
    + labs(x='Interim', y='p(H_1 | x, z) per z draw (amortised, continuous)')
)
tmp = os.path.join(dir_out, f"{file_prefix}_pps_RGEA_p_h1_xz_boxplot.pdf")
p.save(tmp, verbose=False, limitsize=False)
print(f"Saved RGEA p(H_1 | x, z) boxplot to: {tmp}")

# %%

# =============================================================================
# Bar chart: PPS per interim_date. Black = analytic, #17becf = amortised.
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
bs = wide.to_numpy()
n_interim, S = bs.shape
rng2 = np.random.default_rng(seed)
tmp = rng2.integers(0, S, size=(n_interim, bs_B, S))
bs = bs[np.arange(n_interim)[:, None, None], tmp]
bs = (bs > pps_ProbH1_target_lwr_quantile).mean(axis=2)
bs = pd.DataFrame(
    np.quantile(bs, quantiles, axis=1).T, columns=quantile_names,
)
bs['interim_id'] = wide.index.get_level_values('interim_id').to_numpy()
bs['interim_date'] = pd.to_datetime(
    wide.index.get_level_values('interim_date'),
)
bs['method'] = 'Regression endptx amortised'

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
        tmp[['interim_id', 'interim_date', 'pps']].assign(
            method='Regression endptx amortised',
        ),
    ],
    ignore_index=True,
)
tmp['interim_date'] = pd.to_datetime(tmp['interim_date'])
tmp['method'] = pd.Categorical(
    tmp['method'],
    categories=['analytic', 'Regression endptx amortised'], ordered=True,
)
tmp = tmp.merge(
    bs[['interim_id', 'method', 'q025', 'q975']],
    on=['method', 'interim_id'], how='left',
)
tmp.to_pickle(
    os.path.join(dir_out, f"{file_prefix}_pps_RGEA_bootstrap.pkl"),
)

_bar_colours = {'analytic': '#000000', 'Regression endptx amortised': '#17becf'}
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
tmp = os.path.join(dir_out, f"{file_prefix}_pps_RGEA_bars.pdf")
p.save(tmp, verbose=False, limitsize=False)
print(f"Saved RGEA PPS bar chart to: {tmp}")

print("\n✓ Binomial regression-endptx-amortised PPS complete.")
