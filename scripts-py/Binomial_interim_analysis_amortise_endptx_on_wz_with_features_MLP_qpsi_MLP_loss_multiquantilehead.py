#!/usr/bin/env python3
"""
Binomial interim analysis: **amortised** PPS via a DeepSets multi-quantile
network with a **learnable** per-item encoder ``q_tau``
(:class:`Amortiser_PPS_features_MLP_qpsi_MLP_loss_multiquantilehead`).

Parallel to
``Binomial_interim_analysis_amortise_endptx_on_wz_with_features_fixed_qpsi_MLP_loss_multiquantilehead.py``
but with ``q_tau`` learned from raw padded 0/1 item sequences instead of
hardcoded to the identity ``T(x_i) = x_i``. Both scripts share the
training-data cache written by
``Binomial_interim_analyses_make_sim_data.py`` and use the same
``train`` / ``predict_amortised_p_h1_for_one_xz`` / ``save/load``
utilities from :mod:`amortiser_common`.

Outputs (``_RGEB_`` suffix; matches the compare-methods loader).

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/Binomial_interim_analysis_amortise_endptx_on_wz_with_features_MLP_qpsi_MLP_loss_multiquantilehead.py
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
    ggplot, aes, geom_boxplot, geom_col, geom_errorbar, geom_hline,
    scale_fill_manual, scale_x_datetime,
    position_dodge, theme_bw, theme, element_text, labs,
)
import warnings
warnings.filterwarnings('ignore')

from model_binomial import BinomialModel
from amortiser_common import (
    train,
    save_trained_model,
    load_fitted_model,
    predict_amortised_p_h1_for_one_xz,
)
from amortiser_pps_features_MLP_qpsi_MLP_loss_multiquantilehead import (
    Amortiser_PPS_features_MLP_qpsi_MLP_loss_multiquantilehead,
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
PRIOR_A = 1.0
PRIOR_B = 1.0
P_0 = 0.5

pps_H1_min_effect_size_thresh = 0.25
pps_ProbH1_target_lwr_quantile = 0.89
pps_z_total = 200

# Amortised net training config.
NET_Q_TAU_HIDDEN = (32, 32)
NET_EMBED_DIM = 16
NET_HIDDEN = (128, 128, 64)
NET_TAUS = (0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
NET_STEPS = 15000
NET_BATCH = 512            # smaller than fixed variant because
                           # per-sample tensor is (N_max, item_dim).
NET_LR = 1e-3
NET_SEED = seed
NET_N_MAX = N

dir_out = "/Users/or105/sandbox/bIRTistic/py-binomial-interim-amortise-endptx-on-wz-with-features-MLP-qpsi-MLP-loss-multiquantilehead-260711"
os.makedirs(dir_out, exist_ok=True)

file_prefix = "binomial_interim"
NET_CKPT = os.path.join(dir_out, f"{file_prefix}_amortised_pps_net.pkl")

SIM_DATA_DIR = (
    "/Users/or105/sandbox/bIRTistic/py-binomial-interim-make-sim-data-260702"
)
COHORT_PKL = os.path.join(SIM_DATA_DIR, 'binomial_sim_cohort.pkl')

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
    print(f"\nTraining features-MLP amortised net "
          f"({NET_STEPS} steps, batch {NET_BATCH}) ...")

    _amortiser_model = BinomialModel(
        dit=dit,
        dcati=dp.assign(oid=np.arange(1, N + 1, dtype=int)),
        prior_a=PRIOR_A, prior_b=PRIOR_B, p_0=P_0, seed=seed,
    )
    sample_fn = partial(
        _amortiser_model.make_training_data_with_raw_sequences,
        n_max=NET_N_MAX,
    )
    fit = train(
        sample_fn,
        Amortiser_PPS_features_MLP_qpsi_MLP_loss_multiquantilehead,
        net_kwargs={
            'q_tau_hidden_dims': NET_Q_TAU_HIDDEN,
            'embed_dim': NET_EMBED_DIM,
            'hidden_dims': NET_HIDDEN,
        },
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
      f"(measured inside train loop; excludes JIT + data-loading).")

# %%

# =============================================================================
# Per-interim loop:
#   1. fit x-posterior on xi (analytic closed-form)
#   2. build zi (analytic joint draws)
#   3. build padded raw-sequence batch (x, mask_x, z, mask_z, sizes)
#   4. forward-pass net
# =============================================================================


def _padded_bernoulli(k_ones: int, total: int, n_max: int):
    """Return ``x`` shape ``(n_max, 1)`` and ``mask`` shape ``(n_max,)``
    with ``k_ones`` real 1s + ``total - k_ones`` real 0s + zero padding
    to ``n_max``."""
    x = np.zeros((n_max, 1), dtype=np.float32)
    mask = np.zeros((n_max,), dtype=np.float32)
    x[:k_ones, 0] = 1.0
    mask[:total] = 1.0
    return x, mask


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
        print(f"\n[interim {interim_id}] skipped: no missing participants")
        continue
    print(f"\n{'='*70}\nAmortised (features-MLP) PPS for interim {interim_id} "
          f"({interim_date.date()})\n{'='*70}")
    print(f"  n_obs={n_obs:,} | m={interim_m}")

    t0 = time.time()

    # Analytic joint posterior + zi (matches all other Binomial deployment
    # scripts after the fit_closed_form_posterior swap).
    interim_prefix = os.path.join(
        dir_out, f"{file_prefix}_i{interim_id}_analytic",
    )
    model = BinomialModel(
        dit=dit, dcati=dpi,
        prior_a=PRIOR_A, prior_b=PRIOR_B, p_0=P_0, seed=seed,
    )
    model.fit_closed_form_posterior(
        output_file_prefix=interim_prefix,
        save_to_file=True, resume=True, verbose=False,
    )
    zi = model.get_interim_z_from_ypredi(
        f"{interim_prefix}_draws.zarr", interim_m,
        pps_z_total=pps_z_total, seed=seed, keep_order=True,
    )
    ypred_cols = sorted(
        [c for c in zi.columns if c.startswith('ypred_')],
        key=lambda c: int(c.split('_')[1]),
    )
    km_per_draw = zi[ypred_cols].sum(axis=0).to_numpy().astype(np.int32)
    kn = int(dpi['y'].sum())

    # ---- Build the raw-sequences batch (S entries, one per posterior draw).
    # All entries share the same x (observed cohort) but differ in z (per-s
    # km draws from the joint analytic posterior-predictive).
    S = pps_z_total
    x_arr = np.zeros((S, NET_N_MAX, 1), dtype=np.float32)
    mask_x_arr = np.zeros((S, NET_N_MAX), dtype=np.float32)
    z_arr = np.zeros((S, NET_N_MAX, 1), dtype=np.float32)
    mask_z_arr = np.zeros((S, NET_N_MAX), dtype=np.float32)
    sizes_arr = np.tile(
        np.array([n_obs / NET_N_MAX, interim_m / NET_N_MAX], dtype=np.float32),
        (S, 1),
    )
    x_one, x_mask = _padded_bernoulli(kn, n_obs, NET_N_MAX)
    x_arr[:] = x_one[None, :, :]                                # broadcast
    mask_x_arr[:] = x_mask[None, :]
    for s in range(S):
        z_one, z_mask = _padded_bernoulli(int(km_per_draw[s]), interim_m, NET_N_MAX)
        z_arr[s] = z_one
        mask_z_arr[s] = z_mask
    batch = {
        'x':       x_arr,
        'mask_x':  mask_x_arr,
        'z':       z_arr,
        'mask_z':  mask_z_arr,
        'sizes':   sizes_arr,
    }

    p_h1_xz, q_hat, _ = predict_amortised_p_h1_for_one_xz(
        fit, batch, pps_H1_min_effect_size_thresh,
    )
    p_h1_xz_df = pd.DataFrame({
        'item_label':      [dit.iloc[0]['item_label']] * S,
        'item_type':       [dit.iloc[0]['item_type']] * S,
        'item_high_label': [dit.iloc[0]['item_high_label']] * S,
        'p_h1_xz':         p_h1_xz.astype(np.float64),
        'q_hat':           q_hat.astype(np.float64),
        's':               np.arange(1, S + 1, dtype=int),
        'interim_id':      interim_id,
        'interim_date':    interim_date,
        'interim_month_year': interim_month_year,
    })
    p_h1_xz_rows.append(p_h1_xz_df)

    perf_rows.append({
        'item_label':      dit.iloc[0]['item_label'],
        'item_type':       dit.iloc[0]['item_type'],
        'rho':             float('nan'),
        'r2':              float('nan'),
        'interim_id':      interim_id,
        'interim_date':    interim_date,
    })

    mins_interim_id = (time.time() - t0) / 60.0
    pps_timing_rows.append({
        'interim_id': interim_id,
        'mins_interim_id': round(mins_interim_id, 3),
    })
    print(f"  interim {interim_id} amortised (MLP) PPS done in "
          f"{mins_interim_id:.2f} min")

if not p_h1_xz_rows:
    raise RuntimeError("[RGEB-PPS] no interims with missing participants.")

# %%

# =============================================================================
# Aggregate
# =============================================================================

dp_h1_xz = pd.concat(p_h1_xz_rows, ignore_index=True)
dp_h1_xz['pps_H1_min_effect_size_thresh'] = pps_H1_min_effect_size_thresh
dp_h1_xz['pps_ProbH1_target_lwr_quantile'] = pps_ProbH1_target_lwr_quantile
dp_h1_xz['S'] = pps_z_total
perf_all = pd.DataFrame(perf_rows)
mins_total = (time.time() - t_all0) / 60.0
pps_timing = pd.DataFrame(pps_timing_rows)
pps_timing['mins_total'] = round(mins_total, 3)
print(f"\nAll interims RGEB-PPS done in {mins_total:.2f} min")

perf_all.to_csv(
    os.path.join(dir_out, f"{file_prefix}_pps_RGEB_perf.csv"), index=False,
)
pkl_path = os.path.join(dir_out, f"{file_prefix}_pps_RGEB_p_h1_xz.pkl")
dp_h1_xz.to_pickle(pkl_path)
print(f"Saved RGEB P(H_1 | x, z) samples to: {pkl_path}")
pps_timing.to_csv(
    os.path.join(dir_out, f"{file_prefix}_pps_RGEB_timing.csv"), index=False,
)

pps_df = (
    dp_h1_xz.groupby(['interim_id', 'interim_date', 'interim_month_year',
                      'item_label', 'item_type', 'item_high_label'])['p_h1_xz']
    .apply(lambda p: float((p > pps_ProbH1_target_lwr_quantile).mean()))
    .reset_index(name='pps')
)
pps_df['eta'] = pps_ProbH1_target_lwr_quantile
pps_df['S'] = pps_z_total
pps_df.to_csv(os.path.join(dir_out, f"{file_prefix}_pps_RGEB.csv"), index=False)
print(pps_df.to_string(index=False))

tab = pps_df[['interim_id', 'interim_month_year', 'pps']].rename(
    columns={'pps': 'pps_amortised_MLP'},
).merge(
    ppsa[['interim_id', 'pps']].rename(columns={'pps': 'pps_analytic'}),
    on='interim_id', how='left',
)
tab['abs_err'] = (tab['pps_amortised_MLP'] - tab['pps_analytic']).abs()
print("\nAmortised (features-MLP) vs analytic PPS (per interim):")
print(tab.to_string(index=False))
print(f"\nMax abs error: {tab['abs_err'].max():.4f}, mean: {tab['abs_err'].mean():.4f}")

# %%

# =============================================================================
# Perf long-form pkl for compare-methods integration.
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
perf_long['method'] = 'Regression endpt amortised MLP (Strong-Oakley)'
perf_long.to_pickle(
    os.path.join(dir_out, f"{file_prefix}_pps_RGEB_perf_long.pkl"),
)
print("Saved RGEB performance long-form pkl.")

# %%

# =============================================================================
# Boxplot + PPS bar plots.
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
    os.path.join(dir_out, f"{file_prefix}_pps_RGEB_p_h1_xz_boxplot.pkl"),
)

p = (
    ggplot(
        box_stats,
        aes(x='interim_month_year',
            ymin='q025', lower='q25', middle='q50',
            upper='q75', ymax='q975',
            group='interim_month_year'),
    )
    + geom_boxplot(stat='identity', fill='#bcbd22', alpha=0.5,
                   colour='#808080')
    + geom_hline(yintercept=pps_ProbH1_target_lwr_quantile,
                 colour='black', size=1.0)
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        figure_size=(14, 5),
    )
    + labs(x='Interim',
           y='p(H_1 | x, z) per z draw (amortised MLP, continuous)')
)
tmp = os.path.join(dir_out, f"{file_prefix}_pps_RGEB_p_h1_xz_boxplot.pdf")
p.save(tmp, verbose=False, limitsize=False)
print(f"Saved RGEB p(H_1 | x, z) boxplot to: {tmp}")

# %%

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
bs['method'] = 'Regression endptx amortised MLP'

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
            method='Regression endptx amortised MLP',
        ),
    ],
    ignore_index=True,
)
tmp['interim_date'] = pd.to_datetime(tmp['interim_date'])
tmp['method'] = pd.Categorical(
    tmp['method'],
    categories=['analytic', 'Regression endptx amortised MLP'], ordered=True,
)
tmp = tmp.merge(
    bs[['interim_id', 'method', 'q025', 'q975']],
    on=['method', 'interim_id'], how='left',
)
tmp.to_pickle(
    os.path.join(dir_out, f"{file_prefix}_pps_RGEB_bootstrap.pkl"),
)

_bar_colours = {'analytic': '#000000',
                'Regression endptx amortised MLP': '#bcbd22'}
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
tmp = os.path.join(dir_out, f"{file_prefix}_pps_RGEB_bars.pdf")
p.save(tmp, verbose=False, limitsize=False)
print(f"Saved RGEB PPS bar chart to: {tmp}")

print("\n✓ Binomial features-MLP amortised PPS complete.")
