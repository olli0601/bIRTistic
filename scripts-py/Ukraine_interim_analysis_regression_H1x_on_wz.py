#!/usr/bin/env python3
"""
Ukraine interim analysis: regression-based EVSI labels (Strong & Oakley, 2014).

For each interim, draw S=pps_z_total future-data samples and build a per-item
training set for the regression-based label estimator:

  - draw  theta^(s) ~ p(theta | x)               (x-fit posterior draws)
  - draw  z^(s) ~ p(z | theta^(s))               (posterior-predictive ypred)
  - y^(s) := 1{ ratio_item(theta^(s)) > pps_H1_min_effect_size_thresh }   (the H_1 indicator on
            theta^(s); empirically based on x only -- the z-dependence enters
            via the regression covariate W(z) below)
  - W^(s) := per-item summary of z^(s) at baseline (t=0) and endline (t=1):
      * out-of-7  items  : w_t(z) = m^{-1} sum_i z_{it}     (mean response)
      * categorical items: w_t(z) = m^{-1} sum_i 1{ z_{it} >= categorical_threshold }
    plus the direction-aware w_diff and w_ratio (same convention as the endpoint
    ratio: positive means improvement).

The per-item Strong-Oakley label is then estimated by a per-item binomial GLM
``pps_H1_x ~ w_ratio`` (logit link), and the PPS is the Monte-Carlo average
``S^-1 sum_s 1{ pi_hat(W(z^(s))) > pps_ProbH1_target_lwr_quantile }``.

Outputs mirror the other interim scripts (``_RG`` suffix): p_h1_xz pkl, PPS csv,
box-stats pkl, p_h1_xz boxplot, perf long-form pkl (rho / R^2 / time(min)),
timing csv, plus a per-interim diagnostic scatter (w_ratio vs pps_ratio_x with
rho + R^2 annotations).

----------------------------------------------------------------------
Finding: the W summary's signal-to-noise scales with m

For each future sample s and each (item, time) we have m draws from a Bernoulli
or ordinal distribution; the per-(item, time) mean has Var = O(1/m), so as m
shrinks the empirical W(z^(s)) wanders far from its theta-conditional
expectation E[W | theta^(s)]. The direction-aware ratio
``w_ratio = 1 - w_endline / w_baseline`` (lower-is-better) blows up further when
w_baseline approaches 0 with small m. As a result the regression's Pearson
correlation rho(W, pps_ratio_x) collapses for small m, even with the draw-index
alignment fix in place:

| interim | m   | mean |rho| over items | median |rho| |
|---------|-----|------------------------|----------------|
| 1       | 455 | 0.72                   | 0.83           |
| 2       | 329 | 0.76                   | 0.83           |
| 3       | 174 | 0.68                   | 0.75           |
| 4       | 87  | 0.55                   | 0.62           |
| 5       | 30  | 0.40                   | 0.45           |
| 6       | 19  | 0.32                   | 0.36           |
| 7       | 12  | 0.27                   | 0.30           |
| 8       | 2   | 0.09                   | 0.10           |

Conclusion: the Strong-Oakley regression approach becomes prohibitively noisy at
m << 100 because W(z) is dominated by sampling noise; richer summaries (more
sufficient statistics, log/logit scaling to avoid the divide-by-near-zero
explosion) only help to a point. The approach is well-suited to early-to-mid
interims where m is large; for the latest interims one is better off with the
nested-MC / SMC label.

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/Ukraine_interim_analysis_regression.py
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
import statsmodels.api as sm
from plotnine import (
    ggplot, aes, geom_boxplot, geom_hline, geom_jitter,
    geom_line, geom_text, facet_wrap,
    scale_colour_manual, scale_fill_manual, scale_y_continuous,
    theme_bw, theme, element_text, element_rect, labs,
)
import warnings
warnings.filterwarnings('ignore')

from data_loading import read_data_ukraine
from model_pcm import PartialCreditModel
from fit_interim import (
    fit_interim_regress_H1x_on_wz,
)
from utils import _futurama_palette

print("✓ Imports successful")

# %%

# =============================================================================
# Configuration
# =============================================================================

np.random.seed(42)
seed = 123

dir_data = "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/2025_project_Hope_Groups/data"
file_data = os.path.join(dir_data, "Ukraine_Hope_Groups_Baseline_Endline_Wide_Aug6.csv")
dir_out = "/Users/or105/sandbox/bIRTistic/py-ukraine-interim-with-regression-on-H1x-wz-260601"
os.makedirs(dir_out, exist_ok=True)

file_prefix = "pcm_1_interim"
svi_algorithm = 'AutoLowRankMultivariateNormal'
x_formula = "~ time - 1"

pps_z_total = 4000                   # one z-sample per x-posterior draw
pps_H1_min_effect_size_thresh = 0.5
pps_ProbH1_target_lwr_quantile = 0.89
categorical_threshold = 2


print(f"Output dir: {dir_out}")

# %%

# =============================================================================
# Load + preprocess data
# =============================================================================

print("\nLoading Ukraine data...")
raw = read_data_ukraine(file_data)
dp = raw['dp'].copy()
dit = raw['dit'].copy()
dmeta = raw['dmeta'].copy()

tmp = (
    dmeta[['pid', 'time_label', 'displacement_status']]
    .drop_duplicates()
    .rename(columns={'pid': 'pid_label'})
    .dropna(subset=["displacement_status"], how='all')
)
dp = dp.merge(tmp, on=['pid_label', 'time_label'], how='inner', validate='many_to_one')

dp1 = dp[~dp['item_label'].str.contains('agg')].copy()
dp1['y_stan'] = dp1['y'] + 1
dp1 = dp1.merge(dit[['item_label', 'item_type']], on='item_label', how='left')

item_time_df = (
    dp1[['item_type', 'item_label', 'time']]
    .drop_duplicates()
    .sort_values(['item_type', 'time', 'item_label'])
    .reset_index(drop=True)
)
item_time_df['item_time_id'] = item_time_df.groupby('item_type').cumcount() + 1
dp1 = dp1.merge(item_time_df, on=['item_label', 'time', 'item_type'], how='left')
dp1 = dp1.merge(
    dit[['item_type', 'item_type_id']].drop_duplicates(), on='item_type', how='left',
)
dp1 = dp1.sort_values(['item_type_id', 'pid', 'time', 'item_label']).reset_index(drop=True)
dp1['oid'] = range(1, len(dp1) + 1)
dp1['oidt'] = dp1.groupby('item_type').cumcount() + 1

n_full = dp1['pid'].nunique()

_endline_dates = pd.to_datetime(
    dp1.loc[dp1['time_label'] == 'Endline', 'submission_date']
).dropna()
_start = _endline_dates.min().replace(day=1)
_end = (_endline_dates.max() + pd.offsets.MonthEnd(0)).normalize()
month_starts = pd.date_range(_start, _end, freq='MS')
di = pd.DataFrame({
    'interim_id': range(1, len(month_starts) + 1),
    'month_start': month_starts,
    'interim_date': month_starts + pd.offsets.MonthEnd(0),
})
di['interim_month_year'] = di['month_start'].dt.strftime('%Y-%b')
print(f"  n_full={n_full} | looping over {len(di)} interim rounds")
# %%

# =============================================================================
# Per-interim loop: fit (resume) -> z block -> per-item summaries W(z^s)_t,
# per-item logistic GLM (predict P(H_1 | x, z^(s)) by Strong-Oakley regression),
# per-item rho + Gaussian-GLM R^2 on (w_ratio, pps_ratio_x). The diagnostic plot
# is generated in a separate loop below, after all per-interim wa pkls are saved.
# =============================================================================

p_h1_xz_rows = []
perf_rows = []
pps_timing_rows = []
t_all0 = time.time()
for interim_id in di['interim_id']:
    interim_id = int(interim_id)
    interim_date = pd.to_datetime(di.loc[di['interim_id'] == interim_id, 'interim_date'].iloc[0])

    xi = PartialCreditModel.get_interim_data_x(dp1, interim_date)
    if xi.empty or xi['pid'].nunique() < 2:
        print(f"\n[interim {interim_id}] skipped: empty xi"); continue
    interim_m = n_full - xi['pid'].nunique()
    if interim_m <= 0:
        print(f"\n[interim {interim_id}] skipped: no missing participants"); continue
    print(f"\n{'='*70}\nRegression diagnostics for interim {interim_id} ({interim_date.date()})\n{'='*70}")
    print(f"  n_obs={len(xi):,} | n_pid={xi['pid'].nunique()} | m={interim_m}")

    t0 = time.time()

    # fit model on today's data x
    interim_prefix = os.path.join(dir_out, f"{file_prefix}_{interim_id}")
    model = PartialCreditModel(dit=dit, dcati=xi,
                                    x_formula=x_formula, seed=seed)
    fit = model.fit_pyro_svi(
        output_file_prefix=interim_prefix,
        algorithm=svi_algorithm,
        lr=0.01, num_steps=10000, output_samples=pps_z_total,
        resume=True,
        with_core_analyses=True, with_additional_analyses=False, verbose=False,
    )

    # expand to future data keeping draw index linked to posterior draws, and get w(z^s)
    zi = model.get_interim_z_from_ypredi(
        f"{interim_prefix}_draws.zarr", interim_m,
        pps_z_total=pps_z_total, seed=seed, keep_order=True,
    )
    wa = model.get_w(zi, categorical_threshold=categorical_threshold)

    # calculate endpoint on theta | x
    x_ratio = model.get_endpoints_per_draw(
        draws=fit['draws'], categorical_threshold=categorical_threshold,
        endpoint_type='items',
    )

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

    # Save per-interim training pkl (reloaded by the diagnostic-plot loop below).
    pkl_path = os.path.join(dir_out, f"{file_prefix}_i{interim_id}_regression_training.pkl")
    wa.to_pickle(pkl_path)
    print(f"  saved training pkl: {pkl_path}")

    p_h1_xz_interim, perf_interim = fit_interim_regress_H1x_on_wz(wa)
    p_h1_xz_interim['interim_id'] = interim_id
    p_h1_xz_interim['interim_date'] = interim_date
    p_h1_xz_rows.append(p_h1_xz_interim)

    perf_interim['interim_id'] = interim_id
    perf_interim['interim_date'] = interim_date
    perf_rows.append(perf_interim)

    mins_interim_id = (time.time() - t0) / 60.0
    pps_timing_rows.append({'interim_id': interim_id,
                            'mins_interim_id': round(mins_interim_id, 3)})
    print(f"  interim {interim_id} regression done in {mins_interim_id:.2f} min")

if not p_h1_xz_rows:
    raise RuntimeError("[RG-PPS] no interims with missing participants to evaluate.")

# %%

# =============================================================================
# Aggregate across interims; save artifacts mirroring the other interim scripts.
# =============================================================================

p_h1_xz = pd.concat(p_h1_xz_rows, ignore_index=True)
p_h1_xz['pps_H1_min_effect_size_thresh'] = pps_H1_min_effect_size_thresh
p_h1_xz['pps_ProbH1_target_lwr_quantile'] = pps_ProbH1_target_lwr_quantile
p_h1_xz['S'] = pps_z_total
perf_all = pd.concat(perf_rows, ignore_index=True)
mins_total = (time.time() - t_all0) / 60.0
pps_timing = pd.DataFrame(pps_timing_rows)
pps_timing['mins_total'] = round(mins_total, 3)
print(f"\nAll interims RG-PPS done in {mins_total:.2f} min")

perf_all.to_csv(os.path.join(dir_out, f"{file_prefix}_pps_RG_perf.csv"), index=False)
pkl_path = os.path.join(dir_out, f"{file_prefix}_pps_RG_p_h1_xz.pkl")
p_h1_xz.to_pickle(pkl_path)
print(f"Saved RG P(H_1 | x, z) samples to: {pkl_path}")
pps_timing.to_csv(os.path.join(dir_out, f"{file_prefix}_pps_RG_timing.csv"), index=False)

# PPS table (same schema as the other interim scripts).
pps_df = (
    p_h1_xz.groupby(['interim_id', 'interim_date', 'item_label', 'item_type',
                     'item_high_label'])['p_h1_xz']
    .apply(lambda p: float((p > pps_ProbH1_target_lwr_quantile).mean()))
    .reset_index(name='pps')
)
pps_df['eta'] = pps_ProbH1_target_lwr_quantile
pps_df['S'] = pps_z_total
pps_df.to_csv(os.path.join(dir_out, f"{file_prefix}_pps_RG.csv"), index=False)
print(pps_df.head(10).to_string(index=False))

# perf long-form pkl: rho, R^2 (from the Gaussian GLM on w_ratio vs
# pps_ratio_x), and the per-interim wall time (matches the time facet in the
# other scripts' perf_long pkls). No ESS analogue -- the regression has no IS
# weights; rho / R^2 are the natural internal quality diagnostics.
perf = perf_all.merge(di[['interim_id', 'interim_month_year']], on='interim_id', how='left')
perf = perf.merge(pps_timing[['interim_id', 'mins_interim_id']], on='interim_id', how='left')
perf_order = (
    di[di['interim_id'].isin(perf['interim_id'].unique())]
    .sort_values('interim_date')['interim_month_year'].tolist()
)
perf['interim_month_year'] = pd.Categorical(perf['interim_month_year'],
                                            categories=perf_order, ordered=True)
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
perf_long['method'] = 'Regression (Strong-Oakley)'
perf_long.to_pickle(os.path.join(dir_out, f"{file_prefix}_pps_RG_perf_long.pkl"))
print(f"Saved RG performance long-form pkl.")

# %%

# =============================================================================
# Plot: distribution of p(H_1 | x, z) per item across the S samples; save
# box_stats pkl (same schema as the other interim scripts).
# =============================================================================

tmp = p_h1_xz.merge(
    dit[['item_label', 'group_label_long', 'item_label_short']].drop_duplicates(),
    on='item_label', how='left',
)
tmp['item_label_long'] = tmp['group_label_long'] + np.where(
    tmp['item_label_short'].notna(), '\n' + tmp['item_label_short'], ''
)
tmp = tmp.merge(di[['interim_id', 'interim_month_year']], on='interim_id', how='left')
order = (
    di[di['interim_id'].isin(tmp['interim_id'].unique())]
    .sort_values('interim_date')['interim_month_year'].tolist()
)
tmp['interim_month_year'] = pd.Categorical(tmp['interim_month_year'],
                                           categories=order, ordered=True)

box_stats = (
    tmp.groupby(['item_label_long', 'interim_month_year'], observed=True)['p_h1_xz']
    .agg(
        min='min',
        q025=lambda s: s.quantile(0.025),
        q25=lambda s: s.quantile(0.25),
        q50=lambda s: s.quantile(0.50),
        q75=lambda s: s.quantile(0.75),
        q975=lambda s: s.quantile(0.975),
        max='max',
    )
    .reset_index()
)
pkl_path = os.path.join(dir_out, f"{file_prefix}_pps_RG_p_h1_xz_boxplot.pkl")
box_stats.to_pickle(pkl_path)
print(f"Saved RG p(H_1 | x, z) box-stats to: {pkl_path}")

all_items = list(pd.unique(tmp['item_label_long']))
color_dict = dict(zip(all_items, _futurama_palette(len(all_items))))

p = (
    ggplot(box_stats, aes(x='interim_month_year', fill='item_label_long'))
    + geom_boxplot(
        aes(ymin='q025', lower='q25', middle='q50', upper='q75', ymax='q975',
            group='interim_month_year'),
        stat='identity',
    )
    + geom_hline(yintercept=pps_ProbH1_target_lwr_quantile, colour='black', size=1.5)
    + facet_wrap('~ item_label_long', ncol=4)
    + scale_fill_manual(values=color_dict)
    + scale_y_continuous(
        limits=[0, 1],
        breaks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        labels=['0%', '20%', '40%', '60%', '80%', '100%'],
    )
    + theme_bw()
    + theme(
        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
        legend_position='none',
        figure_size=(15, 15),
        strip_text_y=element_text(angle=0),
    )
    + labs(x='Interim', y='p(H_1 | x, z)  [Regression (Strong-Oakley)]')
)
box_pdf = os.path.join(dir_out, f"{file_prefix}_pps_RG_p_h1_xz_boxplot.pdf")
p.save(box_pdf, verbose=False, limitsize=False)
print(f"Saved RG p(H_1 | x, z) boxplot to: {box_pdf}")

# %%

# =============================================================================
# Diagnostic-plot loop: per-interim binomial-GLM illustration. For each item
# scatter (jittered) the binary pps_H1_x indicator against w_ratio and overlay
# the per-item logistic curve from fit_interim_regress_H1x_on_wz, annotated
# with the empirical positive rate and McFadden pseudo-R^2 of the fit. Reloads
# the per-interim wa training pkl saved above.
# =============================================================================

def _per_item_H1x_summary(g: pd.DataFrame) -> pd.Series:
    """Per-item binomial-fit diagnostic: empirical pps_H1_x rate, McFadden
    pseudo-R^2 (1 - deviance/null_deviance from the binomial GLM), and panel-
    relative x positions for an annotation. Returns NaN where degenerate."""
    x = g['w_ratio'].to_numpy()
    y = g['pps_H1_x'].to_numpy()
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 3 or len(np.unique(y)) < 2:
        return pd.Series({'rate': float(y.mean()) if len(y) else np.nan,
                          'mcfadden_r2': np.nan,
                          'x_left': np.nan, 'x_right': np.nan})
    try:
        res = sm.GLM(y, sm.add_constant(x),
                     family=sm.families.Binomial()).fit()
        mcf = float(1.0 - res.deviance / res.null_deviance)
    except Exception:
        mcf = float('nan')
    return pd.Series({'rate': float(y.mean()), 'mcfadden_r2': mcf,
                      'x_left': float(x.min()), 'x_right': float(x.max())})


print("\nBuilding per-interim binomial-fit diagnostic plots...")
for interim_id in di['interim_id']:
    interim_id = int(interim_id)
    wa_pkl = os.path.join(dir_out, f"{file_prefix}_i{interim_id}_regression_training.pkl")
    if not os.path.exists(wa_pkl):
        continue
    wa = pd.read_pickle(wa_pkl)

    _ann = wa.groupby('item_label').apply(_per_item_H1x_summary).reset_index()
    _ann['rate_label'] = _ann['rate'].apply(lambda r: f"P(H1_x)={100*r:.0f}%")
    _ann['r2_label'] = _ann['mcfadden_r2'].apply(
        lambda r: f"McF-R2={100*r:.1f}%" if np.isfinite(r) else "McF-R2=NA"
    )
    _ann['y_top'] = 1.10
    _ann['y_bot'] = -0.10

    # Per-item logistic-fit curve via fit_interim_regress_H1x_on_wz: predicts
    # p_h1_xz at every training-row w_ratio, so plotting it as a line (sorted
    # by w_ratio) traces the fitted sigmoid for each item.
    pred, _ = fit_interim_regress_H1x_on_wz(wa)
    wa_for_merge = wa.assign(s=wa['draw'].astype(int) + 1)[
        ['item_label', 's', 'w_ratio']
    ]
    pred_xy = pred.merge(
        wa_for_merge, on=['item_label', 's'], how='left',
    ).sort_values(['item_label', 'w_ratio']).reset_index(drop=True)

    _items = sorted(wa['item_label'].unique())
    _item_colors = dict(zip(_items, _futurama_palette(len(_items))))

    p = (
        ggplot(wa, aes(x='w_ratio', y='pps_H1_x', colour='item_label'))
        + geom_jitter(alpha=0.25, size=0.5, height=0.05, width=0.0)
        + geom_hline(yintercept=0.5, colour='grey', size=0.3, linetype='dashed')
        + geom_line(pred_xy, aes(x='w_ratio', y='p_h1_xz'),
                    colour='black', size=0.8, inherit_aes=False)
        + geom_text(_ann, aes(x='x_left', y='y_top', label='rate_label'),
                    ha='left', va='top', size=8, inherit_aes=False)
        + geom_text(_ann, aes(x='x_right', y='y_top', label='r2_label'),
                    ha='right', va='top', size=8, inherit_aes=False)
        + facet_wrap('~ item_label', ncol=4, scales='free_x')
        + scale_colour_manual(values=_item_colors)
        + scale_y_continuous(limits=[-0.15, 1.15],
                             breaks=[0.0, 0.25, 0.5, 0.75, 1.0])
        + theme_bw()
        + theme(
            figure_size=(15, 14),
            strip_background=element_rect(fill='none', colour='none'),
            panel_background=element_rect(fill='none', colour='none'),
            legend_position='none',
        )
        + labs(x=f'empirical endpoint summary w(z) at interim id {interim_id}',
               y=f'pps_H1_x (jittered) + fitted sigmoid at interim id {interim_id}')
    )
    pdf_path = os.path.join(dir_out, f"{file_prefix}_i{interim_id}_H1x_vs_wratio.pdf")
    p.save(pdf_path, verbose=False, limitsize=False)
    print(f"  saved diagnostic plot: {pdf_path}")

print("\n✓ regression-H1x PPS + binomial-fit diagnostic plots complete for all interims")
