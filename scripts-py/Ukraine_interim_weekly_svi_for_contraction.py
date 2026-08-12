#!/usr/bin/env python3
"""
Weekly-cadence interim SVI generator for the posterior-contraction
diagnostic (§14.4.9 of dev/amortised_decision_making.md).

The main interim pipeline uses MONTHLY cutoffs (8 interims). For the
contraction diagnostic we want many more (n, width) evaluation points,
so this re-slices the SAME Ukraine trial at WEEKLY cutoffs and produces,
per weekly interim, exactly the artifacts the deepset amortiser needs:

  {prefix}_{k}_data_dp1.csv          observed cohort xi at that cutoff
  {prefix}_{k}_draws.zarr            SVI posterior-predictive z draws
  {prefix}_i{k}_regression_training.pkl   wa (SVI target rho = pps_ratio_x/draw)
  {prefix}_1_data_dit.csv            item table (once)

Same SVI config as the monthly regression run (AutoLowRankMVN, 10k steps,
4000 draws). No regression baselines, no plots.
"""
import os
import sys
import time
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'python'))
import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')
from data_loading import read_data_ukraine
from model_pcm import PartialCreditModel

np.random.seed(42)
seed = 123
dir_data = "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/2025_project_Hope_Groups/data"
file_data = os.path.join(dir_data, "Ukraine_Hope_Groups_Baseline_Endline_Wide_Aug6.csv")
dir_out = "/Users/or105/sandbox/bIRTistic/py-ukraine-interim-weekly-svi-260811"
os.makedirs(dir_out, exist_ok=True)
file_prefix = "pcm_1_interim"
svi_algorithm = 'AutoLowRankMultivariateNormal'
x_formula = "~ time - 1"
pps_z_total = 4000
categorical_threshold = 2
print(f"Output dir: {dir_out}")

# ---- load + preprocess (verbatim from the monthly regression script) ----
raw = read_data_ukraine(file_data)
dp = raw['dp'].copy(); dit = raw['dit'].copy(); dmeta = raw['dmeta'].copy()
tmp = (dmeta[['pid', 'time_label', 'displacement_status']].drop_duplicates()
       .rename(columns={'pid': 'pid_label'}).dropna(subset=["displacement_status"], how='all'))
dp = dp.merge(tmp, on=['pid_label', 'time_label'], how='inner', validate='many_to_one')
dp1 = dp[~dp['item_label'].str.contains('agg')].copy()
dp1['y_stan'] = dp1['y'] + 1
dp1 = dp1.merge(dit[['item_label', 'item_type']], on='item_label', how='left')
item_time_df = (dp1[['item_type', 'item_label', 'time']].drop_duplicates()
                .sort_values(['item_type', 'time', 'item_label']).reset_index(drop=True))
item_time_df['item_time_id'] = item_time_df.groupby('item_type').cumcount() + 1
dp1 = dp1.merge(item_time_df, on=['item_label', 'time', 'item_type'], how='left')
dp1 = dp1.merge(dit[['item_type', 'item_type_id']].drop_duplicates(), on='item_type', how='left')
dp1 = dp1.sort_values(['item_type_id', 'pid', 'time', 'item_label']).reset_index(drop=True)
dp1['oid'] = range(1, len(dp1) + 1)
dp1['oidt'] = dp1.groupby('item_type').cumcount() + 1
n_full = dp1['pid'].nunique()

# save item table once (deepset reads dit['cat_length'])
dit.to_csv(os.path.join(dir_out, f"{file_prefix}_1_data_dit.csv"), index=False)

# ---- WEEKLY cutoffs (was freq='MS' month-starts) ----
_endline_dates = pd.to_datetime(
    dp1.loc[dp1['time_label'] == 'Endline', 'submission_date']).dropna()
_start = _endline_dates.min().normalize()
_end = _endline_dates.max().normalize()
week_ends = pd.date_range(_start, _end, freq='W')
di = pd.DataFrame({'interim_id': range(1, len(week_ends) + 1),
                   'interim_date': week_ends})
print(f"n_full={n_full} | {len(di)} weekly cutoffs {_start.date()}..{_end.date()}")

timing_rows = []
t_all0 = time.time()
for i in range(len(di)):
    interim_id = int(di['interim_id'].iloc[i])
    interim_date = pd.to_datetime(di['interim_date'].iloc[i])
    xi = PartialCreditModel.get_interim_data_x(dp1, interim_date)
    if xi.empty or xi['pid'].nunique() < 2:
        print(f"[week {interim_id} {interim_date.date()}] skip: n<2"); continue
    interim_m = n_full - xi['pid'].nunique()
    if interim_m <= 0:
        print(f"[week {interim_id} {interim_date.date()}] skip: m<=0"); continue
    n_obs = xi['pid'].nunique()
    print(f"\n{'='*60}\nweek {interim_id} ({interim_date.date()}): n={n_obs} m={interim_m}\n{'='*60}")
    t0 = time.time()
    interim_prefix = os.path.join(dir_out, f"{file_prefix}_{interim_id}")
    # observed cohort for the deepset deploy
    xi.to_csv(f"{interim_prefix}_data_dp1.csv", index=False)
    # SVI fit on x -> draws.zarr (posterior-predictive z)
    model = PartialCreditModel(dit=dit, dcati=xi, x_formula=x_formula, seed=seed)
    fit = model.fit_pyro_svi(
        output_file_prefix=interim_prefix, algorithm=svi_algorithm,
        lr=0.01, num_steps=10000, output_samples=pps_z_total, resume=True,
        with_core_analyses=True, with_additional_analyses=False, verbose=False)
    zi = model.get_interim_z_from_ypredi(
        f"{interim_prefix}_draws.zarr", interim_m,
        pps_z_total=pps_z_total, seed=seed, keep_order=True)
    wa = model.get_w(zi, categorical_threshold=categorical_threshold)
    x_ratio = model.get_endpoints_per_draw(
        draws=fit['draws'], categorical_threshold=categorical_threshold,
        endpoint_type='items').rename(columns={'ratio': 'pps_ratio_x'})
    x_ratio['pps_H1_x'] = (x_ratio['pps_ratio_x'] > 0.5).astype(int)
    wa = wa.merge(
        x_ratio[['draw', 'item_label', 'item_type', 'pps_ratio_x', 'pps_H1_x']],
        on=['draw', 'item_label', 'item_type'], how='inner')
    wa.to_pickle(os.path.join(
        dir_out, f"{file_prefix}_i{interim_id}_regression_training.pkl"))
    mins = (time.time() - t0) / 60.0
    timing_rows.append({'interim_id': interim_id, 'interim_date': interim_date,
                        'n_obs': n_obs, 'interim_m': interim_m, 'mins': round(mins, 3)})
    print(f"  week {interim_id} done in {mins:.2f} min")

pd.DataFrame(timing_rows).to_csv(
    os.path.join(dir_out, f"{file_prefix}_weekly_timing.csv"), index=False)
print(f"\nWeekly SVI complete: {len(timing_rows)} interims in "
      f"{(time.time() - t_all0) / 60.0:.1f} min -> {dir_out}")
