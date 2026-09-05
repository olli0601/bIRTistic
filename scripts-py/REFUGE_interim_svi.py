"""Produce the SVI grid for the REFUGE-ED application (§3.9), participant-count interims.
Per interim writes the artifacts the deepset amortiser deploy needs:
  {prefix}_{k}_data_dp1.csv, {prefix}_{k}_draws.zarr, {prefix}_i{k}_regression_training.pkl,
  {prefix}_1_data_dit.csv. MSPSS Likert items -> expected-score endpoint (item_type 'out-of-7')."""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
import numpy as np, pandas as pd
from data_loading import read_data_refuge_ed
from model_pcm import PartialCreditModel

SB = "/Users/or105/sandbox/bIRTistic"
RAW = os.environ.get('REFUGE_XLSX',
    "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/2025_project_Hope_Groups/data/refuge/Youth Baseline & Endline .xlsx")
dir_out = f"{SB}/py-refugee_interim_260831"; os.makedirs(dir_out, exist_ok=True)
file_prefix = "pcm_1_interim"; x_formula = "~ time - 1"; seed = 123
NSTEPS = int(os.environ.get('REFUGE_STEPS', '4000')); S = int(os.environ.get('REFUGE_S', '2000'))
NINT = int(os.environ.get('REFUGE_NINT', '8'))
os.environ.setdefault('PROB_FIT_WIDTH_MULT', '1.4')          # 12 MSPSS items crowd the panel; widen

raw = read_data_refuge_ed(RAW)
dp, dit = raw['dp'].copy(), raw['dit'].copy()
# expected-score endpoint: treat MSPSS Likert as one 'out-of-7' item_type
dit['item_type'] = 'out-of-7'; dit['item_type_id'] = 1
# group_label / group_label_long come from read_data_refuge_ed (per-subscale: Fam/Fri/SO)
dit.to_csv(f"{dir_out}/{file_prefix}_1_data_dit.csv", index=False)

# linked participants only (present at both times)
both = dp.groupby('pid').time.nunique(); linked = set(both[both >= 2].index)
dp1 = dp[dp.pid.isin(linked)].copy()
_R = dp1['y'].astype(int)                                     # raw response 1..7
dp1['y_stan'] = _R                                            # model input (1-indexed)
dp1['y'] = _R - 1                                             # 0-indexed category (diagnostic plots)
dp1['item_type'] = 'out-of-7'; dp1['item_type_id'] = 1
it = (dp1[['item_label', 'time']].drop_duplicates().sort_values(['time', 'item_label']).reset_index(drop=True))
it['item_time_id'] = np.arange(1, len(it) + 1)
dp1 = dp1.merge(it, on=['item_label', 'time'], how='left')
dp1 = dp1.sort_values(['item_type_id', 'pid', 'time', 'item_label']).reset_index(drop=True)
dp1['oid'] = range(1, len(dp1) + 1)
dp1['oidt'] = dp1.groupby('item_type').cumcount() + 1
pids = np.sort(dp1.pid.unique()); n_full = len(pids)
grid = np.unique(np.round(np.linspace(max(30, n_full // NINT), n_full, NINT)).astype(int))
print(f"REFUGE-ED: n_full={n_full} linked, {dp1.item_label.nunique()} items, interims n={grid.tolist()}")

for k, n in enumerate(grid, 1):
    obs = set(pids[:n]); xi = dp1[dp1.pid.isin(obs)].copy()
    xi['oid'] = range(1, len(xi) + 1); xi['oidt'] = xi.groupby('item_type').cumcount() + 1
    m = n_full - n if n < n_full else max(1, n_full // 4)     # future cohort size
    xi.to_csv(f"{dir_out}/{file_prefix}_{k}_data_dp1.csv", index=False)
    pre = f"{dir_out}/{file_prefix}_{k}"
    print(f"\n=== interim {k}: n={n} m={m} ===")
    t0 = time.time()
    model = PartialCreditModel(dit=dit, dcati=xi, x_formula=x_formula, seed=seed)
    fit = model.fit_pyro_svi(output_file_prefix=pre, algorithm='AutoDiagonalNormal',
                             lr=0.01, num_steps=NSTEPS, output_samples=S, resume=True,
                             with_core_analyses=True, with_additional_analyses=False, verbose=False)
    xr = model.get_endpoints_per_draw(draws=fit['draws'], categorical_threshold=2,
                                      endpoint_type='items').rename(columns={'ratio': 'pps_ratio_x'})
    xr['pps_H1_x'] = (xr['pps_ratio_x'] > 0.5).astype(int)
    if 'item_high_label' not in xr.columns:
        xr = xr.merge(dit[['item_label', 'item_high_label']], on='item_label', how='left')
    xr[['draw', 'item_label', 'item_type', 'item_high_label', 'pps_ratio_x', 'pps_H1_x']].to_pickle(
        f"{dir_out}/{file_prefix}_i{k}_regression_training.pkl")
    print(f"  done ({(time.time()-t0)/60:.1f} min)")
print("REFUGE-ED SVI grid complete ->", dir_out)
