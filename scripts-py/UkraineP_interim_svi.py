"""SVI grid for the Ukraine Hope Groups CROSS-ARM, facilitator-matched design
(UkraineP, §15). Instead of the pooled within-participant pre->post change, the
endpoint contrasts, per item, the group-mean item level of the INTERVENTION arm at
ENDLINE (treated + post) against the group-mean of the CONTROL arm at BASELINE
(untreated + pre). The two groups are roughly time-matched and are matched on
facilitator (every facilitator runs both an intervention group and a waitlist-
control group), so the contrast controls for facilitator and secular timing that
the pooled pre->post estimand (~ time - 1 on both arms) confounds.

Design.
  * Two cells of the trial are kept and mapped to the PCM's two 'time' levels:
      time 0  'Baseline'  <-  CONTROL arm, baseline survey   (reference group)
      time 1  'Endline'   <-  INTERVENTION arm, endline survey (treated group)
    The other two cells (intervention-baseline, control-endline) are dropped.
    The cells are DISJOINT sets of people -> unpaired between-group PCM (as in the
    mycelium powder-vs-burger and PISA between-cohort designs; the paired amortiser
    of §14 does NOT apply, so this is an SVI-only reference application).
  * Endpoint rho_j = s_j( w_end,j / w_base,j - 1 ) = direction-aware % change of the
    intervention-endline vs control-baseline group-mean item level (same
    get_endpoints machinery; wb = control-baseline, we = intervention-endline).
  * Interims accrue FACILITATORS (increasingly many group-mean comparisons), ordered
    by the median endline submission date of each facilitator's intervention group.

Per interim writes (same artifacts as the other SVI producers) dp1.csv, draws.zarr,
i{k}_regression_training.pkl, 1_data_dit.csv, and the core fit plots.
"""
# ---- boilerplate ----

import os, sys, time
from pathlib import Path
try:                                   # running as a file
    _root = Path(__file__).resolve().parent.parent
except NameError:                      # pasted / Shift+Enter into the REPL (no __file__)
    _root = Path.cwd()
    if _root.name == 'scripts-py':
        _root = _root.parent
sys.path.insert(0, str(_root / 'python'))
import numpy as np, pandas as pd
from data_loading import read_data_ukraine
from model_pcm import PartialCreditModel

SB = "/Users/or105/sandbox/bIRTistic"
dir_data = ("/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/"
            "OR_Work/2025/2025_project_Hope_Groups/data")
file_data = os.environ.get('UKRAINEP_CSV',
    os.path.join(dir_data, "Ukraine_Hope_Groups_Baseline_Endline_Wide_Aug6.csv"))
dir_out = f"{SB}/py-ukraineP-crossarm-svi-260916"; os.makedirs(dir_out, exist_ok=True)
file_prefix = "pcm_1_interim"; x_formula = "~ time - 1"; seed = 123
NSTEPS = int(os.environ.get('UKRAINEP_STEPS', '10000'))
S = int(os.environ.get('UKRAINEP_S', '4000'))
NINT = int(os.environ.get('UKRAINEP_NINT', '8'))
svi_algorithm = os.environ.get('UKRAINEP_SVI', 'AutoLowRankMultivariateNormal')
os.environ.setdefault('PROB_FIT_WIDTH_MULT', '1.3')

# ---- load + preprocess (mirrors the weekly-svi producer) ----
raw = read_data_ukraine(file_data)
dp = raw['dp'].copy(); dit = raw['dit'].copy()
dp1 = dp[~dp['item_label'].str.contains('agg')].copy()
dp1['y_stan'] = dp1['y'] + 1
dp1 = dp1.merge(dit[['item_label', 'item_type']], on='item_label', how='left')

# ---- restrict to the two cross-arm cells (already carry the right time labels) ----
#   control  + Baseline (time 0) = reference group;  intervention + Endline (time 1) = treated group
cellB = (dp1['treat'] == 0) & (dp1['time'] == 0)     # control-baseline
cellA = (dp1['treat'] == 1) & (dp1['time'] == 1)     # intervention-endline
dp1 = dp1[cellB | cellA].copy()

# item_time_id per (item_type, time); item_type_id (two K-families kept separate -> no K-mixing)
item_time_df = (dp1[['item_type', 'item_label', 'time']].drop_duplicates()
                .sort_values(['item_type', 'time', 'item_label']).reset_index(drop=True))
item_time_df['item_time_id'] = item_time_df.groupby('item_type').cumcount() + 1
dp1 = dp1.merge(item_time_df, on=['item_label', 'time', 'item_type'], how='left')
dp1 = dp1.merge(dit[['item_type', 'item_type_id']].drop_duplicates(), on='item_type', how='left')

dit.to_csv(f"{dir_out}/{file_prefix}_1_data_dit.csv", index=False)

# ---- facilitator accrual order: median endline (intervention) submission date ----
fac_order = (dp1[dp1['time'] == 1].groupby('fid')['submission_date']
             .median().sort_values().index.to_numpy())
n_fac = len(fac_order)
grid = np.unique(np.round(np.linspace(max(4, n_fac // NINT), n_fac, NINT)).astype(int))
print(f"UkraineP cross-arm: {n_fac} facilitators, {dp1.item_label.nunique()} items; "
      f"interims accrue facilitators k={grid.tolist()}")

rows = []
for k, nfac in enumerate(grid, 1):
    keep_fac = set(fac_order[:nfac])
    xi = dp1[dp1['fid'].isin(keep_fac)].copy()
    # sequential pid within the interim (people are disjoint across the two cells)
    xi['pid'] = pd.factorize(xi['pid_label'].astype(str))[0] + 1
    xi = xi.sort_values(['item_type_id', 'pid', 'time', 'item_label']).reset_index(drop=True)
    xi['oid'] = range(1, len(xi) + 1); xi['oidt'] = xi.groupby('item_type').cumcount() + 1
    nb = xi[xi.time == 0].pid.nunique(); ne = xi[xi.time == 1].pid.nunique()
    xi.to_csv(f"{dir_out}/{file_prefix}_{k}_data_dp1.csv", index=False)
    pre = f"{dir_out}/{file_prefix}_{k}"
    print(f"\n=== interim {k}: {nfac} facilitators | ctrl-baseline n={nb}, int-endline n={ne} ===")
    t0 = time.time()
    model = PartialCreditModel(dit=dit, dcati=xi, x_formula=x_formula, seed=seed)
    fit = model.fit_pyro_svi(output_file_prefix=pre, algorithm=svi_algorithm,
                             lr=0.01, num_steps=NSTEPS, output_samples=S, resume=True,
                             with_core_analyses=True, with_additional_analyses=False, verbose=False)
    xr = model.get_endpoints_per_draw(draws=fit['draws'], categorical_threshold=2,
                                      endpoint_type='items').rename(columns={'ratio': 'pps_ratio_x'})
    xr['pps_H1_x'] = (xr['pps_ratio_x'] > 0.5).astype(int)
    if 'item_high_label' not in xr.columns:
        xr = xr.merge(dit[['item_label', 'item_high_label']], on='item_label', how='left')
    xr[['draw', 'item_label', 'item_type', 'item_high_label', 'pps_ratio_x', 'pps_H1_x']].to_pickle(
        f"{dir_out}/{file_prefix}_i{k}_regression_training.pkl")
    rows.append(dict(interim=k, n_fac=nfac, n_ctrl_base=nb, n_int_end=ne,
                     rho_med=float(xr.pps_ratio_x.median())))
    print(f"  rho(int-endline vs ctrl-baseline) median={xr.pps_ratio_x.median():+.3f} "
          f"({(time.time()-t0)/60:.1f} min)")
pd.DataFrame(rows).to_csv(f"{dir_out}/{file_prefix}_interim_index.csv", index=False)
print(f"UkraineP cross-arm SVI grid complete ({len(grid)} interims) ->", dir_out)
