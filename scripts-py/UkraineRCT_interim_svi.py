"""SVI grid for the Ukraine Hope Groups RANDOMISED-TRIAL between-arm design (UkraineRCT).

Corrected estimand (supersedes the parked §14 pre->post pooled change and the §15 UkraineP
cross-arm design). The control arm received NO intervention between baseline and endline, so a
within-participant change is not an intervention effect and even the cross-arm
intervention-endline-vs-control-BASELINE contrast mixes a secular/time component. The scientific
object is the GROUP-MEAN comparison AT ENDLINE, per item j, of the intervention arm vs the control
arm, as a direction-aware relative percent improvement (higher = better):

    group 0  'Control (endline)'       <- control arm, endline survey   (reference)
    group 1  'Intervention (endline)'  <- intervention arm, endline survey (treated)
    rho_j  = direction-aware ( w_intervention,j / w_control,j - 1 )      (same get_endpoints call)

The trial randomises WITHIN facilitator (`fid`): every facilitator runs both arms, and their two
endline surveys are close in calendar time (checked by data_loading.get_data_ukraine_rct: 30/30
facilitators within ~3 weeks), so the pooled between-arm PCM contrast (`~ group - 1`, group = arm)
is unconfounded by facilitator and by secular timing. Disjoint participants across the two cells =>
an unpaired between-group PCM (the paired amortiser of §14 does not apply; SVI-only reference).

Interims accrue every TWO WEEKS: at each biweekly cutoff the cohort is every endline record (both
arms) submitted on/before the cutoff. 26->19 items after dropping the composite `_agg` items; two
item-type families (out-of-7 days K=8, categorical CG-MH K=4) share one fit but keep separate
item_type_id (no K-mixing).

Per interim writes the standard SVI artifacts (dp1.csv, draws.zarr, i{k}_regression_training.pkl,
1_data_dit.csv, and the core fit plots via with_core_analyses=True) plus an interim index and a
per-item rho trajectory CSV = all diagnostics for the SVI run.

    UKRAINERCT_STEPS / _S / _SVI override the fit knobs;   python scripts-py/UkraineRCT_interim_svi.py
"""
import os, sys, time
from pathlib import Path
try:                                   # running as a file
    _root = Path(__file__).resolve().parent.parent
except NameError:                      # pasted into the REPL (no __file__)
    _root = Path.cwd()
    if _root.name == 'scripts-py':
        _root = _root.parent
sys.path.insert(0, str(_root / 'python'))
import numpy as np, pandas as pd
from data_loading import get_data_ukraine_rct
from model_pcm import PartialCreditModel

SB = os.environ.get('UKRAINERCT_SB', "/Users/or105/sandbox/bIRTistic")
dir_data = ("/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/"
            "OR_Work/2025/2025_project_Hope_Groups/data")
file_data = os.environ.get('UKRAINERCT_CSV',
    os.path.join(dir_data, "Ukraine_Hope_Groups_Baseline_Endline_Wide_Aug6.csv"))
dir_out = f"{SB}/py-ukraineRCT-endline-svi-260929"; os.makedirs(dir_out, exist_ok=True)
file_prefix = "pcm_1_interim"; x_formula = "~ group - 1"; seed = 123
NSTEPS = int(os.environ.get('UKRAINERCT_STEPS', '10000'))
S = int(os.environ.get('UKRAINERCT_S', '4000'))
CAT_THR = int(os.environ.get('UKRAINERCT_CATTHR', '2'))          # caseness cut for the CG-MH items: P(y>=2)
H1 = float(os.environ.get('UKRAINERCT_H1', '0.5'))              # success threshold on rho (>50% rel. improvement)
svi_algorithm = os.environ.get('UKRAINERCT_SVI', 'AutoLowRankMultivariateNormal')
os.environ.setdefault('PROB_FIT_WIDTH_MULT', '1.3')

# ---- load: endline both arms, group re-encoded to arm; facilitator date-proximity checked ----
raw = get_data_ukraine_rct(file_data, max_endline_gap_days=21, drop_gap_violators=False)
dp, dit = raw['dp'].copy(), raw['dit'].copy()
raw['date_check'].to_csv(f"{dir_out}/{file_prefix}_facilitator_date_check.csv", index=False)

# drop the composite/aggregate items (not modelled); PCM-prep (y_stan, item_type, ids)
dp1 = dp[~dp['item_label'].str.contains('agg')].copy()
dp1['y_stan'] = dp1['y'] + 1
dp1 = dp1.merge(dit[['item_label', 'item_type']], on='item_label', how='left')

# item_group_id = per-item_type (item x arm) difficulty index; item_type_id keeps the two K-families apart
item_arm_df = (dp1[['item_type', 'item_label', 'group']].drop_duplicates()
               .sort_values(['item_type', 'group', 'item_label']).reset_index(drop=True))
item_arm_df['item_group_id'] = item_arm_df.groupby('item_type').cumcount() + 1
dp1 = dp1.merge(item_arm_df, on=['item_label', 'group', 'item_type'], how='left')
dp1 = dp1.merge(dit[['item_type', 'item_type_id']].drop_duplicates(), on='item_type', how='left')
dit.to_csv(f"{dir_out}/{file_prefix}_1_data_dit.csv", index=False)

# ---- BIWEEKLY interim cutoffs over the endline submission-date span (+ final full cutoff) ----
sd = pd.to_datetime(dp1['submission_date'])
cutoffs = list(pd.date_range(sd.min().normalize(), sd.max().normalize(), freq='2W'))
if not cutoffs or cutoffs[-1] < sd.max().normalize():
    cutoffs.append(sd.max().normalize())                        # capture the final partial fortnight
print(f"UkraineRCT between-arm endline: {dp1.item_label.nunique()} items, "
      f"{dp1['item_type'].nunique()} item types; {len(cutoffs)} biweekly interims "
      f"[{cutoffs[0].date()} .. {cutoffs[-1].date()}]")

rows, rho_rows = [], []
for k, cutoff in enumerate(cutoffs, 1):
    cutoff = pd.to_datetime(cutoff)
    xi = dp1[pd.to_datetime(dp1['submission_date']) <= cutoff].copy()
    nb = xi[xi.group == 0].pid_label.nunique(); ne = xi[xi.group == 1].pid_label.nunique()
    if nb < 2 or ne < 2:                                        # need both arms populated to fit the contrast
        print(f"  [interim {k} {cutoff.date()}] skip: ctrl n={nb}, intv n={ne}")
        continue
    # sequential pid within the interim (people are disjoint across the two arms)
    xi['pid'] = pd.factorize(xi['pid_label'].astype(str))[0] + 1
    xi = xi.sort_values(['item_type_id', 'pid', 'group', 'item_label']).reset_index(drop=True)
    xi['oid'] = range(1, len(xi) + 1); xi['oidt'] = xi.groupby('item_type').cumcount() + 1
    xi.to_csv(f"{dir_out}/{file_prefix}_{k}_data_dp1.csv", index=False)
    pre = f"{dir_out}/{file_prefix}_{k}"
    print(f"\n=== interim {k} ({cutoff.date()}): {xi.fid.nunique()} facilitators | "
          f"control n={nb}, intervention n={ne} ===")
    t0 = time.time()
    model = PartialCreditModel(dit=dit, dcati=xi, x_formula=x_formula, seed=seed)
    fit = model.fit_pyro_svi(output_file_prefix=pre, algorithm=svi_algorithm,
                             lr=0.01, num_steps=NSTEPS, output_samples=S, resume=True,
                             with_core_analyses=True, with_additional_analyses=False, verbose=False)
    xr = model.get_endpoints_per_draw(draws=fit['draws'], categorical_threshold=CAT_THR,
                                      endpoint_type='items').rename(columns={'ratio': 'pps_ratio_x'})
    xr['pps_H1_x'] = (xr['pps_ratio_x'] > H1).astype(int)
    if 'item_high_label' not in xr.columns:
        xr = xr.merge(dit[['item_label', 'item_high_label']], on='item_label', how='left')
    xr[['draw', 'item_label', 'item_type', 'item_high_label', 'pps_ratio_x', 'pps_H1_x']].to_pickle(
        f"{dir_out}/{file_prefix}_i{k}_regression_training.pkl")
    # per-item rho trajectory (median + 95% CrI) for the SVI diagnostics
    for il, g in xr.groupby('item_label'):
        rho_rows.append(dict(interim=k, interim_date=cutoff.date(), item_label=il,
                             item_type=g['item_type'].iloc[0],
                             rho_med=float(g.pps_ratio_x.median()),
                             rho_lo=float(g.pps_ratio_x.quantile(.025)),
                             rho_hi=float(g.pps_ratio_x.quantile(.975)),
                             p_H1=float(g.pps_H1_x.mean())))
    rows.append(dict(interim=k, interim_date=cutoff.date(), n_facilitators=int(xi.fid.nunique()),
                     n_control=nb, n_intervention=ne, rho_med=float(xr.pps_ratio_x.median())))
    print(f"  rho(intervention vs control, endline) median={xr.pps_ratio_x.median():+.3f} "
          f"({(time.time()-t0)/60:.1f} min)")

pd.DataFrame(rows).to_csv(f"{dir_out}/{file_prefix}_interim_index.csv", index=False)
pd.DataFrame(rho_rows).to_csv(f"{dir_out}/{file_prefix}_rho_trajectory.csv", index=False)
print(f"\nUkraineRCT between-arm SVI grid complete ({len(rows)} interims) -> {dir_out}")
