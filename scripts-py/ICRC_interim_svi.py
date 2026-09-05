"""SVI grid for the ICRC community-MHPSS application (Andersen et al. 2022, §3.x),
DASS-21 / DRC arm. Each beneficiary's Depression/Anxiety/Stress subscale total is binned
into the Figure-2 severity levels (Normal..Extremely severe, K=5), giving 3 ordinal items.
PAIRED pre/post (Baseline=pre, Endline=post). Endpoint rho_j = severity reduction
(lower_is_better). Interims accrue every 100 beneficiaries. No calendar time in the data.
Per interim: dp1.csv, draws.zarr, i{k}_regression_training.pkl, 1_data_dit.csv."""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
import numpy as np, pandas as pd
from model_pcm import PartialCreditModel

SB = "/Users/or105/sandbox/bIRTistic"
RAW = os.environ.get('ICRC_XLSX',
    "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/2025_project_Hope_Groups/data/Andersen_Community-Level Mental Health Support Table.XLSX")
dir_out = f"{SB}/py-icrc-dass-drc_260902"; os.makedirs(dir_out, exist_ok=True)
file_prefix = "pcm_1_interim"; x_formula = "~ time - 1"; seed = 123
NSTEPS = int(os.environ.get('ICRC_STEPS', '4000')); S = int(os.environ.get('ICRC_S', '2000'))
STEP = int(os.environ.get('ICRC_STEP', '100'))               # interim every STEP participants
os.environ.setdefault('PROB_FIT_WIDTH_MULT', '1.2')

# DASS-21 (x2 scaled, 0-42) severity bins per Figure 2 -> ordinal 0..4
SEV = ['Normal', 'Mild', 'Moderate', 'Severe', 'Extremely severe']
BINS = {'Depression': [-np.inf, 9, 13, 20, 27, np.inf],
        'Anxiety':    [-np.inf, 7, 9, 14, 19, np.inf],
        'Stress':     [-np.inf, 14, 18, 25, 33, np.inf]}
COL = {'Depression': ('prdep', 'podep'), 'Anxiety': ('pranx', 'poanx'), 'Stress': ('prstr', 'postr')}

raw = pd.read_excel(RAW, sheet_name='Data', header=0)
raw = raw[raw['ctry'] == 2].reset_index(drop=True)           # DRC only (the DASS arm)
need = [c for cc in COL.values() for c in cc]
for c in need:
    raw[c] = pd.to_numeric(raw[c], errors='coerce')
raw = raw.dropna(subset=need).reset_index(drop=True)         # complete pre&post on all 3
raw['pid'] = np.arange(1, len(raw) + 1)
print(f"ICRC DASS/DRC: {len(raw)} beneficiaries with complete pre&post on all 3 subscales")

rows = []
for _, r in raw.iterrows():
    for item, (pc, oc) in COL.items():
        for t, col in ((0, pc), (1, oc)):
            ordv = int(pd.cut([r[col]], BINS[item], labels=[0, 1, 2, 3, 4])[0])
            rows.append(dict(pid=int(r['pid']), time=t,
                             time_label='Baseline' if t == 0 else 'Endline',
                             item_label=f"DASS_{item}", y=ordv))
dp1 = pd.DataFrame(rows)
dp1['pid_label'] = dp1['pid'].astype(str); dp1['fid'] = np.nan; dp1['f_label'] = np.nan
dp1['submission_date'] = pd.NaT; dp1['treat'] = 0.0
dp1['y_stan'] = dp1['y'] + 1                                  # 1..5 model input
dp1['y_label'] = dp1['y'].map({i: f"{i} {s}" for i, s in enumerate(SEV)})  # ordinal-prefixed so x-axis sorts by severity
dp1['item_type'] = 'out-of-7'; dp1['item_type_id'] = 1       # expected-score branch; K set by cat_length

# dit (3 items, K=5, higher DASS = worse -> lower_is_better)
dit = pd.DataFrame({'item_label': [f"DASS_{i}" for i in COL]})
dit['item_type'] = 'out-of-7'; dit['item_type_id'] = 1; dit['cat_length'] = 5
dit['item_label_short'] = [i for i in COL]
dit['group_label'] = 'DASS-21 distress'; dit['group_label_long'] = 'DASS-21 distress'
dit['item_high_label'] = 'lower_is_better'
dit['endpoint_measure'] = 'mean DASS-21 severity level (0-4, pre vs post)'
dit.to_csv(f"{dir_out}/{file_prefix}_1_data_dit.csv", index=False)

it = (dp1[['item_label', 'time']].drop_duplicates().sort_values(['time', 'item_label']).reset_index(drop=True))
it['item_time_id'] = np.arange(1, len(it) + 1)
dp1 = dp1.merge(it, on=['item_label', 'time'], how='left')

pids = np.sort(dp1.pid.unique()); n_full = len(pids)
EARLY = [20, 40, 60, 80]                                      # fine early grid (effect emerges fast)
grid = sorted(set([g for g in EARLY if g < n_full] + list(range(STEP, n_full, STEP)) + [n_full]))
print(f"interims (fine early + every {STEP}): n={grid}")

for k, n in enumerate(grid, 1):
    obs = set(pids[:n]); xi = dp1[dp1.pid.isin(obs)].copy()
    xi = xi.sort_values(['item_type_id', 'pid', 'time', 'item_label']).reset_index(drop=True)
    xi['oid'] = range(1, len(xi) + 1); xi['oidt'] = xi.groupby('item_type').cumcount() + 1
    xi.to_csv(f"{dir_out}/{file_prefix}_{k}_data_dp1.csv", index=False)
    pre = f"{dir_out}/{file_prefix}_{k}"
    print(f"\n=== interim {k}: n={n} ===")
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
print("ICRC DASS/DRC SVI grid complete ->", dir_out)
