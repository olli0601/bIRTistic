"""SVI grid for the mycelium acceptance study (§3.14, Fischer & Hilboesen 2025),
Option A: powder vs burger, pooling all three substrates. BETWEEN-arm design ->
burger = baseline (time 0), powder = endline (time 1); each respondent in one arm.
Endpoint rho_j = relative shift powder-vs-burger. Items: Acceptance (A1-A4),
Disgust (D1-D4), Perceived naturalness (PN); all 7-point -> 'out-of-7' expected score.
Per interim writes dp1.csv, draws.zarr, i{k}_regression_training.pkl, 1_data_dit.csv."""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
import numpy as np, pandas as pd
from model_pcm import PartialCreditModel

SB = "/Users/or105/sandbox/bIRTistic"
RAW = os.environ.get('MYCELIUM_CSV',
    "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/2025_project_Hope_Groups/data/mycelium/Mycelium.csv")
dir_out = f"{SB}/py-mycelium-powdervsburger_260902"; os.makedirs(dir_out, exist_ok=True)
file_prefix = "pcm_1_interim"; x_formula = "~ time - 1"; seed = 123
NSTEPS = int(os.environ.get('MYC_STEPS', '4000')); S = int(os.environ.get('MYC_S', '2000'))
NINT = int(os.environ.get('MYC_NINT', '8'))
os.environ.setdefault('PROB_FIT_WIDTH_MULT', '1.3')

# item -> (paper short code, construct, direction)
ITEMS = {
    'INT1': ('A1', 'Acceptance', 'higher_is_better'),      # would consider consuming
    'ATT1': ('A2', 'Acceptance', 'higher_is_better'),      # would enjoy consuming
    'INT2': ('A3', 'Acceptance', 'higher_is_better'),      # would make me want to consume
    'ATT2': ('A4', 'Acceptance', 'higher_is_better'),      # would make me feel good
    'DISG1': ('D1', 'Disgust', 'lower_is_better'),
    'DISG2': ('D2', 'Disgust', 'lower_is_better'),
    'DISG3': ('D3', 'Disgust', 'lower_is_better'),
    'DISG4': ('D4', 'Disgust', 'lower_is_better'),
    'NATURAL': ('PN', 'Perceived naturalness', 'higher_is_better'),
}
cols = list(ITEMS)

raw = pd.read_csv(RAW, encoding='utf-8-sig', low_memory=False)
raw.columns = [str(c).strip() for c in raw.columns]
raw = raw[raw['Process'].isin([2, 3])].reset_index(drop=True)   # 2=powder, 3=burger; drop cake
raw['pid'] = np.arange(1, len(raw) + 1)
raw['time'] = (raw['Process'] == 2).astype(int)                 # burger=0 (baseline), powder=1 (endline)

long = raw[['pid', 'time'] + cols].melt(id_vars=['pid', 'time'], var_name='item_label', value_name='y')
long['y'] = pd.to_numeric(long['y'], errors='coerce')
long = long.dropna(subset=['y'])
dp1 = pd.DataFrame({
    'time': long['time'].to_numpy(),
    'time_label': long['time'].map({0: 'Baseline', 1: 'Endline'}).to_numpy(),   # burger=Baseline, powder=Endline (get_endpoints keys on these names)
    'pid': long['pid'].to_numpy(), 'pid_label': long['pid'].astype(str).to_numpy(),
    'fid': np.nan, 'f_label': np.nan, 'submission_date': pd.NaT, 'treat': 0.0,
    'item_label': long['item_label'].to_numpy(), 'y': long['y'].astype(int).to_numpy(),
})
_R = dp1['y'].astype(int)
dp1['y_stan'] = _R                                             # 1..7 model input
dp1['y_label'] = _R.astype(str)                               # rating label (1..7)
dp1['y'] = _R - 1                                             # 0-indexed category (plots)
dp1['item_type'] = 'out-of-7'; dp1['item_type_id'] = 1

# dit
dit = pd.DataFrame({'item_label': cols})
dit['item_type'] = 'out-of-7'; dit['item_type_id'] = 1; dit['cat_length'] = 7
dit['item_label_short'] = dit['item_label'].map(lambda c: ITEMS[c][0])
dit['group_label'] = dit['item_label'].map(lambda c: ITEMS[c][1])
dit['group_label_long'] = dit['group_label']
dit['item_high_label'] = dit['item_label'].map(lambda c: ITEMS[c][2])
dit['endpoint_measure'] = 'mean 7-point rating (powder vs burger)'
dit.to_csv(f"{dir_out}/{file_prefix}_1_data_dit.csv", index=False)

it = (dp1[['item_label', 'time']].drop_duplicates().sort_values(['time', 'item_label']).reset_index(drop=True))
it['item_time_id'] = np.arange(1, len(it) + 1)
dp1 = dp1.merge(it, on=['item_label', 'time'], how='left')

# shuffled accrual so both arms are present at every interim
rng = np.random.default_rng(seed)
pids = rng.permutation(np.sort(dp1.pid.unique())); n_full = len(pids)
grid = np.unique(np.round(np.linspace(max(40, n_full // NINT), n_full, NINT)).astype(int))
n_pow = int((raw['time'] == 1).sum()); n_bur = int((raw['time'] == 0).sum())
print(f"mycelium powder-vs-burger: n_full={n_full} (powder={n_pow}, burger={n_bur}), "
      f"{dp1.item_label.nunique()} items, interims n={grid.tolist()}")

for k, n in enumerate(grid, 1):
    obs = set(pids[:n]); xi = dp1[dp1.pid.isin(obs)].copy()
    xi = xi.sort_values(['item_type_id', 'pid', 'time', 'item_label']).reset_index(drop=True)
    xi['oid'] = range(1, len(xi) + 1); xi['oidt'] = xi.groupby('item_type').cumcount() + 1
    xi.to_csv(f"{dir_out}/{file_prefix}_{k}_data_dp1.csv", index=False)
    pre = f"{dir_out}/{file_prefix}_{k}"
    nb = xi[xi.time == 0].pid.nunique(); npw = xi[xi.time == 1].pid.nunique()
    print(f"\n=== interim {k}: n={n} (burger={nb}, powder={npw}) ===")
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
print("mycelium SVI grid complete ->", dir_out)
