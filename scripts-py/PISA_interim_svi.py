"""SVI fits for the PISA Math application (§3.11), BASELINE-ANCHORED over cycles.
Per country the interims are the SUBSEQUENT cycles vs the 2012 baseline: interim =
(country, endline_cycle in {2015,2018,2022}), Baseline=2012, Endline=cycle. Endpoint
rho = 2012->cycle expected-score change per item. Cross-sectional (unpaired cohorts);
the cycle enters as the `time` covariate with shared (anchored) item difficulties = the
PISA trend model. Dichotomous Math items only (mixing K collapses mean-field SVI).
NOTE: 2012 is paper, 2015+ computer -> the 2012->cycle step carries a mode effect.
Per interim writes dp1.csv, draws.zarr, i{k}_regression_training.pkl, {k}_data_dit.csv."""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
import numpy as np, pandas as pd
from data_loading import read_data_pisa
from model_pcm import PartialCreditModel

SB = "/Users/or105/sandbox/bIRTistic"
DATA = os.environ.get('PISA_DIR',
    "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/2025_project_Hope_Groups/data/pisa")
dir_out = f"{SB}/py-pisa-math_260904"; os.makedirs(dir_out, exist_ok=True)
file_prefix = "pcm_1_interim"; x_formula = "~ time - 1"; seed = 123
NSTEPS = int(os.environ.get('PISA_STEPS', '4000')); S = int(os.environ.get('PISA_S', '1000'))
CAP = int(os.environ.get('PISA_CAP', '2000')); BASE_CYCLE = 2012
ENDLINES = [int(x) for x in os.environ.get('PISA_ENDLINES', '2015,2018,2022').split(',')]
os.environ.setdefault('PROB_FIT_WIDTH_MULT', '1.6')

CNAME = {'ARE': 'United Arab Emirates', 'AUS': 'Australia', 'AUT': 'Austria', 'BEL': 'Belgium',
    'BGR': 'Bulgaria', 'BRA': 'Brazil', 'CAN': 'Canada', 'CHE': 'Switzerland', 'CHL': 'Chile',
    'COL': 'Colombia', 'CRI': 'Costa Rica', 'CZE': 'Czechia', 'DEU': 'Germany', 'DNK': 'Denmark',
    'DOM': 'Dominican Republic', 'ESP': 'Spain', 'EST': 'Estonia', 'FIN': 'Finland', 'FRA': 'France',
    'GBR': 'United Kingdom', 'GRC': 'Greece', 'HKG': 'Hong Kong', 'HRV': 'Croatia', 'HUN': 'Hungary',
    'IRL': 'Ireland', 'ISL': 'Iceland', 'ISR': 'Israel', 'ITA': 'Italy', 'JOR': 'Jordan',
    'JPN': 'Japan', 'KOR': 'Korea', 'LTU': 'Lithuania', 'LVA': 'Latvia', 'MAC': 'Macao',
    'MDA': 'Moldova', 'MEX': 'Mexico', 'MKD': 'North Macedonia', 'MNE': 'Montenegro',
    'NLD': 'Netherlands', 'NOR': 'Norway', 'NZL': 'New Zealand', 'PER': 'Peru', 'POL': 'Poland',
    'PRT': 'Portugal', 'QAT': 'Qatar', 'ROU': 'Romania', 'SGP': 'Singapore', 'SVK': 'Slovakia',
    'SVN': 'Slovenia', 'SWE': 'Sweden', 'TAP': 'Chinese Taipei', 'THA': 'Thailand', 'TUR': 'Turkiye',
    'URY': 'Uruguay', 'USA': 'United States', 'ALB': 'Albania', 'GEO': 'Georgia', 'IDN': 'Indonesia',
    'KSV': 'Kosovo', 'MLT': 'Malta', 'IND': 'India'}

import json
comm = json.load(open(f"{DATA}/common_math.json"))
CLIST = os.environ.get('PISA_COUNTRIES', ','.join(comm['countries'])).split(',')
rng = np.random.default_rng(seed)


def cap_students(sub):
    ids = sub['pid'].unique()
    if len(ids) > CAP:
        ids = rng.choice(ids, CAP, replace=False)
    return sub[sub['pid'].isin(ids)]


rows = []; k = 0
for cnt in CLIST:
    raw = read_data_pisa(DATA, cnt, cycles=tuple([BASE_CYCLE] + ENDLINES))
    dp, dit0 = raw['dp'], raw['dit']
    base = dp[dp.cycle == BASE_CYCLE]
    for endc in ENDLINES:
        end = dp[dp.cycle == endc]
        if base.empty or end.empty:
            print(f"skip {cnt} {endc}: missing cycle"); continue
        both = sorted(set(base.item_label) & set(end.item_label))
        kmax = dp[dp.item_label.isin(both)].groupby('item_label').y.max()
        both = [i for i in both if kmax[i] == 1]                 # dichotomous only
        if len(both) < 20:
            print(f"skip {cnt} {endc}: {len(both)} common dichotomous items"); continue
        a = cap_students(base[base.item_label.isin(both)]); b = cap_students(end[end.item_label.isin(both)])
        k += 1
        xi = pd.concat([a.assign(time=0, time_label='Baseline'),
                        b.assign(time=1, time_label='Endline')], ignore_index=True)
        xi['pid'] = pd.factorize(xi['time'].astype(str) + '_' + xi['pid'].astype(str))[0] + 1
        xi['pid_label'] = xi['pid'].astype(str); xi['fid'] = np.nan; xi['f_label'] = np.nan
        xi['submission_date'] = pd.NaT; xi['treat'] = 0.0
        xi['y_stan'] = xi['y'] + 1; xi['y_label'] = xi['y'].astype(str)
        xi['item_type'] = 'out-of-7'; xi['item_type_id'] = 1
        dit = pd.DataFrame({'item_label': both})
        dit['item_type'] = 'out-of-7'; dit['item_type_id'] = 1; dit['cat_length'] = 2
        dit['item_label_short'] = dit['item_label']
        dit['group_label'] = 'PISA Math'
        dit['group_label_long'] = f'{CNAME.get(cnt, cnt)}: 2012 vs {endc}'
        dit['item_high_label'] = 'higher_is_better'
        dit['endpoint_measure'] = 'mean item score (2012 vs cycle)'
        dit.to_csv(f"{dir_out}/{file_prefix}_{k}_data_dit.csv", index=False)
        it = (xi[['item_label', 'time']].drop_duplicates().sort_values(['time', 'item_label']).reset_index(drop=True))
        it['item_time_id'] = np.arange(1, len(it) + 1)
        xi = xi.merge(it, on=['item_label', 'time'], how='left')
        xi = xi.sort_values(['item_type_id', 'pid', 'time', 'item_label']).reset_index(drop=True)
        xi['oid'] = range(1, len(xi) + 1); xi['oidt'] = xi.groupby('item_type_id').cumcount() + 1
        xi.to_csv(f"{dir_out}/{file_prefix}_{k}_data_dp1.csv", index=False)
        print(f"\n=== [{k}] {cnt} 2012 vs {endc}: nB={a.pid.nunique()} nE={b.pid.nunique()} items={len(both)} ===")
        t0 = time.time()
        model = PartialCreditModel(dit=dit, dcati=xi, x_formula=x_formula, seed=seed)
        fit = model.fit_pyro_svi(output_file_prefix=f"{dir_out}/{file_prefix}_{k}",
                                 algorithm='AutoDiagonalNormal', lr=0.01, num_steps=NSTEPS,
                                 output_samples=S, resume=True,
                                 with_core_analyses=os.environ.get('PISA_PLOTS', '0') == '1',
                                 with_additional_analyses=False, verbose=False)
        xr = model.get_endpoints_per_draw(draws=fit['draws'], categorical_threshold=2,
                                          endpoint_type='items').rename(columns={'ratio': 'pps_ratio_x'})
        xr['pps_H1_x'] = (xr['pps_ratio_x'] > 0.5).astype(int); xr['CNT'] = cnt; xr['cycle'] = endc
        if 'item_high_label' not in xr.columns:
            xr = xr.merge(dit[['item_label', 'item_high_label']], on='item_label', how='left')
        xr[['draw', 'item_label', 'item_type', 'item_high_label', 'pps_ratio_x', 'pps_H1_x', 'CNT', 'cycle']].to_pickle(
            f"{dir_out}/{file_prefix}_i{k}_regression_training.pkl")
        rows.append(dict(k=k, CNT=cnt, cycle=endc, items=len(both), rho_med=float(xr.pps_ratio_x.median())))
        print(f"  rho(2012->{endc}) median={xr.pps_ratio_x.median():+.3f}  ({(time.time()-t0)/60:.1f} min)")
pd.DataFrame(rows).to_csv(f"{dir_out}/{file_prefix}_country_index.csv", index=False)
print(f"PISA Math SVI fits complete ({k} interims) ->", dir_out)
