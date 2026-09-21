"""SVI grid for the ImmPort/ImmuneSpace influenza HAI applications (§3.18).
Strains are the ordinal items; day 0 (Baseline) vs the post-vaccination visit
(Endline) are the PAIRED time axis; HAI titre -> ordered category on the 2-fold
dilution ladder (single study-wide K, one item_type_id -> no K-mixing). Endpoint
rho_j = relative rise in mean log2 titre (higher_is_better). No calendar time: a
pseudo-ordering is set by Participant ID, and interims accrue every STEP (=10)
participants. Mirrors ICRC_interim_svi.py; figures as py-icrc-dass-drc_260902.
Per interim: pcm_1_interim_{k}_data_dp1.csv, _draws.zarr, _prob_by_question_fit*,
i{k}_regression_training.pkl, plus 1_data_dit.csv."""
# ---- boilerplate ----

import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
import numpy as np, pandas as pd
from model_pcm import PartialCreditModel
from data_loading import read_data_immport_flu

SB = os.environ.get('IMMPORT_SB', "/Users/or105/sandbox/bIRTistic")
XLSX = os.environ.get('IMMPORT_XLSX',
    "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/"
    "2025_project_Hope_Groups/data/ImmuneSpace_Influenza_Vaccine_Trials_v260918.xlsx")
STUDIES = os.environ.get('IMMPORT_SDY', 'SDY312,SDY314').split(',')
file_prefix = "pcm_1_interim"; x_formula = "~ group - 1"; seed = 123
NSTEPS = int(os.environ.get('IMMPORT_STEPS', '4000')); S = int(os.environ.get('IMMPORT_S', '2000'))
STEP = int(os.environ.get('IMMPORT_STEP', '10'))             # interim every STEP participants
os.environ.setdefault('PROB_FIT_WIDTH_MULT', '1.2')

# Two clinically-standard HAI endpoints per strain (CHMP/CBER; Hobson 1972), computed from
# the fitted category probabilities. y is 0-indexed on the log2 dilution ladder (start 1:5),
# so y>=3 == titre>=1:40 (the ~50% seroprotection threshold).
#   rho 1  SPR   seroprotection RATE at endline, P(titre>=1:40); baseline computed but unused
#                (an absolute endline level, not a change) -> H1 if SPR > 0.70 (CHMP adult)
#   rho 2  GMFR  GMT fold-rise endline/baseline = 2^(mean log2 titre_endline - baseline)
#                -> H1 if fold-rise > 2.5 (CHMP adult)
RHO_SPECS = [
    {'rho_id': 1, 'rho_label': 'SPR: P(titre>=1:40) at endline',
     'reduction': 'threshold', 'threshold': 3, 'compare': 'endline', 'h1_threshold': 0.70},
    {'rho_id': 2, 'rho_label': 'GMT fold-rise (endline/baseline)',
     'reduction': 'mean', 'compare': 'fold_log2', 'h1_threshold': 2.5},
]
H1 = {s['rho_id']: s['h1_threshold'] for s in RHO_SPECS}


def run_study(study):
    dir_out = f"{SB}/py-immport-{study}_260918"; os.makedirs(dir_out, exist_ok=True)
    d = read_data_immport_flu(XLSX, study)
    dp1, dit, K = d['dp'], d['dit'], d['K']
    print(f"\n##### {study}: {dp1.pid.nunique()} paired participants, "
          f"{dp1.item_label.nunique()} strains, K={K}, endline day {d['endline_day']}")
    dit.to_csv(f"{dir_out}/{file_prefix}_1_data_dit.csv", index=False)

    it = (dp1[['item_label', 'group']].drop_duplicates()
          .sort_values(['group', 'item_label']).reset_index(drop=True))
    it['item_group_id'] = np.arange(1, len(it) + 1)
    dp1 = dp1.merge(it, on=['item_label', 'group'], how='left')

    # pseudo-ordering by Participant ID string
    order = dp1[['pid', 'pid_label']].drop_duplicates().sort_values('pid_label')
    pids = order['pid'].to_numpy(); n_full = len(pids)
    grid = sorted(set(list(range(STEP, n_full, STEP)) + [n_full]))
    print(f"interims (every {STEP}, pseudo-order by Participant ID): n={grid}")

    for k, n in enumerate(grid, 1):
        obs = set(pids[:n]); xi = dp1[dp1.pid.isin(obs)].copy()
        xi = xi.sort_values(['item_type_id', 'pid', 'group', 'item_label']).reset_index(drop=True)
        xi['oid'] = range(1, len(xi) + 1); xi['oidt'] = xi.groupby('item_type').cumcount() + 1
        xi.to_csv(f"{dir_out}/{file_prefix}_{k}_data_dp1.csv", index=False)
        pre = f"{dir_out}/{file_prefix}_{k}"
        print(f"\n=== {study} interim {k}: n={n} ===")
        t0 = time.time()
        model = PartialCreditModel(dit=dit, dcati=xi, x_formula=x_formula, seed=seed)
        fit = model.fit_pyro_svi(output_file_prefix=pre, algorithm='AutoDiagonalNormal',
                                 lr=0.01, num_steps=NSTEPS, output_samples=S, resume=True,
                                 with_core_analyses=True, with_additional_analyses=False, verbose=False)
        xr = model.get_endpoints_per_draw(draws=fit['draws'], endpoint_type='items',
                                          rho_specs=RHO_SPECS).rename(columns={'rho': 'pps_rho_x'})
        xr['pps_H1_x'] = (xr['pps_rho_x'] > xr['rho_id'].map(H1)).astype(int)
        xr[['draw', 'item_label', 'item_type', 'item_high_label', 'rho_id', 'rho_label',
            'reduction', 'compare', 'group1', 'group2', 'pps_rho_x', 'pps_H1_x']].to_pickle(
            f"{dir_out}/{file_prefix}_i{k}_regression_training.pkl")
        print(f"  done ({(time.time()-t0)/60:.1f} min)")
    print(f"{study} HAI SVI grid complete ->", dir_out)


if __name__ == "__main__":
    for study in STUDIES:
        run_study(study.strip())
