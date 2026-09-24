"""SVI grid for the ImmPort/ImmuneSpace COVID-19 subpopulation-contrast application (§3.20),
SDY1764 (Distinct antibody responses to SARS-CoV-2 in children and adults). BETWEEN-GROUP
design (MYCELIUM/CAVD-BAMA style): the subpopulation is the group axis — group A -> time 0
('Baseline'), group B -> time 1 ('Endline'); each participant in one group. Item = SARS-CoV-2
serum-neutralisation ID50 titre -> ordered category on the 2-fold ladder. rho = relative
group-B-vs-A shift in mean log2 titre (higher_is_better). Two contrasts via env COVID_GROUP:
'age' (pediatric vs adult) or 'severity' (severe vs mild). Subjects accrue in shuffled order
(both groups present at each interim); SVI at each. Mirrors the folded startme SVI pattern (fit_to_current_data); figures
as py-icrc-dass-drc_260902. Data: ImmuneSpace_COVID19_v260920.xlsx."""
# ---- boilerplate ----

import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
import numpy as np, pandas as pd
from model_pcm import PartialCreditModel
from data_loading import read_data_immport_covid_neut

SB = os.environ.get('COVID_SB', "/Users/or105/sandbox/bIRTistic")
XLSX = os.environ.get('COVID_XLSX',
    "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/"
    "2025_project_Hope_Groups/data/ImmuneSpace_COVID19_v260920.xlsx")
STUDY = os.environ.get('COVID_SDY', 'SDY1764')
GROUPS = os.environ.get('COVID_GROUP', 'age,severity').split(',')
file_prefix = "pcm_1_interim"; x_formula = "~ group - 1"; seed = 123
NSTEPS = int(os.environ.get('COVID_STEPS', '4000')); S = int(os.environ.get('COVID_S', '2000'))
NINT = int(os.environ.get('COVID_NINT', '8'))
os.environ.setdefault('PROB_FIT_WIDTH_MULT', '1.3')


def run_contrast(group):
    dir_out = f"{SB}/py-immport-covid-{STUDY}-{group}_260920"; os.makedirs(dir_out, exist_ok=True)
    d = read_data_immport_covid_neut(XLSX, STUDY, group=group)
    dp1, dit, K = d['dp'], d['dit'], d['K']
    gA, gB = d['groups']
    print(f"\n##### {STUDY} {group}: {dp1.pid.nunique()} subjects, K={K}, contrast {gB} vs {gA}; "
          f"group sizes {dp1.groupby('group_label').pid.nunique().to_dict()}")
    dit.to_csv(f"{dir_out}/{file_prefix}_1_data_dit.csv", index=False)

    it = (dp1[['item_label', 'group']].drop_duplicates()
          .sort_values(['group', 'item_label']).reset_index(drop=True))
    it['item_group_id'] = np.arange(1, len(it) + 1)
    dp1 = dp1.merge(it, on=['item_label', 'group'], how='left')

    rng = np.random.default_rng(seed)
    pids = rng.permutation(np.sort(dp1.pid.unique())); n_full = len(pids)
    grid = np.unique(np.round(np.linspace(max(20, n_full // NINT), n_full, NINT)).astype(int))
    print(f"interims (shuffled accrual, both groups present): n={grid.tolist()}")

    for k, n in enumerate(grid, 1):
        obs = set(pids[:n]); xi = dp1[dp1.pid.isin(obs)].copy()
        xi = xi.sort_values(['item_type_id', 'pid', 'group', 'item_label']).reset_index(drop=True)
        xi['oid'] = range(1, len(xi) + 1); xi['oidt'] = xi.groupby('item_type').cumcount() + 1
        xi.to_csv(f"{dir_out}/{file_prefix}_{k}_data_dp1.csv", index=False)
        pre = f"{dir_out}/{file_prefix}_{k}"
        nb = xi[xi.group == 0].pid.nunique(); ne = xi[xi.group == 1].pid.nunique()
        print(f"\n=== {group} interim {k}: n={n} ({gA}={nb}, {gB}={ne}) ===")
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
    print(f"{STUDY} {group} SVI grid complete ->", dir_out)


if __name__ == "__main__":
    for g in GROUPS:
        run_contrast(g.strip())
