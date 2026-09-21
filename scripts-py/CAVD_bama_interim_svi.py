"""SVI grid for the CAVD DataSpace HIV-vaccine application (§3.19), HVTN 505 (vtn505),
BAMA binding-IgG. HIV vaccine trials have no informative paired baseline (HIV-naive
floor), so the estimand is the single-timepoint (day-196) VACCINE-vs-PLACEBO contrast,
encoded MYCELIUM-style: placebo -> time 0 ('Baseline'), vaccine -> time 1 ('Endline');
each subject in one arm. Antigens are the ordinal items; mfi_delta binned per antigen into
K ordered levels; rho_j = relative vaccine-vs-placebo shift in binding (higher_is_better).
Subjects accrue in shuffled order (so both arms are present at every interim); SVI at each.
Mirrors MYCELIUM_interim_svi.py; figures as py-icrc-dass-drc_260902. Per interim:
pcm_1_interim_{k}_data_dp1.csv, _draws.zarr, _prob_by_question_fit*, i{k}_regression_training.pkl,
plus 1_data_dit.csv.

Data must be pulled first: python python/data_web_extracting.py --cavd "HVTN 505" --out <DATADIR>."""
# ---- boilerplate ----

import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
import numpy as np, pandas as pd
from model_pcm import PartialCreditModel
from data_loading import read_data_cavd_bama_endline

SB = os.environ.get('CAVD_SB', "/Users/or105/sandbox/bIRTistic")
DATADIR = os.environ.get('CAVD_DATADIR',
    "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/"
    "2025_project_Hope_Groups/data/cavd_dataspace")
STUDIES = os.environ.get('CAVD_STUDY', 'vtn505').split(',')
file_prefix = "pcm_1_interim"; x_formula = "~ group - 1"; seed = 123
NSTEPS = int(os.environ.get('CAVD_STEPS', '4000')); S = int(os.environ.get('CAVD_S', '2000'))
NINT = int(os.environ.get('CAVD_NINT', '10')); KCAT = int(os.environ.get('CAVD_KCAT', '3'))
os.environ.setdefault('PROB_FIT_WIDTH_MULT', '1.3')


def run_study(prot):
    dir_out = f"{SB}/py-cavd-{prot}-bama_260918"; os.makedirs(dir_out, exist_ok=True)
    d = read_data_cavd_bama_endline(f"{DATADIR}/cavd_{prot}_BAMA.parquet",
                                    f"{DATADIR}/cavd_{prot}_Demographics.parquet", n_cat=KCAT)
    dp1, dit, K = d['dp'], d['dit'], d['K']
    print(f"\n##### {prot} BAMA between-arm: {dp1.pid.nunique()} subjects, "
          f"{len(d['kept'])} antigens (K={K}, day {d['endline_day']}); dropped {len(d['dropped'])} "
          f"(arm lacked full category coverage)")
    dit.to_csv(f"{dir_out}/{file_prefix}_1_data_dit.csv", index=False)

    it = (dp1[['item_label', 'group']].drop_duplicates()
          .sort_values(['group', 'item_label']).reset_index(drop=True))
    it['item_group_id'] = np.arange(1, len(it) + 1)
    dp1 = dp1.merge(it, on=['item_label', 'group'], how='left')

    # shuffled accrual so both arms (placebo/vaccine) are present at every interim
    rng = np.random.default_rng(seed)
    pids = rng.permutation(np.sort(dp1.pid.unique())); n_full = len(pids)
    grid = np.unique(np.round(np.linspace(max(40, n_full // NINT), n_full, NINT)).astype(int))
    nv = int((dp1.drop_duplicates('pid').group == 1).sum()); npl = dp1.pid.nunique() - nv
    print(f"vaccine={nv}, placebo={npl}; interims n={grid.tolist()}")

    for k, n in enumerate(grid, 1):
        obs = set(pids[:n]); xi = dp1[dp1.pid.isin(obs)].copy()
        xi = xi.sort_values(['item_type_id', 'pid', 'group', 'item_label']).reset_index(drop=True)
        xi['oid'] = range(1, len(xi) + 1); xi['oidt'] = xi.groupby('item_type').cumcount() + 1
        xi.to_csv(f"{dir_out}/{file_prefix}_{k}_data_dp1.csv", index=False)
        pre = f"{dir_out}/{file_prefix}_{k}"
        nvi = xi[xi.group == 1].pid.nunique(); npli = xi[xi.group == 0].pid.nunique()
        print(f"\n=== {prot} interim {k}: n={n} (vaccine={nvi}, placebo={npli}) ===")
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
    print(f"{prot} BAMA between-arm SVI grid complete ->", dir_out)


if __name__ == "__main__":
    for prot in STUDIES:
        run_study(prot.strip())
