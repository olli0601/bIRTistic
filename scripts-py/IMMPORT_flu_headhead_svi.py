"""SVI grid for the ImmPort/ImmuneSpace influenza HEAD-TO-HEAD vaccine application (§3.21),
SDY269 (2008 Systems Biology influenza: LAIV vs TIV) — the JOINT-FIT formulation.

Both vaccine arms are fitted in ONE partial-credit model per interim (not two separate per-arm
fits). The item is the (strain, arm) pair, so each (strain, arm, time) difficulty gets its own
`item_group_id` — the full, non-reduced structure — while the participant abilities are shared
across arms (participants are disjoint people). From this single fit the endpoints are read with
FOUR `get_endpoints_per_draw` calls:
  1. SPR       seroprotection rate P(titre>=1:40) at endline, per (strain, arm)
  2. GMFR      GMT fold-rise endline/baseline, per (strain, arm)
  3. SPR_diff  TIV-minus-LAIV SPR on the SHARED strain (A/Uruguay H3N2), per draw
  4. GMFR_diff TIV-minus-LAIV GMFR on the shared strain, per draw
The two _diff endpoints use the new cross-`arm` capability of `get_endpoints_per_draw`
(`across='arm'`, `across_within='strain'`): because both arms share the fit, the difference is an
EXACT per-draw contrast (no random draw-pairing of independent fits). One joint fit also gives the
amortiser a single target that already carries the between-arm contrast.

Output layout: ONE dir `py-immport-SDY269`. The joint fit files (dp1, draws, prob_by_question_fit)
under prefix `pcm_1_interim_`; the four endpoint posteriors under
`pcm_{SPR,GMFR,SPR_diff,GMFR_diff}_interim_i{k}_regression_training.pkl`. Small arms (~28 paired
each, 56 pooled), interims every STEP=10 pooled participants (arms interleaved so both accrue)."""
# ---- boilerplate ----

import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
import numpy as np, pandas as pd
from model_pcm import PartialCreditModel
from data_loading import read_data_immport_flu_headhead

SB = os.environ.get('IMMPORT_SB', "/Users/or105/sandbox/bIRTistic")
XLSX = os.environ.get('IMMPORT_XLSX',
    "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/"
    "2025_project_Hope_Groups/data/ImmuneSpace_Influenza_Vaccine_Trials_v260918.xlsx")
STUDY = os.environ.get('IMMPORT_SDY', 'SDY269')
ARMS = tuple(os.environ.get('IMMPORT_ARMS', 'LAIV,TIV').split(','))
DIR = f"{SB}/py-immport-{STUDY}"; os.makedirs(DIR, exist_ok=True)
# item_group_id = item_label(strain) x group(arm x phase) carries the fit; the design uses only the
# within-arm paired axis `phase` (0/1) -> fit invariant to the arm-in-group relabelling. Endpoints
# pivot on `phase` and keep `arm` as the stratum.
file_prefix = "pcm_1_interim"; x_formula = "~ phase - 1"; seed = 123
NSTEPS = int(os.environ.get('IMMPORT_STEPS', '4000')); S = int(os.environ.get('IMMPORT_S', '2000'))
STEP = int(os.environ.get('IMMPORT_STEP', '10'))
os.environ.setdefault('PROB_FIT_WIDTH_MULT', '1.2')

# Four endpoints from the one joint fit. y 0-indexed on the log2 ladder (start 1:5) -> y>=3 ==
# titre>=1:40. The _diff specs contrast arm (ref=LAIV, foc=TIV -> TIV-LAIV) within strain, so they
# collapse to the shared strain automatically (non-shared strains have a single arm -> dropped).
REF, FOC = ARMS[0], ARMS[1]                                  # LAIV, TIV
ENDPOINTS = [
    ('SPR',       {'rho_id': 1, 'rho_label': 'SPR: P(titre>=1:40) at endline',
                   'reduction': 'threshold', 'threshold': 3, 'compare': 'endline'}, 0.70),
    ('GMFR',      {'rho_id': 2, 'rho_label': 'GMT fold-rise (endline/baseline)',
                   'reduction': 'mean', 'compare': 'fold_log2'}, 2.5),
    ('SPR_diff',  {'rho_id': 3, 'rho_label': f'SPR difference ({FOC}-{REF}), shared strain',
                   'reduction': 'threshold', 'threshold': 3, 'compare': 'endline',
                   'across': 'arm', 'across_within': 'item_label', 'across_values': [REF, FOC],
                   'across_compare': 'diff'}, 0.0),
    ('GMFR_diff', {'rho_id': 4, 'rho_label': f'GMFR difference ({FOC}-{REF}), shared strain',
                   'reduction': 'mean', 'compare': 'fold_log2',
                   'across': 'arm', 'across_within': 'item_label', 'across_values': [REF, FOC],
                   'across_compare': 'diff'}, 0.0),
]
PKL_COLS = ['draw', 'arm', 'item_label', 'item_type', 'item_high_label', 'rho_id',
            'rho_label', 'reduction', 'compare', 'group1', 'group2', 'pps_rho_x', 'pps_H1_x']


def _interleave(dp1):
    """Pooled accrual order that interleaves the arms (both present at every interim)."""
    seq = []
    per_arm = [dp1[dp1.arm == a][['pid', 'pid_label']].drop_duplicates().sort_values('pid_label')['pid'].tolist()
               for a in ARMS]
    for tup in __import__('itertools').zip_longest(*per_arm):
        seq += [p for p in tup if p is not None]
    return np.array(seq)


def run():
    d = read_data_immport_flu_headhead(XLSX, STUDY, arms=ARMS)
    dp1, dit, K, shared = d['dp'], d['dit'], d['K'], d['shared']
    ngrp = dp1[['item_label', 'group']].drop_duplicates().shape[0]
    print(f"\n##### {STUDY} head-to-head JOINT fit: arms {ARMS}, {dp1.pid.nunique()} pooled paired "
          f"participants, {dit.item_label.nunique()} strain items x {dp1.group.nunique()} conditions "
          f"(arm x phase) = {ngrp} item_group_id, K={K}, endline day {d['endline_day']}\n"
          f"  shared strain(s): {shared}")
    dit.to_csv(f"{DIR}/{file_prefix}_1_data_dit.csv", index=False)

    it = (dp1[['item_label', 'group']].drop_duplicates()
          .sort_values(['group', 'item_label']).reset_index(drop=True))
    it['item_group_id'] = np.arange(1, len(it) + 1)         # (strain, arm, time) -> non-reduced difficulties
    dp1 = dp1.merge(it, on=['item_label', 'group'], how='left')

    pids = _interleave(dp1); n_full = len(pids)
    grid = sorted(set(list(range(STEP, n_full, STEP)) + [n_full]))
    print(f"interims (every {STEP} pooled, arms interleaved): n={grid}")

    for k, n in enumerate(grid, 1):
        obs = set(pids[:n]); xi = dp1[dp1.pid.isin(obs)].copy()
        xi = xi.sort_values(['item_type_id', 'pid', 'group', 'item_label']).reset_index(drop=True)
        xi['oid'] = range(1, len(xi) + 1); xi['oidt'] = xi.groupby('item_type').cumcount() + 1
        xi.to_csv(f"{DIR}/{file_prefix}_{k}_data_dp1.csv", index=False)
        pre = f"{DIR}/{file_prefix}_{k}"
        na = {a: xi[xi.arm == a].pid.nunique() for a in ARMS}
        print(f"\n=== {STUDY} interim {k}: n={n} ({na}) ===")
        t0 = time.time()
        model = PartialCreditModel(dit=dit, dcati=xi, x_formula=x_formula, seed=seed)
        fit = model.fit_pyro_svi(output_file_prefix=pre, algorithm='AutoDiagonalNormal',
                                 lr=0.01, num_steps=NSTEPS, output_samples=S, resume=True,
                                 with_core_analyses=True, with_additional_analyses=False, verbose=False)
        for tag, spec, h1 in ENDPOINTS:                     # four get_endpoints calls off the ONE fit
            xr = model.get_endpoints_per_draw(draws=fit['draws'], endpoint_type='items',
                                              rho_specs=[spec], contrast_col='phase'
                                              ).rename(columns={'rho': 'pps_rho_x'})
            xr['pps_H1_x'] = (xr['pps_rho_x'] > h1).astype(int)
            xr[PKL_COLS].to_pickle(f"{DIR}/pcm_{tag}_interim_i{k}_regression_training.pkl")
        print(f"  done ({(time.time()-t0)/60:.1f} min)")
    print(f"{STUDY} head-to-head JOINT SVI grid -> {DIR}  (fit pcm_1_interim_, endpoints "
          f"pcm_{{SPR,GMFR,SPR_diff,GMFR_diff}}_interim_)")


if __name__ == "__main__":
    run()
