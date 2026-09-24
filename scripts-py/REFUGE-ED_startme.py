"""START HERE — REFUGE-ED (refugee/migrant youth MSPSS social support) amortiser pipeline (§3.9).

The ONE file a human edits for this project. Config split into three separate dictionaries — DATA
(& analysis/estimand), SVI (the reference fit), AMORTISER (the federated deploy + diagnostics) — same
shape as IMMPORT_flu-SDY269_startme.py. Shared machinery imported from fit_to_current_data +
amortiser_common + amortiser_diag_plots.

Paired baseline/endline design: 12 MSPSS perceived-social-support 1-7 Likert items (expected-score
'out-of-7' endpoint); the estimand is the relative GAIN in support (endline vs baseline, higher-is-
better). Single endpoint -> a one-column diagnostic grid over the 12 MSPSS items.

Run stages via env flags (defaults in brackets):
    RUN_SVI=1    (0)  fit the SVI interim grid first (produces py-refugee_interim_260831)
    RUN_DEPLOY=1 (1)  deploy the amortiser against the SVI fit
    RUN_DIAG=1   (1)  build the combined diagnostic figures
    RUN_DEPLOY=0 RUN_DIAG=1 python scripts-py/REFUGE-ED_startme.py   # diagnostics only
"""
import os, sys
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repo root (parent of scripts-py)
sys.path.insert(0, os.path.join(_REPO, 'python'))
from amortiser_common import federated_deploy
from amortiser_diag_plots import FederatedDiagnostics

SB = os.environ.get('REFUGE_SB', "/Users/or105/sandbox/bIRTistic")    # sandbox root for all output dirs
SVI_DIR = 'py-refugee_interim_260831'   # the SVI output dir (shared by SVI['out'] and AMORTISER['svi'])

# =====================================================================================
# 1) DATA & ANALYSIS — the REFUGE-ED workbook + the perceived-support-gain endpoint (the estimand).
# =====================================================================================
DATA = dict(
    xlsx=os.environ.get('REFUGE_XLSX',                              # the REFUGE-ED Youth Baseline & Endline workbook
        "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/"
        "2025_project_Hope_Groups/data/refuge/Youth Baseline & Endline .xlsx"),
)
# ESTIMAND rho: MSPSS support is higher-is-better, so the relative change reads as a GAIN (endline vs
# baseline). h1 = success threshold on rho (>0.5 would be a >50% gain; REFUGE effects are modest, see eta0).
RHO_SPECS = [dict(rho_id=1, rho_label='support_gain', reduction='mean', compare='ratio', h1=0.5,
                  rho_label_long='MSPSS perceived-support relative gain (endline vs baseline)')]

# =====================================================================================
# 2) SVI — the interim schedule + fit settings producing the reference posterior.
# =====================================================================================
SVI = dict(
    out=SVI_DIR,                  # output dir (per interim: dp1 / draws.zarr / prob-plots / endpoint pkl)
    x_formula='~ group - 1',      # PCM design: one difficulty offset per timepoint (group 0=baseline, 1=endline)
    accrual='sorted',             # accrual in natural participant-id (enrolment) order (linked pairs only)
    n_interims=8, floor=30,       # interim grid = linspace(max(floor, n_full//n_interims), n_full, n_interims)
    seed=123,                     # SVI fit seed
    nsteps=4000, samples=2000,    # SVI optimisation steps + posterior draws saved per interim
    prob_width=1.4,               # width multiplier (12 MSPSS items crowd the prob_by_question_fit panel)
)

# =====================================================================================
# 3) AMORTISER — the federated amortiser deploy + combined-diagnostics config. The endpoint (relative
#    change / mean-ratio) is registry S3, deployed with the scalar scale-feat net, warp=none. Passed
#    verbatim to federated_deploy(SB, AMORTISER) + FederatedDiagnostics(SB, AMORTISER).
# =====================================================================================
AMORTISER = dict(
    svi=SVI_DIR,                  # reference SVI dir the amortiser deploys against (== SVI['out'])
    fed='py-refugee-interim-amortise-deepsetXcompAtt-itemamortise-J64-ftheadexpand-bvm-federated_260924',  # parent dir
    title='REFUGE-ED MSPSS',      # figure-title prefix
    item_kind='MSPSS item',       # response-item word (rows = the 12 MSPSS items)
    mixed_units=False,            # single endpoint -> one facet_grid
    calib_prefix='refuge_mspss_amortiser',   # basename of the calibration summary files
    eta0_units='relative gain x100',         # units label for the eta0 success-threshold axes
    ctag_prefix='refuge-ed-fed',             # deploy-cache tag prefix
    endpoints=[dict(
        rho='support_gain',                  # federated subdir + rho_label the deploy filters on
        instance='S3 rel-change',            # registry net-family label
        build='svi',                         # reference build strategy: reuse the whole SVI fit as-is
        net='scale-feat', widetok=0, case_c=2,     # scalar-mean net + token/caseness settings
        warp='none',                         # signed relative change -> no warp
        eta0=0.10, grid='0,0.05,0.10,0.15,0.20')])  # eta0 sweep 0-20% (REFUGE effects are modest)


def _fit_svi():
    """Run the SVI interim grid for this project (DATA + SVI + RHO_SPECS above)."""
    from data_loading import read_data_refuge_ed
    from fit_to_current_data import (fit_interim_grid, assign_item_group_id, sorted_accrual,
                                     linspace_grid, endpoint_frame)
    raw = read_data_refuge_ed(DATA['xlsx'])
    dp, dit = raw['dp'].copy(), raw['dit'].copy()
    dit['item_type'] = 'out-of-7'; dit['item_type_id'] = 1        # MSPSS Likert -> one expected-score item_type
    both = dp.groupby('pid').group.nunique()                     # keep only participants linked at both times
    dp1 = dp[dp.pid.isin(set(both[both >= 2].index))].copy()
    _r = dp1['y'].astype(int)                                    # raw Likert response 1..7
    dp1['y_stan'] = _r; dp1['y'] = _r - 1                        # y_stan = model input (1..K); y = 0-indexed (plots)
    dp1['item_type'] = 'out-of-7'; dp1['item_type_id'] = 1
    dp1 = assign_item_group_id(dp1)                              # item_group_id = (MSPSS item x timepoint) index
    pids = sorted_accrual(dp1)                                   # accrual order over linked participants
    grid = linspace_grid(len(pids), nint=SVI['n_interims'], floor=SVI['floor'])
    print(f"##### REFUGE-ED: {len(pids)} linked participants, {dp1.item_label.nunique()} MSPSS items")
    h1 = {s['rho_id']: s['h1'] for s in RHO_SPECS}
    fit_interim_grid(f"{SB}/{SVI['out']}", dp1, dit, pids, grid, SVI['x_formula'],
                     lambda m, f, k, n, xi: endpoint_frame(m, f, RHO_SPECS, h1=h1),
                     seed=SVI['seed'], nsteps=SVI['nsteps'], output_samples=SVI['samples'],
                     prob_width=SVI['prob_width'], label='REFUGE-ED')


if __name__ == "__main__":
    if os.environ.get('RUN_SVI', '0') == '1':
        _fit_svi()
    if os.environ.get('RUN_DEPLOY', '1') == '1':
        federated_deploy(SB, AMORTISER, run_diagnostics=False)
    if os.environ.get('RUN_DIAG', '1') == '1':
        FederatedDiagnostics(SB, AMORTISER).run()
