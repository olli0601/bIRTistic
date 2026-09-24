"""START HERE — ICRC community-MHPSS (DASS-21 / DRC) amortiser pipeline (§3.10).

The ONE file a human edits for this project. Config split into three separate dictionaries — DATA
(& analysis/estimand), SVI (the reference fit), AMORTISER (the federated deploy + diagnostics) — same
shape as IMMPORT_flu-SDY269_startme.py. Shared machinery imported from fit_to_current_data +
amortiser_common + amortiser_diag_plots.

Paired pre/post design: each beneficiary's DASS-21 Depression/Anxiety/Stress subscale total is binned
into the Figure-2 severity levels (K=5); the endpoint is the relative severity REDUCTION (DASS is
lower-is-better). Single endpoint -> a one-column diagnostic grid over the 3 subscales.

Run stages via env flags (defaults in brackets):
    RUN_SVI=1    (0)  fit the SVI interim grid first (produces py-icrc-dass-drc_260902)
    RUN_DEPLOY=1 (1)  deploy the amortiser against the SVI fit
    RUN_DIAG=1   (1)  build the combined diagnostic figures
    RUN_DEPLOY=0 RUN_DIAG=1 python scripts-py/ICRC-DASS-DRC_startme.py   # diagnostics only
"""
import os, sys
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repo root (parent of scripts-py)
sys.path.insert(0, os.path.join(_REPO, 'python'))
from amortiser_common import federated_deploy
from amortiser_diag_plots import FederatedDiagnostics

SB = os.environ.get('ICRC_SB', "/Users/or105/sandbox/bIRTistic")      # sandbox root for all output dirs
SVI_DIR = 'py-icrc-dass-drc_260902'   # the SVI output dir (shared by SVI['out'] and AMORTISER['svi'])

# =====================================================================================
# 1) DATA & ANALYSIS — the Andersen workbook + the severity-reduction endpoint (the estimand).
# =====================================================================================
DATA = dict(
    xlsx=os.environ.get('ICRC_XLSX',                                # the Andersen community-MHPSS workbook
        "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/"
        "2025_project_Hope_Groups/data/Andersen_Community-Level Mental Health Support Table.XLSX"),
    country=2,                    # ctry code to keep: 2 = DRC (the DASS-21 arm)
)
# ESTIMAND rho: DASS-21 distress is lower-is-better, so the direction-aware relative change reads as a
# severity REDUCTION (post vs pre). h1 = success threshold on rho (here >0.5 = a >50% relative reduction).
RHO_SPECS = [dict(rho_id=1, rho_label='distress_reduction', reduction='mean', compare='ratio', h1=0.5,
                  rho_label_long='DASS-21 distress — relative severity reduction (post vs pre)')]

# =====================================================================================
# 2) SVI — the interim schedule + fit settings producing the reference posterior.
# =====================================================================================
SVI = dict(
    out=SVI_DIR,                  # output dir (per interim: dp1 / draws.zarr / prob-plots / endpoint pkl)
    x_formula='~ group - 1',      # PCM design: one difficulty offset per timepoint (group 0=pre, 1=post)
    accrual='sorted',             # accrual in natural beneficiary-id (enrolment) order
    early=[20, 40, 60, 80],       # extra fine early interims (the effect emerges fast), then...
    step=100,                     # ...every `step` beneficiaries + the full cohort
    seed=123,                     # SVI fit seed
    nsteps=4000, samples=2000,    # SVI optimisation steps + posterior draws saved per interim
    prob_width=1.2,               # width multiplier for the prob_by_question_fit figure
)

# =====================================================================================
# 3) AMORTISER — the federated amortiser deploy + combined-diagnostics config. The endpoint (relative
#    change / mean-ratio) is registry S3, deployed with the scalar scale-feat net, warp=none. Passed
#    verbatim to federated_deploy(SB, AMORTISER) + FederatedDiagnostics(SB, AMORTISER).
# =====================================================================================
AMORTISER = dict(
    svi=SVI_DIR,                  # reference SVI dir the amortiser deploys against (== SVI['out'])
    fed='py-icrc-dass-drc-amortise-deepsetXcompAtt-itemamortise-J64-ftheadexpand-bvm-federated_260924',  # parent dir
    title='ICRC DASS-21 (DRC)',   # figure-title prefix
    item_kind='subscale',         # response-item word (rows = the 3 DASS subscales)
    mixed_units=False,            # single endpoint -> one facet_grid
    calib_prefix='icrc_dass_amortiser',    # basename of the calibration summary files
    eta0_units='relative reduction x100',  # units label for the eta0 success-threshold axes
    ctag_prefix='icrc-dass-fed',           # deploy-cache tag prefix
    endpoints=[dict(
        rho='distress_reduction',          # federated subdir + rho_label the deploy filters on
        instance='S3 rel-change',          # registry net-family label
        build='svi',                       # reference build strategy: reuse the whole SVI fit as-is
        net='scale-feat', widetok=0, case_c=2,   # scalar-mean net + token/caseness settings
        warp='none',                       # signed relative change -> no warp
        eta0=0.50, grid='0,0.50,0.70,0.85,0.95')])   # eta0 sweep (large DASS severity shifts)


def _fit_svi():
    """Run the SVI interim grid for this project (DATA + SVI + RHO_SPECS above)."""
    from data_loading import read_data_icrc_dass
    from fit_to_current_data import (fit_interim_grid, assign_item_group_id, sorted_accrual,
                                     step_grid, endpoint_frame)
    d = read_data_icrc_dass(DATA['xlsx'], country=DATA['country'])
    print(f"##### ICRC DASS/DRC: {d['n']} beneficiaries with complete pre&post on all 3 subscales (K={d['K']})")
    dp1 = assign_item_group_id(d['dp'])                  # item_group_id = (subscale x timepoint) difficulty index
    pids = sorted_accrual(dp1); n_full = len(pids)       # accrual order over beneficiaries
    grid = sorted(set([g for g in SVI['early'] if g < n_full] + step_grid(n_full, step=SVI['step'])))
    h1 = {s['rho_id']: s['h1'] for s in RHO_SPECS}
    fit_interim_grid(f"{SB}/{SVI['out']}", dp1, d['dit'], pids, grid, SVI['x_formula'],
                     lambda m, f, k, n, xi: endpoint_frame(m, f, RHO_SPECS, h1=h1),
                     seed=SVI['seed'], nsteps=SVI['nsteps'], output_samples=SVI['samples'],
                     prob_width=SVI['prob_width'], label='ICRC')


if __name__ == "__main__":
    if os.environ.get('RUN_SVI', '0') == '1':
        _fit_svi()
    if os.environ.get('RUN_DEPLOY', '1') == '1':
        federated_deploy(SB, AMORTISER, run_diagnostics=False)
    if os.environ.get('RUN_DIAG', '1') == '1':
        FederatedDiagnostics(SB, AMORTISER).run()
