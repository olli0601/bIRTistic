"""START HERE — mycelium novel-food acceptance (powder-vs-burger) amortiser pipeline (§3.14 / §21).

The ONE file a human edits for this project. Config split into three separate dictionaries — DATA
(& analysis/estimand), SVI (the reference fit), AMORTISER (the federated deploy + diagnostics) — same
shape as CAVD-hvtn505_startme.py (its nearest analog: a between-arm arm=time contrast). Shared machinery
imported from fit_to_current_data + amortiser_common + amortiser_diag_plots.

Between-arm "Option A" design: the 3x3 survey is collapsed to product powder vs burger (substrates
pooled), encoded MYCELIUM-style with the arm as time — burger = Baseline (group 0), powder = Endline
(group 1); each respondent is in one arm (UNPAIRED). The endpoint is the direction-aware relative mean
shift powder-vs-burger per item (9 items, 3 constructs). Single endpoint -> a one-column diagnostic grid.

Run stages via env flags (defaults in brackets):
    RUN_SVI=1    (0)  fit the SVI interim grid first (produces py-mycelium-powdervsburger_260902)
    RUN_DEPLOY=1 (1)  deploy the amortiser against the SVI fit
    RUN_DIAG=1   (1)  build the combined diagnostic figures
    RUN_SVI=1 python scripts-py/MYCELIUM_startme.py                        # full pipeline from scratch
    RUN_DEPLOY=0 RUN_DIAG=1 python scripts-py/MYCELIUM_startme.py          # diagnostics only
"""
import os, sys
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repo root (parent of scripts-py)
sys.path.insert(0, os.path.join(_REPO, 'python'))
from amortiser_common import federated_deploy
from amortiser_diag_plots import FederatedDiagnostics

SB = os.environ.get('MYCELIUM_SB', "/Users/or105/sandbox/bIRTistic")   # sandbox root for all output dirs
SVI_DIR = 'py-mycelium-powdervsburger_260902'   # the SVI output dir (shared by SVI['out'] and AMORTISER['svi'])

# =====================================================================================
# 1) DATA & ANALYSIS — the mycelium survey CSV + the powder-vs-burger contrast (the estimand).
# =====================================================================================
DATA = dict(
    csv=os.environ.get('MYCELIUM_CSV',                             # the UK Prolific mycelium acceptance survey
        "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/"
        "2025_project_Hope_Groups/data/mycelium/Mycelium.csv"),
)
# ESTIMAND rho: the direction-aware relative mean shift powder-vs-burger, per item (Acceptance/Naturalness
# up = good, Disgust down = good; item_high_label orients each so higher = better). h1 = success threshold
# (>0.5 = a >50% relative gain; the observed effects are modest ~0.1-0.27, see the low eta0 sweep).
RHO_SPECS = [dict(rho_id=1, rho_label='powder_vs_burger', reduction='mean', compare='ratio', h1=0.5,
                  rho_label_long='powder-vs-burger relative effect (direction-aware, per construct)')]

# =====================================================================================
# 2) SVI — the interim schedule + fit settings producing the reference posterior.
# =====================================================================================
SVI = dict(
    out=SVI_DIR,                  # output dir (per interim: dp1 / draws.zarr / prob-plots / endpoint pkl)
    x_formula='~ group - 1',      # PCM design: one difficulty offset per arm (group 0=burger, 1=powder)
    accrual='shuffled',           # shuffled accrual so both arms are present at every interim
    n_interims=8, floor=40,       # interim grid = linspace(max(floor, n_full//n_interims), n_full, n_interims)
    seed=123,                     # RNG for the shuffled accrual + the SVI fit
    nsteps=4000, samples=2000,    # SVI optimisation steps + posterior draws saved per interim
    prob_width=1.3,               # width multiplier for the prob_by_question_fit figure
)

# =====================================================================================
# 3) AMORTISER — the federated amortiser deploy + combined-diagnostics config. The endpoint (relative
#    change / mean-ratio) is registry S3, deployed with the scalar scale-feat net, warp=none. Passed
#    verbatim to federated_deploy(SB, AMORTISER) + FederatedDiagnostics(SB, AMORTISER).
# =====================================================================================
AMORTISER = dict(
    svi=SVI_DIR,                  # reference SVI dir the amortiser deploys against (== SVI['out'])
    fed='py-mycelium-powdervsburger-amortise-deepsetXcompAtt-itemamortise-J64-ftheadexpand-bvm-federated_260924',  # parent dir
    title='Mycelium powder-vs-burger',   # figure-title prefix
    item_kind='construct item',   # response-item word (rows = the 9 acceptance/disgust/naturalness items)
    mixed_units=False,            # single endpoint -> one facet_grid
    calib_prefix='mycelium_amortiser',   # basename of the calibration summary files
    eta0_units='relative change x100',   # units label for the eta0 success-threshold axes
    ctag_prefix='mycelium-pvb-fed',      # deploy-cache tag prefix
    endpoints=[dict(
        rho='powder_vs_burger',          # federated subdir + rho_label the deploy filters on
        instance='S3 rel-change',        # registry net-family label
        build='svi',                     # reference build strategy: reuse the whole SVI fit as-is
        net='scale-feat', widetok=0, case_c=2,     # scalar-mean net + token/caseness settings
        warp='none',                     # signed relative change -> no warp
        eta0=0.10, grid='0,0.10,0.20,0.30')])   # eta0 sweep over the modest observed effect range


def _fit_svi():
    """Run the SVI interim grid for this project (DATA + SVI + RHO_SPECS above)."""
    from data_loading import read_data_mycelium_powdervsburger
    from fit_to_current_data import (fit_interim_grid, assign_item_group_id, shuffled_accrual,
                                     linspace_grid, endpoint_frame)
    d = read_data_mycelium_powdervsburger(DATA['csv'])
    print(f"##### mycelium powder-vs-burger: n={d['n']} (powder={d['n_powder']}, burger={d['n_burger']}), "
          f"{len(d['items'])} items (K={d['K']})")
    dp1 = assign_item_group_id(d['dp'])                   # item_group_id = (item x arm) difficulty index
    pids = shuffled_accrual(dp1, seed=SVI['seed'])        # accrual order over respondents
    grid = linspace_grid(len(pids), nint=SVI['n_interims'], floor=SVI['floor'])
    h1 = {s['rho_id']: s['h1'] for s in RHO_SPECS}
    fit_interim_grid(f"{SB}/{SVI['out']}", dp1, d['dit'], pids, grid, SVI['x_formula'],
                     lambda m, f, k, n, xi: endpoint_frame(m, f, RHO_SPECS, h1=h1),
                     seed=SVI['seed'], nsteps=SVI['nsteps'], output_samples=SVI['samples'],
                     prob_width=SVI['prob_width'], label='mycelium')


if __name__ == "__main__":
    if os.environ.get('RUN_SVI', '0') == '1':
        _fit_svi()
    if os.environ.get('RUN_DEPLOY', '1') == '1':
        federated_deploy(SB, AMORTISER, run_diagnostics=False)
    if os.environ.get('RUN_DIAG', '1') == '1':
        FederatedDiagnostics(SB, AMORTISER).run()
