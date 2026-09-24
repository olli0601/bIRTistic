"""START HERE — CAVD HVTN 505 BAMA amortiser pipeline (§3.19).

The ONE file a human edits for this project. Config is split into three separate dictionaries —
DATA (& analysis/estimand), SVI (the reference fit), AMORTISER (the federated deploy + diagnostics) —
so each concern is edited on its own. The shared machinery is imported: fit_to_current_data (SVI loop),
amortiser_common.federated_deploy (deploy), amortiser_diag_plots.FederatedDiagnostics (diagnostics).

Run stages are selected by env flags (each defaults as shown):
    RUN_SVI=1    (default 0)  fit the SVI interim grid first (slow; usually done once)
    RUN_DEPLOY=1 (default 1)  deploy the amortiser against the SVI fit
    RUN_DIAG=1   (default 1)  build the combined diagnostic figures
e.g.   RUN_SVI=1 python scripts-py/CAVD-hvtn505_startme.py           # full pipeline from scratch
       RUN_DEPLOY=0 RUN_DIAG=1 python scripts-py/CAVD-hvtn505_startme.py   # re-make diagnostics only

Data must be pulled first: python python/data_web_extracting.py --cavd "HVTN 505" --out <CAVD_DATADIR>.
"""
import os, sys
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repo root (parent of scripts-py)
sys.path.insert(0, os.path.join(_REPO, 'python'))
from amortiser_common import federated_deploy
from amortiser_diag_plots import FederatedDiagnostics

SB = os.environ.get('CAVD_SB', "/Users/or105/sandbox/bIRTistic")      # sandbox root: all output dirs live under here
SVI_DIR = 'py-cavd-vtn505-bama_260918'   # the SVI output dir (shared by SVI['out'] and AMORTISER['svi'])

# =====================================================================================
# 1) DATA & ANALYSIS — where the raw data is, how it is loaded/binned, and the ESTIMAND
#    (the rho endpoint every downstream stage targets).
# =====================================================================================
DATA = dict(
    datadir=os.environ.get('CAVD_DATADIR',                           # folder of the pulled DataSpace parquet files
        "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/"
        "2025_project_Hope_Groups/data/cavd_dataspace"),
    protocol='vtn505',            # HVTN study id -> reads cavd_<protocol>_{BAMA,Demographics}.parquet
    n_cat=3,                      # ordinal levels K: mfi_delta binned per antigen into K tertiles (low/interm/high)
    complete_panel=True,          # keep only the rectangular 9-antigen panel measured on ALL subjects (drop sparse)
)
# ESTIMAND rho, declared upfront (its short label = the federated subdir; long label = plot strip). HIV
# vaccine trials have no informative paired baseline, so rho = the single-timepoint vaccine-vs-placebo
# relative mean shift (MYCELIUM-style arm=time). h1 = the success threshold on rho used as the SVI decision.
RHO_SPECS = [dict(rho_id=1, rho_label='vaccine_vs_placebo', reduction='mean', compare='ratio', h1=0.5,
                  rho_label_long='BAMA IgG binding response — vaccine-vs-placebo relative effect')]

# =====================================================================================
# 2) SVI — the interim-analysis schedule + fit settings that PRODUCE the reference posterior.
# =====================================================================================
SVI = dict(
    out=SVI_DIR,                  # output dir (per interim: dp1 / draws.zarr / prob-plots / endpoint pkl)
    x_formula='~ group - 1',      # PCM design: one difficulty offset per arm (group 0=placebo, 1=vaccine)
    accrual='shuffled',           # participant accrual order across interims (both arms present each interim)
    n_interims=10, floor=40,      # interim grid = linspace(max(floor, n_full//n_interims), n_full, n_interims)
    seed=123,                     # RNG for the shuffled accrual + the SVI fit
    nsteps=4000, samples=2000,    # SVI optimisation steps + number of posterior draws saved per interim
    prob_width=1.3,               # width multiplier for the prob_by_question_fit figure
)

# =====================================================================================
# 3) AMORTISER — the federated amortiser deploy + combined-diagnostics config (nets, warps, eta0 sweep).
#    Passed verbatim to federated_deploy(SB, AMORTISER) and FederatedDiagnostics(SB, AMORTISER).
# =====================================================================================
AMORTISER = dict(
    svi=SVI_DIR,                  # the reference SVI dir the amortiser is deployed against (== SVI['out'])
    fed='py-cavd-vtn505-bama-amortise-deepsetXcompAtt-itemamortise-J64-ftheadexpand-bvm-federated_260923',  # parent dir
    title='CAVD HVTN 505 BAMA',   # figure-title prefix
    item_kind='antigen',          # the response-item word (rows of the diagnostic grid)
    mixed_units=False,            # single endpoint -> one facet_grid (no per-endpoint stitching needed)
    calib_prefix='cavd_bama_amortiser',   # basename of the calibration summary files
    eta0_units='relative change x100',     # units label for the eta0 success-threshold axes
    ctag_prefix='cavd-vtn505-fed',         # deploy-cache tag prefix (keeps driver cells distinct)
    endpoints=[dict(                       # one delegated amortiser per endpoint (here: just the one rho)
        rho='vaccine_vs_placebo',          # federated subdir + rho_label the deploy filters on
        instance='S3 rel-change',          # registry net-family label (shown in the calibration plot)
        build='svi',                       # reference build strategy: reuse the whole SVI fit as-is
        net='scale-feat', widetok=0, case_c=2,   # trained net + its token/caseness settings
        warp='none',                       # target warp (none for a signed relative change)
        eta0=0.5, grid='0.0,0.25,0.5,0.75')])    # eta0 success sweep (default anchor + grid, endpoint units)


def _fit_svi():
    """Run the SVI interim grid for this project (the DATA + SVI + RHO_SPECS config above)."""
    from data_loading import read_data_cavd_bama_endline
    from fit_to_current_data import (fit_interim_grid, assign_item_group_id, shuffled_accrual,
                                     linspace_grid, endpoint_frame)
    d = read_data_cavd_bama_endline(f"{DATA['datadir']}/cavd_{DATA['protocol']}_BAMA.parquet",
                                    f"{DATA['datadir']}/cavd_{DATA['protocol']}_Demographics.parquet",
                                    n_cat=DATA['n_cat'], complete_panel=DATA['complete_panel'])
    print(f"##### {DATA['protocol']} BAMA between-arm: {d['dp'].pid.nunique()} subjects, {len(d['kept'])} "
          f"antigens (K={d['K']}, day {d['endline_day']}); dropped {len(d['dropped'])}")
    dp1 = assign_item_group_id(d['dp'])                   # add item_group_id = (antigen x arm) difficulty index
    pids = shuffled_accrual(dp1, seed=SVI['seed'])        # accrual order over subjects
    grid = linspace_grid(len(pids), nint=SVI['n_interims'], floor=SVI['floor'])   # the interim sizes
    h1 = {s['rho_id']: s['h1'] for s in RHO_SPECS}        # per-rho success threshold for pps_H1_x
    fit_interim_grid(f"{SB}/{SVI['out']}", dp1, d['dit'], pids, grid, SVI['x_formula'],
                     lambda m, f, k, n, xi: endpoint_frame(m, f, RHO_SPECS, h1=h1),
                     seed=SVI['seed'], nsteps=SVI['nsteps'], output_samples=SVI['samples'],
                     prob_width=SVI['prob_width'], label=DATA['protocol'])


if __name__ == "__main__":
    if os.environ.get('RUN_SVI', '0') == '1':
        _fit_svi()
    if os.environ.get('RUN_DEPLOY', '1') == '1':
        federated_deploy(SB, AMORTISER, run_diagnostics=False)
    if os.environ.get('RUN_DIAG', '1') == '1':
        FederatedDiagnostics(SB, AMORTISER).run()
