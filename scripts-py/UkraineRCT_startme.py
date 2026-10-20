"""START HERE — Ukraine Hope Groups RANDOMISED-TRIAL amortiser pipeline (UkraineRCT).

The ONE file a human edits for the corrected Ukraine estimand. Same three-dictionary shape as
CAVD-hvtn505_startme.py / UKRAINE_startme.py — DATA (& estimand), SVI (reference fit), AMORTISER
(federated deploy + diagnostics manifest). It supersedes the parked pre->post (§14) and cross-arm
(§15 UkraineP) framings: the control arm received NO intervention between baseline and endline, so
the estimand is the GROUP-MEAN comparison AT ENDLINE of the intervention arm vs the control arm, per
item, as a direction-aware relative percent improvement (higher = better), randomised WITHIN
facilitator (see data_loading.get_data_ukraine_rct + scripts-py/UkraineRCT_interim_svi.py).

Estimand = signed relative change intervention-vs-control at endline. It CAN be negative (a harmful
item), so log2/logit target warps are undefined => warp='none' (the same choice as every other Ukraine
/ mycelium / CAVD between-arm S3 endpoint; §16.3). No rho transformation is applied.

The reference fit is BETWEEN-ARM (disjoint participants, unpaired PCM) with BIWEEKLY interims, produced
by UkraineRCT_interim_svi.py — NOT the weekly paired loop of UKRAINE_startme.py (get_interim_data_x's
both-timepoint pairing does not apply here). So RUN_SVI here delegates to that script.

AMORTISER deploys the amortiser (deepsetScompAtt itemamortise J64) with the BETWEEN-ARM scale-feat net
(scale-feat-between): the same rel-change target rho=we/wb-1 as the paired scale-feat net, but trained on
a BETWEEN-SUBJECTS cohort (one arm slot per participant) to match this disjoint-arm design — vs the
paired scale-feat-paired net used by the §14 within-participant apps. ftheadexpand + STEIN BvM: James-
Stein / empirical-Bayes shrinkage of the per-item power-law exponent p toward the pooled median
(bvm_shrink=True) — borrows strength across the exchangeable items so noisy items no longer over-pin the
conditional width at large n.

Run stages via env flags (defaults in brackets):
    RUN_SVI=1    (0)  fit the biweekly between-arm SVI grid first (delegates to UkraineRCT_interim_svi.py)
    RUN_DEPLOY=1 (1)  deploy the default amortiser (Stein BvM) against the SVI fit
    RUN_DIAG=1   (1)  build the combined diagnostic figures
    RUN_DEPLOY=0 RUN_DIAG=1 python scripts-py/UkraineRCT_startme.py     # diagnostics only
"""
import os, sys, subprocess
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repo root (parent of scripts-py)
sys.path.insert(0, os.path.join(_REPO, 'python'))
from amortiser_common import federated_deploy
from amortiser_diag_plots import FederatedDiagnostics

SB = os.environ.get('UKRAINERCT_SB', "/Users/or105/sandbox/bIRTistic")   # sandbox root for all output dirs
SVI_DIR = 'py-ukraineRCT-endline-svi-260929'   # the between-arm biweekly SVI dir (UkraineRCT_interim_svi.py)

# =====================================================================================
# 1) DATA & ANALYSIS — the RCT between-arm endline estimand (signed relative % improvement).
# =====================================================================================
DATA = dict(
    csv=os.environ.get('UKRAINERCT_CSV',
        "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/"
        "2025_project_Hope_Groups/data/Ukraine_Hope_Groups_Baseline_Endline_Wide_Aug6.csv"),
    categorical_threshold=2,      # caseness cut for the categorical CG-MH items: endpoint P(y>=2)
)
RHO = dict(rho_label='rel_change', h1=0.5,
           rho_label_long='intervention-vs-control relative change at endline (direction-aware per item)')

# =====================================================================================
# 2) SVI — the biweekly between-arm reference fit (delegated to UkraineRCT_interim_svi.py).
# =====================================================================================
SVI = dict(out=SVI_DIR)          # knobs live in UkraineRCT_interim_svi.py (AutoLowRankMVN, 10k, 4000)

# =====================================================================================
# 3) AMORTISER — the DEFAULT federated amortiser deploy + combined-diagnostics MANIFEST.
#    Signed relative change = registry S3, scale-feat net, warp=none, eta0 sweep 0-100%, STEIN BvM.
# =====================================================================================
AMORTISER = dict(
    svi=SVI_DIR,                  # reference SVI dir the amortiser deploys against (== SVI['out'])
    fed=os.environ.get('UKRAINERCT_FED', 'py-ukraineRCT-endline-amortise-widetok-grouprel-case-261020'),  # adopted widetok caseness deploy (UKRAINERCT_FED overrides)
    title='Ukraine Hope Groups RCT',   # figure-title prefix
    item_kind='item',             # response-item word (rows = the Ukraine parenting/MH items)
    mixed_units=False,            # single rho (relative change) -> one facet_grid (rows span both item types)
    calib_prefix='ukraineRCT_amortiser',    # basename of the calibration summary files
    eta0_units='relative change x100',       # units label for the eta0 success-threshold axes
    ctag_prefix='ukraineRCT-fed',            # deploy-cache tag prefix
    bvm_shrink=True,              # STEIN BvM: James-Stein/EB shrink of the per-item power-law exponent p
    endpoints=[dict(              # MANIFEST: one delegated amortiser per endpoint (here the single rho)
        rho='rel_change',                    # federated subdir + rho_label
        instance='S3 rel-change',            # registry net-family label (the DEFAULT Ukraine instance)
        build='svi',                         # reference build strategy: reuse the whole SVI fit as-is
        # The estimand is a CASENESS (thr=2) group rel% — a threshold functional — so the token-correct net is
        # the WIDE (cumulative-exceedance) group-rel% caseness net (not the scalar mean token, §14.4.18 wall).
        # Measured: widetok-grouprel-case PIT 0.127/marg 0.189 vs scalar scale-feat-paired 0.134/0.201. The
        # residual gap to the other between-arm apps (~0.09) is the between-subjects design noise (half effective
        # n per arm), NOT the token. scale-feat-paired (scalar) stays registered as the fallback.
        net=os.environ.get('UKRAINERCT_NET', 'widetok-grouprel-case'), widetok=1, case_c=2,  # token-correct caseness net
        warp='none',                         # signed relative change -> NO warp (log2/logit undefined)
        eta0=0.2, grid='0,0.1,0.2,0.3,0.4')])      # eta0 success sweep 0-40% (anchor 20% rel. improvement)


def _fit_svi():
    """Delegate to the between-arm biweekly SVI producer (single source of truth for the RCT reference)."""
    subprocess.run([sys.executable, os.path.join(_REPO, 'scripts-py', 'UkraineRCT_interim_svi.py')],
                   env=dict(os.environ, UKRAINERCT_SB=SB, UKRAINERCT_CSV=DATA['csv'],
                            UKRAINERCT_H1='0.2'),   # SVI success threshold aligned to the deploy eta0 anchor
                   cwd=_REPO, check=True)


if __name__ == "__main__":
    if os.environ.get('RUN_SVI', '0') == '1':
        _fit_svi()
    if os.environ.get('RUN_DEPLOY', '1') == '1':
        federated_deploy(SB, AMORTISER, run_diagnostics=False)
    if os.environ.get('RUN_DIAG', '1') == '1':
        FederatedDiagnostics(SB, AMORTISER).run()
