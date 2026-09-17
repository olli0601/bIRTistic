#!/usr/bin/env bash
# =============================================================================
# ICRC DASS-21 / DRC amortiser deployment (§3.10). Records the exact RAGD_*
# environment used to deploy the item-general J64 amortiser on the ICRC SVI grid
# via the SHARED ragged deploy driver (same driver + classes as Ukraine/REFUGE;
# only the environment differs). Reproducible record.
#
#   Step 1  reference SVI grid  ->  py-icrc-dass-drc_260902/
#   Step 2  amortiser deploy    ->  py-icrc-dass-drc-amortise-...-ftheadexpand_260903/
#
# Endpoint: DASS-21 subscale totals binned to Figure-2 severity (K=5), paired
# pre/post; interims every 100 beneficiaries (+ early 20/40/60/80), up to n_full.
# Encoder/net: itemamortise J64; calibration: expanding head-ft + affine (no BvM).
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
SB=/Users/or105/sandbox/bIRTistic
DRV=scripts-py/Ukraine_interim_analysis_amortise_endpt_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead_contraction_bvm.py

# ---- Step 1: SVI reference grid (writes RAGD_RGE inputs). Run once, then comment out. ----
pixi run python scripts-py/ICRC_interim_svi.py

# ---- Step 2: amortiser deploy ----
RAGD_BASE=$SB/py-ukraine-interim-amortise-deepsetXcompAtt-itemamortise-J64-scale-feat-260831 \
RAGD_RGE=$SB/py-icrc-dass-drc_260902 \
RAGD_OUT=$SB/py-icrc-dass-drc-amortise-deepsetXcompAtt-itemamortise-J64-scale-feat-ftheadexpand_260903 \
RAGD_HEADFT=expand \
RAGD_NREF=1700 \
RAG_S=200 \
RAGD_ETA0=0.50 \
RAGD_ETAH=0.89 \
RAGD_ETA0GRID=1 \
RAGD_ETA0GRID_VALS=0,0.50,0.70,0.85,0.95 \
  pixi run python "$DRV"

# Notes:
#   RAGD_NREF=1700 is the trial total N; it must exceed the largest interim n
#     (n_full=1669) so the final interim keeps a non-empty future cohort (m>0).
#     If the analysis N changes, set RAGD_NREF to that total (> max interim n).
#   RAGD_ETA0=0.50 + the 50-95% sweep grid reflect the large DASS severity shifts.
#   RAGD_BVM left unset (=0): expanding head-ft + affine only (matches the
#     '-ftheadexpand' output dir, no '-bvm' suffix).
