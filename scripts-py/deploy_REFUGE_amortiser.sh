#!/usr/bin/env bash
# =============================================================================
# REFUGE-ED amortiser deployment (§3.9). Records the exact RAGD_* environment
# used to deploy the item-general J64 amortiser on the REFUGE SVI grid via the
# SHARED ragged deploy driver (same driver + classes as Ukraine/ICRC; only the
# environment differs). Reproducible record — re-running writes the same dir.
#
#   Step 1  reference SVI grid  ->  py-refugee_interim_260831/
#   Step 2  amortiser deploy    ->  py-refugee-interim-amortise-...-ftheadexpand_260831/
#
# Encoder/net: itemamortise J64 (scale + item features); head_mode from the net.
# Calibration: expanding head fine-tune + affine median-shift (RAGD_BVM unset = 0).
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
SB=/Users/or105/sandbox/bIRTistic
DRV=scripts-py/Ukraine_interim_analysis_amortise_endpt_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead_contraction_bvm.py

# ---- Step 1: SVI reference grid (writes RAGD_RGE inputs). Run once, then comment out. ----
pixi run python scripts-py/REFUGE_interim_svi.py

# ---- Step 2: amortiser deploy ----
RAGD_BASE=$SB/py-ukraine-interim-amortise-deepsetXcompAtt-itemamortise-J64-scale-feat-260831 \
RAGD_RGE=$SB/py-refugee_interim_260831 \
RAGD_OUT=$SB/py-refugee-interim-amortise-deepsetXcompAtt-itemamortise-J64-scale-feat-ftheadexpand_260831 \
RAGD_HEADFT=expand \
RAGD_NREF=400 \
RAG_S=200 \
RAGD_ETA0=0.10 \
RAGD_ETAH=0.89 \
RAGD_ETA0GRID=1 \
RAGD_ETA0GRID_VALS=0,0.05,0.10,0.15,0.20 \
  pixi run python "$DRV"

# Notes:
#   RAGD_NREF=400 is the trial total N; it must exceed the largest interim n
#     (324 here) so the final interim keeps a non-empty future cohort (m>0).
#   RAGD_ETA0=0.10 is the decision effect threshold (REFUGE effects are modest);
#     RAGD_ETA0GRID_VALS is the eta0 sweep grid (0-20% endpoint change).
#   RAGD_BVM is left unset (=0): expanding head-ft + affine only (matches the
#     '-ftheadexpand' output dir, no '-bvm' suffix).
