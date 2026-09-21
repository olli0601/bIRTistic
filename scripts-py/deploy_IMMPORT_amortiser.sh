#!/usr/bin/env bash
# =============================================================================
# ImmPort influenza HAI amortiser deployment (§3.18). Records the exact RAGD_*
# environment used to deploy the item-general J64 amortiser (+ head BvM) on the
# ImmPort SVI grids via the SHARED ragged deploy driver (same driver + classes as
# Ukraine/REFUGE/ICRC; only the environment differs). Reproducible record.
#
#   Step 1  reference SVI grid   -> py-immport-SDY{312,314}_260918/   (IMMPORT_flu_interim_svi.py)
#   Step 2  GMFR ref grid + deploy -> the SINGLE per-study amortiser dir
#           py-immport-SDY{312,314}-amortise-deepsetXcompAtt-itemamortise-J64-Cwarp-ftheadexpand-bvm_260920/
#           under prefix pcm_gmfr_interim_* (SPR outputs share the SAME dir under pcm_spr_interim_*,
#           written by deploy_IMMPORT_amortiser_spr.sh). The base SVI dir stays SVI-only.
#
# Endpoint deployed: rho_GMFR = GMT fold-rise (endline/baseline) = 2^(mean log2 titre diff)
#   (rho_id 2). Encoder/net: item-amortised J64 scalar-mean token (deepsetXcompAtt);
#   calibration: expanding head-ft + affine + §14.4.26 head BvM between-first (RAGD_BVM=1).
# Success sweep: eta0 = GMFR fold thresholds {2, 2.5(CHMP), 3, 3.5}; etaH=0.89.
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
SB=/Users/or105/sandbox/bIRTistic
DRV=scripts-py/Ukraine_interim_analysis_amortise_endpt_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead_contraction_bvm.py
BASE=$SB/py-ukraine-interim-amortise-deepsetXcompAtt-itemamortise-J64-scale-feat-260831

# ---- Step 1: SVI reference grids (writes the endpoint pkls). Run once, then comment out. ----
# pixi run python scripts-py/IMMPORT_flu_interim_svi.py

# ---- Step 2: GMFR ref grid (into the amortiser dir) + amortiser deploy, per study ----
#   N_REF = enrolled total, > max interim n. SDY312: 79 paired (84 enrolled); SDY314: 89 (92).
deploy () {   # $1=SDY  $2=N_REF  $3=ctag
  local A=$SB/py-immport-$1-amortise-deepsetXcompAtt-itemamortise-J64-Cwarp-ftheadexpand-bvm_260920
  IMMPORT_SDY=$1 IMMPORT_RHO_ID=2 IMMPORT_AMORTDIR=$A \
    pixi run python scripts-py/IMMPORT_flu_gmfr_refgrid.py
  RAGD_BASE=$BASE \
  RAGD_RGE=$A \
  RAGD_OUT=$A \
  RAGD_FILEPREFIX=pcm_gmfr_interim \
  RAGD_CTAG=immport-${3} \
  RAGD_WARP=log2 \
  RAGD_HEADFT=expand \
  RAGD_BVM=1 \
  RAGD_NREF=$2 \
  RAG_S=200 \
  RAGD_ETA0=2.5 \
  RAGD_ETAH=0.89 \
  RAGD_ETA0GRID=1 \
  RAGD_ETA0GRID_VALS=2.0,2.5,3.0,3.5 \
  RAGD_PDFS=1 \
    pixi run python "$DRV"
}

deploy SDY312 84 sdy312
deploy SDY314 92 sdy314

# Result (calibrated head-ft-expand+bvm vs SVI GMFR reference):
#   SDY312  PIT-KS 0.089  marg-KS 0.114   (baseline 0.318 / 0.433)
#   SDY314  PIT-KS 0.082  marg-KS 0.090   (baseline 0.333 / 0.512)
# Notes:
#   RAGD_BVM=1 adds the §14.4.26 head BvM between-first correction (per-item power-law SVI
#     targets) on top of expanding head-ft + affine.
#   eta0 is on the GMFR (fold) scale; 250% == the CHMP >2.5-fold criterion. At that bar the
#     high-baseline Immune Signatures cohorts mostly fail (futility); a >=2-fold rise (200%)
#     is met only by A/Uruguay. The eta0 sweep captures the full go/no-go picture.
