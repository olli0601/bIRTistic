#!/usr/bin/env bash
# =============================================================================
# ImmPort influenza HAI amortiser deployment — rho_SPR (seroprotection rate, §3.18).
# Deploys the WIDE-TOKEN item-general J64 amortiser (cumulative-exceedance token =
# categorically sufficient for the threshold functional, §14.4.18) + head BvM, on the
# ImmPort SVI grids for the SPR endpoint. Companion to deploy_IMMPORT_amortiser.sh (GMFR).
#
#   Step 0  retrain the SPR net (once): the item-general J64 trainer with
#           PP_WIDETOK=1 (cumulative-exceedance token), PP_CASEENDLINE=1 (target = endline
#           caseness Pr(y_e>=c) = SPR), PP_ALLCASE=1 (caseness-only, uniform [0,1] scale),
#           PP_META4=1 (caseness threshold c/Kmax in metadata).
#   Step 1  SPR-only reference grids  -> py-immport-SDY{312,314}_260918-spr/
#           (IMMPORT_RHO_ID=1 -> rho_id==1 as pps_ratio_x, symlink dp1/draws/dit).
#   Step 2  amortiser deploy          -> py-immport-SDY{312,314}-amortise-...-spr-widetok-...-bvm_260919/
#
# Endpoint: rho_SPR = endline seroprotection rate P(titre>=1:40) = P(y>=3) (baseline
#   computed but unused). eta0 is on the [0,1] SPR scale; 0.70 = the CHMP >70% criterion.
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
SB=/Users/or105/sandbox/bIRTistic
DRV=scripts-py/Ukraine_interim_analysis_amortise_endpt_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead_contraction_bvm.py
BASE=$SB/py-ukraine-interim-amortise-deepsetXcompAtt-itemamortise-J64-widetok-spr-260919

# ---- Step 0: retrain the SPR-capable net (once; ~100 min). Comment out after. ----
# PP_TAG=itemamortise-J64-widetok-spr-260919 PP_WIDETOK=1 PP_CASEENDLINE=1 PP_ALLCASE=1 PP_META4=1 \
# PP_STEPS=6000 PP_B=24 \
#   pixi run python scripts-py/Ukraine_interim_analysis_amortise_endpt_deepsetXcompAtt_ragged_itemamortise_explore.py

# ---- Step 1+2: SPR ref grid (into the SHARED amortiser dir, prefix pcm_spr_interim_*) + deploy ----
#   Writes into the SAME per-study amortiser dir as GMFR (deploy_IMMPORT_amortiser.sh); the two
#   endpoints are distinguished by file prefix, not by directory. Base SVI dir stays SVI-only.
deploy () {   # $1=SDY  $2=N_REF  $3=ctag
  local A=$SB/py-immport-$1-amortise-deepsetXcompAtt-itemamortise-J64-Cwarp-ftheadexpand-bvm_260920
  IMMPORT_SDY=$1 IMMPORT_RHO_ID=1 IMMPORT_AMORTDIR=$A \
    pixi run python scripts-py/IMMPORT_flu_gmfr_refgrid.py
  RAGD_BASE=$BASE \
  RAGD_RGE=$A \
  RAGD_OUT=$A \
  RAGD_FILEPREFIX=pcm_spr_interim \
  RAGD_CTAG=immport-${3}-spr \
  RAGD_WARP=logit \
  RAGD_WIDETOK=1 \
  RAGD_KMAX=10 \
  RAGD_CASE_C=3 \
  RAGD_HEADFT=expand \
  RAGD_BVM=1 \
  RAGD_NREF=$2 \
  RAG_S=200 \
  RAGD_ETA0=0.70 \
  RAGD_ETAH=0.89 \
  RAGD_ETA0GRID=1 \
  RAGD_ETA0GRID_VALS=0.5,0.6,0.7,0.8 \
  RAGD_PDFS=1 \
    pixi run python "$DRV"
}

deploy SDY312 84 sdy312
deploy SDY314 92 sdy314

# Endpoint scale is [0,1] (a rate), so eta0 = 0.70 is the CHMP >70% seroprotection bar.
# The wide cumulative-exceedance token makes the pooled statistic the empirical category
# CDF (sufficient for the threshold), which the scalar-mean GMFR net could not carry.
