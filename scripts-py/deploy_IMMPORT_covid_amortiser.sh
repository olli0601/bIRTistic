#!/usr/bin/env bash
# =============================================================================
# ImmPort COVID-19 subpopulation-contrast amortiser deployment (§3.20), SDY1764.
# Deploys the item-general J64 amortiser (scalar-mean token) + head BvM on the BETWEEN-GROUP
# SVI grids for the neutralisation-titre contrast. The endpoint rho is a SIGNED relative-change
# (group2-vs-group1 mean log2 titre) -> no warp (log2/logit are for folds/rates). The COVID SVI
# pkls already carry `pps_ratio_x`, so the SVI grid dir IS the reference grid (RAGD_RGE=grid dir;
# no rho-id transform needed). eta0 sweeps the effect-size threshold; 0 = "group2 > group1".
#   grids: py-immport-covid-SDY1764-{age,severity}_260920/  (IMMPORT_covid_interim_svi.py)
#   deploy -> py-immport-covid-SDY1764-{age,severity}-amortise-...-ftheadexpand-bvm_260920/
# Calibration to SVI: age PIT-KS 0.135 / marg 0.115; severity 0.122 / 0.104.
# PPS: severity severe>mild PPS 0.94 at eta0=0 (grades to 0 by +40%); age group2(pediatric) is
# confidently LOWER so P(pediatric>adult)=0 (the negative effect is in the p_rho plot; flip
# gA/gB in the loader to score adult>pediatric).
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
SB=/Users/or105/sandbox/bIRTistic
DRV=scripts-py/Ukraine_interim_analysis_amortise_endpt_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead_contraction_bvm.py
BASE=$SB/py-ukraine-interim-amortise-deepsetXcompAtt-itemamortise-J64-scale-feat-260831

# Step 1 (once): SVI grids -> COVID_GROUP=age,severity python scripts-py/IMMPORT_covid_interim_svi.py
# Step 2: amortiser deploy per contrast (N_REF > max interim n = 79; eta0 on the rho scale)
deploy () {   # $1=contrast
  local RGE=$SB/py-immport-covid-SDY1764-$1_260920
  local OUT=$SB/py-immport-covid-SDY1764-$1-amortise-deepsetXcompAtt-itemamortise-J64-ftheadexpand-bvm_260920
  RAGD_BASE=$BASE \
  RAGD_RGE=$RGE \
  RAGD_OUT=$OUT \
  RAGD_FILEPREFIX=pcm_1_interim \
  RAGD_CTAG=covid-$1 \
  RAGD_HEADFT=expand \
  RAGD_BVM=1 \
  RAGD_NREF=100 \
  RAG_S=200 \
  RAGD_ETA0=0 \
  RAGD_ETAH=0.89 \
  RAGD_ETA0GRID=1 \
  RAGD_ETA0GRID_VALS=-0.4,-0.2,0,0.2,0.4 \
  RAGD_PDFS=1 \
    pixi run python "$DRV"
}

deploy age
deploy severity
