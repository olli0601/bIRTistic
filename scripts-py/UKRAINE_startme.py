"""START HERE — Ukraine Hope Groups parenting trial amortiser pipeline (§14, §16).

The ONE file a human edits for this project. Config split into three separate dictionaries — DATA
(& analysis/estimand), SVI (the reference fit), AMORTISER (the federated deploy + diagnostics) — same
shape as CAVD-hvtn505_startme.py. Shared machinery imported from fit_to_current_data + amortiser_common
+ amortiser_diag_plots. Ukraine is the ORIGINAL application and the net-training source; here it is wired
onto the same startme/federated pattern as the newer apps, deploying the DEFAULT amortiser.

Paired baseline/endline parenting trial with a real accrual CALENDAR, so interims accrue by DATE at
WEEKLY cadence (not participant count). Two item types share one fit — the days-in-week practices
(`out-of-7`, K=8, endpoint = mean days) and the caseness items (`categorical`, K=4, endpoint =
P(y>=2)) — and the single endpoint rho is the direction-aware SIGNED relative change endline-vs-
baseline (some items lower-is-better, e.g. violence; the ratio is oriented so higher = better). Signed
rho => NO target warp (log2/logit undefined; §16.3). One rho -> a one-column diagnostic grid over items.

This writes only to NEW output dirs (SVI_DIR / AMORTISER['fed'] below dated 260924); it does NOT touch
any existing py-ukraine-* analysis directory.

Run stages via env flags (defaults in brackets):
    RUN_SVI=1    (0)  fit the weekly SVI interim grid first (~13 min; produces py-ukraine-interim-weekly-svi-260924)
    RUN_DEPLOY=1 (1)  deploy the default amortiser against the SVI fit
    RUN_DIAG=1   (1)  build the combined diagnostic figures
    RUN_SVI=1 python scripts-py/UKRAINE_startme.py                     # full pipeline from scratch
    RUN_DEPLOY=0 RUN_DIAG=1 python scripts-py/UKRAINE_startme.py       # diagnostics only
"""
import os, sys
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repo root (parent of scripts-py)
sys.path.insert(0, os.path.join(_REPO, 'python'))
from amortiser_common import federated_deploy
from amortiser_diag_plots import FederatedDiagnostics

SB = os.environ.get('UKRAINE_SB', "/Users/or105/sandbox/bIRTistic")   # sandbox root for all output dirs
SVI_DIR = 'py-ukraine-interim-weekly-svi-260924'   # NEW weekly-SVI dir (leaves the existing -260811 untouched)

# =====================================================================================
# 1) DATA & ANALYSIS — the Ukraine Hope Groups wide CSV + the signed relative-change endpoint (estimand).
# =====================================================================================
DATA = dict(
    csv=os.environ.get('UKRAINE_CSV',                              # the Ukraine Hope Groups baseline/endline wide CSV
        "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/"
        "2025_project_Hope_Groups/data/Ukraine_Hope_Groups_Baseline_Endline_Wide_Aug6.csv"),
    categorical_threshold=2,      # caseness cut for the categorical (CG-MH) items: their endpoint is P(y>=2)
)
# ESTIMAND rho: the direction-aware SIGNED relative change endline-vs-baseline, per item. For lower-is-
# better items (violence, MH, depression, child behaviour) the ratio is 1 - g2/g1; else g2/g1 - 1 — so
# higher rho always = better. h1 = success threshold (>0.5 = a >50% relative improvement).
RHO = dict(rho_label='rel_change', h1=0.5,
           rho_label_long='signed relative change, endline vs baseline (direction-aware per item)')

# =====================================================================================
# 2) SVI — the WEEKLY interim schedule + fit settings producing the reference posterior. The knobs match
#    the published Ukraine weekly reference (AutoLowRankMVN, 10k steps, 4000 draws), not the count-loop default.
# =====================================================================================
SVI = dict(
    out=SVI_DIR,                  # output dir (per weekly interim: dp1 / draws.zarr / prob-plots / endpoint pkl)
    x_formula='~ group - 1',      # PCM design: one difficulty offset per timepoint (group 0=baseline, 1=endline)
    cadence='weekly',             # interims accrue by DATE (submission_date), one cutoff per week-ending
    seed=123,                     # SVI fit seed
    nsteps=10000, samples=4000,   # SVI optimisation steps + posterior draws saved per interim (Ukraine reference)
    algorithm='AutoLowRankMultivariateNormal',   # the Ukraine reference guide family (low-rank MVN)
)

# =====================================================================================
# 3) AMORTISER — the DEFAULT federated amortiser deploy + combined-diagnostics config. The signed relative
#    change is registry S3, deployed with the scalar scale-feat net (the Ukraine-trained default), warp=none,
#    eta0 sweep over the full 0-100% range. Passed verbatim to federated_deploy + FederatedDiagnostics.
# =====================================================================================
AMORTISER = dict(
    svi=SVI_DIR,                  # reference SVI dir the amortiser deploys against (== SVI['out'])
    fed='py-ukraine-interim-amortise-deepsetXcompAtt-itemamortise-J64-ftheadexpand-bvm-federated_260924',  # parent dir
    title='Ukraine Hope Groups',  # figure-title prefix
    item_kind='item',             # response-item word (rows = the Ukraine parenting/MH items)
    mixed_units=False,            # single rho (relative change) -> one facet_grid (rows span both item types)
    calib_prefix='ukraine_amortiser',      # basename of the calibration summary files
    eta0_units='relative change x100',      # units label for the eta0 success-threshold axes
    ctag_prefix='ukraine-fed',              # deploy-cache tag prefix
    endpoints=[dict(
        rho='rel_change',                   # federated subdir + rho_label
        instance='S3 rel-change',           # registry net-family label (the DEFAULT Ukraine instance)
        build='svi',                        # reference build strategy: reuse the whole SVI fit as-is
        net='scale-feat', widetok=0, case_c=2,     # the default Ukraine-trained scalar-mean net
        warp='none',                        # signed relative change -> no warp
        eta0=0.5, grid='0,0.25,0.5,0.75,1.0')])   # eta0 sweep over the full 0-100% range


def _fit_svi():
    """Run the WEEKLY SVI interim grid for Ukraine (DATA + SVI + RHO above)."""
    from data_loading import read_data_ukraine
    from fit_to_current_data import fit_interim_grid_weekly, weekly_dates
    raw = read_data_ukraine(DATA['csv'])
    dp, dit, dmeta = raw['dp'].copy(), raw['dit'].copy(), raw['dmeta'].copy()
    dit = dit[~dit['item_label'].astype(str).str.contains('agg')].reset_index(drop=True)   # drop aggregate items (not deployed)
    # attach displacement_status (a dmeta covariate carried through), then the Ukraine preprocessing:
    tmp = (dmeta[['pid', 'group_label', 'displacement_status']].drop_duplicates()
           .rename(columns={'pid': 'pid_label'}).dropna(subset=['displacement_status'], how='all'))
    dp = dp.merge(tmp, on=['pid_label', 'group_label'], how='inner', validate='many_to_one')
    dp1 = dp[~dp['item_label'].str.contains('agg')].copy()          # drop the aggregate/composite items
    dp1['y_stan'] = dp1['y'] + 1                                    # model input is 1..K
    dp1 = dp1.merge(dit[['item_label', 'item_type']], on='item_label', how='left')
    # item_group_id = per-item_type (item x timepoint) difficulty index (categorical + out-of-7 kept separate)
    it = (dp1[['item_type', 'item_label', 'group']].drop_duplicates()
          .sort_values(['item_type', 'group', 'item_label']).reset_index(drop=True))
    it['item_group_id'] = it.groupby('item_type').cumcount() + 1
    dp1 = dp1.merge(it, on=['item_label', 'group', 'item_type'], how='left')
    dp1 = dp1.merge(dit[['item_type', 'item_type_id']].drop_duplicates(), on='item_type', how='left')
    dp1 = dp1.sort_values(['item_type_id', 'pid', 'group', 'item_label']).reset_index(drop=True)
    dp1['oid'] = range(1, len(dp1) + 1); dp1['oidt'] = dp1.groupby('item_type').cumcount() + 1
    dates = weekly_dates(dp1)                                       # weekly cutoffs over the endline date span
    cat_thr, h1 = DATA['categorical_threshold'], RHO['h1']
    print(f"##### Ukraine weekly: {dp1.pid.nunique()} paired participants, {dp1.item_label.nunique()} items, "
          f"{dp1.item_type.nunique()} item types; {len(dates)} weekly cutoffs")

    def on_fit(m, f, k, n, xi):   # legacy direction-aware endpoint path (categorical_threshold per item type)
        xr = m.get_endpoints_per_draw(draws=f['draws'], categorical_threshold=cat_thr,
                                      endpoint_type='items').rename(columns={'ratio': 'pps_ratio_x'})
        xr['pps_H1_x'] = (xr['pps_ratio_x'] > h1).astype(int)
        xr['rho_label'] = RHO['rho_label']; xr['rho_label_long'] = RHO['rho_label_long']
        return xr

    fit_interim_grid_weekly(f"{SB}/{SVI['out']}", dp1, dit, dates, SVI['x_formula'], on_fit,
                            seed=SVI['seed'], nsteps=SVI['nsteps'], output_samples=SVI['samples'],
                            algorithm=SVI['algorithm'], label='Ukraine')


if __name__ == "__main__":
    if os.environ.get('RUN_SVI', '0') == '1':
        _fit_svi()
    if os.environ.get('RUN_DEPLOY', '1') == '1':
        federated_deploy(SB, AMORTISER, run_diagnostics=False)
    if os.environ.get('RUN_DIAG', '1') == '1':
        FederatedDiagnostics(SB, AMORTISER).run()
