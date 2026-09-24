"""START HERE — SDY269 head-to-head (LAIV vs TIV) amortiser pipeline (§17).

The ONE file a human edits for this project. Config is split into three separate dictionaries —
DATA (& analysis), SVI (the reference joint fit), AMORTISER (the 5-endpoint federated deploy +
diagnostics). Shared machinery imported from fit_to_current_data + amortiser_common + amortiser_diag_plots.

ONE joint partial-credit fit per interim (item = strain, group = arm x phase); endpoints are read
one-by-one off that fit (the two cross-arm diffs are EXACT per-draw contrasts). The rho set itself is
DATA-DERIVED (which strains are shared between arms), so it is declared by the loader (d['rho_specs'])
rather than hardcoded here.

Run stages via env flags (defaults in brackets):
    RUN_SVI=1    (0)  fit the joint SVI grid first (produces py-immport-SDY269)
    RUN_DEPLOY=1 (1)  deploy the 5 delegated amortisers
    RUN_DIAG=1   (1)  build the combined diagnostic figures
    RUN_DEPLOY=0 RUN_DIAG=1 python scripts-py/IMMPORT_flu-SDY269_startme.py   # diagnostics only
"""
import os, sys
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repo root (parent of scripts-py)
sys.path.insert(0, os.path.join(_REPO, 'python'))
from amortiser_common import federated_deploy
from amortiser_diag_plots import FederatedDiagnostics

SB = os.environ.get('IMMPORT_SB', "/Users/or105/sandbox/bIRTistic")   # sandbox root for all output dirs
SVI_DIR = 'py-immport-SDY269'    # the joint-fit output dir (shared by SVI['out'] and every AMORTISER['svi'])

# =====================================================================================
# 1) DATA & ANALYSIS — the ImmuneSpace export, the study + arms. The rho endpoints (SPR/GMFR per arm +
#    the two cross-arm diffs, each with its own h1) are declared by the loader since they depend on which
#    strain is shared between the arms; see read_data_immport_flu_headhead(...)['rho_specs'].
# =====================================================================================
DATA = dict(
    xlsx=os.environ.get('IMMPORT_XLSX',                              # the ImmuneSpace influenza export workbook
        "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/"
        "2025_project_Hope_Groups/data/ImmuneSpace_Influenza_Vaccine_Trials_v260918.xlsx"),
    study='SDY269',               # the head-to-head study id
    arms=('LAIV', 'TIV'),         # the two vaccine arms compared (reference, focal)
)

# =====================================================================================
# 2) SVI — the interim schedule + fit settings producing the joint reference posterior.
# =====================================================================================
SVI = dict(
    out=SVI_DIR,                  # output dir (per interim: joint dp1 / draws.zarr / prob-plots / joint rho pkl)
    x_formula='~ phase - 1',      # PCM design on the within-arm PAIRED axis phase (0=baseline, 1=endline)
    accrual='interleaved',        # accrual interleaves the two arms so both accrue together
    step=10,                      # interim grid = every `step` pooled participants + the full cohort
    seed=123,                     # SVI fit seed
    nsteps=4000, samples=2000,    # SVI optimisation steps + posterior draws saved per interim
    prob_width=1.2,               # width multiplier for the prob_by_question_fit figure
)

# =====================================================================================
# 3) AMORTISER — the 5 delegated amortisers (SPR + GMFR per arm on the shared strain, and the standardised
#    between-arm SPR difference). Mixed-unit columns (rate / fold / difference) so the contraction &
#    rho-vs-eta0 figures stitch one own-y-scale column per endpoint. Passed to federated_deploy + diagnostics.
# =====================================================================================
_spr = dict(net='widetok-spr', widetok=1, case_c=3, eta0=0.70, grid='0.5,0.6,0.7,0.8')   # SPR net + eta0 sweep
_gmfr = dict(net='scale-feat', widetok=0, case_c=2, eta0=2.5, grid='2.0,2.5,3.0,3.5')     # GMFR net + eta0 sweep
AMORTISER = dict(
    svi=SVI_DIR,                  # reference SVI dir the amortisers deploy against (== SVI['out'])
    fed=f'py-immport-{DATA["study"]}-amortise-deepsetXcompAtt-itemamortise-J64-ftheadexpand-bvm-federated_260922',
    title='SDY269 head-to-head', item_kind='strain', mixed_units=True, calib_prefix='headhead_amortiser',
    eta0_units='endpoint units x100', ctag_prefix='sdy269-fed', nref=64,   # nref fixed (small paired arms)
    endpoints=[  # build strategies: 'level' = per-arm phase subset; 'diffstd' = between-arm standardised diff
        dict(rho='LAIV_spr', instance='S1 rate', build='level', arm='LAIV', joint='LAIV_spr', **_spr, warp='logit'),
        dict(rho='TIV_spr', instance='S1 rate', build='level', arm='TIV', joint='TIV_spr', **_spr, warp='logit'),
        dict(rho='LAIV_gmfr', instance='S2 fold', build='level', arm='LAIV', joint='LAIV_gmfr', **_gmfr, warp='log2'),
        dict(rho='TIV_gmfr', instance='S2 fold', build='level', arm='TIV', joint='TIV_gmfr', **_gmfr, warp='log2'),
        dict(rho='TIV-LAIV_spr_diff', instance='S5 diff', build='diffstd', joint=('LAIV_spr', 'TIV_spr'),
             net='groupdiff', widetok=1, case_c=3, warp='none', eta0=0.5, grid='0.0,0.2,0.5,0.8')])


def _fit_svi():
    """Run the joint SVI interim grid (DATA + SVI above; rho set from the loader)."""
    import pandas as pd
    from data_loading import read_data_immport_flu_headhead
    from model_pcm import PartialCreditModel
    from fit_to_current_data import (fit_interim_grid, assign_item_group_id, interleaved_accrual, step_grid)
    d = read_data_immport_flu_headhead(DATA['xlsx'], DATA['study'], arms=DATA['arms'])
    dp1 = assign_item_group_id(d['dp'])                  # item_group_id = strain x (arm x phase) difficulty index
    rho = d['rho_specs']                                 # the 6 loader-declared rhos (short/long label, kind, arm, h1)
    pids = interleaved_accrual(dp1, DATA['arms'])         # arm-interleaved accrual order
    grid = step_grid(len(pids), step=SVI['step'])
    print(f"##### {DATA['study']} JOINT fit: {dp1.pid.nunique()} pooled participants, "
          f"{dp1.item_label.nunique()} strains, shared {d['shared']}, {len(rho)} rho")
    out = f"{SB}/{SVI['out']}"; any_rows = []
    hh_cols = ['draw', 'arm', 'item_label', 'item_type', 'item_high_label', 'rho_id', 'rho_label',
               'rho_label_long', 'reduction', 'compare', 'group1', 'group2', 'pps_rho_x', 'pps_H1_x']

    def _one_rho(m, f, spec):                            # one endpoint off the joint fit, tagged with its labels
        xr = m.get_endpoints_per_draw(draws=f['draws'], endpoint_type='items', rho_specs=[spec],
                                      contrast_col='phase').rename(columns={'rho': 'pps_rho_x'})
        if spec.get('kind') == 'level':                  # a per-arm level -> restrict to its arm
            xr = xr[xr['arm'] == spec['arm']]
        xr['rho_label'] = spec['rho_label']; xr['rho_label_long'] = spec['rho_label_long']
        xr['pps_H1_x'] = (xr['pps_rho_x'] > spec['h1']).astype(int)
        return xr[hh_cols]

    def on_fit(m, f, k, n, xi):                          # concat all rhos (shared draw index) + composite '>=1 met'
        joint = pd.concat([_one_rho(m, f, s) for s in rho], ignore_index=True)
        lvl = joint[joint.rho_label.isin([s['rho_label'] for s in rho if s['is_level']])]
        _, pps = PartialCreditModel.joint_any_met(lvl, group_keys=['item_label', 'arm'])
        pps['interim'] = k; any_rows.append(pps)
        return joint
    fit_interim_grid(out, dp1, d['dit'], pids, grid, SVI['x_formula'], on_fit,
                     seed=SVI['seed'], nsteps=SVI['nsteps'], output_samples=SVI['samples'],
                     prob_width=SVI['prob_width'], label=DATA['study'])
    pd.concat(any_rows, ignore_index=True).to_csv(f"{out}/headhead_any_met.csv", index=False)


if __name__ == "__main__":
    if os.environ.get('RUN_SVI', '0') == '1':
        _fit_svi()
    if os.environ.get('RUN_DEPLOY', '1') == '1':
        federated_deploy(SB, AMORTISER, run_diagnostics=False)
    if os.environ.get('RUN_DIAG', '1') == '1':
        FederatedDiagnostics(SB, AMORTISER).run()
