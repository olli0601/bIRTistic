"""START HERE — SDY312 influenza HAI amortiser pipeline (§3.18 / §16.3).

The ONE file a human edits for this project. Config is split into three separate dictionaries —
DATA (& analysis/estimand), SVI (the reference fit), AMORTISER (the federated deploy + diagnostics;
here a LIST of the three §16.3 head/target variants that share one SVI fit). Shared machinery imported
from fit_to_current_data + amortiser_common + amortiser_diag_plots.

Run stages via env flags (defaults in brackets):
    RUN_SVI=1    (0)  fit the SVI grid first (produces py-immport-SDY312_260918)
    RUN_DEPLOY=1 (1)  deploy each variant amortiser against the SVI fit
    RUN_DIAG=1   (1)  build the combined diagnostic figures per variant
    RUN_SVI=1 python scripts-py/IMMPORT_flu-SDY312_startme.py            # full pipeline (all 3 variants)
    RUN_DEPLOY=0 RUN_DIAG=1 python scripts-py/IMMPORT_flu-SDY312_startme.py   # diagnostics only
"""
import os, sys
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repo root (parent of scripts-py)
sys.path.insert(0, os.path.join(_REPO, 'python'))
from amortiser_common import federated_deploy
from amortiser_diag_plots import FederatedDiagnostics

SB = os.environ.get('IMMPORT_SB', "/Users/or105/sandbox/bIRTistic")   # sandbox root for all output dirs

# =====================================================================================
# 1) DATA & ANALYSIS — the ImmuneSpace export + the two CHMP endpoints (the estimands).
# =====================================================================================
DATA = dict(
    xlsx=os.environ.get('IMMPORT_XLSX',                              # the ImmuneSpace influenza export workbook
        "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/"
        "2025_project_Hope_Groups/data/ImmuneSpace_Influenza_Vaccine_Trials_v260918.xlsx"),
    studies=['SDY312'],           # study id(s) to fit; add 'SDY314' to also produce that (non-federated) fit
)
# The two endpoints, declared upfront (rho_id -> federated subdir SPR/GMFR). h1_threshold = the CHMP bar:
# SPR (endline seroprotection rate P(titre>=1:40)) > 0.70; GMFR (GMT fold-rise endline/baseline) > 2.5-fold.
RHO_SPECS = [dict(rho_id=1, rho_label='SPR: P(titre>=1:40) at endline', red_fn=6, thr=3,
                  h1_threshold=0.70, rho_label_long='seroprotection rate  P(titre ≥ 1:40) at endline'),
             dict(rho_id=2, rho_label='GMT fold-rise (endline/baseline)', red_fn=1,
                  h1_threshold=2.5, rho_label_long='GMT fold-rise (endline / baseline)'),
             # SCR = third CHMP criterion, per-participant seroconversion (red_fn 2: mean_u P(>=4-fold rise),
             # via the paired joint G=p1'p2). No thr (d>=2 baked in). CHMP bar > 0.40 (adults 18-60).
             dict(rho_id=3, rho_label='SCR: P(>=4-fold rise)', red_fn=2,
                  h1_threshold=0.40, rho_label_long='seroconversion rate  P(>=4-fold rise or neg->>=1:40)')]

# =====================================================================================
# 2) SVI — the interim schedule + fit settings producing the reference posterior per study.
# =====================================================================================
SVI = dict(
    out='py-immport-{study}_260918',   # output dir pattern per study (the amortiser's reference dir)
    x_formula='~ group - 1',           # PCM design: one difficulty offset per timepoint (group 0=baseline, 1=endline)
    accrual='ordered',                 # accrual pseudo-ordered by participant-id string (paired within-subject design)
    step=10,                           # interim grid = every `step` participants + the full cohort
    seed=123,                          # SVI fit seed
    nsteps=4000, samples=2000,         # SVI optimisation steps + posterior draws saved per interim
)

# =====================================================================================
# 3) AMORTISER — the three §16.3 head/target variants (same SVI fit + endpoints; deploy-only difference:
#    base = symmetric head/no warp; cwarp = target warp logit(SPR)/log2(GMFR) [the default]; bfreeq =
#    free-quantile skew head). Each entry is passed to federated_deploy + FederatedDiagnostics.
# =====================================================================================
_spr = dict(net='widetok-spr', widetok=1, case_c=3, eta0=0.70, grid='0.5,0.6,0.7,0.8')   # SPR net + eta0 sweep
_gmfr = dict(net='scale-feat-paired', widetok=0, case_c=2, eta0=2.5, grid='2.0,2.5,3.0,3.5')     # GMFR net + eta0 sweep
# SCR net (S6, per-participant seroconversion): widetok + pre-pool SAB + mean-pool, deploy warp=none (SCR is
# interior ~0.19-0.51, not boundary-skewed like SPR; logit is a cheap deploy fallback). SCR_NET switches the
# post-pool variant: 'widetok-scr' = post-OFF (expected best, per-pp pattern §19.4); 'widetok-scr-postSAB' = post-on.
_scr = dict(net=os.environ.get('SCR_NET', 'widetok-scr'), widetok=1, case_c=3, eta0=0.40, grid='0.3,0.4,0.5,0.6')


def _variant(variant, tag, warp_spr, warp_gmfr, headmode='', with_scr=False):
    """Build one variant's amortiser config; `tag` -> federated parent dir, warps/headmode = the §16.3 remedy.
    `with_scr` appends the SCR (S6) endpoint — only the adopted (cwarp) variant deploys it, since SCR has its
    own net + warp and does not share the SPR/GMFR head experiment."""
    eps = [dict(rho='SPR', instance='S1 rate', rho_id=1, build='rho_id', warp=warp_spr, **_spr),
           dict(rho='GMFR', instance='S2 fold', rho_id=2, build='rho_id', warp=warp_gmfr, **_gmfr)]
    if with_scr:
        eps.append(dict(rho='SCR', instance='S6 seroconv', rho_id=3, build='rho_id', warp='none', **_scr))
    return dict(
        svi='py-immport-SDY312_260918',    # shared reference SVI dir (all 3 variants deploy against it)
        fed=f'py-immport-SDY312-amortise-deepsetScompAtt-itemamortise-J64-{tag}-federated_260924',
        title='SDY312 flu HAI', item_kind='strain', mixed_units=False, calib_prefix='flu_amortiser',
        eta0_units='endpoint units x100', headmode=headmode, ctag_prefix=f'flu-SDY312-{variant}',
        endpoints=eps)


AMORTISER = [_variant('base', 'ftheadexpand-bvm_260919', 'none', 'none'),
             _variant('cwarp', 'Cwarp-ftheadexpand-bvm_260920', 'logit', 'log2', with_scr=True),  # §16.3 adopted default; SCR deployed here
             _variant('bfreeq', 'Bfreeq-ftheadexpand-bvm_260920', 'none', 'none', headmode='freeq')]


def _fit_svi(study):
    """Run the SVI interim grid for one study (DATA + SVI + RHO_SPECS above)."""
    import pandas as pd
    from data_loading import read_data_immport_flu
    from model_pcm import PartialCreditModel
    from fit_to_current_data import (fit_interim_grid, assign_item_group_id, ordered_accrual,
                                     step_grid, endpoint_frame, FILE_PREFIX)
    d = read_data_immport_flu(DATA['xlsx'], study)
    dp1 = assign_item_group_id(d['dp'])                  # item_group_id = (strain x timepoint) difficulty index
    pids = ordered_accrual(dp1)                          # accrual order over participants
    grid = step_grid(len(pids), step=SVI['step'])        # interim sizes
    print(f"##### {study}: {dp1.pid.nunique()} paired participants, {dp1.item_label.nunique()} strains, K={d['K']}")
    h1 = {s['rho_id']: s['h1_threshold'] for s in RHO_SPECS}; any_rows = []
    reductions = [dict(red_fn=s['red_fn'], rho=s['rho_label'], thr=s.get('thr')) for s in RHO_SPECS]

    def on_fit(m, f, k, n, xi):                          # SPR+GMFR+SCR frame + BOTH composites (>=1 met / all met)
        xr = endpoint_frame(m, f, RHO_SPECS, h1=h1)
        _, pps_any = PartialCreditModel.joint_any_met(xr, group_keys=['item_label'])   # CHMP: SCR v SPR v GMFR
        _, pps_all = PartialCreditModel.joint_all_met(xr, group_keys=['item_label'])   # strict: SCR ^ SPR ^ GMFR
        pps = pps_any.merge(pps_all, on='item_label'); pps['interim'] = k; pps['n'] = n
        any_rows.append(pps)
        return xr
    out = f"{SB}/{SVI['out'].format(study=study)}"
    fit_interim_grid(out, dp1, d['dit'], pids, grid, SVI['x_formula'], on_fit,
                     seed=SVI['seed'], nsteps=SVI['nsteps'], output_samples=SVI['samples'], label=study,
                     reductions=reductions)
    pd.concat(any_rows, ignore_index=True).to_csv(f"{out}/{FILE_PREFIX}_any_met.csv", index=False)


if __name__ == "__main__":
    if os.environ.get('RUN_SVI', '0') == '1':
        for study in DATA['studies']:
            _fit_svi(study)
    for cfg in AMORTISER:
        if os.environ.get('RUN_DEPLOY', '1') == '1':
            federated_deploy(SB, cfg, run_diagnostics=False)
        if os.environ.get('RUN_DIAG', '1') == '1':
            FederatedDiagnostics(SB, cfg).run()
