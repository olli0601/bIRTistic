"""START HERE — ImmPort/ImmuneSpace COVID-19 SDY1764 subpopulation-contrast amortiser pipeline (§3.20 / §22).

The ONE file a human edits for this project. Config split into three separate dictionaries — DATA
(& analysis/estimand), SVI (the reference fits), AMORTISER (a LIST of the two contrast deploys) — same
shape as IMMPORT_flu-SDY312_startme.py (which also carries a list). Shared machinery imported from
fit_to_current_data + amortiser_common + amortiser_diag_plots.

SDY1764 (Distinct antibody responses to SARS-CoV-2 in children vs adults). BETWEEN-GROUP design
(HVTN 505 / mycelium style): the subpopulation is the group axis, encoded arm=time (group A -> time 0
'Baseline', group B -> time 1 'Endline'); each subject in one group. Item = SARS-CoV-2 serum-
neutralisation ID50 titre -> ordered category on the 2-fold ladder. rho = the direction-aware SIGNED
relative shift in mean log2 titre (group-B vs group-A, higher_is_better) -> registry S3, warp NONE.
TWO contrasts share the recipe: 'age' (pediatric vs adult) and 'severity' (severe vs mild).

Run stages via env flags (defaults in brackets):
    RUN_SVI=1    (0)  fit the two SVI interim grids first (produces py-immport-covid-SDY1764-{age,severity}_260920)
    RUN_DEPLOY=1 (1)  deploy each contrast's amortiser against its SVI fit
    RUN_DIAG=1   (1)  build the combined diagnostic figures per contrast
    RUN_DEPLOY=0 RUN_DIAG=1 python scripts-py/IMMPORT-covid-SDY1764_startme.py   # diagnostics only
"""
import os, sys
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repo root (parent of scripts-py)
sys.path.insert(0, os.path.join(_REPO, 'python'))
from amortiser_common import federated_deploy
from amortiser_diag_plots import FederatedDiagnostics

SB = os.environ.get('COVID_SB', "/Users/or105/sandbox/bIRTistic")     # sandbox root for all output dirs

# =====================================================================================
# 1) DATA & ANALYSIS — the ImmuneSpace COVID export, the study, and the two subpopulation contrasts.
# =====================================================================================
DATA = dict(
    xlsx=os.environ.get('COVID_XLSX',                              # the ImmuneSpace COVID-19 neutralisation export
        "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/"
        "2025_project_Hope_Groups/data/ImmuneSpace_COVID19_v260920.xlsx"),
    study='SDY1764',              # the one COVID SDY with an ordered-titre assay usable for the PCM
    contrasts=['age', 'severity'],   # subpopulation axes: age = pediatric vs adult; severity = severe vs mild
    categorical_threshold=2,      # caseness cut passed to the endpoint extraction (ignored for the out-of-7 titre)
)
# ESTIMAND rho: the direction-aware SIGNED relative shift in mean log2 titre, group-B vs group-A. h1 = the
# success threshold (>0.5 = a >50% relative titre gain); the eta0 sweep below spans the signed effect range.
RHO = dict(rho_label='neut_contrast', h1=0.5)

# =====================================================================================
# 2) SVI — the interim schedule + fit settings producing each contrast's reference posterior. Reuses the
#    existing 260920 grids (RUN_SVI defaults off); regenerating writes the same dir names.
# =====================================================================================
SVI = dict(
    out='py-immport-covid-SDY1764-{contrast}_260920',   # per-contrast SVI dir (the amortiser's reference)
    x_formula='~ group - 1',      # PCM design: one difficulty offset per subpopulation (group 0=A, 1=B)
    accrual='shuffled',           # shuffled accrual so both subpopulations are present at every interim
    n_interims=8, floor=20,       # interim grid = linspace(max(floor, n_full//n_interims), n_full, n_interims)
    seed=123,                     # RNG for the shuffled accrual + the SVI fit
    nsteps=4000, samples=2000,    # SVI optimisation steps + posterior draws saved per interim
    prob_width=1.3,               # width multiplier for the prob_by_question_fit figure
)

# =====================================================================================
# 3) AMORTISER — a LIST of the two contrast deploys (each its own SVI ref + federated dir). The signed
#    relative change is registry S3, deployed with the scalar scale-feat net, warp=none; the eta0 sweep is
#    SIGNED (0 = "group-B > group-A"). Each entry -> federated_deploy(SB, cfg) + FederatedDiagnostics(SB, cfg).
# =====================================================================================
def _contrast(name, long):
    """One contrast's amortiser config; `name` -> SVI ref + federated parent dir, `long` -> plot strip label."""
    return dict(
        svi=SVI['out'].format(contrast=name),   # reuse this contrast's existing SVI grid as the reference
        fed=f'py-immport-covid-SDY1764-{name}-amortise-deepsetXcompAtt-itemamortise-J64-ftheadexpand-bvm-federated_260924',
        title=f'COVID SDY1764 {name}', item_kind='neut readout', mixed_units=False,
        calib_prefix=f'covid_{name}_amortiser', eta0_units='relative change x100', ctag_prefix=f'covid-{name}-fed',
        endpoints=[dict(
            rho='neut_contrast',                # federated subdir + rho_label
            long=long,                          # explicit facet-strip label (the covid pkl carries no rho_label_long)
            instance='S3 rel-change',           # registry net-family label
            build='svi',                        # reference build strategy: reuse the whole SVI fit as-is
            net='scale-feat', widetok=0, case_c=2,     # the default scalar-mean net + token/caseness settings
            warp='none',                        # signed relative change -> no warp
            eta0=0.0, grid='-0.4,-0.2,0,0.2,0.4')])   # signed eta0 sweep (0 = group-B > group-A)


AMORTISER = [
    _contrast('age', 'SARS-CoV-2 neutralisation — pediatric vs adult (relative log2-titre shift)'),
    _contrast('severity', 'SARS-CoV-2 neutralisation — severe vs mild (relative log2-titre shift)')]


def _fit_svi(contrast):
    """Run the SVI interim grid for one contrast (DATA + SVI + RHO above)."""
    from data_loading import read_data_immport_covid_neut
    from fit_to_current_data import (fit_interim_grid, assign_item_group_id, shuffled_accrual, linspace_grid)
    d = read_data_immport_covid_neut(DATA['xlsx'], DATA['study'], group=contrast)
    dp1, dit, K = d['dp'], d['dit'], d['K']; gA, gB = d['groups']
    print(f"##### {DATA['study']} {contrast}: {dp1.pid.nunique()} subjects, K={K}, contrast {gB} vs {gA}; "
          f"group sizes {dp1.groupby('group_label').pid.nunique().to_dict()}")
    dp1 = assign_item_group_id(dp1)                       # item_group_id = (item x subpopulation) difficulty index
    pids = shuffled_accrual(dp1, seed=SVI['seed'])        # accrual order over subjects
    grid = linspace_grid(len(pids), nint=SVI['n_interims'], floor=SVI['floor'])
    cat_thr, h1 = DATA['categorical_threshold'], RHO['h1']

    def on_fit(m, f, k, n, xi):   # legacy direction-aware endpoint path (matches the original covid producer)
        xr = m.get_endpoints_per_draw(draws=f['draws'], categorical_threshold=cat_thr,
                                      endpoint_type='items').rename(columns={'ratio': 'pps_ratio_x'})
        xr['pps_H1_x'] = (xr['pps_ratio_x'] > h1).astype(int)
        if 'item_high_label' not in xr.columns:
            xr = xr.merge(dit[['item_label', 'item_high_label']], on='item_label', how='left')
        return xr[['draw', 'item_label', 'item_type', 'item_high_label', 'pps_ratio_x', 'pps_H1_x']]

    fit_interim_grid(f"{SB}/{SVI['out'].format(contrast=contrast)}", dp1, dit, pids, grid, SVI['x_formula'],
                     on_fit, seed=SVI['seed'], nsteps=SVI['nsteps'], output_samples=SVI['samples'],
                     prob_width=SVI['prob_width'], label=f"{DATA['study']}-{contrast}")


if __name__ == "__main__":
    if os.environ.get('RUN_SVI', '0') == '1':
        for c in DATA['contrasts']:
            _fit_svi(c)
    for cfg in AMORTISER:
        if os.environ.get('RUN_DEPLOY', '1') == '1':
            federated_deploy(SB, cfg, run_diagnostics=False)
        if os.environ.get('RUN_DIAG', '1') == '1':
            FederatedDiagnostics(SB, cfg).run()
