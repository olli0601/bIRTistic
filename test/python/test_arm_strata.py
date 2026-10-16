"""Arm-stratified effect sizes (SDY269 head-to-head): define_effect_sizes + reduction_stan_data +
the rho_by_cat_qu reader's cross-arm `across` step. Verifies that a strain shared between two arms
yields ONE comparison per arm (no item_group_id collision), and that the reader emits per-arm
SPR/GMFR plus an EXACT per-draw cross-arm difference on the shared strain only."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import arviz as az
import pytest

_repo_root = Path(__file__).resolve().parents[2]
_python_dir = _repo_root / 'python'
if str(_python_dir) not in sys.path:
    sys.path.insert(0, str(_python_dir))

from utils import define_effect_sizes, reduction_stan_data
from model_irt import IRTModel


@pytest.fixture
def headhead():
    """Two arms (LAIV/TIV), strains S1 (shared), S2 (LAIV only), S3 (TIV only), phases Baseline/Endline."""
    igid = {('LAIV', 'S1', 'Baseline'): 1, ('LAIV', 'S1', 'Endline'): 2,
            ('LAIV', 'S2', 'Baseline'): 3, ('LAIV', 'S2', 'Endline'): 4,
            ('TIV', 'S1', 'Baseline'): 5, ('TIV', 'S1', 'Endline'): 6,
            ('TIV', 'S3', 'Baseline'): 7, ('TIV', 'S3', 'Endline'): 8}
    dcati = pd.DataFrame([
        dict(item_type_id=1, item_label=strain, item_group_id=gid, phase_label=phase, arm=arm,
             group_label=f'{arm}_{phase.lower()}', item_type='categorical')
        for (arm, strain, phase), gid in igid.items()])
    dit = pd.DataFrame({'item_type_id': [1, 1, 1], 'item_label': ['S1', 'S2', 'S3'],
                        'item_type': ['categorical'] * 3, 'construct': ['HAI'] * 3,
                        'item_high_label': ['higher_is_better'] * 3})
    reductions = [dict(red_fn=6, rho='spr', thr=3), dict(red_fn=1, rho='gmfr')]
    drho = define_effect_sizes(dcati, dit, reductions=reductions, contrast_col='phase_label',
                               contrast_values=('Baseline', 'Endline'), strata_col='arm')
    return dcati, dit, drho


def test_comparisons_are_arm_stratified(headhead):
    _, _, drho = headhead
    assert drho.attrs['strata_col'] == 'arm'
    item = drho[drho.scope == 'item'][['item_label', 'stratum']].drop_duplicates()
    assert set(map(tuple, item.values)) == {('S1', 'LAIV'), ('S2', 'LAIV'), ('S1', 'TIV'), ('S3', 'TIV')}


def test_no_item_group_id_collision_across_arms(headhead):
    dcati, _, drho = headhead
    sd = reduction_stan_data(drho, dcati)
    ucmp = drho[['jc', 'scope', 'item_label', 'stratum']].drop_duplicates('jc')
    jc_l = ucmp[(ucmp.item_label == 'S1') & (ucmp.stratum == 'LAIV') & (ucmp.scope == 'item')].jc.iloc[0]
    jc_t = ucmp[(ucmp.item_label == 'S1') & (ucmp.stratum == 'TIV') & (ucmp.scope == 'item')].jc.iloc[0]
    # the shared strain resolves to its OWN baseline item_group_id in each arm
    assert {sd['question_cmp'][jc_l - 1][0][0], sd['question_cmp'][jc_t - 1][0][0]} == {1, 5}


def test_reader_per_arm_and_cross_arm_diff(headhead):
    dcati, _, drho = headhead
    sd = reduction_stan_data(drho, dcati)
    ucmp = drho[['jc', 'scope', 'item_label', 'stratum', 'r']].drop_duplicates('jc').sort_values('jc')
    jcmap = {(r.item_label, r.stratum): r.jc for r in ucmp[ucmp.scope == 'item'].itertuples()}
    spr = {('S1', 'LAIV'): 0.8, ('S2', 'LAIV'): 0.7, ('S1', 'TIV'): 0.6, ('S3', 'TIV'): 0.5}
    gmfr = {('S1', 'LAIV'): 4.0, ('S2', 'LAIV'): 3.0, ('S1', 'TIV'): 2.0, ('S3', 'TIV'): 1.5}
    R, Qcmp = sd['R'], sd['Qcmp']
    r_spr = int(drho[drho.rho == 'spr'].r.iloc[0]); r_gmfr = int(drho[drho.rho == 'gmfr'].r.iloc[0])
    arr = np.zeros((1, 2, R, Qcmp))
    for (lab, st), jc in jcmap.items():
        arr[0, :, r_spr - 1, jc - 1] = spr[(lab, st)]
        arr[0, :, r_gmfr - 1, jc - 1] = gmfr[(lab, st)]
    idata = az.from_dict(posterior={'rho_by_cat_qu_pr': arr})

    class _Stub:
        pass
    stub = _Stub(); stub.drho = drho; stub._resolve_draws = IRTModel._resolve_draws
    reader = IRTModel._rho_by_cat_qu_endpoints_per_draw.__get__(stub, _Stub)
    across = [dict(rho_label='spr', out_label='TIV-LAIV_spr_diff', ref='LAIV', foc='TIV', compare='diff'),
              dict(rho_label='gmfr', out_label='TIV-LAIV_gmfr_diff', ref='LAIV', foc='TIV', compare='diff')]
    out = reader(idata, None, 'rho_by_cat_qu_pr', 'items', verbose=False, across=across)
    d0 = out[out.draw == 0]

    def val(arm, lab, col):
        return float(d0[(d0.arm == arm) & (d0.item_label == lab)][col].iloc[0])

    assert abs(val('LAIV', 'S1', 'spr') - 0.8) < 1e-9
    assert abs(val('TIV', 'S1', 'gmfr') - 2.0) < 1e-9
    assert abs(val('LAIV', 'S2', 'gmfr') - 3.0) < 1e-9
    diff_rows = d0[d0.arm == 'TIV-LAIV']
    assert set(diff_rows.item_label) == {'S1'}                       # cross-arm diff only on the SHARED strain
    assert len(diff_rows) == 1                                       # both diff columns on ONE row
    assert abs(val('TIV-LAIV', 'S1', 'TIV-LAIV_spr_diff') - (0.6 - 0.8)) < 1e-9
    assert abs(val('TIV-LAIV', 'S1', 'TIV-LAIV_gmfr_diff') - (2.0 - 4.0)) < 1e-9
