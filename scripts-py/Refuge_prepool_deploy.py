"""§19 apple-to-apple: deploy the PRE-POOL any-J BvM amortiser on REFUGE-ED through the exact same
federated manifest as the post-pool default (scale-feat-paired), everything identical except the
net family (pre-pool SAB, post-pool attention off). Produces the calibration csv comparable to §19.2.

Run: .pixi/envs/default/bin/python scripts-py/Refuge_prepool_deploy.py
"""
import os, sys
sys.path.insert(0, 'python')
from amortiser_common import federated_deploy
from amortiser_diag_plots import FederatedDiagnostics

SB = "/Users/or105/sandbox/bIRTistic"
SVI_DIR = 'py-refugee_interim_260831'
AMORTISER = dict(
    svi=SVI_DIR,
    fed='py-refugee-interim-amortise-deepsetXcompAtt-itemamortise-J64-scale-feat-prepool-bvm-federated_261006',
    title='REFUGE-ED MSPSS (pre-pool)', item_kind='MSPSS item', mixed_units=False,
    calib_prefix='refuge_mspss_prepool_amortiser', eta0_units='relative gain x100',
    ctag_prefix='refuge-ed-prepool',
    endpoints=[dict(rho='support_gain', instance='S3 rel-change (pre-pool)', build='svi',
                    net='scale-feat-paired-prepool', widetok=0, case_c=2, warp='none',
                    eta0=0.10, grid='0,0.05,0.10,0.15,0.20')])

if __name__ == "__main__":
    federated_deploy(SB, AMORTISER, run_diagnostics=False)
    FederatedDiagnostics(SB, AMORTISER).run()
