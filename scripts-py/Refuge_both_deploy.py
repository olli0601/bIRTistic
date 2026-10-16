"""§19.3 apple-to-apple: deploy the BOTH (pre-pool SAB + post-pool attention) any-J BvM amortiser on
REFUGE-ED through the same federated manifest as the post-pool default, everything identical except the
net family. Produces the calibration csv comparable to §19.2 / §19.3.

Run: .pixi/envs/default/bin/python scripts-py/Refuge_both_deploy.py
"""
import os, sys
sys.path.insert(0, 'python')
from amortiser_common import federated_deploy
from amortiser_diag_plots import FederatedDiagnostics

SB = "/Users/or105/sandbox/bIRTistic"
SVI_DIR = 'py-refugee_interim_260831'
AMORTISER = dict(
    svi=SVI_DIR,
    fed='py-refugee-interim-amortise-deepsetScompAtt-itemamortise-J64-scale-feat-both-bvm-federated_261006',
    title='REFUGE-ED MSPSS (both)', item_kind='MSPSS item', mixed_units=False,
    calib_prefix='refuge_mspss_both_amortiser', eta0_units='relative gain x100',
    ctag_prefix='refuge-ed-both',
    endpoints=[dict(rho='support_gain', instance='S3 rel-change (both)', build='svi',
                    net='scale-feat-paired-both', widetok=0, case_c=2, warp='none',
                    eta0=0.10, grid='0,0.05,0.10,0.15,0.20')])

if __name__ == "__main__":
    federated_deploy(SB, AMORTISER, run_diagnostics=False)
    FederatedDiagnostics(SB, AMORTISER).run()
