"""Generic FEDERATED amortiser deploy for the SDY269 head-to-head (§17.5). Each declared endpoint
(rho_spec) resolves to its own trained network, warp and token; this script builds the per-rho
reference by SUBSETTING the existing joint SVI fit (no re-fit) and deploys the matching net, so the
head-to-head's diagnostics come from several federated deploy directories — one per rho.

Reference construction per rho_spec (reuse of the joint fit, cf. §17.4 S5):
  - dp1: the joint fit's rows for the shared strain, filtered to the rho's arm; `group` = the paired
    phase (0/1) for a within-arm level (SPR/GMFR); item = the strain. oid/item_group_id re-indexed.
  - draws.zarr: the joint fit's posterior `ypred` subset to those observations (the future cohort).
  - endpoint pkl: the rho's per-draw value taken straight from the joint frame (pps_rho_x -> pps_ratio_x).
Then the ragged BvM driver deploys BASE=<the rho's net> with its token/warp. The cross-arm diff (S5)
is deployed separately (standardised target) by IMMPORT_headhead_S5 note; see project memory."""
# ---- boilerplate ----

import os, sys, glob, re, shutil
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
import numpy as np, pandas as pd, xarray as xr

SB = os.environ.get('IMMPORT_SB', "/Users/or105/sandbox/bIRTistic")
SRC = f"{SB}/py-immport-SDY269"; pfx = "pcm_1_interim"
URU = 'A/Uruguay/716/2007(H3N2)'
NET = {  # registry instance -> trained net dir + token/warp
    'widetok-spr':  ("py-ukraine-interim-amortise-deepsetXcompAtt-itemamortise-J64-widetok-spr-260919", 1, 'logit', 3),
    'scale-feat':   ("py-ukraine-interim-amortise-deepsetXcompAtt-itemamortise-J64-scale-feat-260831", 0, 'log2', 2),
}
# the head-to-head level manifest on the shared strain: (rho_label, arm, endpoint-net, endpoint rho_label in the joint frame)
MANIFEST = [
    ('LAIV_spr',  'LAIV', 'widetok-spr', 'LAIV_spr'),
    ('TIV_spr',   'TIV',  'widetok-spr', 'TIV_spr'),
    ('LAIV_gmfr', 'LAIV', 'scale-feat',  'LAIV_gmfr'),
    ('TIV_gmfr',  'TIV',  'scale-feat',  'TIV_gmfr'),
]
DRV = "scripts-py/Ukraine_interim_analysis_amortise_endpt_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead_contraction_bvm.py"


def build_ref(rho_label, arm, joint_rho):
    """Subset the joint fit -> a per-(arm, shared-strain) paired reference (group=phase)."""
    ref = f"{SB}/py-immport-SDY269-{rho_label}ref_260922"; os.makedirs(ref, exist_ok=True)
    ditj = pd.read_csv(f"{SRC}/{pfx}_1_data_dit.csv"); K = int(ditj[ditj.item_label == URU].cat_length.iloc[0])
    it = 'categorical' if 'spr' in rho_label else 'out-of-7'
    pd.DataFrame([dict(item_label=URU, item_type=it, item_type_id=1, cat_length=K,
                       item_label_short=np.nan, construct='HAI titre', construct_long='A/Uruguay/07 (H3N2)',
                       item_high_label='higher_is_better', endpoint_measure=rho_label,
                       cat_labels=ditj[ditj.item_label == URU].cat_labels.iloc[0])]
                 ).to_csv(f"{ref}/{pfx}_1_data_dit.csv", index=False)
    interims = sorted(int(re.search(r'_i(\d+)_', f).group(1)) for f in glob.glob(f"{SRC}/{pfx}_i*_regression_training.pkl"))
    for k in interims:
        dp = pd.read_csv(f"{SRC}/{pfx}_{k}_data_dp1.csv")
        e = dp[(dp.item_label == URU) & dp.group_label.str.startswith(arm)].copy()   # this arm, both phases
        e['group'] = (e.group_label.str.contains('endline')).astype(int)             # phase 0/1
        e = e.sort_values(['group', 'pid']).reset_index(drop=True)
        orig_oid = e['oid'].to_numpy()
        e['group_label'] = np.where(e.group == 0, 'Baseline', 'Endline')
        e['item_group_id'] = e['group'] + 1; e['item_type'] = it; e['item_type_id'] = 1
        e['oid'] = np.arange(1, len(e) + 1); e['oidt'] = np.arange(1, len(e) + 1)
        cols = ['pid', 'pid_label', 'group', 'group_label', 'item_label', 'y', 'y_stan',
                'item_type', 'item_type_id', 'item_group_id', 'oid', 'oidt']
        e[cols].to_csv(f"{ref}/{pfx}_{k}_data_dp1.csv", index=False)
        d = xr.open_zarr(f"{SRC}/{pfx}_{k}_draws.zarr", group='posterior')
        yps = d['ypred'].values[:, :, orig_oid - 1]
        zp = f"{ref}/{pfx}_{k}_draws.zarr"; shutil.rmtree(zp, ignore_errors=True)
        xr.Dataset({'ypred': (('chain', 'draw', 'ypred_dim_0'), yps)}).to_zarr(zp, group='posterior', mode='w')
        x = pd.read_pickle(f"{SRC}/{pfx}_i{k}_regression_training.pkl")
        r = x[(x.rho_label == joint_rho) & (x.item_label == URU)][['draw', 'pps_rho_x', 'pps_H1_x']].copy()
        r['item_label'] = URU; r['item_type'] = it; r['item_high_label'] = 'higher_is_better'
        r.rename(columns={'pps_rho_x': 'pps_ratio_x'}).to_pickle(f"{ref}/{pfx}_i{k}_regression_training.pkl")
    return ref, len(interims)


if __name__ == "__main__":
    for rho_label, arm, netkey, joint_rho in MANIFEST:
        netdir, widetok, warp, case_c = NET[netkey]
        ref, nint = build_ref(rho_label, arm, joint_rho)
        out = f"{SB}/py-immport-SDY269-{rho_label}deploy_260922"
        ctag = f"sdy269-{rho_label}"
        shutil.rmtree(f"/private/tmp/claude-501/-Users-or105-git-bIRTistic/49372e22-d11b-4882-a442-d1c60bcbdfb0/scratchpad/ragged_deploy_cells_{ctag}", ignore_errors=True)
        env = dict(os.environ, RAGD_BASE=f"{SB}/{netdir}", RAGD_RGE=ref, RAGD_OUT=out,
                   RAGD_WARP=warp, RAGD_WIDETOK=str(widetok), RAGD_KMAX='10', RAGD_CASE_C=str(case_c),
                   RAGD_BVM='1', RAGD_HEADFT='expand', RAGD_NREF='64', RAGD_CTAG=ctag, RAGD_PDFS='0')
        print(f"\n===== deploy {rho_label} ({arm}, net={netkey}, warp={warp}) -> {out}")
        import subprocess
        subprocess.run([sys.executable, DRV], env=env, cwd=os.path.join(os.path.dirname(__file__), '..'))
    print("\nfederated head-to-head level deploys done.")
