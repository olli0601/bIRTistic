"""Build a single-endpoint reference grid for the amortiser deploy, written into the shared
amortiser dir (IMMPORT_AMORTDIR) under an endpoint-specific prefix pcm_<rho>_interim, so BOTH
endpoints' amortiser inputs+outputs live in ONE dir and the base SVI dir stays SVI-only. The
ragged deploy driver pivots the reference pkls on a single `pps_ratio_x`; our ImmPort endpoint
pkls carry TWO rho (SPR + GMFR) in long form. For the chosen rho we (a) write
pcm_<rho>_interim_i{k}_regression_training.pkl (rho as pps_ratio_x) into the amortiser dir, and
(b) symlink the SVI dp1.csv / draws.zarr / dit.csv under the same prefix. Deploy with
RAGD_RGE=RAGD_OUT=<amortiser dir>, RAGD_FILEPREFIX=pcm_<rho>_interim. Idempotent."""
# ---- boilerplate ----

import os, sys, glob, re
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
import pandas as pd

SB = os.environ.get('IMMPORT_SB', "/Users/or105/sandbox/bIRTistic")
STUDIES = os.environ.get('IMMPORT_SDY', 'SDY312,SDY314').split(',')
RHO_ID = int(os.environ.get('IMMPORT_RHO_ID', '2'))       # 1 = SPR, 2 = GMFR
RHO_SUFFIX = {1: 'spr', 2: 'gmfr'}.get(RHO_ID, f'rho{RHO_ID}')
SRC_PREFIX = "pcm_1_interim"                              # source SVI grid prefix
OUT_PREFIX = f"pcm_{RHO_SUFFIX}_interim"                  # endpoint prefix in the amortiser dir


def _link(src, link):
    if os.path.islink(link):
        os.remove(link)
    if not os.path.exists(link):
        os.symlink(src, link)


def build_refgrid(study):
    src = f"{SB}/py-immport-{study}_260918"               # SVI dir (read)
    dst = os.environ.get('IMMPORT_AMORTDIR',              # amortiser dir (write); default = SVI dir
                         f"{SB}/py-immport-{study}-amortise-deepsetXcompAtt-itemamortise-J64-ftheadexpand-bvm_260919")
    os.makedirs(dst, exist_ok=True)
    ks = sorted(int(re.search(r'_i(\d+)_', p).group(1))
                for p in glob.glob(f"{src}/{SRC_PREFIX}_i*_regression_training.pkl"))
    _link(f"{src}/{SRC_PREFIX}_1_data_dit.csv", f"{dst}/{OUT_PREFIX}_1_data_dit.csv")
    n = 0
    for k in ks:
        _link(f"{src}/{SRC_PREFIX}_{k}_data_dp1.csv", f"{dst}/{OUT_PREFIX}_{k}_data_dp1.csv")
        _link(f"{src}/{SRC_PREFIX}_{k}_draws.zarr", f"{dst}/{OUT_PREFIX}_{k}_draws.zarr")
        x = pd.read_pickle(f"{src}/{SRC_PREFIX}_i{k}_regression_training.pkl")
        g = x[x['rho_id'] == RHO_ID].rename(columns={'pps_rho_x': 'pps_ratio_x'})
        g = g[['draw', 'item_label', 'item_type', 'item_high_label', 'pps_ratio_x', 'pps_H1_x']]
        g.to_pickle(f"{dst}/{OUT_PREFIX}_i{k}_regression_training.pkl")
        n += 1
    print(f"{study}: {n} {RHO_SUFFIX.upper()} pkls + data symlinks under {OUT_PREFIX}_* -> {dst}")
    return dst


if __name__ == "__main__":
    for study in STUDIES:
        build_refgrid(study.strip())
