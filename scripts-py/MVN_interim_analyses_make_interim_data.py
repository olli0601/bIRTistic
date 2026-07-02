#!/usr/bin/env python3
"""
Standalone builder of per-(J, interim) future-data block ``zi`` from the
simulated x cohort in ``mvn_sim_data.pkl``.

Loads the simulated x cohorts from the simulations directory (the pickle
produced by ``MVN_interim_analyses_make_sim_data.py``), fits the outer
x-cohort HMC posterior per (J, interim) (``resume=True`` on cached
zarrs), expands ``ypred`` into a long-form ``zi`` of ``interim_m`` new
participants via :meth:`MVNModel.get_interim_z_from_ypredi`, and stores
both ``zi`` and the first ``PPS_Z_TOTAL`` posterior ``mu`` draws (aligned
on the same draw index) in a single per-J pkl in the simulations
directory.

Each ``mvn_J{J}_interim_data.pkl`` is a dict keyed by ``interim_id`` ->
``{'zi', 'mu_draws', 'n_obs', 'interim_m', 'interim_date',
'interim_month_year', 'dpi'}``. Downstream nested-MC / regression
scripts load these pkls and never need to touch the zarrs themselves.

Skips J for which the pkl already exists.

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/MVN_interim_analyses_make_interim_data.py
"""

# %%

import os
import sys
import warnings
from pathlib import Path

try:
    project_root = Path(__file__).resolve().parent.parent
except NameError:
    project_root = Path.cwd()
    if project_root.name == 'scripts-py':
        project_root = project_root.parent

sys.path.insert(0, str(project_root / 'python'))

import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')

from model_mvn import MVNModel

print("Imports successful")

# %%

# =============================================================================
# Configuration. ``DIR_OUT`` is the simulations directory shared with
# ``MVN_interim_analyses_make_sim_data.py``; interim_data pkls + outer
# HMC zarrs live there.
# =============================================================================

DIR_OUT = os.path.join(
    "/Users/or105/sandbox/bIRTistic", "py-mvn-interim-simulations-260609",
)
SIM_PKL = os.path.join(DIR_OUT, 'mvn_sim_data.pkl')

PPS_Z_TOTAL = 4000
HMC_CHAINS = 4
HMC_ITER_WARMUP = 500
HMC_ITER_SAMPLING = 1000

# %%

# =============================================================================
# Load simulated x cohorts.
# =============================================================================

if not os.path.exists(SIM_PKL):
    raise FileNotFoundError(
        f"{SIM_PKL} missing; "
        f"run MVN_interim_analyses_make_sim_data.py first."
    )
sim_data = pd.read_pickle(SIM_PKL)
simu_params = sim_data['simu_params']
di = sim_data['interim_grid']
print(f"Loaded sim_data ({len(simu_params['J_grid'])} J cells) from {SIM_PKL}")

# %%

# =============================================================================
# Build zi + mu_draws per (J, interim).
# =============================================================================

for J in simu_params['J_grid']:
    pkl_path = os.path.join(DIR_OUT, f'mvn_J{J}_interim_data.pkl')
    if os.path.exists(pkl_path):
        print(f"J={J}: cached interim data at {pkl_path}; skipping.")
        continue

    cell = sim_data['cells'][J]
    dit, dp, R_chol = cell['dit'], cell['dp'], cell['R_chol']
    model_init_kwargs = dict(
        J=J, K_chol=R_chol, sigma=simu_params['sigma'],
        prior_tau=simu_params['prior_tau'],
        mu_0_baseline=simu_params['mu_0_baseline'],
    )
    dir_J = os.path.join(DIR_OUT, f"nested_mc_J{J}")
    os.makedirs(dir_J, exist_ok=True)

    zi_by_interim = {}
    for i in range(len(di)):
        interim_id = int(di.iloc[i]['interim_id'])
        interim_date = di.iloc[i]['interim_date']
        interim_month_year = di.iloc[i]['interim_month_year']
        dpi = MVNModel.get_interim_data_x(
            dp[dp['submission_date'] <= interim_date],
        )
        n_obs = int(dpi['pid'].nunique())
        interim_m = simu_params['N_full'] - n_obs
        if interim_m <= 0:
            print(f"  J={J} interim {interim_id}: m=0, skipping.")
            continue
        interim_prefix = os.path.join(dir_J, f"mvn_pps_i{interim_id}")
        model_x = MVNModel(
            dit=dit, dcati=dpi, seed=simu_params['seed'],
            **model_init_kwargs,
        )
        fit_x = model_x.fit_pyro_hmc(
            output_file_prefix=f"{interim_prefix}_x",
            chains=HMC_CHAINS, iter_warmup=HMC_ITER_WARMUP,
            iter_sampling=HMC_ITER_SAMPLING,
            save_to_file=True, resume=True, verbose=False,
        )
        zi = model_x.get_interim_z_from_ypredi(
            f"{interim_prefix}_x_draws.zarr", interim_m,
            pps_z_total=PPS_Z_TOTAL, seed=simu_params['seed'],
            keep_order=True,
        )
        mu_draws_all = fit_x['posterior_samples']['mu']      # (S_total, J)
        mu_draws = np.asarray(mu_draws_all[:PPS_Z_TOTAL, :], dtype=np.float64)
        zi_by_interim[interim_id] = {
            'zi':                 zi,
            'mu_draws':           mu_draws,
            'n_obs':              n_obs,
            'interim_m':          interim_m,
            'interim_date':       interim_date,
            'interim_month_year': interim_month_year,
            'dpi':                dpi,
        }
        print(f"  J={J} interim {interim_id} ({interim_date.date()}):"
              f" zi has {len(zi)} rows, m={interim_m},"
              f" S={PPS_Z_TOTAL}; mu_draws {mu_draws.shape}.")

    pd.to_pickle(zi_by_interim, pkl_path)
    print(f"J={J}: saved {len(zi_by_interim)} interim zi blocks to {pkl_path}")
