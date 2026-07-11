#!/usr/bin/env python3
"""
Cache the training-data cohort + amortiser training samples for the
Binomial §3.1 / §11.6 prototype.

Two artifacts are produced:

- ``binomial_sim_cohort.pkl``: the fixed simulated Bernoulli cohort of
  size ``N = 500`` spaced over one year, plus the ``dp`` /
  ``dit`` / ``di`` (monthly interim grid) tables. Consumed by the
  deployment and comparison scripts so they share a single ground-truth
  cohort.
- ``binomial_amortiser_training_data.pkl``: a stack of
  ``(features, rho)`` batches drawn by
  :meth:`BinomialModel.make_training_data_with_features`. Consumed by
  the training loop in
  ``Binomial_interim_analysis_amortise_endptx_on_wz_with_features_fixed_qpsi_MLP_loss_multiquantilehead.py``
  when ``resume=True``.

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/Binomial_interim_analyses_make_sim_data.py
"""

# %%

import os
import sys
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

from model_binomial import BinomialModel

print("✓ Imports successful")

# %%

# =============================================================================
# Configuration (must stay in sync with the deployment / test scripts)
# =============================================================================

seed = 123
N = 500
TRUE_P = 0.4
PRIOR_A = 1.0
PRIOR_B = 1.0
P_0 = 0.5

# Amortiser training-data budget. Batches are (2048 samples wide) so the
# training loop can consume them via a simple index slice.
TRAIN_TOTAL = 65_536         # >= NET_BATCH * NET_STEPS / re-use factor
TRAIN_BATCH_STORE = 4096
NET_N_MAX = N

dir_out = "/Users/or105/sandbox/bIRTistic/py-binomial-interim-make-sim-data-260702"
os.makedirs(dir_out, exist_ok=True)
print(f"Output dir: {dir_out}")

# %%

# =============================================================================
# Simulate the fixed Bernoulli cohort of size N = 500 over one year and
# build the monthly interim grid.
# =============================================================================

start = pd.Timestamp('2025-01-01')
end = pd.Timestamp('2025-12-31')
dates = pd.date_range(start, end, periods=N)
rng = np.random.default_rng(seed)
y = rng.binomial(1, TRUE_P, size=N)

dp = pd.DataFrame({
    'pid': np.arange(N),
    'submission_date': dates,
    'y': y.astype(int),
})

dit = pd.DataFrame({
    'item_label': ['intervention'],
    'item_type': ['binomial'],
    'item_high_label': ['higher_is_better'],
})

month_starts = pd.date_range(start, end, freq='MS')
di = pd.DataFrame({
    'interim_id': range(1, len(month_starts) + 1),
    'month_start': month_starts,
    'interim_date': month_starts + pd.offsets.MonthEnd(0),
})
di['interim_month_year'] = di['month_start'].dt.strftime('%Y-%b')

cohort = {
    'seed': seed,
    'N': N,
    'TRUE_P': TRUE_P,
    'PRIOR_A': PRIOR_A,
    'PRIOR_B': PRIOR_B,
    'P_0': P_0,
    'dp': dp,
    'dit': dit,
    'di': di,
}
cohort_path = os.path.join(dir_out, 'binomial_sim_cohort.pkl')
pd.to_pickle(cohort, cohort_path)
print(f"  saved cohort to {cohort_path}")

# %%

# =============================================================================
# Amortiser training data. Uses the model instance's prior_a, prior_b, p_0.
# =============================================================================

model = BinomialModel(
    dit=dit,
    dcati=dp.assign(oid=np.arange(1, N + 1, dtype=int)),
    prior_a=PRIOR_A, prior_b=PRIOR_B, p_0=P_0, seed=seed,
)

n_batches = TRAIN_TOTAL // TRAIN_BATCH_STORE

# ---- Features-fixed cache -------------------------------------------------
rng_train = np.random.default_rng(seed + 1)
fixed_batches = []
fixed_rho = []
for b in range(n_batches):
    batch, rho = model.make_training_data_with_features(
        rng_train, TRAIN_BATCH_STORE, n_max=NET_N_MAX,
    )
    fixed_batches.append(batch['features'])
    fixed_rho.append(rho)

features_fixed_all = np.concatenate(fixed_batches, axis=0)
rho_fixed_all = np.concatenate(fixed_rho, axis=0)
print(f"  features-fixed: features {features_fixed_all.shape}, "
      f"rho {rho_fixed_all.shape}")

training_fixed = {
    'batch': {'features': features_fixed_all},
    'rho': rho_fixed_all,
    'seed': seed + 1,
    'n_max': NET_N_MAX,
    'fixed_total': True,
    'prior_a': PRIOR_A,
    'prior_b': PRIOR_B,
    'p_0': P_0,
}
fixed_path = os.path.join(
    dir_out, 'binomial_amortiser_training_data_features_fixed.pkl',
)
pd.to_pickle(training_fixed, fixed_path)
print(f"  saved features-fixed training data to {fixed_path}")

# ---- Features-MLP cache (raw padded item sequences) ------------------------
rng_train_mlp = np.random.default_rng(seed + 2)
x_all = []
mask_x_all = []
z_all = []
mask_z_all = []
sizes_all = []
rho_mlp_all = []
for b in range(n_batches):
    batch, rho = model.make_training_data_with_raw_sequences(
        rng_train_mlp, TRAIN_BATCH_STORE, n_max=NET_N_MAX,
    )
    x_all.append(batch['x'])
    mask_x_all.append(batch['mask_x'])
    z_all.append(batch['z'])
    mask_z_all.append(batch['mask_z'])
    sizes_all.append(batch['sizes'])
    rho_mlp_all.append(rho)

training_mlp = {
    'batch': {
        'x':       np.concatenate(x_all, axis=0),
        'mask_x':  np.concatenate(mask_x_all, axis=0),
        'z':       np.concatenate(z_all, axis=0),
        'mask_z':  np.concatenate(mask_z_all, axis=0),
        'sizes':   np.concatenate(sizes_all, axis=0),
    },
    'rho': np.concatenate(rho_mlp_all, axis=0),
    'seed': seed + 2,
    'n_max': NET_N_MAX,
    'fixed_total': True,
    'prior_a': PRIOR_A,
    'prior_b': PRIOR_B,
    'p_0': P_0,
}
print(
    "  features-MLP: x {}, mask_x {}, z {}, mask_z {}, sizes {}, rho {}".format(
        training_mlp['batch']['x'].shape,
        training_mlp['batch']['mask_x'].shape,
        training_mlp['batch']['z'].shape,
        training_mlp['batch']['mask_z'].shape,
        training_mlp['batch']['sizes'].shape,
        training_mlp['rho'].shape,
    )
)
mlp_path = os.path.join(
    dir_out, 'binomial_amortiser_training_data_features_MLP.pkl',
)
pd.to_pickle(training_mlp, mlp_path)
print(f"  saved features-MLP training data to {mlp_path}")

print("\n✓ Binomial sim-data cache complete.")
