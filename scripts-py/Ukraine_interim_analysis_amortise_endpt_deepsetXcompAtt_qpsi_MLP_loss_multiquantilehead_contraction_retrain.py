#!/usr/bin/env python3
"""
Contraction-pressure retrain of the deepset PPS amortiser (§14.4.9 next
step). Adds the two cheap mechanisms proposed in the memory note
`deepset-contraction-calibration`:

  (1) BvM aux features: aux = (n/N, m/N, 1/sqrt(n), 1/sqrt(m)) so the head
      has the exact posterior-SD rate to scale the predicted spread.
  (3) n-stratified mini-batches: each batch draws n uniformly across
      log-spaced strata of [2, N) (vs iid U), so the n-dependence gets a
      clean gradient every step.

Trains to a NEW dir, then runs the prior-predictive contraction probe on
the fresh net: does the predicted 90% width now scale with n on the
TRAINING distribution? (Fast checkpoint before the full deployment test.)

The baseline net's width was flat on its own prior-predictive
distribution (slope -0.02), and prediction error SATURATES for n >~ 15
because training conditions on BOTH cohorts (n + m = N always identifies
the full-population endpoint). So this run also tells us whether the
n-signal exists to be learned at all.
"""
import os
import sys
import time
import warnings
from functools import partial
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'python'))
import numpy as np
import pandas as pd
import jax
import jax.numpy as jnp

warnings.filterwarnings('ignore')
from amortiser_common import train as amortiser_train, save_trained_model, load_fitted_model
from amortiser_pps_features_deepsetXcompAtt_qpsi_MLP_loss_multiquantilehead import (
    Amortiser_PPS_features_deepsetXcompAtt_qpsi_MLP_loss_multiquantilehead as DSNet)
from model_pcm import PartialCreditModel

SB = "/Users/or105/sandbox/bIRTistic"
RGE = f"{SB}/py-ukraine-interim-with-regression-on-endptx-wz-260601"
dir_out = f"{SB}/py-ukraine-interim-amortise-endptx-on-wz-with-features-deepsetXcompAtt-qpsi-MLP-loss-multiquantilehead-contraction-260811"
os.makedirs(dir_out, exist_ok=True)
file_prefix = "pcm_1_interim"
N_FULL = 503
NET_STEPS = int(os.environ.get('UKR_CONTRACT_STEPS', 8000))
NET_BATCH = 32
NET_QUERIES = 4
NET_TAUS = (0.05, 0.25, 0.5, 0.75, 0.95)

dit = pd.read_csv(f"{RGE}/{file_prefix}_1_data_dit.csv")
wa1 = pd.read_pickle(f"{RGE}/{file_prefix}_i1_regression_training.pkl")
items = (wa1[['item_label', 'item_type', 'item_high_label']].drop_duplicates()
         .sort_values(['item_type', 'item_label']).reset_index(drop=True))
J = len(items)
itype = items.item_type.map({'out-of-7': 0., 'categorical': 1.}).to_numpy(np.float32)
ihigh = items.item_high_label.map({'higher_is_better': 1., 'lower_is_better': 0.}).to_numpy(np.float32)
items = items.merge(dit[['item_label', 'cat_length']].drop_duplicates(), on='item_label', how='left')
K_O7 = int(items.loc[items.item_type == 'out-of-7', 'cat_length'].iloc[0])
K_CAT = int(items.loc[items.item_type == 'categorical', 'cat_length'].iloc[0])

_base = partial(
    PartialCreditModel.make_training_data_with_participant_tokens_prior,
    item_type_idx=itype, item_high_idx=ihigh, n_max=N_FULL,
    K_out_of_7=K_O7, K_categorical=K_CAT, categorical_threshold=2,
    queries_per_sample=NET_QUERIES, target_clip=20.0, threshold_scale=3.5)

# log-spaced n strata for mechanism (3)
_STRATA = np.geomspace(2, N_FULL - 1, 12)


def sample_contraction(rng, S):
    """Wrap the base sampler: (3) n-stratified draw, (1) BvM aux features."""
    batch, rho = _base(rng, S)
    # (1) append 1/sqrt(n), 1/sqrt(m) to aux (recover n, m from n/N, m/N)
    n = np.maximum(batch['aux'][:, 0] * N_FULL, 1.0)
    m = np.maximum(batch['aux'][:, 1] * N_FULL, 1.0)
    batch['aux'] = np.concatenate(
        [batch['aux'], (1.0 / np.sqrt(n))[:, None].astype(np.float32),
         (1.0 / np.sqrt(m))[:, None].astype(np.float32)], axis=-1)
    return batch, rho


# ---- train (or load) ----
CKPT = f"{dir_out}/{file_prefix}_amortised_pps_net.pkl"
if os.path.exists(CKPT):
    print(f"loading cached contraction net {CKPT}")
    fit = load_fitted_model(CKPT)
else:
    print(f"training contraction-pressure deepset ({NET_STEPS} steps, aux-dim=4) ...")
    fit = amortiser_train(
        sample_contraction, DSNet, net_kwargs={},
        pps_ProbH1_lwr_quantiles_mesh=NET_TAUS, num_steps=NET_STEPS,
        batch_size=NET_BATCH, lr=1e-3, seed=123, verbose=True)
    save_trained_model(fit, CKPT)
    print(f"saved -> {CKPT}  ({fit.get('training_mins', float('nan')):.1f} min)")

# ---- prior-predictive contraction probe on the fresh net ----
nkw = dict(fit['net_kwargs']); nkw.setdefault('num_quantiles', 5); net = DSNet(**nkw)


@jax.jit
def fwd(p, b):
    return net.apply(p, b)


rng = np.random.default_rng(7)
batch, rho = sample_contraction(rng, 3000)
B = batch['aux'].shape[0]; qs = np.empty((B, 5))
for i in range(0, B, 500):
    sl = slice(i, i + 500)
    bb = {k: jnp.asarray(v[sl]) for k, v in batch.items()}
    qs[sl] = np.maximum.accumulate(np.asarray(fwd(fit['params'], bb)), 1)
n = np.rint(batch['aux'][:, 0] * N_FULL).astype(int)
width = qs[:, 4] - qs[:, 0]; err = np.abs(rho - qs[:, 2])
df = pd.DataFrame(dict(n=n, width=width, err=err))
d = df[(df.width > 0) & (df.n > 0)]
sl_w = float(np.polyfit(np.log(d.n), np.log(d.width), 1)[0])
sl_e = float(np.polyfit(np.log(d[d.err > 0].n), np.log(d[d.err > 0].err), 1)[0])
bins = np.geomspace(2, N_FULL, 9); df['nb'] = pd.cut(df.n, bins)
g = df.groupby('nb', observed=True).agg(nmid=('n', 'median'), width=('width', 'median'),
                                        err=('err', 'median'), c=('n', 'size'))
print("\nPRIOR-PREDICTIVE probe on CONTRACTION-RETRAINED net:")
print(g.to_string())
print(f"\nwidth log-log slope: {sl_w:.3f}  (baseline was -0.02; target ~ SVI -0.59)")
print(f"error log-log slope: {sl_e:.3f}")
print(f"median width n<80: {df[df.n<80].width.median():.3f}  n>400: {df[df.n>400].width.median():.3f}")
pd.DataFrame(dict(slope_width=[sl_w], slope_err=[sl_e])).to_csv(
    f"{dir_out}/{file_prefix}_prior_contraction_probe.csv", index=False)
print(f"\nprobe done -> {dir_out}")
