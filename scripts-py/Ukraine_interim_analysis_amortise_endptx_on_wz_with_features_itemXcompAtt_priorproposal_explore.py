"""
§14.4.19 proposal ablation for the item-token `itemXcompAtt` amortiser.

Question: what exactly erodes the sharp z-response of the §14.2.5 hand-token
encoder? We hold the ARCHITECTURE fixed (itemXcompAtt) and swap only the
TRAINING PROPOSAL, then measure the encoder's z-tracking with one fixed probe
identical across cells.

Proposal knobs (env):
  IXC_TAG       output-dir suffix / cell name
  IXC_DECOUPLE  1 -> draw m independently of n (variable N, small m possible)
  IXC_DR        domain-randomise z time-effect  beta_z = beta + N(0, DR)
  IXC_NLO/NHI   fixed-total band  N ~ U[NLO, NHI]  (coupled m = N - n)
  IXC_MLO/MHI   z-cohort size band when decoupled
  IXC_STEPS     training steps (default 6000)

  cell recipes:
    fixed503-noDR   DECOUPLE=0 NLO=NHI=503 DR=0     (= the §14.2.5 recipe)
    fixed503-DR     DECOUPLE=0 NLO=NHI=503 DR=0.8   (isolate DR)
    decouple-noDR   DECOUPLE=1 MLO=2 MHI=503 DR=0   (isolate decoupling)
    decouple-DR     DECOUPLE=1 MLO=2 MHI=503 DR=0.8 (full deepset-style)
    Nband-noDR      DECOUPLE=0 NLO=500 NHI=3000 DR=0 (#2 practical: vary N, big m)

z-tracking probe (fixed seed, same for every cell): pick R reference scenarios;
per scenario fix an x-cohort, then draw S future cohorts z^(s) from a posterior-
like jitter of the time effect (rho^(s) varies through beta^(s), z^(s) reflects
it, x is held fixed); forward the net; corr(median_s, rho^(s)) per item-type.
Because x is fixed, the median varies only through the encoder's z-response, so
the correlation is exactly the z-tracking measured in §14.4.18.
"""

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path

project_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(project_root / 'python'))

from model_pcm import PartialCreditModel
from amortiser_common import train as amortiser_train, save_trained_model, _jax_pytree

# ---------------- config ----------------
DIR_RGE = "/Users/or105/sandbox/bIRTistic/py-ukraine-interim-with-regression-on-endptx-wz-260601"
SB = "/Users/or105/sandbox/bIRTistic"
TAG = os.environ.get('IXC_TAG', 'fixed503-noDR')
DECOUPLE = os.environ.get('IXC_DECOUPLE', '0') == '1'
DR = float(os.environ.get('IXC_DR', '0'))
N_LO = int(os.environ.get('IXC_NLO', '503')); N_HI = int(os.environ.get('IXC_NHI', '503'))
M_LO = int(os.environ.get('IXC_MLO', '2')); M_HI = int(os.environ.get('IXC_MHI', '503'))
STEPS = int(os.environ.get('IXC_STEPS', '6000'))
BATCH = int(os.environ.get('IXC_BATCH', '64'))
AUX_REF = 503.0
dir_out = f"{SB}/py-ukraine-interim-amortise-endptx-itemXcompAtt-priorproposal-{TAG}"
os.makedirs(dir_out, exist_ok=True)
file_prefix = "pcm_1_interim"
NET_TAUS = (0.05, 0.25, 0.5, 0.75, 0.95)
NET_SEED = 123
print(f"CONFIG tag={TAG} decouple={DECOUPLE} DR={DR} N=[{N_LO},{N_HI}] "
      f"m=[{M_LO},{M_HI}] steps={STEPS}")

# ---------------- item metadata (mirror the itemXcompAtt script) ----------------
wa0 = pd.read_pickle(f"{DIR_RGE}/{file_prefix}_i1_regression_training.pkl")
items_df = (wa0[['item_label', 'item_type', 'item_high_label']].drop_duplicates()
            .sort_values(['item_type', 'item_label']).reset_index(drop=True))
J = len(items_df)
item_type_idx = items_df['item_type'].map({'out-of-7': 0.0, 'categorical': 1.0}).to_numpy(np.float32)
item_high_idx = items_df['item_high_label'].map(
    {'higher_is_better': 1.0, 'lower_is_better': 0.0}).to_numpy(np.float32)
dit_df = pd.read_csv(f"{DIR_RGE}/{file_prefix}_1_data_dit.csv")
items_df = items_df.merge(dit_df[['item_label', 'cat_length']].drop_duplicates(),
                          on='item_label', how='left')
item_klevels = items_df['cat_length'].to_numpy(np.float32)
K_O7 = int(item_klevels[item_type_idx == 0.0][0])
K_CAT = int(item_klevels[item_type_idx == 1.0][0])
typ = items_df['item_type'].tolist()

FLEX_KW = dict(item_type_idx=item_type_idx, item_high_idx=item_high_idx,
               N_lo=N_LO, N_hi=N_HI, decouple=DECOUPLE, m_lo=M_LO, m_hi=M_HI,
               dr=DR, aux_ref=AUX_REF, K_out_of_7=K_O7, K_categorical=K_CAT,
               categorical_threshold=2, target_clip=20.0, threshold_scale=3.5,
               M_ref=300)

# ---------------- per-item target std (pilot) ----------------
pilot_rng = np.random.default_rng(NET_SEED + 1)
_pb, _prho = PartialCreditModel.make_training_data_with_item_tokens_prior_flex(
    pilot_rng, 2000, **FLEX_KW)
item_std = _prho.reshape(2000, J).std(0).astype(np.float32)
item_std = np.where(item_std > 1e-6, item_std, 1.0)
np.save(f"{dir_out}/{file_prefix}_item_std.npy", item_std)


def sample_fn(rng, S):
    batch, rho = PartialCreditModel.make_training_data_with_item_tokens_prior_flex(
        rng, S, **FLEX_KW)
    return batch, rho / item_std[batch['query_idx']]


# ---------------- train ----------------
from amortiser_pps_features_itemXcompAtt_qpsi_MLP_loss_multiquantilehead import (
    Amortiser_PPS_features_itemXcompAtt_qpsi_MLP_loss_multiquantilehead as Net)
print(f"Training itemXcompAtt ({STEPS} steps) ...")
fit = amortiser_train(sample_fn, Net,
                      net_kwargs={'q_tok_hidden': (32, 32), 'embed_dim': 32, 'hidden_dims': (64, 64)},
                      pps_ProbH1_lwr_quantiles_mesh=NET_TAUS,
                      num_steps=STEPS, batch_size=BATCH, lr=1e-3, seed=NET_SEED, verbose=True)
save_trained_model(fit, f"{dir_out}/{file_prefix}_amortised_pps_net.pkl")
apply_fn, params = fit['apply_fn'], fit['params']

# ---------------- z-tracking probe (fixed seed, identical across cells) ----------------
PR = np.random.default_rng(20260824)     # SAME across all cells
R_SCEN, S_DRAW = 25, 45
N_X, M_Z = 132, 371                       # deployment-like sizes (n=132, m=N-n at N=503)
JIT = 0.5                                 # posterior-like jitter of the time effect
type_groups = [(np.flatnonzero(item_type_idx == 0.0), K_O7),
               (np.flatnonzero(item_type_idx == 1.0), K_CAT)]
meta = np.stack([item_type_idx, item_high_idx], -1).astype(np.float32)
eye = np.eye(J, dtype=np.float32)
theta_ref = PR.standard_normal(400)


def _sim(theta, beta, tau, lam, js, K):
    latent = theta[:, None, None] + beta[None, None, :]
    inc = lam[None, :, :, None] * (latent[:, :, :, None] - tau[None, :, :, :])
    eta = np.concatenate([np.zeros((theta.shape[0], js.size, 2, 1)), np.cumsum(inc, -1)], -1)
    eta -= eta.max(-1, keepdims=True); p = np.exp(eta); p /= p.sum(-1, keepdims=True)
    cdf = np.cumsum(p, -1); u = PR.random((theta.shape[0], js.size, 2))
    k = (cdf < u[..., None]).sum(-1)
    return (k.astype(np.float32) / (K - 1)), p


def _endpoint(p_bar, js, K):
    if K == K_O7:
        ks = np.arange(1, K + 1, dtype=np.float64); w = (p_bar * ks[None, None, :]).sum(-1)
    else:
        w = p_bar[:, :, 1:].sum(-1)
    wb = np.maximum(w[:, 0], 1e-3); we = w[:, 1]
    high = item_high_idx[js] > 0.5
    return np.where(high, we / wb - 1.0, 1.0 - we / wb)


def _tokens(w_x, w_z, chg, n):
    wb = np.maximum(w_z[:, 0], 1e-3); we = w_z[:, 1]
    rat = np.clip(np.where(item_high_idx > 0.5, we / wb - 1.0, 1.0 - we / wb), -20, 20)
    tb = np.concatenate([w_x, w_z, rat[:, None], meta], -1)
    rho_m = np.nan_to_num(pd.DataFrame(chg).corr(method='spearman').to_numpy(), nan=0.0)
    np.fill_diagonal(rho_m, 1.0); lam_s = n / (n + 50.0)
    kh = lam_s * rho_m + (1 - lam_s) * eye
    tok = np.empty((J, J, 8), np.float32)
    tok[..., :7] = tb[None]; tok[..., 7] = kh
    return tok


def _fwd_batch(tok_stack):
    """tok_stack (S, J, J, 8) -> medians (S, J). ONE forward call."""
    Sd = tok_stack.shape[0]; B = Sd * J
    batch = dict(tokens=tok_stack.reshape(B, J, 8), mask=np.ones((B, J), np.float32),
                 query_idx=np.tile(np.arange(J, dtype=np.int32), Sd),
                 aux=np.tile(np.array([[N_X / AUX_REF, M_Z / AUX_REF]], np.float32), (B, 1)))
    preds = np.asarray(apply_fn(params, _jax_pytree(batch)))     # (S*J, K_tau)
    return preds[:, 2].reshape(Sd, J)                            # median (tau=0.5)


corr_rows = []
for r in range(R_SCEN):
    beta0 = PR.standard_normal(2)
    tpar = {}
    for js, K in type_groups:
        tau = 3.5 * PR.standard_normal((js.size, 2, K - 1))
        lam = np.abs(PR.standard_t(3.0, size=(js.size, 2))); lam[0, 0] = 1.0
        tpar[K] = (js, tau, lam)
    # fixed x-cohort from beta0
    theta_x = PR.standard_normal(N_X)
    w_x = np.zeros((J, 2), np.float32); chg = np.zeros((N_X, J))
    for K, (js, tau, lam) in tpar.items():
        Yx, _ = _sim(theta_x, beta0, tau, lam, js, K)
        w_x[js] = Yx.mean(0); chg[:, js] = Yx[:, :, 1] - Yx[:, :, 0]
    rhos = np.zeros((S_DRAW, J)); tok_stack = np.zeros((S_DRAW, J, J, 8), np.float32)
    for s in range(S_DRAW):
        beta_s = beta0 + PR.normal(0, JIT, 2)                     # posterior-like draw
        theta_z = PR.standard_normal(M_Z)
        w_z = np.zeros((J, 2), np.float32)
        for K, (js, tau, lam) in tpar.items():
            _, p_ref = _sim(theta_ref, beta_s, tau, lam, js, K)
            rhos[s, js] = np.clip(_endpoint(p_ref.mean(0), js, K), -20, 20)
            Yz, _ = _sim(theta_z, beta_s, tau, lam, js, K)
            w_z[js] = Yz.mean(0)
        tok_stack[s] = _tokens(w_x, w_z, chg, N_X)
    meds = _fwd_batch(tok_stack)                                  # (S_DRAW, J)
    for j in range(J):
        if np.std(meds[:, j]) < 1e-9 or np.std(rhos[:, j]) < 1e-9:
            continue
        corr_rows.append(dict(typ=typ[j], c=float(np.corrcoef(meds[:, j], rhos[:, j])[0, 1])))

d = pd.DataFrame(corr_rows); g = d.groupby('typ')['c'].mean()
print(f"\nZ-TRACKING [{TAG}]  o7={g.get('out-of-7', float('nan')):.3f}  "
      f"cat={g.get('categorical', float('nan')):.3f}  (n_pairs={len(d)})")
pd.DataFrame([dict(tag=TAG, decouple=DECOUPLE, dr=DR, n_lo=N_LO, n_hi=N_HI,
                   o7=float(g.get('out-of-7', np.nan)),
                   cat=float(g.get('categorical', np.nan)))]).to_csv(
    f"{dir_out}/{file_prefix}_ztracking.csv", index=False)
print("DONE", TAG)
