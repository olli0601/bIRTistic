"""Item-amortised deepset PPS amortiser (variable J = 2..JMAX items).

Extends the ragged prior-proposal recipe (§14.4.2-4) to amortise over the ITEM SET:
each training batch draws one item count J_b ~ U{2..JMAX} (fixed-J-per-batch -> dense,
NO padding, NO mask) and B synthetic items each with its own category count K, endpoint
type (expected-score | caseness), direction and PCM parameters. The net is unchanged
(reads J = item_metadata.shape[1], generic metadata width M=3); target is standardised by
a single GLOBAL sigma, and the power-law head learns per-item scale C_j. head_mode=powerlaw.

Deploy the trained net on any application's items (2..JMAX) via the §14.4.6 head-ft + affine
and §14.4.7 BvM correction re-fit to that application's SVI.
"""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
import numpy as np, jax, jax.numpy as jnp, optax
from amortiser_common import save_trained_model
from amortiser_pps_features_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead import (
    Amortiser_PPS_features_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead as Net)

SB = "/Users/or105/sandbox/bIRTistic"
OUT_TAG = os.environ.get('PP_TAG', f"itemamortise-J{os.environ.get('PP_JMAX','64')}-260831")
dir_out = f"{SB}/py-ukraine-interim-amortise-deepsetXcompAtt-{OUT_TAG}"
os.makedirs(dir_out, exist_ok=True)
file_prefix = "pcm_1_interim"
THRESH, M_REF, N_REF, CLIP = 3.5, 2000, 503, 20.0
JMAX = int(os.environ.get('PP_JMAX', 64)); KMAX = 10
KVALS = np.array([2, 4, 5, 7, 8, 9, 10]); KP = np.array([.05, .2, .2, .2, .1, .15, .1]); KP /= KP.sum()
NET_STEPS = int(os.environ.get('PP_STEPS', 6000))
B = int(os.environ.get('PP_B', 24)); Q = 4
EMBED_DIM = int(os.environ.get('PP_EMBED', 32))
HEAD_HIDDEN = tuple(int(x) for x in os.environ.get('PP_HEAD', '64,64').split(','))
HEAD_SCALE = os.environ.get('PP_HEADSCALE', '0') == '1'   # #2: data-driven per-item scale at the head
HEAD_RAW = os.environ.get('PP_HEADRAW', '0') == '1'       # #3: raw level/ratio fingerprint at the head
META4 = os.environ.get('PP_META4', '0') == '1'            # #3: add caseness threshold c/K_max to metadata
META_DIM = 4 if META4 else 3
NET_KW = dict(num_quantiles=5, embed_dim=EMBED_DIM,
              q_tau_hidden=(EMBED_DIM, EMBED_DIM), q_tok_hidden=(EMBED_DIM, EMBED_DIM),
              q_query_hidden=(EMBED_DIM, EMBED_DIM), hidden_dims=HEAD_HIDDEN, head_mode='powerlaw',
              head_scale=HEAD_SCALE, head_raw=HEAD_RAW)
net = Net(**NET_KW)
TH_REF = np.random.default_rng(0).standard_normal(M_REF).astype(np.float64)   # fixed reference abilities
print(f"ITEM-AMORTISE tag={OUT_TAG} JMAX={JMAX} KMAX={KMAX} B={B} steps={NET_STEPS}")


def _probs(theta, beta, tau, lam, K):
    """theta (N,), beta (2,), tau (nK,2,K-1), lam (nK,2) -> p (N, nK, 2, K)."""
    latent = theta[:, None, None] + beta[None, None, :]                       # (N,1,2)
    inc = lam[None, :, :, None] * (latent[:, :, :, None] - tau[None, :, :, :])  # (N,nK,2,K-1)
    eta = np.concatenate([np.zeros((*inc.shape[:3], 1)), np.cumsum(inc, -1)], -1)
    eta -= eta.max(-1, keepdims=True); p = np.exp(eta); p /= p.sum(-1, keepdims=True)
    return p


def gen_example(rng, Jb):
    """One synthetic item set of size Jb + PCM params."""
    K = rng.choice(KVALS, size=Jb, p=KP)
    typ = (rng.random(Jb) < 0.5).astype(int)          # 0 expected-score, 1 caseness
    dr = (rng.random(Jb) < 0.5).astype(int)           # 1 higher-is-better
    cc = np.maximum(2, (K + 1) // 2 + 1)              # caseness threshold per item (2..K)
    beta = rng.standard_normal(2)
    taus = [(THRESH * rng.standard_normal((2, int(Kj) - 1))).astype(np.float64) for Kj in K]
    lam = np.abs(rng.standard_t(3.0, size=(Jb, 2))).astype(np.float64); lam[0, 0] = 1.0
    return dict(K=K.astype(int), typ=typ, dr=dr, cc=cc.astype(int), beta=beta, taus=taus, lam=lam)


def _wbe(p, Kval, idx_in_bucket, ex):
    """endpoint level w (2,) for one item, by its type."""
    a = idx_in_bucket
    if ex['typ_bucket'][a] == 0:                       # expected score E[k]
        ks = np.arange(1, Kval + 1)
        return (p[:, a] * ks[None, None, :]).sum(-1).mean(0)          # (2,)
    c = ex['cc_bucket'][a]                             # caseness Pr(y>=c)
    return p[:, a, :, c - 1:].sum(-1).mean(0)


def _bucket(ex, Kval):
    idx = np.flatnonzero(ex['K'] == Kval)
    tau = np.stack([ex['taus'][i] for i in idx]); lam = ex['lam'][idx]
    return idx, tau, lam


def rho_ex(ex):
    rho = np.zeros(len(ex['K']))
    for Kval in np.unique(ex['K']):
        idx, tau, lam = _bucket(ex, Kval)
        p = _probs(TH_REF, ex['beta'], tau, lam, Kval)                # (M_REF, nK, 2, K)
        ex['typ_bucket'] = ex['typ'][idx]; ex['cc_bucket'] = ex['cc'][idx]
        for a, i in enumerate(idx):
            w = _wbe(p, Kval, a, ex); wb = max(w[0], 1e-3); we = w[1]
            rho[i] = np.clip((we / wb - 1) if ex['dr'][i] else (1 - we / wb), -CLIP, CLIP)
    return rho


def sim_cohort(ex, nab, rng):
    Jb = len(ex['K']); Y = np.zeros((nab, Jb, 2), np.float32); th = rng.standard_normal(nab)
    for Kval in np.unique(ex['K']):
        idx, tau, lam = _bucket(ex, Kval)
        p = _probs(th, ex['beta'], tau, lam, Kval)                    # (nab,nK,2,K)
        cdf = np.cumsum(p, -1); u = rng.random((nab, len(idx), 2))
        k = (cdf < u[..., None]).sum(-1)
        Y[:, idx, :] = (k.astype(np.float32) / (Kval - 1))
    return Y


# ---- global sigma pilot (single scalar; power-law head learns per-item C_j) ----
_pr = np.random.default_rng(1)
_pilot = [rho_ex(gen_example(_pr, int(_pr.integers(2, JMAX + 1)))) for _ in range(600)]
SIG_GLOB = float(np.maximum(np.std(np.concatenate([r for r in _pilot])), 1e-3))
np.save(f"{dir_out}/{file_prefix}_item_std.npy", np.array([SIG_GLOB], np.float32))
np.save(f"{dir_out}/{file_prefix}_meta_dim.npy", np.array([META_DIM], np.int32))   # deploy reads this
print(f"global sigma = {SIG_GLOB:.3f} (target standardised by this scalar); META_DIM={META_DIM} head_scale={HEAD_SCALE} head_raw={HEAD_RAW}")

T_X, T_Z = B * 120, B * 120                                            # fixed ragged totals


def _sizes(rng, total):
    n = np.geomspace(4, N_REF - 3, B) * rng.uniform(0.7, 1.4, B)
    n = np.clip(np.round(n * total / n.sum()), 2, N_REF - 3).astype(int)
    while n.sum() != total:
        j = rng.integers(B); s = 1 if n.sum() < total else -1
        if 2 <= n[j] + s <= N_REF - 3:
            n[j] += s
    return n


def sample_batch(rng):
    Jb = int(rng.integers(2, JMAX + 1)); Qb = min(Q, Jb)
    exs = [gen_example(rng, Jb) for _ in range(B)]
    rhos = np.stack([rho_ex(ex) for ex in exs])                        # (B, Jb)
    nb = _sizes(rng, T_X); mb = _sizes(rng, T_Z)
    xl = [sim_cohort(exs[b], int(nb[b]), rng) for b in range(B)]
    zl = [sim_cohort(exs[b], int(mb[b]), rng) for b in range(B)]
    x_flat = np.concatenate(xl, 0); z_flat = np.concatenate(zl, 0)
    x_seg = np.repeat(np.arange(B), nb).astype(np.int32); z_seg = np.repeat(np.arange(B), mb).astype(np.int32)
    def _mc(ex):
        cols = [ex['typ'], ex['dr'], ex['K'] / KMAX]
        if META4:                                   # caseness threshold (0 for expected-score items)
            cols.append(np.where(ex['typ'] == 1, ex['cc'], 0) / KMAX)
        return np.stack(cols, -1)
    META = np.stack([_mc(ex) for ex in exs]).astype(np.float32)     # (B, Jb, META_DIM)
    qidx = np.stack([rng.integers(0, Jb, Qb) for _ in range(B)]).astype(np.int32)
    aux = np.stack([1/np.sqrt(nb), 1/np.sqrt(mb), nb/N_REF, mb/N_REF], -1).astype(np.float32)
    y = (rhos[np.arange(B)[:, None], qidx] / SIG_GLOB).astype(np.float32)
    return (dict(x_flat=jnp.asarray(x_flat), x_seg=jnp.asarray(x_seg),
                 z_flat=jnp.asarray(z_flat), z_seg=jnp.asarray(z_seg),
                 item_metadata=jnp.asarray(META), query_idx=jnp.asarray(qidx), aux=jnp.asarray(aux)),
            jnp.asarray(y))


b0, y0 = sample_batch(np.random.default_rng(7))
params = net.init(jax.random.PRNGKey(0), b0)


def interval_score(pr, y):
    pr = jnp.maximum.accumulate(pr, -1)
    l05, l25, med, u75, u95 = pr[..., 0], pr[..., 1], pr[..., 2], pr[..., 3], pr[..., 4]
    is10 = (u95 - l05) + 20.0*jnp.maximum(l05 - y, 0) + 20.0*jnp.maximum(y - u95, 0)
    is50 = (u75 - l25) + 4.0*jnp.maximum(l25 - y, 0) + 4.0*jnp.maximum(y - u75, 0)
    return jnp.mean(is10 + is50 + jnp.abs(y - med))


sched = optax.warmup_cosine_decay_schedule(1e-4, 1e-3, NET_STEPS // 10, NET_STEPS, 1e-5)
opt = optax.adam(sched); ost = opt.init(params)


def _loss(params, batch, y):
    return interval_score(net.apply(params, batch), y)


@jax.jit
def step(params, ost, batch, y):
    loss, g = jax.value_and_grad(_loss)(params, batch, y)
    up, ost = opt.update(g, ost, params)
    return optax.apply_updates(params, up), ost, loss


rng = np.random.default_rng(2)
print(f"training item-amortised net ({NET_STEPS} steps, B={B}, J~U[2,{JMAX}]) ...")
t0 = time.time()
for it in range(NET_STEPS):
    batch, y = sample_batch(rng)
    params, ost, loss = step(params, ost, batch, y)
    if (it + 1) % max(1, NET_STEPS // 20) == 0:
        print(f"  step {it+1}/{NET_STEPS}  IS={float(loss):.4f}")
mins = (time.time() - t0) / 60.0
fit = {'params': params, 'apply_fn': net.apply, 'net_class_name': type(net).__name__,
       'net_class_module': type(net).__module__, 'net_kwargs': dict(NET_KW),
       'training_mins': mins, 'global_sigma': SIG_GLOB, 'JMAX': JMAX, 'KMAX': KMAX,
       'item_amortised': True, 'pps_ProbH1_lwr_quantiles_mesh': (0.05, 0.25, 0.5, 0.75, 0.95)}
save_trained_model(fit, f"{dir_out}/{file_prefix}_amortised_pps_net.pkl")
print(f"saved ({mins:.1f} min) -> {dir_out}")
