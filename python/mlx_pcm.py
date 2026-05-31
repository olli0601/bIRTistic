"""
MLX-Metal port of the partial-credit ncats log-likelihood and prior.

Only the move-kernel-relevant pieces of
``src/numpyro/partial_credit_model_ncats_v260413.pyro`` are ported here — the
loglik (`get_log_likelihood_of_partial_credit_model_ncats`) and the prior
(`get_prior_of_partial_credit_model_ncats`). The generated quantities
(``ordered_prob_by_cat_qu_fit``) stay on JAX/CPU because they use scatter; the
SMC move kernel only needs loglik + prior, both of which are pure gather +
matmul + cumsum + logsumexp ops and run cleanly on the Metal GPU.

Used by :mod:`fit_interim_mlx` to drive the batched-vmap SMC sampler on Apple
GPU; the rest of the codebase (CPU/JAX) is unaffected.
"""

import math

import numpy as np
import mlx.core as mx


# FoldedStudentT(df=3, loc=0, scale=1) log-pdf for x > 0:
# log(2) + lgamma((nu+1)/2) - lgamma(nu/2) - 0.5*log(nu*pi) - (nu+1)/2 * log(1 + x^2/nu)
# With nu=3: (nu+1)/2 = 2, nu/2 = 1.5, (nu+1)/2 = 2.
_LGAMMA_2 = 0.0                                              # lgamma(2) = log(1!) = 0
_LGAMMA_1p5 = math.lgamma(1.5)
_FOLDED_T3_CONST = math.log(2.0) + _LGAMMA_2 - _LGAMMA_1p5 - 0.5 * math.log(3.0 * math.pi)


def prepare_pcm_data_mlx(stan_data: dict) -> dict:
    """Convert a partial-credit stan_data dict to a dict of MLX arrays + static
    Python-int shapes used by :func:`loglik_pcm_mlx`. ``y`` is stored as an MLX
    int32 array (1-based, matching the model) and can be overridden per call so
    the SMC move can swap in the future sample's outcomes without rebuilding the
    rest of the design.
    """
    C = int(stan_data['C'])
    Q = np.asarray(stan_data['Q'], dtype=int).tolist()
    K = np.asarray(stan_data['K'], dtype=int).tolist()
    cat_type = np.asarray(stan_data['cat_type'], dtype=int)

    starts, ends = [], []
    s = 0
    for i in range(1, len(cat_type)):
        if cat_type[i] != cat_type[i - 1]:
            starts.append(s); ends.append(i); s = i
    starts.append(s); ends.append(len(cat_type))

    skill_off = [0] * C
    load_off = [0] * C
    for c in range(1, C):
        skill_off[c] = skill_off[c - 1] + Q[c - 1] * (K[c - 1] - 1)
        load_off[c] = load_off[c - 1] + (Q[c - 1] - 1)

    return dict(
        C=C, Q=Q, K=K, starts=starts, ends=ends,
        skill_off=skill_off, load_off=load_off,
        unit_of_obs_0=mx.array(np.asarray(stan_data['unit_of_obs'], dtype=np.int32) - 1),
        question_of_obs_0=mx.array(np.asarray(stan_data['question_of_obs'], dtype=np.int32) - 1),
        X=mx.array(np.asarray(stan_data['X'], dtype=np.float32)),
        y=mx.array(np.asarray(stan_data['y'], dtype=np.int32)),
    )


def _pcm_get_etas(load_c, skill_c, latent_obs_c, q_c0, cumsum_c):
    """PCM logits for one category block; mirrors ``pcm_get_etas`` in the .pyro."""
    inc = load_c[q_c0][:, None] * (latent_obs_c[:, None] - skill_c[q_c0, :])
    eta_cum = inc @ cumsum_c
    return mx.concatenate([mx.zeros((latent_obs_c.shape[0], 1)), eta_cum], axis=1)


def _cat_lpmf(y_c, eta_c):
    """Categorical-logit log-pmf; ``y_c`` is 1-based (matches stan_data['y'])."""
    y_idx = y_c - 1
    logp = eta_c - mx.logsumexp(eta_c, axis=1, keepdims=True)
    return mx.take_along_axis(logp, y_idx[:, None], axis=1)[:, 0]


def loglik_pcm_mlx(data, params, y=None):
    """Sum of per-obs log-likelihood for a single parameter set; vmappable.

    Pass ``y`` to override ``data['y']`` (used by the batched SMC where ``y``
    varies across future samples but the rest of the design is shared).
    """
    if y is None:
        y = data['y']
    lfu = params['latent_factor_unit']
    bp = params['latent_factor_beta']
    skill_all = params['skill_thresholds']
    load_all = params['loadings_questions_m1']
    X = data['X']
    u_obs = data['unit_of_obs_0']
    q_obs = data['question_of_obs_0']

    total = mx.zeros(())
    for c in range(data['C']):
        Qc = data['Q'][c]
        Kc = data['K'][c]
        Km1c = Kc - 1
        s_obs = data['starts'][c]
        e_obs = data['ends'][c]

        y_c = y[s_obs:e_obs]
        u_c = u_obs[s_obs:e_obs]
        q_c0 = q_obs[s_obs:e_obs]

        s_skill = data['skill_off'][c]
        e_skill = s_skill + Qc * Km1c
        skill_c = skill_all[s_skill:e_skill].reshape(Qc, Km1c)

        s_load = data['load_off'][c]
        e_load = s_load + (Qc - 1)
        load_c = mx.concatenate([mx.array([1.0]), load_all[s_load:e_load]])

        cumsum_c = mx.triu(mx.ones((Km1c, Km1c)))
        latent_obs_c = lfu[u_c] + X[s_obs:e_obs, :] @ bp
        eta_c = _pcm_get_etas(load_c, skill_c, latent_obs_c, q_c0, cumsum_c)
        total = total + _cat_lpmf(y_c, eta_c).sum()
    return total


def logprior_pcm_mlx(params):
    """log p(theta) for the partial-credit ncats priors:
    latent_factor_unit ~ N(0, 1/sqrt(1-1/U)) (ZeroSumNormal marginal),
    latent_factor_beta ~ N(0,1), skill_thresholds ~ N(0,3.5),
    loadings_questions_m1 ~ FoldedStudentT(3,0,1). ``U`` inferred from the
    latent_factor_unit length. Single (1-D) parameter set."""
    lfu = params['latent_factor_unit']
    bp = params['latent_factor_beta']
    skill = params['skill_thresholds']
    load = params['loadings_questions_m1']

    U = lfu.shape[-1]
    s2z2 = 1.0 / (1.0 - 1.0 / U)                                # variance = 1/(1-1/U)
    inv_2s2 = 0.5 / s2z2
    log_norm_lfu = -0.5 * math.log(2.0 * math.pi * s2z2)
    log_norm_beta = -0.5 * math.log(2.0 * math.pi)
    log_norm_skill = -0.5 * math.log(2.0 * math.pi * 3.5 * 3.5)

    lp = mx.sum(-inv_2s2 * lfu * lfu + log_norm_lfu)
    lp = lp + mx.sum(-0.5 * bp * bp + log_norm_beta)
    lp = lp + mx.sum(-0.5 / (3.5 * 3.5) * skill * skill + log_norm_skill)
    lp = lp + mx.sum(_FOLDED_T3_CONST - 2.0 * mx.log(1.0 + load * load / 3.0))
    return lp
