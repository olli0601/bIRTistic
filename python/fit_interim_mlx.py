"""
MLX-Metal port of the batched-SMC interim sampler.

Drives the Case-A SMC PPS over all S future-data samples at once on the Apple
Metal GPU: the MALA move kernel + (S, K, D) particle tensor + the partial-credit
log-likelihood / log-prior all live in MLX (see :mod:`mlx_pcm`); the
orchestration (β bisection, systematic resampling, endpoint-based label) stays
on CPU/numpy/JAX so the heavy per-tempering compute is one fused Metal kernel.

Per-sample β masking — each s's β advances independently (no shared
worst-case schedule), so easy samples finish their tempering in a few steps and
the move is exercised only where it's needed. This is the structural fix for
the over-tempering that hurt the JAX-CPU batched prototype.
"""

from typing import Optional
import time

import arviz as az
import numpy as np
import pandas as pd
from scipy.special import softmax

import mlx.core as mx

from mlx_pcm import prepare_pcm_data_mlx, loglik_pcm_mlx, logprior_pcm_mlx
from fit_interim import (
    _stack_posterior_theta,
    _interim_make_x_stan,
    _interim_make_z_stan_batched,
    _ratio_per_draw_from_params,
)


def _unpack(u, U, P, L, Ld):
    return dict(
        latent_factor_unit=u[:U],
        latent_factor_beta=u[U:U + P],
        skill_thresholds=u[U + P:U + P + L],
        loadings_questions_m1=mx.exp(u[U + P + L:U + P + L + Ld]),
    )


def _make_logpost(x_data, z_data, U, P, L, Ld):
    """Closure over the (shared) x-cohort design and (shared) z-cohort design;
    ``y`` is passed in per call so each future sample can use its own outcomes."""
    sl_l_start = U + P + L

    def logpost(u, y, beta):
        pr = _unpack(u, U, P, L, Ld)
        return (logprior_pcm_mlx(pr)
                + loglik_pcm_mlx(x_data, pr)
                + beta * loglik_pcm_mlx(z_data, pr, y=y)
                + mx.sum(u[sl_l_start:sl_l_start + Ld]))            # log-Jacobian

    def llz_only(u, y):
        pr = _unpack(u, U, P, L, Ld)
        return loglik_pcm_mlx(z_data, pr, y=y)

    return logpost, llz_only


def fit_interim_SMC_batched_GPU_PPS_of_posterior_xz_from_x(
    xi: pd.DataFrame,
    zi: pd.DataFrame,
    dit: pd.DataFrame,
    fitting_method_args: dict,
    draws=None,
    draws_file: Optional[str] = None,
    pps_z_total: int = 10,
    pps_H1_def: float = 0.5,
    pps_ProbH1_thresh: float = 0.89,
    categorical_threshold: int = 3,
    output_file_prefix: Optional[str] = None,
    save_to_file: bool = True,
    verbose: bool = True,
):
    """Batched SMC PPS on Apple Metal via MLX.

    All ``pps_z_total`` future samples are evolved together as one ``(S, K, D)``
    particle tensor on the GPU; each sample has its OWN β_s that advances
    independently by adaptive bisection (per-sample masking), so easy samples
    finish in a handful of temperatures and the move kernel never over-tempers.

    Parameters
    ----------
    fitting_method_args : dict
        ``n_particles`` (K, default 128), ``ess_frac_target`` (default 0.5),
        ``n_move_steps`` (default 20), ``init_step_size`` (default 0.02),
        ``target_accept`` (default 0.574), ``max_temps`` (default 300),
        ``x_formula`` (default ``"~ time - 1"``), ``seed`` (default 123).

    Returns
    -------
    (p_h1_xz, schedule) : tuple of pd.DataFrame
        ``p_h1_xz``: per-(item, sample) probability that the improvement ratio
        exceeds ``pps_H1_def``. ``schedule``: one row per shared tempering step
        with ``n_active`` (samples not yet at β=1), ``beta_min``/``beta_mean`` /
        ``beta_max``, ``ess_frac_min`` / ``mean``, ``accept``, ``step_size``,
        ``move_secs``; plus ``mins_total``, ``n_particles``.
    """
    if draws is None and draws_file is None:
        raise ValueError("Provide either draws or draws_file.")
    if draws is None:
        draws = az.from_zarr(draws_file)
    if 'src_pid' not in zi.columns:
        raise ValueError("zi must carry 'src_pid' (use get_interim_z_from_ypredi).")
    vprint = print if verbose else (lambda *a, **k: None)

    fma = {
        'n_particles': 128, 'ess_frac_target': 0.5, 'n_move_steps': 20,
        'init_step_size': 0.02, 'target_accept': 0.574, 'max_temps': 300,
        'x_formula': "~ time - 1", 'seed': 123,
        **fitting_method_args,
    }
    K = int(fma['n_particles'])
    ess_frac_target = float(fma['ess_frac_target'])
    n_move_steps = int(fma['n_move_steps'])
    target_accept = float(fma['target_accept'])
    max_temps = int(fma['max_temps'])
    x_formula = fma['x_formula']
    seed = int(fma['seed'])
    S = int(pps_z_total)

    theta = _stack_posterior_theta(draws)
    n_draw = int(theta['latent_factor_beta'].shape[0])
    U = int(theta['latent_factor_unit'].shape[1])
    P = int(theta['latent_factor_beta'].shape[1])
    L = int(theta['skill_thresholds'].shape[1])
    Ld = int(theta['loadings_questions_m1'].shape[1])
    D = U + P + L + Ld

    # CPU-side data prep.
    x_stan = _interim_make_x_stan(xi, dit, x_formula)
    z_stan, Y_jax = _interim_make_z_stan_batched(zi, dit, S, x_formula)
    Y_np = np.asarray(Y_jax, dtype=np.int32)            # (S, N_total_z)

    # MLX data (lives on GPU).
    x_data = prepare_pcm_data_mlx(x_stan)
    z_data = prepare_pcm_data_mlx(z_stan)
    Y_mlx = mx.array(Y_np)
    logpost_fn, llz_only_fn = _make_logpost(x_data, z_data, U, P, L, Ld)
    grad_lp = mx.grad(logpost_fn, argnums=0)

    # Double vmap over (S samples x K particles); per-s y and beta are broadcast
    # across the K particles of that sample.
    def per_s_logpost(u_K_D, y_N, beta_scalar):
        return mx.vmap(lambda u: logpost_fn(u, y_N, beta_scalar))(u_K_D)

    def per_s_grad(u_K_D, y_N, beta_scalar):
        return mx.vmap(lambda u: grad_lp(u, y_N, beta_scalar))(u_K_D)

    def per_s_llz(u_K_D, y_N):
        return mx.vmap(lambda u: llz_only_fn(u, y_N))(u_K_D)

    batched_logpost = mx.vmap(per_s_logpost, in_axes=(0, 0, 0))
    batched_grad = mx.vmap(per_s_grad, in_axes=(0, 0, 0))
    batched_llz = mx.vmap(per_s_llz, in_axes=(0, 0))

    # Init particles (K base x-posterior draws broadcast to S samples).
    rng = np.random.default_rng(seed)
    idx0 = rng.choice(n_draw, size=K, replace=False)
    sub = {k: np.asarray(theta[k])[idx0] for k in
           ('latent_factor_unit', 'latent_factor_beta', 'skill_thresholds', 'loadings_questions_m1')}
    u_single = np.concatenate([
        sub['latent_factor_unit'], sub['latent_factor_beta'],
        sub['skill_thresholds'], np.log(sub['loadings_questions_m1']),
    ], axis=1).astype(np.float32)                                          # (K, D)
    u_np = np.broadcast_to(u_single[None], (S, K, D)).copy()               # (S, K, D)
    u = mx.array(u_np)

    # Initial llz per (s, particle).
    llz = np.asarray(batched_llz(u, Y_mlx), dtype=np.float64)              # (S, K)
    llz = np.where(np.isfinite(llz), llz, -1e30)

    betas = np.zeros(S, dtype=np.float64)
    eps = float(fma['init_step_size'])
    mx.random.seed(seed)
    rows = []
    t_all = time.time()
    vprint(f"[SMC-mlx] S={S} samples, K={K} particles per s, D={D}, GPU={mx.default_device()}")

    for t_step in range(max_temps):
        active = betas < 1.0 - 1e-9
        if not active.any():
            break

        # Per-s adaptive Δβ via vectorized numpy bisection.
        cap = 1.0 - betas                                                  # (S,)
        lo = np.zeros(S); hi = cap.copy()

        def ess_at(db_vec):
            w = softmax(db_vec[:, None] * llz, axis=1)
            return 1.0 / (K * np.sum(w ** 2, axis=1))

        ok_hi = ess_at(hi) >= ess_frac_target
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            ok = ess_at(mid) >= ess_frac_target
            lo = np.where(ok, mid, lo); hi = np.where(ok, hi, mid)
        db = np.where(ok_hi, cap, lo)
        db = np.where(active, db, 0.0)                                     # frozen samples don't advance
        betas_new = betas + db
        wn = softmax(db[:, None] * llz, axis=1)
        ess_s = 1.0 / (K * np.sum(wn ** 2, axis=1))

        # Per-s systematic resample (active samples only; frozen reuse current u).
        anc = np.tile(np.arange(K), (S, 1))
        for s in range(S):
            if not active[s]:
                continue
            positions = (rng.random() + np.arange(K)) / K
            anc[s] = np.searchsorted(np.cumsum(wn[s]), positions)
        u_np = np.asarray(u)
        u_np = u_np[np.arange(S)[:, None], anc]                            # (S, K, D)
        u = mx.array(u_np)

        # MALA move (n_move_steps GPU steps), per-particle beta = betas_new[s].
        t0 = time.time()
        betas_mx = mx.array(betas_new.astype(np.float32))
        n_acc = 0.0
        for _ in range(n_move_steps):
            g = batched_grad(u, Y_mlx, betas_mx)
            noise = mx.random.normal(shape=u.shape)
            prop = u + 0.5 * eps ** 2 * g + eps * noise
            gp = batched_grad(prop, Y_mlx, betas_mx)
            lp_u = batched_logpost(u, Y_mlx, betas_mx)
            lp_p = batched_logpost(prop, Y_mlx, betas_mx)
            sqd_fwd = mx.sum((prop - u - 0.5 * eps ** 2 * g) ** 2, axis=-1)
            sqd_bwd = mx.sum((u - prop - 0.5 * eps ** 2 * gp) ** 2, axis=-1)
            log_acc = lp_p - lp_u + (sqd_fwd - sqd_bwd) / (2 * eps ** 2)
            log_unif = mx.log(mx.random.uniform(shape=lp_u.shape))
            accept = log_unif < log_acc
            u = mx.where(accept[..., None], prop, u)
            mx.eval(u)                                          # keep lazy graph shallow
            n_acc += float(mx.mean(accept))
        acc = n_acc / n_move_steps
        move_secs = time.time() - t0

        eps = float(np.clip(eps * np.exp(0.5 * (acc - target_accept)), 1e-4, 1.0))

        # Refresh llz on the moved particles.
        llz = np.asarray(batched_llz(u, Y_mlx), dtype=np.float64)
        llz = np.where(np.isfinite(llz), llz, -1e30)

        rows.append({
            'temp': t_step + 1, 'n_active': int(active.sum()),
            'beta_min': round(float(betas_new.min()), 4),
            'beta_mean': round(float(betas_new.mean()), 4),
            'beta_max': round(float(betas_new.max()), 4),
            'ess_frac_min': round(float(ess_s.min()), 3),
            'ess_frac_mean': round(float(ess_s.mean()), 3),
            'accept': round(acc, 3), 'step_size': round(eps, 5),
            'move_secs': round(move_secs, 2),
        })
        vprint(f"  temp {t_step+1}: active={int(active.sum())}/{S} "
               f"beta [{betas_new.min():.3f}, {betas_new.mean():.3f}, {betas_new.max():.3f}] "
               f"| ESS/K min={ess_s.min():.2f} mean={ess_s.mean():.2f} "
               f"| accept={acc:.2f} | move {move_secs:.1f}s")
        betas = betas_new

    mins_total = (time.time() - t_all) / 60.0
    schedule = pd.DataFrame(rows)
    schedule['n_particles'] = K
    schedule['mins_total'] = round(mins_total, 4)
    vprint(f"\n[SMC-mlx] all S={S} reached beta=1 in {len(schedule)} steps, {mins_total:.2f} min total")

    # Final per-sample labels (uniform post-move weights) via JAX/CPU ratio.
    sl_u, sl_b = slice(0, U), slice(U, U + P)
    sl_s, sl_l = slice(U + P, U + P + L), slice(U + P + L, U + P + L + Ld)
    u_final = np.asarray(u)                                                # (S, K, D)

    import jax.numpy as jnp
    p_rows = []
    for s in range(S):
        particles = {
            'latent_factor_unit': jnp.asarray(u_final[s, :, sl_u]),
            'latent_factor_beta': jnp.asarray(u_final[s, :, sl_b]),
            'skill_thresholds': jnp.asarray(u_final[s, :, sl_s]),
            'loadings_questions_m1': jnp.asarray(np.exp(u_final[s, :, sl_l])),
        }
        ratio = _ratio_per_draw_from_params(particles, x_stan, xi, dit, categorical_threshold)
        sample_p = (
            ratio.assign(ind=(ratio['ratio'] > pps_H1_def).astype(float))
            .groupby(['item_label', 'item_type', 'item_high_label'])['ind']
            .mean().reset_index(name='p_h1_xz')
        )
        sample_p['s'] = s + 1
        p_rows.append(sample_p)

    p_h1_xz = pd.concat(p_rows, ignore_index=True)
    p_h1_xz['pps_H1_def'] = pps_H1_def
    p_h1_xz['pps_ProbH1_thresh'] = pps_ProbH1_thresh
    p_h1_xz['S'] = S

    if save_to_file and output_file_prefix is not None:
        p_h1_xz.to_pickle(f"{output_file_prefix}_p_h1_xz_SMCmlx.pkl")
        schedule.to_csv(f"{output_file_prefix}_smcmlx_schedule.csv", index=False)
        vprint(f"\nSaved MLX-SMC P(H_1 | x, z) to: {output_file_prefix}_p_h1_xz_SMCmlx.pkl")
    return p_h1_xz, schedule
