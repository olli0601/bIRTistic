"""
Interim-analysis helpers for predictive probability of success (PPS).

- ``get_interim_z``: build the 'missing' future data block z from the interim
  cohort xi, the final cohort xf, and posterior-predictive draws ypred.
- ``fit_interim_posterior_xz_with_nested_monte_carlo``: Monte-Carlo estimate
  of p(H_1 | x, z_s) across S hypothetical future datasets, refitting the
  model to (xi + z_s) for each s.
"""

from typing import Optional
import time

import arviz as az
import numpy as np
import pandas as pd
import jax
import jax.numpy as jnp
import jax.scipy.stats as jstats
from jax import random
from scipy.special import softmax
import statsmodels.api as sm

# IRT-specific endpoint helpers + cohort builder live on IRTModel
# (python/model_irt.py). Algorithm callers route through a Model instance.


def _fit_interim_posterior_xz_with_nested_monte_carlo_one_sample(
    s_idx,
    xi,
    zi,
    dit,
    model_cls,
    interim_method_args,
    fit_method_args,
):
    """
    Refit the model on (xi + z_{s_idx}) and return the per-item p(H_1 | x, z_s)
    frame for a single Monte-Carlo sample. Module-level so it is picklable for
    process-based parallelism. Inherits the two arg dicts from
    :func:`fit_interim_posterior_xz_with_nested_monte_carlo`.
    """
    ima = interim_method_args
    fma = dict(fit_method_args)
    x_formula = fma.pop('x_formula')   # KeyError if caller omits

    vprint = print if ima['verbose'] else (lambda *args, **kwargs: None)
    s_label = s_idx + 1
    vprint(f"\n--- PPS sample {s_label} ---")

    # This sample's predicted outcomes become y_stan for the new participants;
    # concat onto the interim data and rebuild the cohort via the model-class
    # hook (IRT rebuilds pid/oid/oidt; Binomial coalesces y_stan -> y).
    ycol = f'ypred_{s_idx}'
    tmp = zi.drop(
        columns=['y', 'y_stan'] + [c for c in zi.columns if c.startswith('ypred_') and c != ycol],
        errors='ignore',
    ).rename(columns={ycol: 'y_stan'})

    xzi = model_cls.get_interim_data_x(
        pd.concat([xi, tmp], ignore_index=True),
    )

    s_model = model_cls(dit=dit, dcati=xzi, x_formula=x_formula,
                        seed=ima['seed'] + s_label)
    # SVI-only post-fit IRT analyses are off by default in the driver.
    # HMC fit signatures absorb the with_* kwargs via **method_kwargs.
    out_prefix = ima['output_file_prefix']
    fit = getattr(s_model, ima['fit_method'])(
        output_file_prefix=(f"{out_prefix}_s{s_label}" if out_prefix else None),
        save_to_file=False,
        resume=False,
        with_core_analyses=False,
        with_additional_analyses=False,
        verbose=ima['verbose'],
        **fma,
    )

    # Per-item P(H_1 | x, z_s). Each model class supplies its own
    # ratio definition (IRT: 1 - p_end/p_base for lower_is_better;
    # Binomial: 1 - p/p_0) and shares the same get_p_h1 aggregator
    # (right-tail r > pps_H1_def).
    ratio_xz = s_model.get_endpoints_per_draw(
        draws=fit['draws'], endpoint_type='items',
    )
    sample_p = (
        s_model.get_p_h1(ratio_xz, pps_H1_def=ima['pps_H1_def'])
        .rename(columns={'p_h1': 'p_h1_xz'})
    )
    sample_p['s'] = s_label
    return sample_p


def fit_interim_posterior_xz_with_nested_monte_carlo(
    model,
    zi: pd.DataFrame,
    interim_method_args: dict,
    fit_method_args: dict,
) -> pd.DataFrame:
    """
    Monte-Carlo estimate of p(H_1 | x, z_s) over S hypothetical future datasets.

    For each sample ``s`` the ``ypred_s`` column of ``zi`` becomes the outcome
    ``y_stan`` for the new participants; this is concatenated onto
    ``model.dcati`` and rebuilt via ``type(model).get_interim_data_x``. A fresh
    model is refit in memory (``save_to_file=False``) and the per-item
    probability that the directional improvement ratio exceeds ``pps_H1_def``
    is recorded.

    Parameters
    ----------
    model : :class:`model.Model`
        Carries ``dcati`` (= xi), ``dit``, ``x_formula`` and the subclass used
        to build per-sample refits via ``type(model)(...)``.
    zi : pd.DataFrame
        Future-data block (must hold ``ypred_0 .. ypred_{S-1}`` columns).
    interim_method_args : dict
        PPS-loop tuning. ALL keys are required (no defaults):
        ``pps_z_total``, ``pps_H1_def``, ``pps_ProbH1_thresh``,
        ``fit_method``, ``seed``, ``save_to_file``, ``verbose``,
        ``cpu_n``, ``output_file_prefix``.
    fit_method_args : dict
        Kwargs forwarded to ``getattr(model, fit_method)`` in each per-s refit
        (e.g. ``algorithm``, ``lr``, ``num_steps``, ``output_samples``,
        ``chains``, ``iter_warmup``, ``iter_sampling``). Must include
        ``x_formula`` (popped before forwarding); no implicit default.

    Returns
    -------
    pd.DataFrame
        p_h1_xz: long form, one row per (item, sample) with columns
        ``item_label``, ``item_type``, ``item_high_label``, ``p_h1_xz``, ``s``,
        plus the recorded ``pps_H1_def``, ``pps_ProbH1_thresh`` and ``S``.
    """
    ima = interim_method_args
    fma = fit_method_args

    vprint = print if ima['verbose'] else (lambda *args, **kwargs: None)
    work = dict(
        xi=model.dcati, zi=zi, dit=model.dit,
        model_cls=type(model),
        interim_method_args=ima, fit_method_args=fma,
    )

    if ima['cpu_n'] == 1:
        rows = [
            _fit_interim_posterior_xz_with_nested_monte_carlo_one_sample(s_idx, **work)
            for s_idx in range(ima['pps_z_total'])
        ]
    else:
        import multiprocessing as _mp
        from concurrent.futures import ProcessPoolExecutor
        vprint(f"Running {ima['pps_z_total']} PPS refits over {ima['cpu_n']} worker processes...")
        ctx = _mp.get_context('spawn')
        with ProcessPoolExecutor(max_workers=ima['cpu_n'], mp_context=ctx) as ex:
            futures = [
                ex.submit(
                    _fit_interim_posterior_xz_with_nested_monte_carlo_one_sample,
                    s_idx, **work,
                )
                for s_idx in range(ima['pps_z_total'])
            ]
            rows = [f.result() for f in futures]

    p_h1_xz = pd.concat(rows, ignore_index=True)
    p_h1_xz['pps_H1_def'] = ima['pps_H1_def']
    p_h1_xz['pps_ProbH1_thresh'] = ima['pps_ProbH1_thresh']
    p_h1_xz['S'] = ima['pps_z_total']

    if ima['save_to_file'] and ima['output_file_prefix'] is not None:
        pkl_path = f"{ima['output_file_prefix']}_p_h1_xz.pkl"
        p_h1_xz.to_pickle(pkl_path)
        vprint(f"\nSaved P(H_1 | x, z) samples to: {pkl_path}")

    return p_h1_xz


def fit_interim_posterior_xz_from_z_with_IS_reweight(
    model,
    zi: pd.DataFrame,
    interim_method_args: dict,
    draws=None,
    draws_file: Optional[str] = None,
):
    """
    Importance-sampling estimate of p(H_1 | x, z_s) for fixed x. Reweights
    the x-posterior draws by ``softmax_k log p(z_s | theta_k)``; no per-s refit.

    Parameters
    ----------
    model : :class:`model.Model`
        Carries ``dit`` / ``dcati`` / ``x_formula`` and exposes
        ``eval_loglik`` / ``get_endpoints_per_draw`` /
        ``get_stacked_posterior`` / ``make_stan_data``.
    zi : pd.DataFrame
        Future-data block from
        :meth:`Model.get_interim_z_from_ypredi`; must carry ``src_pid`` +
        ``ypred_0 .. ypred_{S-1}`` columns.
    interim_method_args : dict
        ALL keys required (no defaults): ``pps_z_total``, ``pps_H1_def``,
        ``pps_ProbH1_thresh``, ``output_file_prefix``,
        ``save_to_file``, ``verbose``. No ``fit_method_args`` -- IS does
        not refit.
    draws, draws_file
        Supply one. The x-fit posterior arviz object or its zarr.

    Returns
    -------
    (p_h1_xz, is_perf) : tuple of pd.DataFrame
    """
    if draws is None and draws_file is None:
        raise ValueError("Provide either draws or draws_file.")
    if draws is None:
        draws = az.from_zarr(draws_file)
    if 'src_pid' not in zi.columns:
        raise ValueError("zi must carry 'src_pid' (use get_interim_z_from_ypredi).")

    ima = interim_method_args
    vprint = print if ima['verbose'] else (lambda *args, **kwargs: None)

    # Per-draw H_1 ratios from the x-posterior; 'draw' aligns with theta_k
    # below. param_name + categorical_threshold come from model defaults.
    x_ratio = model.get_endpoints_per_draw(
        draws=draws, endpoint_type='items',
    )
    x_ratio = x_ratio.assign(
        ind=(x_ratio['ratio'] > ima['pps_H1_def']).astype(float),
    )

    theta = model.get_stacked_posterior(draws)
    n_draw = int(next(iter(theta.values())).shape[0])

    p_h1_xz = []
    perf_rows = []
    for s_idx in range(ima['pps_z_total']):
        s_label = s_idx + 1
        vprint(f"\n--- IS sample {s_label}/{ima['pps_z_total']} ---")

        # z_s stan_data via the model hook (IRT does sort + oid/oidt rebuild
        # + make_stan_data; Binomial overrides per its data shape).
        z_stan = model.make_stan_data_from_zi(zi, s_idx)

        # log p(z_s | theta_k) for every draw via vmap over the stacked params.
        loglik = jax.vmap(lambda p: model.eval_loglik(z_stan, p))(theta)
        w = np.asarray(jax.nn.softmax(loglik.sum(axis=1)))

        sum_w2 = float(np.sum(w ** 2))
        perf_rows.append({
            's': s_label,
            'N': n_draw,
            'ess': float(1.0 / sum_w2),
            'ess_over_n': float(1.0 / (n_draw * sum_w2)),
            'ew2': float(np.mean(w ** 2)),
        })

        wdf = pd.DataFrame({'draw': np.arange(n_draw), 'w': w})
        m = x_ratio.merge(wdf, on='draw')
        sample_p = (
            m.groupby(['item_label', 'item_type', 'item_high_label'])
            .apply(lambda g: float((g['w'] * g['ind']).sum() / g['w'].sum()))
            .reset_index(name='p_h1_xz')
        )
        sample_p['s'] = s_label
        p_h1_xz.append(sample_p)

    p_h1_xz = pd.concat(p_h1_xz, ignore_index=True)
    p_h1_xz['pps_H1_def'] = ima['pps_H1_def']
    p_h1_xz['pps_ProbH1_thresh'] = ima['pps_ProbH1_thresh']
    p_h1_xz['S'] = ima['pps_z_total']
    is_perf = pd.DataFrame(perf_rows)

    if ima['save_to_file'] and ima['output_file_prefix'] is not None:
        pkl_path = f"{ima['output_file_prefix']}_p_h1_xz_IS.pkl"
        p_h1_xz.to_pickle(pkl_path)
        is_perf.to_csv(f"{ima['output_file_prefix']}_is_perf.csv", index=False)
        vprint(f"\nSaved IS P(H_1 | x, z) samples to: {pkl_path}")

    return p_h1_xz, is_perf


# =============================================================================
# Mitigations for IS weight degeneracy at early interims (large m = N_full - n_x)
# =============================================================================
# At early interims p(theta | x) is diffuse and the future block z is highly
# informative, so the self-normalised weights of
# :func:`fit_interim_importance_sampling_of_posterior_xz_from_x` collapse onto a
# single draw (ESS/N -> 1/N). The two functions below implement the standard
# remedies, both specialised to the partial credit ncats model:
#   - moment-matching IS (no refit; cannot rescue a fully-collapsed base),
#   - SMC with resample-move (relocates particles; the only method that crosses
#     a large x -> x,z gap, at the cost of many tempering steps).


def fit_interim_posterior_xz_with_IS_moment_matching(
    model,
    zi: pd.DataFrame,
    interim_method_args: dict,
    draws=None,
    draws_file: Optional[str] = None,
):
    """
    Mean-match IS for p(H_1 | x, z) (Paananen 2021).

    Reads ``eval_loglik`` / ``eval_log_prior`` / ``positive_params`` /
    ``x_formula`` from ``model``.

    Parameters
    ----------
    model : :class:`model.Model`
    zi : pd.DataFrame
        Future-data block (must carry ``src_pid`` + ``ypred_*`` columns).
    interim_method_args : dict
        ALL keys required (no defaults): ``pps_z_total``, ``pps_H1_def``,
        ``pps_ProbH1_thresh``, ``output_file_prefix``,
        ``save_to_file``, ``verbose``.
    draws, draws_file
        Supply one. The x-fit posterior or its zarr.

    Returns
    -------
    (p_h1_xz, mm) : tuple of pd.DataFrame
    """
    if draws is None and draws_file is None:
        raise ValueError("Provide either draws or draws_file.")
    if draws is None:
        draws = az.from_zarr(draws_file)
    if 'src_pid' not in zi.columns:
        raise ValueError("zi must carry 'src_pid' (use get_interim_z_from_ypredi).")
    ima = interim_method_args
    vprint = print if ima['verbose'] else (lambda *args, **kwargs: None)

    eval_loglik = model.eval_loglik
    eval_log_prior = model.eval_log_prior
    positive_params = model.positive_params

    theta = model.get_stacked_posterior(draws)
    n_draw = int(next(iter(theta.values())).shape[0])
    x_stan = model.make_stan_data_from_xi()

    # Build the match-space matrix Theta_u (K, D) and its block layout.
    layout, blocks = [], []
    cursor = 0
    for name in theta.keys():
        arr = np.asarray(theta[name])
        is_pos = name in positive_params
        block = np.log(arr) if is_pos else arr
        d = block.shape[1]
        layout.append((name, cursor, cursor + d, is_pos))
        cursor += d
        blocks.append(block)
    Theta_u = np.concatenate(blocks, axis=1)
    mu_p = Theta_u.mean(axis=0)
    sd_p = Theta_u.std(axis=0) + 1e-8
    pos_cols = np.concatenate(
        [np.arange(a, b) for (_, a, b, is_pos) in layout if is_pos]
    ) if any(is_pos for *_, is_pos in layout) else np.array([], dtype=int)

    def _u_to_params(u):
        return {
            name: (jnp.exp(u[a:b]) if is_pos else u[a:b])
            for (name, a, b, is_pos) in layout
        }

    mu_p_j, sd_p_j = jnp.asarray(mu_p), jnp.asarray(sd_p)

    # Compile logw and llz_exact ONCE per interim: only z's `y` changes across
    # s_idx (the rest of z_stan -- C/N/Q/K/cat_type/X/unit_of_obs/question_of_obs
    # -- is invariant within an interim). Passing z_y as the only traced arg and
    # closing over the rest avoids the 200x JIT-recompile that previously
    # dominated MM wall time at ~10 min/interim.
    z_stan_template = model.make_stan_data_from_zi(zi, 0)

    def _logw_inner(u, z_y):
        z_stan_local = {**z_stan_template, 'y': z_y}
        pr = _u_to_params(u)
        jac = u[pos_cols].sum() if pos_cols.size else 0.0
        lt = (eval_log_prior(pr)
              + eval_loglik(x_stan, pr).sum()
              + eval_loglik(z_stan_local, pr).sum()
              + jac)
        return lt - jstats.norm.logpdf(u, mu_p_j, sd_p_j).sum()

    logw_batch = jax.jit(jax.vmap(_logw_inner, in_axes=(0, None)))

    def _llz_inner(p, z_y):
        z_stan_local = {**z_stan_template, 'y': z_y}
        return eval_loglik(z_stan_local, p).sum()

    llz_batch = jax.jit(jax.vmap(_llz_inner, in_axes=(0, None)))

    def _u_batch_to_params(u_batch):
        return {
            name: (jnp.exp(u_batch[:, a:b]) if is_pos else u_batch[:, a:b])
            for (name, a, b, is_pos) in layout
        }

    Theta_u_j = jnp.asarray(Theta_u)
    rows = []
    p_h1_xz = []
    t_all = time.time()
    for s_idx in range(ima['pps_z_total']):
        z_stan = model.make_stan_data_from_zi(zi, s_idx)
        z_y = jnp.asarray(z_stan['y']).astype(jnp.int32)

        logw0 = np.asarray(logw_batch(Theta_u_j, z_y))
        w0 = softmax(logw0)
        ess0 = float(1.0 / (n_draw * np.sum(w0 ** 2)))
        llz_exact = np.asarray(llz_batch(theta, z_y))
        corr = float(np.corrcoef(logw0, llz_exact)[0, 1])

        mu_w = (w0[:, None] * Theta_u).sum(axis=0)
        Theta_u_star = Theta_u + (mu_w - mu_p)[None, :]
        logw1 = np.asarray(logw_batch(jnp.asarray(Theta_u_star), z_y))
        w1 = softmax(logw1)
        ess1 = float(1.0 / (n_draw * np.sum(w1 ** 2)))

        rows.append({
            's': s_idx + 1, 'N': n_draw,
            'ess_over_n_base': ess0, 'ess_over_n_meanmatch': ess1,
            'basew_vs_exact_corr': corr,
        })
        vprint(f"  [MM] s={s_idx+1}: ESS/N base={ess0:.4f} -> mean-match={ess1:.4f} "
               f"(base-weight vs exact corr={corr:.2f})")

        # p(H_1 | x, z_s): improvement ratio re-evaluated on the SHIFTED draws,
        # then averaged with the mean-match weights w1. Delegated to the model
        # so the param-name + posterior-shape stays model-aware.
        theta_batch = _u_batch_to_params(jnp.asarray(Theta_u_star))
        ratio = model.get_endpoints_per_draw_from_theta_batch(
            theta_batch, x_stan, endpoint_type='items',
        )
        ratio = ratio.assign(
            ind=(ratio['ratio'] > ima['pps_H1_def']).astype(float),
        )
        wdf = pd.DataFrame({'draw': np.arange(n_draw), 'w': w1})
        m = ratio.merge(wdf, on='draw')
        sample_p = (
            m.groupby(['item_label', 'item_type', 'item_high_label'])
            .apply(lambda g: float((g['w'] * g['ind']).sum() / g['w'].sum()))
            .reset_index(name='p_h1_xz')
        )
        sample_p['s'] = s_idx + 1
        p_h1_xz.append(sample_p)

    mm = pd.DataFrame(rows)
    mm['mins'] = round((time.time() - t_all) / 60.0, 4)
    p_h1_xz = pd.concat(p_h1_xz, ignore_index=True)
    p_h1_xz['pps_H1_def'] = ima['pps_H1_def']
    p_h1_xz['pps_ProbH1_thresh'] = ima['pps_ProbH1_thresh']
    p_h1_xz['S'] = ima['pps_z_total']

    out_prefix = ima['output_file_prefix']
    if ima['save_to_file'] and out_prefix is not None:
        mm.to_csv(f"{out_prefix}_momentmatch.csv", index=False)
        p_h1_xz.to_pickle(f"{out_prefix}_p_h1_xz_MM.pkl")
        vprint(f"\nSaved moment-matching diagnostics to: {out_prefix}_momentmatch.csv")
    return p_h1_xz, mm


def fit_interim_posterior_xz_from_x_with_SMC_resampling(
    model,
    zi: pd.DataFrame,
    interim_method_args: dict,
    fit_method_args: dict,
    draws=None,
    draws_file: Optional[str] = None,
):
    """
    SMC sampler with resample-move for p(theta | x, z_s).

    Reads ``eval_loglik`` / ``eval_loglik_annealed`` / ``eval_log_prior`` /
    ``get_stacked_posterior`` from ``model``. The PCM-specific param layout
    (latent_factor_unit / latent_factor_beta / skill_thresholds /
    loadings_questions_m1) is currently still hard-coded here.

    Parameters
    ----------
    model : :class:`model.Model`
    zi : pd.DataFrame
    interim_method_args : dict
        ALL keys required (no defaults): ``output_file_prefix``,
        ``save_to_file``, ``verbose``.
    fit_method_args : dict
        ALL keys required: ``s_idx``, ``n_particles``, ``ess_frac_target``,
        ``n_move_steps``, ``init_step_size``, ``target_accept``,
        ``max_temps``, ``seed``.
    draws, draws_file
        Supply one.

    Returns
    -------
    (schedule, particles) : tuple
    """
    if draws is None and draws_file is None:
        raise ValueError("Provide either draws or draws_file.")
    if draws is None:
        draws = az.from_zarr(draws_file)
    if 'src_pid' not in zi.columns:
        raise ValueError("zi must carry 'src_pid' (use get_interim_z_from_ypredi).")
    ima = interim_method_args
    fma = fit_method_args
    vprint = print if ima['verbose'] else (lambda *args, **kwargs: None)

    s_idx = fma['s_idx']
    n_particles = fma['n_particles']
    ess_frac_target = fma['ess_frac_target']
    n_move_steps = fma['n_move_steps']
    init_step_size = fma['init_step_size']
    target_accept = fma['target_accept']
    max_temps = fma['max_temps']
    seed = fma['seed']

    eval_loglik = model.eval_loglik
    eval_loglik_annealed = model.eval_loglik_annealed
    eval_log_prior = model.eval_log_prior

    theta = model.get_stacked_posterior(draws)
    n_draw = int(theta['latent_factor_beta'].shape[0])
    U = int(theta['latent_factor_unit'].shape[1])
    P = int(theta['latent_factor_beta'].shape[1])
    L = int(theta['skill_thresholds'].shape[1])
    Ld = int(theta['loadings_questions_m1'].shape[1])

    x_stan = model.make_stan_data_from_xi()
    z_stan = model.make_stan_data_from_zi(zi, s_idx)
    pcm_keys = ('latent_factor_unit', 'latent_factor_beta',
                'skill_thresholds', 'loadings_questions_m1')

    # Match-space layout u = [latent_factor_unit | latent_factor_beta |
    # skill_thresholds | log(loadings_questions_m1)]. Loadings are the only
    # positive param, log-transformed so MALA moves in unconstrained space.
    sl_u, sl_b = slice(0, U), slice(U, U + P)
    sl_s, sl_l = slice(U + P, U + P + L), slice(U + P + L, U + P + L + Ld)

    def _u_to_params(u):           # single particle (1-D u)
        return {
            'latent_factor_unit': u[sl_u],
            'latent_factor_beta': u[sl_b],
            'skill_thresholds': u[sl_s],
            'loadings_questions_m1': jnp.exp(u[sl_l]),
        }

    def _logpost(u, beta):
        pr = _u_to_params(u)
        return (eval_log_prior(pr)
                + eval_loglik(x_stan, pr).sum()
                + eval_loglik_annealed(z_stan, {**pr, 'temperature': beta}).sum()
                + u[sl_l].sum())          # Jacobian of the log-transform

    _grad = jax.grad(_logpost, argnums=0)

    def _mala_sweep(u, beta, eps, key):
        """n_move_steps MALA steps, vectorised over particles; beta/eps traced."""
        def step(carry, k):
            u, n_acc = carry
            g = jax.vmap(lambda ui: _grad(ui, beta))(u)
            k1, k2 = random.split(k)
            prop = u + 0.5 * eps ** 2 * g + eps * random.normal(k1, u.shape)
            gp = jax.vmap(lambda ui: _grad(ui, beta))(prop)
            lp_u = jax.vmap(lambda ui: _logpost(ui, beta))(u)
            lp_p = jax.vmap(lambda ui: _logpost(ui, beta))(prop)
            logq_fwd = -jnp.sum((prop - u - 0.5 * eps ** 2 * g) ** 2, axis=1) / (2 * eps ** 2)
            logq_bwd = -jnp.sum((u - prop - 0.5 * eps ** 2 * gp) ** 2, axis=1) / (2 * eps ** 2)
            log_acc = lp_p - lp_u + logq_bwd - logq_fwd
            accept = jnp.log(random.uniform(k2, (u.shape[0],))) < log_acc
            u = jnp.where(accept[:, None], prop, u)
            return (u, n_acc + accept.mean()), None
        (u, n_acc), _ = jax.lax.scan(step, (u, 0.0), random.split(key, n_move_steps))
        return u, n_acc / n_move_steps

    _mala_jit = jax.jit(_mala_sweep)

    def _llz_at(u):
        v = np.asarray(jax.vmap(lambda ui: eval_loglik(z_stan, _u_to_params(ui)).sum())(u))
        return np.where(np.isfinite(v), v, -1e30)

    def _systematic_resample(w, key):
        positions = (random.uniform(key) + np.arange(len(w))) / len(w)
        return np.searchsorted(np.cumsum(w), positions)

    # Init particles (subsample x-posterior draws) into match space.
    rng = np.random.default_rng(seed)
    idx0 = rng.choice(n_draw, size=n_particles, replace=False)
    sub = {k: np.asarray(theta[k])[idx0] for k in pcm_keys}
    u = jnp.asarray(np.concatenate([
        sub['latent_factor_unit'], sub['latent_factor_beta'],
        sub['skill_thresholds'], np.log(sub['loadings_questions_m1']),
    ], axis=1))
    llz = _llz_at(u)

    key = random.PRNGKey(seed)
    eps = float(init_step_size)
    beta = 0.0
    rows = []
    t_all = time.time()
    vprint(f"[SMC] resample-move (MALA) on s={s_idx+1}, K={n_particles} particles")
    for t in range(max_temps):
        def ess_frac(db):
            wn = softmax(db * llz)
            return 1.0 / (n_particles * np.sum(wn ** 2))
        lo, hi = 0.0, 1.0 - beta
        if ess_frac(hi) >= ess_frac_target:
            db = hi
        else:
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                if ess_frac(mid) >= ess_frac_target:
                    lo = mid
                else:
                    hi = mid
            db = lo
        beta_new = beta + db
        wn = softmax(db * llz)
        ess_temper = float(1.0 / (n_particles * np.sum(wn ** 2)))

        key, ksub = random.split(key)
        u = u[_systematic_resample(wn, ksub)]

        t0 = time.time()
        key, kmove = random.split(key)
        u, acc = _mala_jit(u, jnp.asarray(beta_new), jnp.asarray(eps), kmove)
        acc = float(acc)
        move_secs = time.time() - t0
        # Adapt the MALA step size toward target_accept (MALA optimum ~0.574).
        eps = float(np.clip(eps * np.exp(0.5 * (acc - target_accept)), 1e-4, 1.0))
        llz = _llz_at(u)

        rows.append({
            'temp': t + 1, 'beta': round(beta_new, 5), 'd_beta': round(db, 5),
            'ess_frac_temper': round(ess_temper, 3), 'accept': round(acc, 3),
            'step_size': round(eps, 5), 'move_secs': round(move_secs, 2),
        })
        vprint(f"  temp {t+1}: beta {beta:.4f}->{beta_new:.4f} (dbeta={db:.4f}) "
               f"| tempering ESS/K={ess_temper:.2f} | accept={acc:.2f} | move {move_secs:.1f}s")
        beta = beta_new
        if beta >= 1.0 - 1e-9:
            break

    mins_total = (time.time() - t_all) / 60.0
    schedule = pd.DataFrame(rows)
    schedule['n_particles'] = n_particles
    schedule['mins_total'] = round(mins_total, 4)
    u = jnp.asarray(u)
    particles = {
        'latent_factor_unit': u[:, sl_u],
        'latent_factor_beta': u[:, sl_b],
        'skill_thresholds': u[:, sl_s],
        'loadings_questions_m1': jnp.exp(u[:, sl_l]),
    }
    vprint(f"\n[SMC] reached beta={beta:.4f} in {len(schedule)} steps, {mins_total:.2f} min total")
    out_prefix = ima['output_file_prefix']
    if ima['save_to_file'] and out_prefix is not None:
        schedule.to_csv(f"{out_prefix}_smc_schedule.csv", index=False)
        vprint(f"Saved SMC schedule to: {out_prefix}_smc_schedule.csv")
    return schedule, particles


def _fit_interim_SMC_one_sample(
    s_idx,
    xi,
    zi,
    dit,
    model_cls,
    draws_file,
    interim_method_args,
    fit_method_args,
):
    """Run the SMC sampler for a single future-data sample and return its
    per-item p(H_1 | x, z_s) plus a one-row schedule summary. Module-level so
    it is picklable for process-based parallelism. The post-move particles are
    uniformly weighted, so p(H_1 | x, z_s) is the fraction with
    ``ratio > pps_H1_def``. Inherits the two arg dicts from
    :func:`fit_interim_SMC_PPS`.
    """
    ima = interim_method_args
    fma = {**fit_method_args, 's_idx': s_idx}
    x_formula = fma.pop('x_formula')   # KeyError if caller omits

    model = model_cls(dit=dit, dcati=xi, x_formula=x_formula)
    schedule, particles = fit_interim_posterior_xz_from_x_with_SMC_resampling(
        model=model, zi=zi, draws_file=draws_file,
        interim_method_args={
            'verbose': ima['verbose'],
            'save_to_file': False,
            'output_file_prefix': None,
        },
        fit_method_args=fma,
    )
    x_stan = model.make_stan_data_from_xi()
    ratio = model.get_endpoints_per_draw_from_theta_batch(
        particles, x_stan, endpoint_type='items',
    )
    sample_p = (
        ratio.assign(ind=(ratio['ratio'] > ima['pps_H1_def']).astype(float))
        .groupby(['item_label', 'item_type', 'item_high_label'])['ind']
        .mean().reset_index(name='p_h1_xz')
    )
    sample_p['s'] = s_idx + 1
    summary = {
        's': s_idx + 1,
        'n_particles': int(schedule['n_particles'].iloc[0]),
        'n_steps': int(len(schedule)),
        'beta_final': float(schedule['beta'].iloc[-1]),
        'mins': float(schedule['mins_total'].iloc[0]),
    }
    return sample_p, summary


def _fit_interim_posterior_xz_from_x_with_SMC_resampling_per_sample(
    model,
    zi: pd.DataFrame,
    draws_file: str,
    interim_method_args: dict,
    fit_method_args: dict,
):
    """
    SMC resample-move PPS over S future datasets at a fixed x (Case A).

    Runs :func:`fit_interim_posterior_xz_from_x_with_SMC_resampling` for each
    future sample z_s and scores p(H_1 | x, z_s) as the fraction of the moved
    particles whose per-item improvement ratio exceeds ``pps_H1_def``. The S
    SMC runs are independent and parallelised over ``cpu_n`` ``spawn``
    workers (each re-imports JAX and reloads ``draws_file``).

    Parameters
    ----------
    model : :class:`model.Model`
    zi : pd.DataFrame
        Future-data block (must carry ``src_pid`` + ``ypred_*``).
    draws_file : str
        Path to the x-fit zarr (passed instead of in-memory draws so workers
        can reload it).
    interim_method_args : dict
        ALL keys required (no defaults): ``pps_z_total``, ``pps_H1_def``,
        ``pps_ProbH1_thresh``, ``cpu_n``, ``output_file_prefix``,
        ``save_to_file``, ``verbose``.
    fit_method_args : dict
        SMC tuning forwarded to
        :func:`fit_interim_posterior_xz_from_x_with_SMC_resampling`.
        ALL keys required: ``n_particles``, ``ess_frac_target``,
        ``n_move_steps``, ``init_step_size``, ``target_accept``,
        ``max_temps``, ``seed``, ``x_formula``. ``s_idx`` is set per sample.

    Returns
    -------
    (p_h1_xz, smc_summary) : tuple of pd.DataFrame
    """
    ima = interim_method_args
    fma = fit_method_args

    vprint = print if ima['verbose'] else (lambda *args, **kwargs: None)
    work = dict(
        xi=model.dcati, zi=zi, dit=model.dit,
        model_cls=type(model),
        draws_file=draws_file,
        interim_method_args=ima, fit_method_args=fma,
    )
    if ima['cpu_n'] == 1:
        results = [_fit_interim_SMC_one_sample(s, **work) for s in range(ima['pps_z_total'])]
    else:
        import multiprocessing as _mp
        from concurrent.futures import ProcessPoolExecutor
        vprint(f"Running {ima['pps_z_total']} SMC samples over {ima['cpu_n']} worker processes...")
        ctx = _mp.get_context('spawn')
        with ProcessPoolExecutor(max_workers=ima['cpu_n'], mp_context=ctx) as ex:
            futures = [ex.submit(_fit_interim_SMC_one_sample, s, **work) for s in range(ima['pps_z_total'])]
            results = [f.result() for f in futures]

    p_h1_xz = pd.concat([r[0] for r in results], ignore_index=True)
    p_h1_xz['pps_H1_def'] = ima['pps_H1_def']
    p_h1_xz['pps_ProbH1_thresh'] = ima['pps_ProbH1_thresh']
    p_h1_xz['S'] = ima['pps_z_total']
    smc_summary = pd.DataFrame([r[1] for r in results])

    out_prefix = ima['output_file_prefix']
    if ima['save_to_file'] and out_prefix is not None:
        p_h1_xz.to_pickle(f"{out_prefix}_p_h1_xz_SMC.pkl")
        smc_summary.to_csv(f"{out_prefix}_smc_summary.csv", index=False)
        vprint(f"\nSaved SMC P(H_1 | x, z) samples to: {out_prefix}_p_h1_xz_SMC.pkl")
    return p_h1_xz, smc_summary


# =============================================================================
# Strong-Oakley regression-based label estimator (Case A)
# =============================================================================
# 1. Callers build the per-item training set ``wa`` inline via
#    ``model.get_interim_z_from_ypredi`` (keep_order=True) + ``model.get_w``
#    + ``model.get_endpoints_per_draw`` + ``model.get_p_h1``.
# 2. ``fit_interim_regress_H1x_on_wz`` fits a per-item binomial GLM
#    ``pps_H1_x ~ w_ratio`` and predicts ``p_h1_xz`` = pi_hat(W(z^(s))) for
#    every draw, plus per-item Pearson rho and Gaussian-GLM R^2 diagnostics.
# 3. ``fit_interim_regress_H1x_on_wz_per_item_summary`` is the per-item helper
#    used both inside (2) and by the diagnostic plot scripts.


def fit_interim_regress_H1x_on_wz_per_item_summary(g: pd.DataFrame) -> pd.Series:
    """Per-item diagnostic: Pearson rho on ``(w_ratio, pps_ratio_x)``,
    Gaussian-GLM R^2 (= 1 - deviance / null_deviance, exact for the identity-
    link Gaussian fit), and panel-relative positions for an annotation plot.
    Returns NaN when there are too few finite rows to fit."""
    x = g['w_ratio'].to_numpy()
    y = g['pps_ratio_x'].to_numpy()
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 3:
        return pd.Series({'rho': np.nan, 'r2': np.nan,
                          'x_left': np.nan, 'x_right': np.nan, 'y_top': np.nan})
    rho = float(np.corrcoef(x, y)[0, 1])
    res = sm.GLM(y, sm.add_constant(x), family=sm.families.Gaussian()).fit()
    r2 = float(1.0 - res.deviance / res.null_deviance)
    y_lo, y_hi = float(y.min()), float(y.max())
    y_top = y_hi + 0.04 * (y_hi - y_lo if y_hi > y_lo else 1.0)
    return pd.Series({'rho': rho, 'r2': r2,
                      'x_left': float(x.min()), 'x_right': float(x.max()),
                      'y_top': y_top})


def fit_interim_regress_H1x_on_wz(wa: pd.DataFrame):
    """
    Per-item Strong-Oakley regression label: fit a binomial GLM
    ``pps_H1_x ~ w_ratio`` (logit link) on each item's S training rows and
    predict ``p_h1_xz = pi_hat(W(z^(s)))`` for every draw. Also returns the per-
    item diagnostic (rho, R^2) used by the regression diagnostic plot. When
    ``pps_H1_x`` is constant for an item or the GLM fails to converge, fall back
    to the empirical mean (a proper constant pi_hat).

    Returns
    -------
    (p_h1_xz, perf) : tuple of pd.DataFrame
        - ``p_h1_xz``: one row per (item, draw) with ``item_label``,
          ``item_type``, ``item_high_label``, ``p_h1_xz`` (predicted label),
          ``s`` (= draw + 1).
        - ``perf``: one row per item with ``item_label``, ``item_type``,
          ``rho`` and ``r2``.
    """
    p_rows = []
    perf_rows = []
    for item in wa['item_label'].unique():
        g = wa[wa['item_label'] == item].copy()
        mask = np.isfinite(g['w_ratio']) & np.isfinite(g['pps_H1_x'])
        g_ok = g[mask]
        if len(g_ok) < 3 or g_ok['pps_H1_x'].nunique() < 2:
            pi_hat = np.full(len(g),
                             float(g['pps_H1_x'].mean()) if len(g) else np.nan,
                             dtype=float)
        else:
            X_fit = sm.add_constant(g_ok['w_ratio'].to_numpy())
            try:
                res = sm.GLM(g_ok['pps_H1_x'].to_numpy(), X_fit,
                             family=sm.families.Binomial()).fit()
                w_fill = float(np.nanmedian(g['w_ratio']))
                X_all = sm.add_constant(g['w_ratio'].fillna(w_fill).to_numpy())
                pi_hat = np.asarray(res.predict(X_all), dtype=float)
            except Exception:
                pi_hat = np.full(len(g),
                                 float(g_ok['pps_H1_x'].mean()), dtype=float)
        gout = g[['item_label', 'item_type', 'item_high_label', 'draw']].copy()
        gout['p_h1_xz'] = pi_hat
        gout['s'] = gout['draw'].astype(int) + 1
        p_rows.append(gout.drop(columns='draw'))

        s = fit_interim_regress_H1x_on_wz_per_item_summary(g)
        perf_rows.append({
            'item_label': item, 'item_type': g['item_type'].iloc[0],
            'rho': float(s['rho']), 'r2': float(s['r2']),
        })

    p_h1_xz = pd.concat(p_rows, ignore_index=True)
    perf = pd.DataFrame(perf_rows)
    return p_h1_xz, perf


def fit_interim_regress_endptx_on_wz(wa: pd.DataFrame, pps_H1_def: float = 0.5):
    """
    Per-item Strong-Oakley regression label using the continuous endpoint ratio
    as the target (more discriminative than binarising to ``pps_H1_x``).

    For each item we fit a Gaussian GLM (identity link)
    ``pps_ratio_x ~ w_ratio`` on the S training rows, then convert the
    predictive distribution of ratio | W(z) into a label probability via
    ``p_h1_xz = P(ratio > pps_H1_def | W) = 1 - Phi((pps_H1_def - mu_hat) / sigma_hat)``,
    where ``sigma_hat = sqrt(res.scale)`` is the fitted residual SD.

    When ``w_ratio`` is degenerate or the GLM fails, falls back to a constant
    ``mu_hat = mean(pps_ratio_x)`` with the empirical SD; if the SD is not
    positive/finite, returns the indicator ``1{mu_hat > pps_H1_def}``.

    Returns
    -------
    (p_h1_xz, perf) : tuple of pd.DataFrame
        - ``p_h1_xz``: one row per (item, draw) with ``item_label``,
          ``item_type``, ``item_high_label``, ``p_h1_xz``, ``mu_hat``,
          ``sigma_hat``, ``s`` (= draw + 1).
        - ``perf``: one row per item with ``item_label``, ``item_type``,
          ``rho`` (Pearson on ``w_ratio`` vs ``pps_ratio_x``) and ``r2``
          (Gaussian-GLM deviance R^2).
    """
    from scipy.stats import norm

    p_rows = []
    perf_rows = []
    for item in wa['item_label'].unique():
        g = wa[wa['item_label'] == item].copy()
        mask = np.isfinite(g['w_ratio']) & np.isfinite(g['pps_ratio_x'])
        g_ok = g[mask]
        if len(g_ok) < 3:
            mu_const = float(g['pps_ratio_x'].mean()) if len(g) else np.nan
            mu_hat = np.full(len(g), mu_const, dtype=float)
            sigma = (float(g['pps_ratio_x'].std(ddof=1))
                     if len(g) >= 2 else float('nan'))
        else:
            X_fit = sm.add_constant(g_ok['w_ratio'].to_numpy())
            try:
                res = sm.GLM(g_ok['pps_ratio_x'].to_numpy(), X_fit,
                             family=sm.families.Gaussian()).fit()
                w_fill = float(np.nanmedian(g['w_ratio']))
                X_all = sm.add_constant(g['w_ratio'].fillna(w_fill).to_numpy())
                mu_hat = np.asarray(res.predict(X_all), dtype=float)
                sigma = float(np.sqrt(res.scale)) if res.scale > 0 else float('nan')
            except Exception:
                mu_hat = np.full(len(g), float(g_ok['pps_ratio_x'].mean()),
                                 dtype=float)
                sigma = float(g_ok['pps_ratio_x'].std(ddof=1))

        if not np.isfinite(sigma) or sigma <= 0:
            pi_hat = (mu_hat > pps_H1_def).astype(float)
        else:
            pi_hat = 1.0 - norm.cdf((pps_H1_def - mu_hat) / sigma)

        gout = g[['item_label', 'item_type', 'item_high_label', 'draw']].copy()
        gout['p_h1_xz'] = pi_hat
        gout['mu_hat'] = mu_hat
        gout['sigma_hat'] = sigma
        gout['s'] = gout['draw'].astype(int) + 1
        p_rows.append(gout.drop(columns='draw'))

        s = fit_interim_regress_H1x_on_wz_per_item_summary(g)
        perf_rows.append({
            'item_label': item, 'item_type': g['item_type'].iloc[0],
            'rho': float(s['rho']), 'r2': float(s['r2']),
        })

    p_h1_xz = pd.concat(p_rows, ignore_index=True)
    perf = pd.DataFrame(perf_rows)
    return p_h1_xz, perf
