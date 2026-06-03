"""
Interim-analysis helpers for predictive probability of success (PPS).

- ``get_interim_z``: build the 'missing' future data block z from the interim
  cohort xi, the final cohort xf, and posterior-predictive draws ypred.
- ``fit_interim_MC_of_posterior_xz``: Monte-Carlo estimate of p(H_1 | x, z_s)
  across S hypothetical future datasets, refitting the partial credit model to
  (xi + z_s) for each s.
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

# IRT-specific endpoint helpers now live on IRTModel
# (python/model_irt.py). Algorithm callers route through a Model instance.
from model_pcm import PartialCreditModelNCats
# Pure-pandas helpers live in interim_helpers.py (OO-port step 1). Re-exported
# here for the convenience of the interim-analysis scripts, which already had
# imports against this module before the refactor.
from interim_helpers import (
    _load_ypred,
    get_interim_x,
    get_interim_z_from_ypredf,
    get_interim_z_from_ypredi,
)


def _fit_interim_MC_of_posterior_xz_one_sample(
    s_idx, xi, zi, dit, output_file_prefix, model_cls,
    fitting_method_args, categorical_threshold, pps_H1_def, seed, verbose,
):
    """
    Refit the model on (xi + z_{s_idx}) and return the per-item p(H_1 | x, z_s)
    frame for a single Monte-Carlo sample. Module-level so it is picklable for
    process-based parallelism.
    """
    vprint = print if verbose else (lambda *args, **kwargs: None)
    s_label = s_idx + 1
    vprint(f"\n--- PPS sample {s_label} ---")

    # This sample's predicted outcomes become y_stan for the new participants;
    # concat onto the interim data and rebuild the cohort.
    ycol = f'ypred_{s_idx}'
    tmp = zi.drop(
        columns=['y', 'y_stan'] + [c for c in zi.columns if c.startswith('ypred_') and c != ycol],
        errors='ignore',
    ).rename(columns={ycol: 'y_stan'})
    xzi = get_interim_x(pd.concat([xi, tmp], ignore_index=True))

    # Refit the model to the augmented data (x, z_s), in memory only.
    fma = dict(fitting_method_args)
    x_formula = fma.pop('x_formula', "~ time - 1")
    model = model_cls(dit=dit, dcati=xzi, x_formula=x_formula,
                      seed=seed + s_label)
    fit = model.fit_pyro_svi(
        output_file_prefix=f"{output_file_prefix}_s{s_label}",
        save_to_file=False,
        resume=False,
        with_core_analyses=False,
        with_additional_analyses=False,
        verbose=verbose,
        **fma,
    )

    # Per-item P(H_1 | x, z_s) = fraction of draws with ratio > pps_H1_def.
    ratio_xz = model.get_endpoints_per_draw(
        draws=fit['draws'],
        categorical_threshold=categorical_threshold,
        endpoint_type='items',
        param_name='ordered_prob_by_cat_qu_pr',
        verbose=verbose,
    )
    sample_p = (
        ratio_xz.groupby(['item_label', 'item_type', 'item_high_label'])['ratio']
        .apply(lambda r: float((r > pps_H1_def).mean()))
        .reset_index(name='p_h1_xz')
    )
    sample_p['s'] = s_label
    return sample_p


def fit_interim_MC_of_posterior_xz(
    xi: pd.DataFrame,
    zi: pd.DataFrame,
    dit: pd.DataFrame,
    output_file_prefix: str,
    fitting_method_args: dict,
    pps_z_total: int = 10,
    pps_H1_def: float = 0.5,
    pps_ProbH1_thresh: float = 0.89,
    categorical_threshold: int = 3,
    model_cls=PartialCreditModelNCats,
    seed: int = 123,
    save_to_file: Optional[bool] = True,
    verbose: bool = True,
    cpu_n: int = 1,
) -> pd.DataFrame:
    """
    Monte-Carlo estimate of p(H_1 | x, z_s) over S hypothetical future datasets.

    For each sample ``s`` (0 .. ``pps_z_total`` - 1) the column ``ypred_s`` of
    ``zi`` becomes the outcome ``y_stan`` for the new participants; this is
    concatenated onto the interim data ``xi`` and rebuilt into a cohort frame via
    :func:`get_interim_x`. The partial credit model is refit in memory
    (``save_to_file=False``) and the per-item probability that the directional
    improvement ratio exceeds ``pps_H1_def`` is recorded.

    Parameters
    ----------
    xi : pd.DataFrame
        Interim cohort.
    zi : pd.DataFrame
        Future-data block from :func:`get_interim_z` (must hold ``ypred_0 ..
        ypred_{pps_z_total-1}`` columns).
    dit : pd.DataFrame
        Item metadata.
    output_file_prefix : str
        Path prefix; the result is written to ``{output_file_prefix}_p_h1_xz.pkl``.
    pps_z_total : int, default 10
        Number of Monte-Carlo samples (ypred draws) to use.
    pps_H1_def : float, default 0.5
        H_1 threshold on the per-draw improvement ratio (ratio > pps_H1_def).
    pps_ProbH1_thresh : float, default 0.89
        Decision threshold on p(H_1 | data); recorded on the output for the
        downstream PPS aggregation.
    categorical_threshold : int, default 3
        Categorical aggregation threshold passed to
        :func:`get_endpoints_per_draw`.
    model_cls : type, default :class:`model_pcm.PartialCreditModelNCats`
        Model subclass instantiated per sample as
        ``model_cls(dit=dit, dcati=xzi, x_formula=..., seed=...)``. The worker
        calls ``model.fit_pyro_svi(output_file_prefix=...,
        **fitting_method_args_without_x_formula)`` and consumes ``model
        .get_endpoints_per_draw`` to score the resulting cohort. Spawn workers
        pickle the class reference; supply an importable top-level class.
    fitting_method_args : dict
        Keyword arguments forwarded to ``model.fit_pyro_svi`` (e.g.
        ``algorithm``, ``lr``, ``num_steps``, ``output_samples``). May include
        ``x_formula`` -- popped before forwarding and used to construct the
        model.  Required (pass ``{}`` for the fit defaults). Must not include
        any of the keys this function sets itself.
    seed : int, default 123
        Base seed; sample ``s`` uses ``seed + s``.
    cpu_n : int, default 1
        Number of worker processes for the S independent refits. ``1`` runs the
        loop serially. ``> 1`` uses a ``spawn`` ``ProcessPoolExecutor`` (each
        worker re-imports JAX/NumPyro and fits fresh).

    Returns
    -------
    pd.DataFrame
        p_h1_xz: long form, one row per (item, sample) with columns
        ``item_label``, ``item_type``, ``item_high_label``, ``p_h1_xz``, ``s``,
        plus the recorded ``pps_H1_def``, ``pps_ProbH1_thresh`` and ``S``.
    """
    vprint = print if verbose else (lambda *args, **kwargs: None)

    work = dict(
        xi=xi, zi=zi, dit=dit, output_file_prefix=output_file_prefix,
        model_cls=model_cls, fitting_method_args=fitting_method_args,
        categorical_threshold=categorical_threshold, pps_H1_def=pps_H1_def,
        seed=seed, verbose=verbose,
    )

    if cpu_n == 1:
        rows = [_fit_interim_MC_of_posterior_xz_one_sample(s_idx, **work) for s_idx in range(pps_z_total)]
    else:
        import multiprocessing as _mp
        from concurrent.futures import ProcessPoolExecutor
        vprint(f"Running {pps_z_total} PPS refits over {cpu_n} worker processes...")
        ctx = _mp.get_context('spawn')
        with ProcessPoolExecutor(max_workers=cpu_n, mp_context=ctx) as ex:
            futures = [ex.submit(_fit_interim_MC_of_posterior_xz_one_sample, s_idx, **work) for s_idx in range(pps_z_total)]
            rows = [f.result() for f in futures]

    p_h1_xz = pd.concat(rows, ignore_index=True)
    p_h1_xz['pps_H1_def'] = pps_H1_def
    p_h1_xz['pps_ProbH1_thresh'] = pps_ProbH1_thresh
    p_h1_xz['S'] = pps_z_total

    if save_to_file:
        pkl_path = f"{output_file_prefix}_p_h1_xz.pkl"
        p_h1_xz.to_pickle(pkl_path)
        vprint(f"\nSaved P(H_1 | x, z) samples to: {pkl_path}")

    return p_h1_xz


def fit_interim_IS_reweight(
    model,
    zi: pd.DataFrame,
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
    """
    Importance-sampling estimate of p(H_1 | x, z_s) for fixed x. Model-aware
    rewrite of :func:`fit_interim_importance_sampling_of_posterior_xz_from_x`
    that consumes a :class:`Model` instance (which holds dit / dcati / x_formula
    and exposes ``eval_loglik`` / ``endpoints_per_draw`` /
    ``stack_posterior_theta`` / ``make_stan_data``).

    Returns the same ``(p_h1_xz, is_perf)`` tuple as the legacy free function.
    """
    if draws is None and draws_file is None:
        raise ValueError("Provide either draws or draws_file.")
    if draws is None:
        draws = az.from_zarr(draws_file)
    if 'src_pid' not in zi.columns:
        raise ValueError("zi must carry 'src_pid' (use get_interim_z_from_ypredi).")
    vprint = print if verbose else (lambda *args, **kwargs: None)

    # Per-draw H_1 ratios from the x-posterior; 'draw' aligns with theta_k below.
    x_ratio = model.get_endpoints_per_draw(
        draws=draws,
        categorical_threshold=categorical_threshold,
        endpoint_type='items',
        param_name='ordered_prob_by_cat_qu_fit',
    )
    x_ratio = x_ratio.assign(ind=(x_ratio['ratio'] > pps_H1_def).astype(float))

    theta = model.stack_posterior_theta(draws)
    # First param array's leading dim is K (number of posterior draws).
    n_draw = int(next(iter(theta.values())).shape[0])

    p_h1_xz = []
    perf_rows = []
    for s_idx in range(pps_z_total):
        s_label = s_idx + 1
        vprint(f"\n--- IS sample {s_label}/{pps_z_total} ---")

        # z_s: resampled x participants (unit = src_pid) carrying this draw's
        # predicted outcomes; rebuild oid/oidt for the stan_data builder.
        zcol = f'ypred_{s_idx}'
        z_dcati = zi.assign(pid=zi['src_pid'], y_stan=zi[zcol].astype(int))
        z_dcati['y'] = z_dcati['y_stan'] - 1
        z_dcati = z_dcati.sort_values(
            ['item_type_id', 'pid', 'time', 'item_label']
        ).reset_index(drop=True)
        z_dcati['oid'] = np.arange(1, len(z_dcati) + 1)
        z_dcati['oidt'] = z_dcati.groupby('item_type').cumcount() + 1
        z_stan = model.make_stan_data(z_dcati, model.x_formula)

        # log p(z_s | theta_k) for every draw via vmap over the stacked params.
        loglik = jax.vmap(lambda p: model.eval_loglik(z_stan, p))(theta)
        # Numerically-stable softmax over total-loglik per draw.
        w = np.asarray(jax.nn.softmax(loglik.sum(axis=1)))  # (K,) sums to 1

        sum_w2 = float(np.sum(w ** 2))
        perf_rows.append({
            's': s_label,
            'N': n_draw,
            'ess': float(1.0 / sum_w2),
            'ess_over_n': float(1.0 / (n_draw * sum_w2)),
            'ew2': float(np.mean(w ** 2)),
        })

        # Self-normalised IS estimate of P(H_1 | x, z_s) per item.
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
    p_h1_xz['pps_H1_def'] = pps_H1_def
    p_h1_xz['pps_ProbH1_thresh'] = pps_ProbH1_thresh
    p_h1_xz['S'] = pps_z_total
    is_perf = pd.DataFrame(perf_rows)

    if save_to_file and output_file_prefix is not None:
        pkl_path = f"{output_file_prefix}_p_h1_xz_IS.pkl"
        p_h1_xz.to_pickle(pkl_path)
        is_perf.to_csv(f"{output_file_prefix}_is_perf.csv", index=False)
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


def _stack_posterior_theta(draws) -> dict:
    """Flatten the (chain, draw) axes of an arviz posterior and return the
    parameter arrays any partial-credit-family log-likelihood consumes (only the
    params present are returned, so credit / ordered-logit draws also work)."""
    post = draws.posterior
    names = (
        'latent_factor_unit', 'latent_factor_beta',
        'skill_thresholds', 'skill_thresholds_1', 'skill_thresholds_incs',
        'loadings_questions_m1',
    )
    return {
        name: jnp.asarray(np.asarray(post[name].values).reshape(-1, *post[name].shape[2:]))
        for name in names if name in post
    }


def _interim_make_x_stan(model) -> dict:
    """stan_data for the x cohort (the full p(x | theta) factor)."""
    x_dcati = model.dcati.copy()
    if 'y_stan' not in x_dcati:
        x_dcati['y_stan'] = x_dcati['y'] + 1
    x_dcati = x_dcati.sort_values(
        ['item_type_id', 'pid', 'time', 'item_label']
    ).reset_index(drop=True)
    x_dcati['oid'] = np.arange(1, len(x_dcati) + 1)
    x_dcati['oidt'] = x_dcati.groupby('item_type').cumcount() + 1
    return model.make_stan_data(x_dcati, model.x_formula)


def _interim_make_z_stan(model, zi: pd.DataFrame, s_idx: int) -> dict:
    """stan_data for future-data sample s (mirrors the IS function)."""
    zcol = f'ypred_{s_idx}'
    z_dcati = zi.assign(pid=zi['src_pid'], y_stan=zi[zcol].astype(int))
    z_dcati['y'] = z_dcati['y_stan'] - 1
    z_dcati = z_dcati.sort_values(
        ['item_type_id', 'pid', 'time', 'item_label']
    ).reset_index(drop=True)
    z_dcati['oid'] = np.arange(1, len(z_dcati) + 1)
    z_dcati['oidt'] = z_dcati.groupby('item_type').cumcount() + 1
    return model.make_stan_data(z_dcati, model.x_formula)


def _ratio_per_draw_from_params(theta_batch, x_stan, model, categorical_threshold,
                                ordered_prob_eval):
    """Per-(draw, item) improvement ratio for an arbitrary batch of parameter sets.

    Evaluates ``ordered_prob_by_cat_qu_fit`` for every draw in ``theta_batch``
    (dict of ``(K, ...)`` jnp arrays) on the x-cohort design ``x_stan``, wraps the
    result as an arviz posterior and runs ``model.get_endpoints_per_draw``. Lets
    reweighted (moment-matched) draws or moved (SMC) particles be scored for
    p(H_1 | x, z) exactly like a refit, rather than reusing the frozen x-ratio.
    """
    ordprob = np.asarray(jax.vmap(lambda p: ordered_prob_eval(x_stan, p))(theta_batch))  # (K, L)
    idata = az.from_dict(posterior={'ordered_prob_by_cat_qu_fit': ordprob[None, ...]})
    return model.get_endpoints_per_draw(
        draws=idata,
        categorical_threshold=categorical_threshold,
        endpoint_type='items', param_name='ordered_prob_by_cat_qu_fit',
        verbose=False,
    )


def fit_interim_IS_moment_matching(
    model,
    zi: pd.DataFrame,
    fitting_method_args: Optional[dict] = None,
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
    """
    Mean-match IS for p(H_1 | x, z) (Paananen 2021), model-aware rewrite of
    :func:`fit_interim_IS_moment_matching_of_posterior_xz_from_x`. Reads
    ``eval_loglik`` / ``eval_log_prior`` / ``positive_params`` from ``model``; the
    legacy callable overrides in ``fitting_method_args`` are ignored. The
    only honoured ``fitting_method_args`` key is ``x_formula`` (carried on
    ``model.x_formula`` -- kept for shim compatibility).

    Returns the same ``(p_h1_xz, mm)`` tuple as the legacy function.
    """
    if draws is None and draws_file is None:
        raise ValueError("Provide either draws or draws_file.")
    if draws is None:
        draws = az.from_zarr(draws_file)
    if 'src_pid' not in zi.columns:
        raise ValueError("zi must carry 'src_pid' (use get_interim_z_from_ypredi).")
    vprint = print if verbose else (lambda *args, **kwargs: None)
    if fitting_method_args is None:
        fitting_method_args = {}

    eval_loglik = model.eval_loglik
    eval_log_prior = model.eval_log_prior
    positive_params = model.positive_params

    theta = model.stack_posterior_theta(draws)
    n_draw = int(next(iter(theta.values())).shape[0])
    x_stan = _interim_make_x_stan(model)

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
    z_stan_template = _interim_make_z_stan(model, zi, 0)

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
    for s_idx in range(pps_z_total):
        z_stan = _interim_make_z_stan(model, zi, s_idx)
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
        # then averaged with the mean-match weights w1.
        ratio = _ratio_per_draw_from_params(
            _u_batch_to_params(jnp.asarray(Theta_u_star)),
            x_stan, model, categorical_threshold,
            ordered_prob_eval=model.eval_outcome_for_endpoint,
        )
        ratio = ratio.assign(ind=(ratio['ratio'] > pps_H1_def).astype(float))
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
    p_h1_xz['pps_H1_def'] = pps_H1_def
    p_h1_xz['pps_ProbH1_thresh'] = pps_ProbH1_thresh
    p_h1_xz['S'] = pps_z_total

    if save_to_file and output_file_prefix is not None:
        mm.to_csv(f"{output_file_prefix}_momentmatch.csv", index=False)
        p_h1_xz.to_pickle(f"{output_file_prefix}_p_h1_xz_MM.pkl")
        vprint(f"\nSaved moment-matching diagnostics to: {output_file_prefix}_momentmatch.csv")
    return p_h1_xz, mm


def fit_interim_SMC_resample(
    model,
    zi: pd.DataFrame,
    fitting_method_args: dict,
    draws=None,
    draws_file: Optional[str] = None,
    output_file_prefix: Optional[str] = None,
    save_to_file: bool = True,
    verbose: bool = True,
):
    """
    SMC sampler with resample-move for p(theta | x, z_s), model-aware rewrite
    of :func:`fit_interim_SMC_resample_of_posterior_xz_from_x`. Reads
    ``eval_loglik`` / ``eval_loglik_annealed`` / ``eval_log_prior`` /
    ``stack_posterior_theta`` from ``model``; the legacy callable overrides in
    ``fitting_method_args`` are ignored. The PCM-specific param layout
    (latent_factor_unit / latent_factor_beta / skill_thresholds /
    loadings_questions_m1) is currently still hard-coded here -- generalising
    that layout to credit / ordered-logit / Binomial is OO-port step 6+.

    Returns the same ``(schedule, particles)`` tuple as the legacy function.
    """
    if draws is None and draws_file is None:
        raise ValueError("Provide either draws or draws_file.")
    if draws is None:
        draws = az.from_zarr(draws_file)
    if 'src_pid' not in zi.columns:
        raise ValueError("zi must carry 'src_pid' (use get_interim_z_from_ypredi).")
    vprint = print if verbose else (lambda *args, **kwargs: None)

    if fitting_method_args is None:
        fitting_method_args = {}
    fitting_method_args = {
        's_idx': 0,
        'n_particles': 128,
        'ess_frac_target': 0.5,
        'n_move_steps': 20,
        'init_step_size': 0.02,
        'target_accept': 0.574,
        'max_temps': 150,
        'seed': 123,
        **fitting_method_args,
    }
    s_idx = fitting_method_args['s_idx']
    n_particles = fitting_method_args['n_particles']
    ess_frac_target = fitting_method_args['ess_frac_target']
    n_move_steps = fitting_method_args['n_move_steps']
    init_step_size = fitting_method_args['init_step_size']
    target_accept = fitting_method_args['target_accept']
    max_temps = fitting_method_args['max_temps']
    seed = fitting_method_args['seed']

    eval_loglik = model.eval_loglik
    eval_loglik_annealed = model.eval_loglik_annealed
    eval_log_prior = model.eval_log_prior

    theta = model.stack_posterior_theta(draws)
    n_draw = int(theta['latent_factor_beta'].shape[0])
    U = int(theta['latent_factor_unit'].shape[1])
    P = int(theta['latent_factor_beta'].shape[1])
    L = int(theta['skill_thresholds'].shape[1])
    Ld = int(theta['loadings_questions_m1'].shape[1])

    x_stan = _interim_make_x_stan(model)
    z_stan = _interim_make_z_stan(model, zi, s_idx)
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
    if save_to_file and output_file_prefix is not None:
        schedule.to_csv(f"{output_file_prefix}_smc_schedule.csv", index=False)
        vprint(f"Saved SMC schedule to: {output_file_prefix}_smc_schedule.csv")
    return schedule, particles


def _fit_interim_SMC_one_sample(
    s_idx, xi, zi, dit, draws_file, fitting_method_args,
    pps_H1_def, categorical_threshold, verbose,
):
    """Run the SMC sampler for a single future-data sample and return its
    per-item p(H_1 | x, z_s) plus a one-row schedule summary. Module-level so it
    is picklable for process-based parallelism. The post-move particles are
    uniformly weighted, so p(H_1 | x, z_s) is the fraction with ratio > pps_H1_def.
    """
    fma = {**fitting_method_args, 's_idx': s_idx}
    # Build a PCM model from (xi, dit). Spawn workers re-import every module
    # anyway, so the lazy import keeps the helper picklable.
    from model_pcm import PartialCreditModelNCats
    model = PartialCreditModelNCats(
        dit=dit, dcati=xi,
        x_formula=fma.get('x_formula', "~ time - 1"),
    )
    schedule, particles = fit_interim_SMC_resample(
        model=model, zi=zi, draws_file=draws_file,
        fitting_method_args=fma, save_to_file=False, verbose=verbose,
    )
    x_stan = _interim_make_x_stan(model)
    ratio = _ratio_per_draw_from_params(
        particles, x_stan, model, categorical_threshold,
        ordered_prob_eval=model.eval_outcome_for_endpoint,
    )
    sample_p = (
        ratio.assign(ind=(ratio['ratio'] > pps_H1_def).astype(float))
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


def fit_interim_SMC_PPS_of_posterior_xz_from_x(
    xi: pd.DataFrame,
    zi: pd.DataFrame,
    dit: pd.DataFrame,
    fitting_method_args: dict,
    draws_file: str,
    pps_z_total: int = 10,
    pps_H1_def: float = 0.5,
    pps_ProbH1_thresh: float = 0.89,
    categorical_threshold: int = 3,
    cpu_n: int = 1,
    output_file_prefix: Optional[str] = None,
    save_to_file: bool = True,
    verbose: bool = True,
):
    """
    SMC resample-move PPS over S future datasets at a fixed x (Case A).

    Runs :func:`fit_interim_SMC_resample_of_posterior_xz_from_x` for each future
    sample z_s and scores p(H_1 | x, z_s) as the fraction of the moved particles
    whose per-item improvement ratio exceeds ``pps_H1_def`` (the particles are
    uniformly weighted). The S SMC runs are independent and parallelised over
    ``cpu_n`` ``spawn`` workers (each re-imports JAX and reloads ``draws_file``),
    mirroring :func:`fit_interim_MC_of_posterior_xz`.

    Parameters
    ----------
    fitting_method_args : dict
        SMC tuning forwarded to the per-sample sampler (``n_particles``,
        ``ess_frac_target``, ``n_move_steps``, ``init_step_size``,
        ``target_accept``, ``max_temps``, ``x_formula``, ``seed`` ...); ``s_idx``
        is set per sample. Must be picklable (no in-memory callables) when
        ``cpu_n > 1``.
    draws_file : str
        Path to the x-fit zarr (passed instead of in-memory draws so workers can
        reload it).
    cpu_n : int, default 1
        Worker processes for the S SMC runs.

    Returns
    -------
    (p_h1_xz, smc_summary) : tuple of pd.DataFrame
        - ``p_h1_xz``: one row per (item, sample s) with ``p_h1_xz``, ``s`` and
          the recorded ``pps_H1_def`` / ``pps_ProbH1_thresh`` / ``S``.
        - ``smc_summary``: one row per s with ``n_particles``, ``n_steps``,
          ``beta_final`` and ``mins`` (per-sample SMC wall time).
    """
    vprint = print if verbose else (lambda *args, **kwargs: None)
    work = dict(
        xi=xi, zi=zi, dit=dit, draws_file=draws_file,
        fitting_method_args=fitting_method_args, pps_H1_def=pps_H1_def,
        categorical_threshold=categorical_threshold, verbose=False,
    )
    if cpu_n == 1:
        results = [_fit_interim_SMC_one_sample(s, **work) for s in range(pps_z_total)]
    else:
        import multiprocessing as _mp
        from concurrent.futures import ProcessPoolExecutor
        vprint(f"Running {pps_z_total} SMC samples over {cpu_n} worker processes...")
        ctx = _mp.get_context('spawn')
        with ProcessPoolExecutor(max_workers=cpu_n, mp_context=ctx) as ex:
            futures = [ex.submit(_fit_interim_SMC_one_sample, s, **work) for s in range(pps_z_total)]
            results = [f.result() for f in futures]

    p_h1_xz = pd.concat([r[0] for r in results], ignore_index=True)
    p_h1_xz['pps_H1_def'] = pps_H1_def
    p_h1_xz['pps_ProbH1_thresh'] = pps_ProbH1_thresh
    p_h1_xz['S'] = pps_z_total
    smc_summary = pd.DataFrame([r[1] for r in results])

    if save_to_file and output_file_prefix is not None:
        p_h1_xz.to_pickle(f"{output_file_prefix}_p_h1_xz_SMC.pkl")
        smc_summary.to_csv(f"{output_file_prefix}_smc_summary.csv", index=False)
        vprint(f"\nSaved SMC P(H_1 | x, z) samples to: {output_file_prefix}_p_h1_xz_SMC.pkl")
    return p_h1_xz, smc_summary


def fit_interim_SMC_PPS(
    model,
    zi: pd.DataFrame,
    fitting_method_args: dict,
    draws_file: str,
    pps_z_total: int = 10,
    pps_H1_def: float = 0.5,
    pps_ProbH1_thresh: float = 0.89,
    categorical_threshold: int = 3,
    cpu_n: int = 1,
    output_file_prefix: Optional[str] = None,
    save_to_file: bool = True,
    verbose: bool = True,
):
    """
    Model-aware wrapper around
    :func:`fit_interim_SMC_PPS_of_posterior_xz_from_x`. Unpacks
    ``model.dcati`` / ``model.dit`` / ``model.x_formula`` into the legacy
    free-function signature. The multiprocessing worker
    (:func:`_fit_interim_SMC_one_sample`) still takes the picklable
    ``(xi, zi, dit)`` triple -- a deeper refactor that hands a Model into
    spawn workers is deferred to OO-port step 6+ once we know what state the
    non-PCM subclasses carry.
    """
    fma = {**(fitting_method_args or {}), 'x_formula': model.x_formula}
    return fit_interim_SMC_PPS_of_posterior_xz_from_x(
        xi=model.dcati, zi=zi, dit=model.dit,
        fitting_method_args=fma, draws_file=draws_file,
        pps_z_total=pps_z_total,
        pps_H1_def=pps_H1_def, pps_ProbH1_thresh=pps_ProbH1_thresh,
        categorical_threshold=categorical_threshold,
        cpu_n=cpu_n,
        output_file_prefix=output_file_prefix,
        save_to_file=save_to_file, verbose=verbose,
    )


# =============================================================================
# Strong-Oakley regression-based label estimator (Case A)
# =============================================================================
# 1. ``get_interim_endpt_and_w_from_poi`` builds the per-item training set
#    ``wa`` from the x-posterior draws: per-(item, draw) the actual endpoint
#    ratio ``pps_ratio_x`` (and its H_1 indicator), and per-(item, draw) the
#    posterior-predictive summary ``W(z^(s))_t`` at baseline / endline plus the
#    direction-aware ``w_diff`` / ``w_ratio``.
# 2. ``fit_interim_regress_H1x_on_wz`` fits a per-item binomial GLM
#    ``pps_H1_x ~ w_ratio`` and predicts ``p_h1_xz`` = pi_hat(W(z^(s))) for
#    every draw, plus per-item Pearson rho and Gaussian-GLM R^2 diagnostics.
# 3. ``fit_interim_regress_H1x_on_wz_per_item_summary`` is the per-item helper
#    used both inside (2) and by the diagnostic plot scripts.


def get_interim_endpt_and_w_from_poi(
    model,
    draws,
    draws_file: str,
    interim_m: int,
    pps_z_total: int,
    pps_H1_def: float = 0.5,
    pps_ProbH1_thresh: float = 0.89,
    categorical_threshold: int = 3,
    seed: int = 123,
) -> pd.DataFrame:
    """
    Per-(item, draw) endpoint ratio from the x-posterior + per-(item, draw)
    summary ``W(z^(s))_t`` for the Strong-Oakley regression training set.

    For each draw s = 0..pps_z_total-1 (in posterior order; see
    :func:`_load_ypred` with ``keep_order=True``):

    - load ``ypred[s]`` (the posterior-predictive at every xi observation) and
      build a future-data block z by resampling ``interim_m`` xi participants
      with replacement (seeded), inheriting their ypred values;
    - compute per-(item, time) summary ``W(z^(s))_t``:
        * ``out-of-7``  items -> mean of y_stan (1..7),
        * ``categorical`` items -> proportion of responses
          ``>= categorical_threshold``;
    - pivot to ``w_baseline`` / ``w_endline`` and derive the direction-aware
      ``w_diff`` and ``w_ratio`` (positive => improvement, matching the
      endpoint-ratio convention).

    The actual per-draw endpoint ratio ``pps_ratio_x`` and the H_1 indicator
    ``pps_H1_x = 1{pps_ratio_x > pps_H1_def}`` come from
    ``model.endpoints_per_draw(...)`` on the x-fit posterior; per-item
    summaries ``pps_ProbH1_x`` and the decision indicator ``pps_H1_yes`` are
    also added.

    Returns
    -------
    wa : pd.DataFrame
        One row per (item, draw) with columns
        ``item_label, item_type, item_high_label, draw,
        w_baseline, w_endline, w_diff, w_ratio,
        pps_ratio_x, pps_H1_x``.
    """
    xi = model.dcati
    dit = model.dit

    # Load ypred in posterior-draw order so row s == posterior draw s.
    ypred = _load_ypred(draws_file, pps_z_total, np.random.default_rng(seed), keep_order=True)
    S = ypred.shape[0]

    # Resample participants (with replacement, seeded) for the shadow z block.
    rng = np.random.default_rng(seed)
    tmp = pd.DataFrame({
        'src_pid': rng.choice(xi['pid'].unique(), size=interim_m, replace=True),
        'pid': np.arange(1, interim_m + 1) + int(xi['pid'].max()),
    })
    zi = (
        tmp.merge(xi.rename(columns={'pid': 'src_pid'}), on='src_pid')
        .sort_values(['pid', 'oid']).reset_index(drop=True)
    )
    ypred_cols = [f'ypred_{s}' for s in range(S)]
    zi[ypred_cols] = ypred[:, zi['oid'].to_numpy() - 1].T

    # Per-(item, time) summaries W(z^(s))_t, type-specific.
    zi_o7 = zi[zi['item_type'] == 'out-of-7']
    zi_cat = zi[zi['item_type'] == 'categorical']
    w_o7 = (
        zi_o7.groupby(['item_label', 'item_type', 'time'])[ypred_cols]
        .mean().reset_index()
    )
    zi_cat_bin = zi_cat[ypred_cols].ge(categorical_threshold).astype(float)
    zi_cat_bin = pd.concat(
        [zi_cat[['item_label', 'item_type', 'time']].reset_index(drop=True),
         zi_cat_bin.reset_index(drop=True)], axis=1,
    )
    w_cat = (
        zi_cat_bin.groupby(['item_label', 'item_type', 'time'])[ypred_cols]
        .mean().reset_index()
    )
    wa = pd.concat([w_o7, w_cat], ignore_index=True).melt(
        id_vars=['item_label', 'item_type', 'time'],
        value_vars=ypred_cols, var_name='s_col', value_name='w',
    )
    wa['draw'] = wa['s_col'].str.replace('ypred_', '').astype(int)
    wa = wa.drop(columns='s_col')

    t_min, t_max = wa['time'].min(), wa['time'].max()
    wa = (
        wa.pivot_table(index=['item_label', 'item_type', 'draw'],
                       columns='time', values='w')
        .reset_index()
    )
    wa.columns.name = None
    wa = wa.rename(columns={t_min: 'w_baseline', t_max: 'w_endline'})

    # Direction-aware diff / ratio (mirrors the endpoint-ratio computation).
    wa = wa.merge(
        dit[['item_label', 'item_high_label']].drop_duplicates(),
        on='item_label', how='left',
    )
    wa['w_diff'] = np.nan
    wa['w_ratio'] = np.nan
    tmp = wa['item_high_label'] == 'lower_is_better'
    wa.loc[tmp, 'w_diff'] = wa.loc[tmp, 'w_baseline'] - wa.loc[tmp, 'w_endline']
    wa.loc[tmp, 'w_ratio'] = 1 - wa.loc[tmp, 'w_endline'] / wa.loc[tmp, 'w_baseline']
    tmp = wa['item_high_label'] == 'higher_is_better'
    wa.loc[tmp, 'w_diff'] = wa.loc[tmp, 'w_endline'] - wa.loc[tmp, 'w_baseline']
    wa.loc[tmp, 'w_ratio'] = wa.loc[tmp, 'w_endline'] / wa.loc[tmp, 'w_baseline'] - 1

    # Actual per-draw endpoint ratio from the x-fit posterior.
    x_ratio = model.get_endpoints_per_draw(
        draws=draws,
        categorical_threshold=categorical_threshold,
        endpoint_type='items',
        param_name='ordered_prob_by_cat_qu_fit',
    )
    x_ratio = x_ratio.rename(columns={'ratio': 'pps_ratio_x'})
    x_ratio['pps_H1_x'] = (x_ratio['pps_ratio_x'] > pps_H1_def).astype(int)
    tmp = (
        x_ratio.groupby(['item_type', 'item_label'])
        .agg(pps_ProbH1_x=('pps_H1_x', 'mean'))
        .reset_index()
    )
    tmp['pps_H1_yes'] = (tmp['pps_ProbH1_x'] > pps_ProbH1_thresh).astype(int)
    x_ratio = x_ratio.merge(
        tmp[['item_label', 'pps_ProbH1_x', 'pps_H1_yes']], on='item_label', how='left',
    )

    wa = wa.merge(
        x_ratio[['draw', 'item_label', 'item_type', 'pps_ratio_x', 'pps_H1_x']],
        on=['draw', 'item_label', 'item_type'], how='inner',
    )
    return wa


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
