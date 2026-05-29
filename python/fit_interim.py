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

from get_endpoints import get_endpoints_per_draw
from fit_partial_credit_model import (
    fit_partial_credit_model_ncats_pyrosvi,
    _fit_partial_credit_make_stan_data,
    eval_loglik_partial_credit_model_ncats,
    eval_loglik_partial_credit_model_ncats_with_annealing,
    get_prior_of_partial_credit_model_ncats,
    get_ordered_prob_of_partial_credit_model_ncats,
)


def _load_ypred(draws_file: str, pps_z_total: int, rng) -> np.ndarray:
    """
    Load the ``ypred`` posterior-predictive draws from a zarr, flatten the
    (chain, draw) dims, and return ``pps_z_total`` randomly selected draws as a
    ``(pps_z_total, N_total)`` array.
    """
    ypred = az.from_zarr(draws_file).posterior['ypred'].values  # (chain, draw, N_total)
    ypred = ypred.reshape(-1, ypred.shape[-1])                  # (n_draw, N_total)
    return ypred[rng.choice(ypred.shape[0], size=pps_z_total, replace=False)]

#%%
def get_interim_x(dp1: pd.DataFrame,
                  interim_date: Optional[str] = None
                  ) -> pd.DataFrame:
    """
    Slice dp1 to the cohort observed on/before ``interim_date``, keeping only
    (participant, item) responses present at BOTH baseline and endline, then
    re-index pids / oids / oidt.

    Participants who are missing a timepoint for some items are retained for
    their complete items; only the incomplete (pid, item) responses are
    dropped (not the whole participant).

    Used by the interim-analysis scripts to slice the longitudinal panel into
    monthly cohorts ready to be passed to the ncats fit helpers.
    """
    dcati = dp1

    if interim_date is not None:
        dcati = dcati[dcati['submission_date'] <= interim_date]

    if dcati.empty:
        return dcati

    # Keep (pid, item_label) responses observed at both timepoints; drop only
    # the incomplete item-responses, not the whole participant.
    n_per_pid_item = (
        dcati.groupby(['pid', 'item_label'])['time']
        .nunique().reset_index(name='n_times')
    )
    keep = n_per_pid_item.loc[n_per_pid_item['n_times'] == 2, ['pid', 'item_label']]
    dcati = dcati.merge(keep, on=['pid', 'item_label'], how='inner')
    if dcati.empty:
        return dcati

    pid_map = pd.DataFrame({'pid_orig': sorted(dcati['pid'].unique())})
    pid_map['pid_new'] = range(1, len(pid_map) + 1)
    dcati = dcati.merge(pid_map, left_on='pid', right_on='pid_orig')
    dcati['pid'] = dcati['pid_new']
    dcati = dcati.drop(columns=['pid_orig', 'pid_new'])

    dcati = dcati.sort_values(['item_type_id', 'pid', 'time', 'item_label']).reset_index(drop=True)
    dcati['oid'] = range(1, len(dcati) + 1)
    dcati['oidt'] = dcati.groupby('item_type').cumcount() + 1
    return dcati


def get_interim_z_from_ypredf(xi: pd.DataFrame, xf: pd.DataFrame, draws_file: str,
                              pps_z_total: int = 10, seed: int = 123) -> pd.DataFrame:
    """
    Build the future-data block z from the FINAL-cohort fit's predictions.

    z consists of the participants present in the final cohort ``xf`` but not in
    the interim cohort ``xi``, drawn at random (without replacement). Each of
    ``pps_z_total`` posterior-predictive draws read from ``draws_file`` (the
    final-cohort fit) is attached as a column ``ypred_0 .. ypred_{pps_z_total-1}``
    holding that draw's predicted outcome at the corresponding observation
    (matched via ``oid`` - 1, since ``oid`` is 1-indexed).

    Parameters
    ----------
    xi : pd.DataFrame
        Interim cohort (long form; must include ``pid`` and ``oid``).
    xf : pd.DataFrame
        Final cohort; its participant set is a superset of ``xi``.
    draws_file : str
        Path to the final-cohort fit zarr (holding ``ypred``).
    pps_z_total : int, default 10
        Number of posterior-predictive draws (Monte-Carlo samples) to attach.
    seed : int, default 123
        RNG seed for sampling new participants and draws.

    Returns
    -------
    pd.DataFrame
        zi: the ``xf`` rows for the sampled new participants, sorted by ``oid``,
        with one ``ypred_s`` column per selected draw.
    """
    interim_m = xf['pid'].nunique() - xi['pid'].nunique()
    if interim_m <= 0:
        raise ValueError(
            f"[get_interim_z_from_ypredf] no new participants: xf has "
            f"{xf['pid'].nunique()}, xi has {xi['pid'].nunique()}."
        )

    rng = np.random.default_rng(seed)
    ypred = _load_ypred(draws_file, pps_z_total, rng)  # (pps_z_total, N_total_xf)

    new_pids = np.setdiff1d(xf['pid'].unique(), xi['pid'].unique())
    zi_pids = rng.choice(new_pids, size=interim_m, replace=False)
    zi_pids.sort()

    zi = xf.merge(pd.DataFrame({'pid': zi_pids}), on='pid', how='inner')
    zi.sort_values('oid', inplace=True)
    zi.reset_index(drop=True, inplace=True)

    zi[[f'ypred_{s}' for s in range(pps_z_total)]] = ypred[:, zi['oid'].to_numpy() - 1].T
    return zi


def get_interim_z_from_ypredi(xi: pd.DataFrame, draws_file: str, interim_m: int,
                              pps_z_total: int = 10, seed: int = 123) -> pd.DataFrame:
    """
    Build the future-data block z from the INTERIM fit's predictions.

    Unlike :func:`get_interim_z_from_ypredf` there is no final cohort: the
    interim fit (``draws_file``) only contains the ``xi`` participants, and the
    number of new participants required (``interim_m``) typically exceeds the
    number of distinct ``xi`` participants. New participants are therefore
    resampled WITH replacement from ``xi``, each assigned a fresh ``pid`` offset
    above ``xi``'s maximum so they don't collide on a later concat, and their
    predicted outcomes are taken from the interim fit's ``ypred`` at the original
    ``xi`` observation positions.

    Parameters
    ----------
    xi : pd.DataFrame
        Interim cohort (long form; must include ``pid`` and ``oid``).
    draws_file : str
        Path to the interim-cohort fit zarr (holding ``ypred`` for ``xi``).
    interim_m : int
        Number of new (shadow) participants to create.
    pps_z_total : int, default 10
        Number of posterior-predictive draws (Monte-Carlo samples) to attach.
    seed : int, default 123
        RNG seed for sampling participants and draws.

    Returns
    -------
    pd.DataFrame
        zi: ``interim_m`` resampled ``xi`` participants (fresh pids), sorted by
        ``pid`` then ``oid``, with one ``ypred_s`` column per selected draw.
    """
    rng = np.random.default_rng(seed)
    ypred = _load_ypred(draws_file, pps_z_total, rng)  # (pps_z_total, N_total_xi)

    src_pids = rng.choice(xi['pid'].unique(), size=interim_m, replace=True)
    tmp = pd.DataFrame({
        'src_pid': src_pids,
        'pid': np.arange(1, interim_m + 1) + int(xi['pid'].max()),
    })
    # Keep ``src_pid`` (the original xi unit id) so downstream importance
    # sampling can index that unit's latent factor in the x-posterior.
    zi = tmp.merge(xi.rename(columns={'pid': 'src_pid'}), on='src_pid')
    zi.sort_values(['pid', 'oid'], inplace=True)
    zi.reset_index(drop=True, inplace=True)

    # ypred is indexed by the ORIGINAL xi observation (oid - 1); each resampled
    # row keeps its source oid, so the same lookup applies.
    zi[[f'ypred_{s}' for s in range(pps_z_total)]] = ypred[:, zi['oid'].to_numpy() - 1].T
    return zi


def _fit_interim_MC_of_posterior_xz_one_sample(
    s_idx, xi, zi, dit, output_file_prefix, fitting_method,
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
    fit = fitting_method(
        dit,
        xzi,
        output_file_prefix=f"{output_file_prefix}_s{s_label}",
        seed=seed + s_label,
        save_to_file=False,
        resume=False,
        with_core_analyses=False,
        with_additional_analyses=False,
        verbose=verbose,
        **fitting_method_args,
    )

    # Per-item P(H_1 | x, z_s) = fraction of draws with ratio > pps_H1_def.
    ratio_xz = get_endpoints_per_draw(
        dcati=xzi,
        dit=dit,
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
    fitting_method=fit_partial_credit_model_ncats_pyrosvi,
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
    fitting_method : callable, default ``fit_partial_credit_model_ncats_pyrosvi``
        Model-fitting function called per sample. It must accept ``(dit, xzi)``
        positionally plus the keyword args ``output_file_prefix``, ``seed``,
        ``save_to_file``, ``resume``, ``with_core_analyses``,
        ``with_additional_analyses`` (all set by this function), and return a
        dict with a ``'draws'`` idata.
    fitting_method_args : dict
        Keyword arguments forwarded to ``fitting_method`` (e.g. ``algorithm``,
        ``x_formula``, ``lr``, ``num_steps``, ``output_samples``); required,
        specify at the call site (pass ``{}`` for the fit defaults). Must not
        include any of the keys this function sets itself.
    seed : int, default 123
        Base seed; sample ``s`` uses ``seed + s``.
    cpu_n : int, default 1
        Number of worker processes for the S independent refits. ``1`` runs the
        loop serially. ``> 1`` uses a ``spawn`` ``ProcessPoolExecutor`` (each
        worker re-imports JAX/NumPyro and fits fresh, so ``fitting_method`` must
        be an importable top-level function, not a lambda/closure).

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
        fitting_method=fitting_method, fitting_method_args=fitting_method_args,
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


def fit_interim_importance_sampling_of_posterior_xz_from_x(
    xi: pd.DataFrame,
    zi: pd.DataFrame,
    dit: pd.DataFrame,
    draws=None,
    draws_file: Optional[str] = None,
    pps_z_total: int = 10,
    pps_H1_def: float = 0.5,
    pps_ProbH1_thresh: float = 0.89,
    categorical_threshold: int = 3,
    x_formula: str = "~ time - 1",
    eval_loglik=eval_loglik_partial_credit_model_ncats,
    output_file_prefix: Optional[str] = None,
    save_to_file: bool = True,
    verbose: bool = True,
):
    """
    Importance-sampling estimate of p(H_1 | x, z_s) for fixed x (Case A).

    Instead of refitting the model for each future dataset z_s (cf.
    :func:`fit_interim_MC_of_posterior_xz`), reuse the existing x-posterior
    draws ``theta_k ~ p(theta | x)`` and reweight them by the likelihood of the
    future data (dev/amortised_decision_making.md, Step 7 Case A):

        w_k^(s)  ∝ p(z_s | theta_k)
        p(H_1 | x, z_s)_item ≈ (sum_k w_k^(s) 1[ratio_{k,item} > pps_H1_def])
                               / (sum_k w_k^(s))

    The per-draw improvement ratio (the H_1 basis) comes from
    :func:`get_endpoints_per_draw` on the x-fit; the IS weights come from
    ``eval_loglik``. Each z_s participant is a resampled x participant
    (``zi['src_pid']``), so its latent factor exists in ``theta_k``.

    Parameters
    ----------
    xi : pd.DataFrame
        Interim cohort the x-posterior was fit on.
    zi : pd.DataFrame
        Future-data block from :func:`get_interim_z_from_ypredi` (must carry the
        ``src_pid`` column + ``ypred_0 .. ypred_{pps_z_total-1}``).
    dit : pd.DataFrame
        Item metadata.
    draws : arviz.InferenceData, optional
        In-memory x-fit posterior (``fit['draws']``). Provide this or ``draws_file``.
    draws_file : str, optional
        Path to the x-fit zarr. Provide this or ``draws``.
    pps_z_total, pps_H1_def, pps_ProbH1_thresh, categorical_threshold, x_formula
        As in :func:`fit_interim_MC_of_posterior_xz`.
    eval_loglik : callable, default ``eval_loglik_partial_credit_model_ncats``
        ``(stan_data, params) -> array[N_total]`` pointwise log-likelihood used
        for the IS weights (swap for the credit / ordered-logit equivalents).
    output_file_prefix : str, optional
        If given and ``save_to_file``, writes ``{prefix}_p_h1_xz_IS.pkl`` and
        ``{prefix}_is_perf.csv``.

    Returns
    -------
    (p_h1_xz, is_perf) : tuple of pd.DataFrame
        - ``p_h1_xz``: one row per (item, sample s) with ``p_h1_xz``, ``s`` and
          the recorded ``pps_H1_def`` / ``pps_ProbH1_thresh`` / ``S``.
        - ``is_perf``: one row per sample s with IS diagnostics ``N`` (number of
          particles), ``ess`` (effective sample size = 1/sum(w^2)),
          ``ess_over_n`` (ESS / N) and ``ew2`` (the second-order weight moment
          E(w^2) = mean(w^2) — the weight-degeneracy diagnostic).
    """
    if draws is None and draws_file is None:
        raise ValueError("Provide either draws or draws_file.")
    if draws is None:
        draws = az.from_zarr(draws_file)
    if 'src_pid' not in zi.columns:
        raise ValueError("zi must carry 'src_pid' (use get_interim_z_from_ypredi).")
    vprint = print if verbose else (lambda *args, **kwargs: None)

    # Per-draw H_1 ratios from the x-posterior; 'draw' aligns with theta_k below
    # (both flatten the (chain, draw) axes in C order).
    x_ratio = get_endpoints_per_draw(
        dcati=xi, dit=dit, draws=draws,
        categorical_threshold=categorical_threshold,
        endpoint_type='items', param_name='ordered_prob_by_cat_qu_fit',
        verbose=verbose,
    )
    x_ratio = x_ratio.assign(ind=(x_ratio['ratio'] > pps_H1_def).astype(float))

    # Stack the posterior draws each loglik back-end may consume; include only
    # the params present so partial-credit/credit (skill_thresholds) and
    # ordered-logit (skill_thresholds_1 / skill_thresholds_incs) all work.
    post = draws.posterior
    theta_names = (
        'latent_factor_unit', 'latent_factor_beta',
        'skill_thresholds', 'skill_thresholds_1', 'skill_thresholds_incs',
        'loadings_questions_m1',
    )
    theta = {
        name: jnp.asarray(np.asarray(post[name].values).reshape(-1, *post[name].shape[2:]))
        for name in theta_names if name in post
    }
    n_draw = int(theta['latent_factor_beta'].shape[0])

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
        z_stan = _fit_partial_credit_make_stan_data(
            dit=dit, dcati=z_dcati, x_formula=x_formula, verbose=False,
        )

        # log p(z_s | theta_k) for every draw via vmap over the stacked params.
        loglik = jax.vmap(lambda p: eval_loglik(z_stan, p))(theta)  # (K, N_total_z)
        # Normalised weights via the numerically-stable softmax (= exp(logw -
        # logsumexp(logw))); keeps the numerator in log-space throughout.
        w = np.asarray(jax.nn.softmax(loglik.sum(axis=1)))  # (K,) sums to 1

        # IS performance diagnostics for this sample. ess = 1/sum(w^2);
        # ew2 = E(w^2) = mean(w^2), the second-order moment of the weights.
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


def _interim_make_x_stan(xi: pd.DataFrame, dit: pd.DataFrame, x_formula: str) -> dict:
    """stan_data for the x cohort (the full p(x | theta) factor)."""
    x_dcati = xi.copy()
    if 'y_stan' not in x_dcati:
        x_dcati['y_stan'] = x_dcati['y'] + 1
    x_dcati = x_dcati.sort_values(
        ['item_type_id', 'pid', 'time', 'item_label']
    ).reset_index(drop=True)
    x_dcati['oid'] = np.arange(1, len(x_dcati) + 1)
    x_dcati['oidt'] = x_dcati.groupby('item_type').cumcount() + 1
    return _fit_partial_credit_make_stan_data(
        dit=dit, dcati=x_dcati, x_formula=x_formula, verbose=False,
    )


def _interim_make_z_stan(zi: pd.DataFrame, dit: pd.DataFrame, s_idx: int, x_formula: str) -> dict:
    """stan_data for future-data sample s (mirrors the IS function)."""
    zcol = f'ypred_{s_idx}'
    z_dcati = zi.assign(pid=zi['src_pid'], y_stan=zi[zcol].astype(int))
    z_dcati['y'] = z_dcati['y_stan'] - 1
    z_dcati = z_dcati.sort_values(
        ['item_type_id', 'pid', 'time', 'item_label']
    ).reset_index(drop=True)
    z_dcati['oid'] = np.arange(1, len(z_dcati) + 1)
    z_dcati['oidt'] = z_dcati.groupby('item_type').cumcount() + 1
    return _fit_partial_credit_make_stan_data(
        dit=dit, dcati=z_dcati, x_formula=x_formula, verbose=False,
    )


def _ratio_per_draw_from_params(theta_batch, x_stan, xi, dit, categorical_threshold,
                                ordered_prob_eval=get_ordered_prob_of_partial_credit_model_ncats):
    """Per-(draw, item) improvement ratio for an arbitrary batch of parameter sets.

    Evaluates ``ordered_prob_by_cat_qu_fit`` for every draw in ``theta_batch``
    (dict of ``(K, ...)`` jnp arrays) on the x-cohort design ``x_stan``, wraps the
    result as an arviz posterior and runs :func:`get_endpoints_per_draw`. Lets
    reweighted (moment-matched) draws or moved (SMC) particles be scored for
    p(H_1 | x, z) exactly like a refit, rather than reusing the frozen x-ratio.
    """
    ordprob = np.asarray(jax.vmap(lambda p: ordered_prob_eval(x_stan, p))(theta_batch))  # (K, L)
    idata = az.from_dict(posterior={'ordered_prob_by_cat_qu_fit': ordprob[None, ...]})
    return get_endpoints_per_draw(
        dcati=xi, dit=dit, draws=idata,
        categorical_threshold=categorical_threshold,
        endpoint_type='items', param_name='ordered_prob_by_cat_qu_fit',
        verbose=False,
    )


def fit_interim_IS_moment_matching_of_posterior_xz_from_x(
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
    """
    Moment-matching importance sampling (Paananen et al. 2021, mean-match step)
    for p(H_1 | x, z) at a fixed x (Case A), as a *mitigation experiment* for the
    weight degeneracy of
    :func:`fit_interim_importance_sampling_of_posterior_xz_from_x`.

    Method. Work in a "match space" ``u`` where the ``positive_params`` are
    log-transformed and everything else is identity. Approximate the proposal
    p(theta | x) by a diagonal Gaussian ``q(u) = N(u; mu_p, diag sd_p^2)`` fit to
    the x-posterior draws (this Gaussian approximation is intrinsic to applying
    MMIS without the exact variational density). The match-space weights are

        log w(u) = logprior(theta(u)) + loglik_x(theta(u)) + loglik_z(theta(u))
                   + sum(u_positive)          # Jacobian of the log-transform
                   - logq(u)

    The mean-match step shifts every draw by ``(mu_weighted - mu_proposal)`` and
    recomputes the weights. **Finding:** when the base weights already collapse
    onto a single draw (early interims), ``mu_weighted`` *is* that draw, so the
    shift only piles the particles onto it and ESS does not recover.

    Parameters
    ----------
    fitting_method_args : dict
        Method tuning (required; specify at the call site). Defaults are applied
        for any missing keys:
        ``x_formula`` (design formula, default ``"~ time - 1"``),
        ``eval_loglik`` (pointwise log-lik callable, default
        ``eval_loglik_partial_credit_model_ncats``),
        ``logprior`` (callable ``params -> scalar``; default the partial credit
        prior ``get_prior_of_partial_credit_model_ncats``),
        ``positive_params`` (params log-transformed into match space, default
        ``('loadings_questions_m1',)``).

    Returns
    -------
    (p_h1_xz, mm) : tuple of pd.DataFrame
        - ``p_h1_xz``: one row per (item, sample s) with the mean-match estimate
          ``p_h1_xz`` = sum_k w1_k 1[ratio(theta*_k) > pps_H1_def] / sum_k w1_k,
          where ``theta*`` are the shifted draws (so the improvement ratio is
          re-evaluated on the shifted params, not the frozen x-ratio), plus
          ``s`` / ``pps_H1_def`` / ``pps_ProbH1_thresh`` / ``S``.
        - ``mm``: one row per s with ``N``, ``ess_over_n_base``,
          ``ess_over_n_meanmatch``, ``basew_vs_exact_corr`` (Pearson corr of the
          diagonal-Gaussian base log-weights against the exact IS log-weights
          ``loglik_z`` — a check that the Gaussian proposal approximation is sound)
          and ``mins`` (wall-clock minutes for the whole call). ESS and E(w^2)
          derive from ``ess_over_n_* `` and ``N``: ess = ess_over_n*N,
          E(w^2) = 1/(N^2 * ess_over_n).
    """
    if draws is None and draws_file is None:
        raise ValueError("Provide either draws or draws_file.")
    if draws is None:
        draws = az.from_zarr(draws_file)
    if 'src_pid' not in zi.columns:
        raise ValueError("zi must carry 'src_pid' (use get_interim_z_from_ypredi).")
    vprint = print if verbose else (lambda *args, **kwargs: None)

    fitting_method_args = {
        'x_formula': "~ time - 1",
        'eval_loglik': eval_loglik_partial_credit_model_ncats,
        'logprior': get_prior_of_partial_credit_model_ncats,
        'positive_params': ('loadings_questions_m1',),
        **fitting_method_args,
    }
    x_formula = fitting_method_args['x_formula']
    eval_loglik = fitting_method_args['eval_loglik']
    logprior = fitting_method_args['logprior']
    positive_params = fitting_method_args['positive_params']

    theta = _stack_posterior_theta(draws)
    n_draw = int(theta['latent_factor_beta'].shape[0])
    x_stan = _interim_make_x_stan(xi, dit, x_formula)

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

    def _make_logw(z_stan):
        def logw(u):
            pr = _u_to_params(u)
            jac = u[pos_cols].sum() if pos_cols.size else 0.0
            lt = (logprior(pr)
                  + eval_loglik(x_stan, pr).sum()
                  + eval_loglik(z_stan, pr).sum()
                  + jac)
            return lt - jstats.norm.logpdf(u, mu_p_j, sd_p_j).sum()
        return jax.jit(jax.vmap(logw))

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
        z_stan = _interim_make_z_stan(zi, dit, s_idx, x_formula)
        logw_fn = _make_logw(z_stan)

        logw0 = np.asarray(logw_fn(Theta_u_j))
        w0 = softmax(logw0)
        ess0 = float(1.0 / (n_draw * np.sum(w0 ** 2)))
        llz_exact = np.asarray(jax.vmap(lambda p: eval_loglik(z_stan, p).sum())(theta))
        corr = float(np.corrcoef(logw0, llz_exact)[0, 1])

        mu_w = (w0[:, None] * Theta_u).sum(axis=0)
        Theta_u_star = Theta_u + (mu_w - mu_p)[None, :]
        logw1 = np.asarray(logw_fn(jnp.asarray(Theta_u_star)))
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
            x_stan, xi, dit, categorical_threshold,
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


def fit_interim_SMC_resample_of_posterior_xz_from_x(
    xi: pd.DataFrame,
    zi: pd.DataFrame,
    dit: pd.DataFrame,
    fitting_method_args: dict,
    draws=None,
    draws_file: Optional[str] = None,
    output_file_prefix: Optional[str] = None,
    save_to_file: bool = True,
    verbose: bool = True,
):
    """
    SMC sampler with resample-move for p(theta | x, z_s) at a fixed x (Case A),
    the mitigation that actually crosses a large proposal/target gap.

    Bridges ``pi_beta ∝ p(theta|x) p(z_s|theta)^beta`` from beta=0 (the proposal,
    = the x-posterior draws) to beta=1 (the target). Each tempering step:

      1. adaptively pick the next beta so the tempering ESS ≈ ``ess_frac_target``
         * ``n_particles`` (bisection on the incremental weights),
      2. systematic-resample the particles by the incremental weights,
      3. MOVE them with ``n_move_steps`` MALA (Metropolis-adjusted Langevin) steps
         invariant to ``pi_beta`` — the relocation that plain reweighting (and
         moment matching) cannot do. The step size adapts toward ``target_accept``.

    Compile-once kernel. The move is a hand-rolled vectorised MALA in pure JAX
    with ``beta`` (and the step size) as *traced* arguments, so the gradient/step
    function compiles once and is reused for every temperature (numpyro's
    multi-chain NUTS instead vmaps model args inconsistently between init and
    sample, which a per-temperature beta cannot satisfy). The annealed z-term is
    evaluated through ``eval_loglik_annealed`` with the temperature in ``params``.

    Notes. Specialised to the partial credit ncats model (priors hard-coded). The
    move works in a partly-unconstrained space (loadings via ``log``); the
    ``latent_factor_unit`` prior uses the plain ``Normal(0, 1/sqrt(1-1/U))``
    ZeroSumNormal marginal, dropping only the (here immaterial) sum-to-zero
    identifiability.

    Parameters
    ----------
    fitting_method_args : dict
        Method tuning (required; specify at the call site). Defaults are applied
        for any missing keys: ``s_idx`` (future-data
        sample, default 0), ``n_particles`` (default 128), ``ess_frac_target``
        (default 0.5), ``n_move_steps`` (MALA steps per temperature, default 20),
        ``init_step_size`` (default 0.02), ``target_accept`` (default 0.574),
        ``max_temps`` (default 150), ``x_formula`` (default ``"~ time - 1"``),
        ``eval_loglik`` (default ``eval_loglik_partial_credit_model_ncats``),
        ``eval_loglik_annealed`` (default
        ``eval_loglik_partial_credit_model_ncats_with_annealing``),
        ``seed`` (default 123).

    Returns
    -------
    (schedule, particles) : (pd.DataFrame, dict)
        - ``schedule``: one row per tempering step with ``temp``, ``beta``,
          ``d_beta``, ``ess_frac_temper`` and ``move_secs``.
        - ``particles``: the final move-step particles (dict of jnp arrays,
          shape ``(n_particles, *event)``), approximately ~ pi_beta at the
          last beta reached.
    """
    if draws is None and draws_file is None:
        raise ValueError("Provide either draws or draws_file.")
    if draws is None:
        draws = az.from_zarr(draws_file)
    if 'src_pid' not in zi.columns:
        raise ValueError("zi must carry 'src_pid' (use get_interim_z_from_ypredi).")
    vprint = print if verbose else (lambda *args, **kwargs: None)

    fitting_method_args = {
        's_idx': 0,
        'n_particles': 128,
        'ess_frac_target': 0.5,
        'n_move_steps': 20,
        'init_step_size': 0.02,
        'target_accept': 0.574,
        'max_temps': 150,
        'x_formula': "~ time - 1",
        'eval_loglik': eval_loglik_partial_credit_model_ncats,
        'eval_loglik_annealed': eval_loglik_partial_credit_model_ncats_with_annealing,
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
    x_formula = fitting_method_args['x_formula']
    eval_loglik = fitting_method_args['eval_loglik']
    eval_loglik_annealed = fitting_method_args['eval_loglik_annealed']
    seed = fitting_method_args['seed']

    theta = _stack_posterior_theta(draws)
    n_draw = int(theta['latent_factor_beta'].shape[0])
    U = int(theta['latent_factor_unit'].shape[1])
    P = int(theta['latent_factor_beta'].shape[1])
    L = int(theta['skill_thresholds'].shape[1])
    Ld = int(theta['loadings_questions_m1'].shape[1])

    x_stan = _interim_make_x_stan(xi, dit, x_formula)
    z_stan = _interim_make_z_stan(zi, dit, s_idx, x_formula)
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
        return (get_prior_of_partial_credit_model_ncats(pr)
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
    schedule, particles = fit_interim_SMC_resample_of_posterior_xz_from_x(
        xi=xi, zi=zi, dit=dit, draws_file=draws_file,
        fitting_method_args=fma, save_to_file=False, verbose=verbose,
    )
    x_stan = _interim_make_x_stan(xi, dit, fma.get('x_formula', "~ time - 1"))
    ratio = _ratio_per_draw_from_params(particles, x_stan, xi, dit, categorical_threshold)
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


def _interim_make_z_stan_batched(zi: pd.DataFrame, dit: pd.DataFrame, pps_z_total: int, x_formula: str):
    """One z stan_data (shared design) + a ``(S, N_total)`` matrix ``Y`` of the
    per-sample outcomes, all in the same observation order. Across future samples
    s the unit / question / design are identical (same resampled participants and
    items); only the outcome ``ypred_s`` differs, so the SMC move's per-obs linear
    predictor ``eta`` is shared and only the categorical target changes with s."""
    z_dcati = zi.assign(pid=zi['src_pid'])
    z_dcati = z_dcati.sort_values(
        ['item_type_id', 'pid', 'time', 'item_label']
    ).reset_index(drop=True)
    z_dcati['oid'] = np.arange(1, len(z_dcati) + 1)
    z_dcati['oidt'] = z_dcati.groupby('item_type').cumcount() + 1
    z0 = z_dcati.assign(y_stan=z_dcati['ypred_0'].astype(int))
    z0['y'] = z0['y_stan'] - 1
    z_stan = _fit_partial_credit_make_stan_data(dit=dit, dcati=z0, x_formula=x_formula, verbose=False)
    Y = np.stack(
        [z_dcati[f'ypred_{s}'].astype(int).to_numpy() for s in range(pps_z_total)], axis=0
    )  # (S, N_total) of y_stan (1-based, matching stan_data['y'])
    return z_stan, jnp.asarray(Y)


def fit_interim_SMC_batched_PPS_of_posterior_xz_from_x(
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
    """
    Batched SMC resample-move PPS over all S future datasets at once (Case A).

    Vectorised alternative to :func:`fit_interim_SMC_PPS_of_posterior_xz_from_x`:
    instead of S independent SMC runs (each its own adaptive schedule + JIT
    compile), all S future samples share **one** tempering schedule and are moved
    together as a single ``(S, K, D)`` particle tensor under a vmapped MALA kernel
    that compiles once. The shared schedule advances by the most conservative
    sample (``min_s`` tempering ESS hits the target), so easy samples are
    over-tempered slightly; in exchange there is no per-sample recompile and the
    move is a single fused kernel.

    Particles live in the unconstrained match space
    ``u = [latent_factor_unit | latent_factor_beta | skill_thresholds | log loadings]``;
    the move target for sample s is ``log p(theta|x) + beta * log p(z_s|theta)``
    (prior + full x-likelihood + tempered z-likelihood, with the log-Jacobian).
    Final per-sample particles are scored for p(H_1 | x, z_s) as the fraction with
    item improvement ratio > pps_H1_def (uniform post-move weights).

    Parameters
    ----------
    fitting_method_args : dict
        Tuning (defaults applied for missing keys): ``n_particles`` (K, default
        128), ``ess_frac_target`` (default 0.5), ``n_move_steps`` (MALA steps per
        temperature, default 20), ``init_step_size`` (default 0.02),
        ``target_accept`` (default 0.574), ``max_temps`` (default 300),
        ``x_formula`` (default ``"~ time - 1"``), ``seed`` (default 123).

    Returns
    -------
    (p_h1_xz, schedule) : tuple of pd.DataFrame
        - ``p_h1_xz``: one row per (item, sample s) with ``p_h1_xz``, ``s`` and the
          recorded ``pps_H1_def`` / ``pps_ProbH1_thresh`` / ``S``.
        - ``schedule``: one row per shared tempering step with ``temp``, ``beta``,
          ``d_beta``, ``ess_frac_min`` / ``ess_frac_mean`` (across samples),
          ``accept``, ``step_size`` and ``move_secs``; plus ``n_particles`` and
          ``mins_total`` columns.
    """
    if draws is None and draws_file is None:
        raise ValueError("Provide either draws or draws_file.")
    if draws is None:
        draws = az.from_zarr(draws_file)
    if 'src_pid' not in zi.columns:
        raise ValueError("zi must carry 'src_pid' (use get_interim_z_from_ypredi).")
    vprint = print if verbose else (lambda *args, **kwargs: None)

    fma = {
        'n_particles': 128, 'ess_frac_target': 0.5, 'n_move_steps': 20,
        'init_step_size': 0.02, 'target_accept': 0.574, 'max_temps': 300,
        'x_formula': "~ time - 1", 'seed': 123,
        **fitting_method_args,
    }
    K = fma['n_particles']
    ess_frac_target = fma['ess_frac_target']
    n_move_steps = fma['n_move_steps']
    target_accept = fma['target_accept']
    max_temps = fma['max_temps']
    x_formula = fma['x_formula']
    seed = fma['seed']
    S = pps_z_total

    theta = _stack_posterior_theta(draws)
    n_draw = int(theta['latent_factor_beta'].shape[0])
    U = int(theta['latent_factor_unit'].shape[1])
    P = int(theta['latent_factor_beta'].shape[1])
    L = int(theta['skill_thresholds'].shape[1])
    Ld = int(theta['loadings_questions_m1'].shape[1])

    x_stan = _interim_make_x_stan(xi, dit, x_formula)
    z_stan, Y = _interim_make_z_stan_batched(zi, dit, S, x_formula)

    sl_u, sl_b = slice(0, U), slice(U, U + P)
    sl_s, sl_l = slice(U + P, U + P + L), slice(U + P + L, U + P + L + Ld)

    def _u_to_params(u):
        return {
            'latent_factor_unit': u[sl_u],
            'latent_factor_beta': u[sl_b],
            'skill_thresholds': u[sl_s],
            'loadings_questions_m1': jnp.exp(u[sl_l]),
        }

    def _logbase(u):
        pr = _u_to_params(u)
        return (get_prior_of_partial_credit_model_ncats(pr)
                + eval_loglik_partial_credit_model_ncats(x_stan, pr).sum()
                + u[sl_l].sum())          # Jacobian of the log-transform

    def _llz(u, y):
        pr = _u_to_params(u)
        return eval_loglik_partial_credit_model_ncats({**z_stan, 'y': y}, pr).sum()

    def _logpost(u, y, beta):
        return _logbase(u) + beta * _llz(u, y)

    _grad = jax.grad(_logpost, argnums=0)

    # Batched over (S samples, K particles); y is shared within a sample.
    def _b_llz(u, Y):
        return jax.vmap(jax.vmap(_llz, (0, None)), (0, 0))(u, Y)            # (S, K)

    def _b_logpost(u, Y, beta):
        return jax.vmap(jax.vmap(_logpost, (0, None, None)), (0, 0, None))(u, Y, beta)

    def _b_grad(u, Y, beta):
        return jax.vmap(jax.vmap(_grad, (0, None, None)), (0, 0, None))(u, Y, beta)

    def _mala_sweep(u, Y, beta, eps, key):
        def step(carry, k):
            u, n_acc = carry
            g = _b_grad(u, Y, beta)
            k1, k2 = random.split(k)
            prop = u + 0.5 * eps ** 2 * g + eps * random.normal(k1, u.shape)
            gp = _b_grad(prop, Y, beta)
            lp_u = _b_logpost(u, Y, beta)
            lp_p = _b_logpost(prop, Y, beta)
            logq_fwd = -jnp.sum((prop - u - 0.5 * eps ** 2 * g) ** 2, axis=-1) / (2 * eps ** 2)
            logq_bwd = -jnp.sum((u - prop - 0.5 * eps ** 2 * gp) ** 2, axis=-1) / (2 * eps ** 2)
            log_acc = lp_p - lp_u + logq_bwd - logq_fwd
            accept = jnp.log(random.uniform(k2, lp_u.shape)) < log_acc
            u = jnp.where(accept[..., None], prop, u)
            return (u, n_acc + accept.mean()), None
        (u, n_acc), _ = jax.lax.scan(step, (u, 0.0), random.split(key, n_move_steps))
        return u, n_acc / n_move_steps

    _mala_jit = jax.jit(_mala_sweep)

    def _systematic_resample(w, key):
        positions = (random.uniform(key) + np.arange(len(w))) / len(w)
        return np.searchsorted(np.cumsum(w), positions)

    # Init: same K base draws for every sample, broadcast to (S, K, D).
    rng = np.random.default_rng(seed)
    idx0 = rng.choice(n_draw, size=K, replace=False)
    sub = {k: np.asarray(theta[k])[idx0] for k in
           ('latent_factor_unit', 'latent_factor_beta', 'skill_thresholds', 'loadings_questions_m1')}
    u_single = np.concatenate([
        sub['latent_factor_unit'], sub['latent_factor_beta'],
        sub['skill_thresholds'], np.log(sub['loadings_questions_m1']),
    ], axis=1)                                                              # (K, D)
    u = jnp.asarray(np.broadcast_to(u_single[None], (S, *u_single.shape)).copy())  # (S, K, D)
    llz = np.asarray(_b_llz(u, Y))                                          # (S, K)
    llz = np.where(np.isfinite(llz), llz, -1e30)

    key = random.PRNGKey(seed)
    eps = float(fma['init_step_size'])
    beta = 0.0
    rows = []
    t_all = time.time()
    vprint(f"[SMC-batched] all S={S} samples, K={K} particles, one shared schedule")
    for t in range(max_temps):
        # Shared d_beta: most conservative sample (min ESS) hits the target.
        def worst_ess(db):
            wn = softmax(db * llz, axis=1)
            return float(np.min(1.0 / (K * np.sum(wn ** 2, axis=1))))
        lo, hi = 0.0, 1.0 - beta
        if worst_ess(hi) >= ess_frac_target:
            db = hi
        else:
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                if worst_ess(mid) >= ess_frac_target:
                    lo = mid
                else:
                    hi = mid
            db = lo
        beta_new = beta + db
        wn = softmax(db * llz, axis=1)                                     # (S, K)
        ess_s = 1.0 / (K * np.sum(wn ** 2, axis=1))                        # (S,)

        # Per-sample systematic resampling.
        key, ksub = random.split(key)
        sub_keys = random.split(ksub, S)
        anc = np.stack([_systematic_resample(wn[s], sub_keys[s]) for s in range(S)])  # (S, K)
        u = jnp.asarray(np.asarray(u)[np.arange(S)[:, None], anc])

        t0 = time.time()
        key, kmove = random.split(key)
        u, acc = _mala_jit(u, Y, jnp.asarray(beta_new), jnp.asarray(eps), kmove)
        acc = float(acc)
        move_secs = time.time() - t0
        eps = float(np.clip(eps * np.exp(0.5 * (acc - target_accept)), 1e-4, 1.0))
        llz = np.asarray(_b_llz(u, Y))
        llz = np.where(np.isfinite(llz), llz, -1e30)

        rows.append({
            'temp': t + 1, 'beta': round(beta_new, 5), 'd_beta': round(db, 5),
            'ess_frac_min': round(float(ess_s.min()), 3),
            'ess_frac_mean': round(float(ess_s.mean()), 3),
            'accept': round(acc, 3), 'step_size': round(eps, 5),
            'move_secs': round(move_secs, 2),
        })
        vprint(f"  temp {t+1}: beta {beta:.4f}->{beta_new:.4f} (dbeta={db:.4f}) "
               f"| ESS/K min={ess_s.min():.2f} mean={ess_s.mean():.2f} "
               f"| accept={acc:.2f} | move {move_secs:.1f}s")
        beta = beta_new
        if beta >= 1.0 - 1e-9:
            break

    mins_total = (time.time() - t_all) / 60.0
    schedule = pd.DataFrame(rows)
    schedule['n_particles'] = K
    schedule['mins_total'] = round(mins_total, 4)
    vprint(f"\n[SMC-batched] reached beta={beta:.4f} in {len(schedule)} steps, {mins_total:.2f} min total")

    # Per-sample labels from the moved particles (uniform weights).
    u_np = jnp.asarray(u)
    p_rows = []
    for s in range(S):
        particles = {
            'latent_factor_unit': u_np[s, :, sl_u],
            'latent_factor_beta': u_np[s, :, sl_b],
            'skill_thresholds': u_np[s, :, sl_s],
            'loadings_questions_m1': jnp.exp(u_np[s, :, sl_l]),
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
        p_h1_xz.to_pickle(f"{output_file_prefix}_p_h1_xz_SMCbatched.pkl")
        schedule.to_csv(f"{output_file_prefix}_smcbatched_schedule.csv", index=False)
        vprint(f"\nSaved batched-SMC P(H_1 | x, z) to: {output_file_prefix}_p_h1_xz_SMCbatched.pkl")
    return p_h1_xz, schedule
