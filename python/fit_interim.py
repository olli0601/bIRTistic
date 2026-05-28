"""
Interim-analysis helpers for predictive probability of success (PPS).

- ``get_interim_z``: build the 'missing' future data block z from the interim
  cohort xi, the final cohort xf, and posterior-predictive draws ypred.
- ``fit_interim_MC_of_posterior_xz``: Monte-Carlo estimate of p(H_1 | x, z_s)
  across S hypothetical future datasets, refitting the partial credit model to
  (xi + z_s) for each s.
"""

from typing import Optional

import arviz as az
import numpy as np
import pandas as pd

from get_endpoints import get_endpoints_per_draw
from fit_partial_credit_model import fit_partial_credit_model_ncats_pyrosvi


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
    zi = tmp.merge(xi.rename(columns={'pid': 'src_pid'}), on='src_pid').drop(columns='src_pid')
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
    pps_z_total: int = 10,
    pps_H1_def: float = 0.5,
    pps_ProbH1_thresh: float = 0.89,
    categorical_threshold: int = 3,
    fitting_method=fit_partial_credit_model_ncats_pyrosvi,
    fitting_method_args: Optional[dict] = None,
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
    fitting_method_args : dict, optional
        Extra keyword arguments forwarded to ``fitting_method`` (e.g.
        ``algorithm``, ``x_formula``, ``lr``, ``num_steps``, ``output_samples``).
        Must not include any of the keys this function sets itself.
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
    fitting_method_args = fitting_method_args or {}

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
