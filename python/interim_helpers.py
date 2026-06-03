"""
Pure-pandas / arviz / numpy helpers shared across the interim-analysis
algorithms. These functions have no model-specific math: any generative
model whose draws zarr carries a ``posterior['ypred']`` of shape
``(chain, draw, N_total)`` can reuse them as-is.

Migrated out of ``fit_interim.py`` in step 1 of the OO-port refactor
(see ``dev/interim_model_blueprint_plan.md``).
"""

from typing import Optional

import arviz as az
import numpy as np
import pandas as pd


def _load_ypred(draws_file: str, pps_z_total: int, rng,
                keep_order: bool = False) -> np.ndarray:
    """
    Load the ``ypred`` posterior-predictive draws from a zarr, flatten the
    (chain, draw) dims, and return ``pps_z_total`` draws as a
    ``(pps_z_total, N_total)`` array.

    ``keep_order=False`` (default) randomly selects ``pps_z_total`` draws via
    ``rng.choice`` (used by the IS/MC scripts where a permuted subsample is
    fine). ``keep_order=True`` returns the first ``pps_z_total`` draws in
    posterior order (no shuffling), so row ``s`` of the returned matrix
    corresponds to posterior draw ``s`` -- required when downstream code merges
    on the posterior-draw index (e.g. against ``get_endpoints_per_draw`` output).
    """
    ypred = az.from_zarr(draws_file).posterior['ypred'].values  # (chain, draw, N_total)
    ypred = ypred.reshape(-1, ypred.shape[-1])                  # (n_draw, N_total)
    if keep_order:
        return ypred[:pps_z_total]
    return ypred[rng.choice(ypred.shape[0], size=pps_z_total, replace=False)]


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
