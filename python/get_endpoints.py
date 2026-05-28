"""
Analysis functions for Bayesian IRT models.

This module provides functions for analyzing posterior draws from IRT models,
specifically endpoint extraction and effect size computation.
"""

from typing import Literal, Optional
import os
import pandas as pd
import numpy as np
import arviz as az

from utils import _map_cq_id_to_item_structure


def _resolve_draws(draws, draws_file: Optional[str], param_name: str, verbose: bool = True) -> np.ndarray:
    """
    Return the ``(chain, draw, cq_id)`` array for ``param_name`` from either an
    in-memory ArviZ InferenceData (``draws``) or a zarr path (``draws_file``).
    Exactly one of the two must be provided.
    """
    vprint = print if verbose else (lambda *args, **kwargs: None)
    if draws is None and draws_file is None:
        raise ValueError("Provide either draws or draws_file.")
    if draws is None:
        if not os.path.exists(draws_file):
            raise FileNotFoundError(f"Draws file not found: {draws_file}")
        vprint(f"Loading draws from: {draws_file}")
        draws = az.from_zarr(draws_file)
    if param_name not in draws.posterior.data_vars:
        raise ValueError(
            f"No data_vars '{param_name}'. Available: {list(draws.posterior.data_vars)}"
        )
    return draws.posterior[param_name].values


def _make_po(po_arr: np.ndarray, dp1: pd.DataFrame, dit: pd.DataFrame) -> pd.DataFrame:
    """
    Reshape a posterior probability array of shape ``(chain, draw, cq_id)`` into
    a long DataFrame with columns ``['.draw', 'cq_id', 'prob']`` and merge in
    the item structure (``item_type_id``, ``item_time_id``, ``y``) recovered
    from ``cq_id`` via :func:`_map_cq_id_to_item_structure`.
    """
    po = po_arr.reshape(-1, po_arr.shape[-1])
    po = pd.DataFrame.from_records(po)
    po = (
        po.melt(var_name='cq_id', value_name='prob', ignore_index=False)
        .reset_index()
        .rename(columns={'index': '.draw'})
    )
    cq_map = _map_cq_id_to_item_structure(dp1, dit)
    return po.merge(cq_map, on='cq_id')


def _get_endpoints_per_draw(
    po: pd.DataFrame,
    dp1: pd.DataFrame,
    dit: pd.DataFrame,
    categorical_threshold: int = 3,
    endpoint_type: Literal["items", "item_groups"] = "items",
) -> pd.DataFrame:
    """
    Compute per-draw endpoint scalars + directional ``diff`` and ``ratio``
    for both categorical and out-of-7 items.

    Parameters
    ----------
    po : pd.DataFrame
        Long-form posterior probabilities as returned by :func:`_make_po`.
        Required columns: ``.draw``, ``cq_id``, ``prob``, ``item_type_id``,
        ``item_time_id``, ``y``.
    dp1 : pd.DataFrame
        Pre-processed data table; provides (``item_type_id``, ``item_time_id``)
        → (``item_label``, ``item_type``, ``time_label``) mapping.
    dit : pd.DataFrame
        Item metadata; provides ``group_label`` and ``item_high_label``.
    categorical_threshold : int, default 3
        For ``item_type == 'categorical'`` items the scalar per draw is
        ``sum_{y >= categorical_threshold} prob``.
    endpoint_type : {"items", "item_groups"}, default "items"
        ``items`` keeps per-item rows; ``item_groups`` averages per draw across
        items within ``group_label``.

    Returns
    -------
    pd.DataFrame
        po per-draw frame with one row per (``.draw``, item-or-group,
        ``time_label`` pair) and columns ``Baseline``, ``Endline``, ``diff``,
        ``ratio``. ``diff`` / ``ratio`` follow the sign convention from
        ``item_high_label``: improvement is positive.
    """
    parts = []
    for item_type in ('categorical', 'out-of-7'):
        tmp = dp1[dp1['item_type'] == item_type][
            ['item_type_id', 'item_label', 'item_time_id', 'item_type', 'time_label']
        ].drop_duplicates()
        sub = po.merge(tmp, on=['item_type_id', 'item_time_id'])
        if sub.empty:
            continue
        if item_type == 'categorical':
            sub = sub[sub['y'] >= categorical_threshold]
            sub = sub.assign(_w=sub['prob'])
        else:
            sub = sub.assign(_w=sub['y'] * sub['prob'])
        sub = sub.groupby(
            ['.draw', 'item_type_id', 'item_label', 'item_time_id', 'item_type', 'time_label']
        ).agg(value=('_w', 'sum')).reset_index()
        parts.append(sub)
    po = pd.concat(parts, ignore_index=True)

    po = po.merge(
        dit[['item_type', 'item_label', 'group_label']],
        on=['item_type', 'item_label'],
    )

    if endpoint_type == 'item_groups':
        id_vars = ['item_type', 'group_label']
        po = po.groupby(['.draw', 'item_type', 'time_label', 'group_label']).agg(
            value=('value', 'mean')
        ).reset_index()
    else:
        id_vars = ['item_type_id', 'item_type', 'item_label', 'group_label']

    po = po.pivot_table(
        index=['.draw'] + id_vars,
        columns='time_label',
        values='value',
    ).reset_index()
    po = po.dropna(subset=['Baseline', 'Endline'])
    po = po.merge(
        dit[['item_type', 'group_label', 'item_high_label']].drop_duplicates(),
        on=['item_type', 'group_label'],
    )

    # Compute differences and ratios (direction depends on item_high_label)
    po['diff'] = np.nan
    po['ratio'] = np.nan
    tmp = po['item_high_label'] == 'lower_is_better'
    po.loc[tmp, 'diff'] = po.loc[tmp, 'Baseline'] - po.loc[tmp, 'Endline']
    po.loc[tmp, 'ratio'] = 1 - po.loc[tmp, 'Endline'] / po.loc[tmp, 'Baseline']
    tmp = po['item_high_label'] == 'higher_is_better'
    po.loc[tmp, 'diff'] = po.loc[tmp, 'Endline'] - po.loc[tmp, 'Baseline']
    po.loc[tmp, 'ratio'] = po.loc[tmp, 'Endline'] / po.loc[tmp, 'Baseline'] - 1

    return po


def get_endpoints(
    dp1: pd.DataFrame,
    dit: pd.DataFrame,
    draws_file: Optional[str] = None,
    draws=None,
    categorical_threshold: int = 3,
    endpoint_type: Literal["items", "item_groups"] = "items",
    param_name: str = "ordered_prob_by_cat_qu_pr",
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Get endpoints from ordered logit model draws.

    Loads posterior probability draws, computes per-(draw, item-or-group)
    Baseline/Endline scalars + directional ``diff``/``ratio`` via
    :func:`_get_endpoints_per_draw`, then summarises across draws with
    quantile summaries (2.5%, 25%, 50%, 75%, 97.5%).

    Parameters
    ----------
    dp1 : pd.DataFrame
        Input data with columns: item_type_id, item_label, item_time_id, item_type, time_label, y.
    dit : pd.DataFrame
        Item definitions with columns: item_type_id, item_label, cat_length, item_type,
        group_label, item_label_short, group_label_long, item_high_label.
    draws_file : str, optional
        Path to zarr directory holding the posterior of an ordered_prob parameter.
        Provide this or ``draws``.
    draws : arviz.InferenceData, optional
        In-memory posterior (e.g. ``fit['draws']``) to use instead of reading
        ``draws_file``. Provide this or ``draws_file``.
    categorical_threshold : int, default 3
        Threshold for aggregating categorical responses. Categories >= this value
        are summed.
    endpoint_type : {"items", "item_groups"}, default "items"
        ``items`` keeps individual items; ``item_groups`` averages across items
        within ``group_label``.
    param_name : {"ordered_prob_by_cat_qu_pr", "ordered_prob_by_cat_qu_fit"}
        Which ordered probability parameter to load.

    Returns
    -------
    pd.DataFrame
        Endpoint summaries with columns:
        - item_type, group_label, group_label_long
        - item_label, item_label_short (if endpoint_type="items")
        - item_high_label
        - variable: one of Baseline, Endline, diff, ratio
        - q_lower, iqr_lower, median, iqr_upper, q_upper
    """
    vprint = print if verbose else (lambda *args, **kwargs: None)
    if endpoint_type not in ("items", "item_groups"):
        raise ValueError("endpoint_type must be either 'items' or 'item_groups'")
    if param_name not in ("ordered_prob_by_cat_qu_pr", "ordered_prob_by_cat_qu_fit"):
        raise ValueError(
            "param_name must be either 'ordered_prob_by_cat_qu_pr' or 'ordered_prob_by_cat_qu_fit'"
        )

    vprint("Computing per-draw endpoints...")
    po = _make_po(_resolve_draws(draws, draws_file, param_name, verbose=verbose), dp1, dit)
    po = _get_endpoints_per_draw(
        po, dp1, dit,
        categorical_threshold=categorical_threshold,
        endpoint_type=endpoint_type,
    )

    if endpoint_type == "item_groups":
        id_vars = ['item_type', 'group_label']
    else:
        id_vars = ['item_type_id', 'item_type', 'item_label', 'group_label']

    # Reshape to long format for summarisation
    po = po.melt(
        id_vars=['.draw'] + id_vars,
        value_vars=['diff', 'ratio', 'Baseline', 'Endline'],
        var_name='variable',
        value_name='value',
    )

    # Quantile summaries
    quantiles = [0.025, 0.25, 0.5, 0.75, 0.975]
    quantile_names = ['q_lower', 'iqr_lower', 'median', 'iqr_upper', 'q_upper']
    pos = po.groupby(id_vars + ['variable'])['value'].quantile(quantiles).unstack()
    pos.columns = quantile_names
    pos = pos.reset_index()

    # Merge with item metadata
    if endpoint_type == "item_groups":
        tmp = dit[['item_type', 'group_label', 'group_label_long', 'item_high_label']].drop_duplicates()
        pos = pos.merge(tmp, on=['item_type', 'group_label'])
    else:
        tmp = dit[
            ['item_type', 'item_label', 'item_label_short', 'group_label',
             'group_label_long', 'item_high_label']
        ].drop_duplicates()
        pos = pos.merge(tmp, on=['item_type', 'item_label', 'group_label'])

    vprint(f"Computed endpoints for {len(pos)} item-variable combinations")
    return pos


def get_endpoints_per_draw(
    dcati: pd.DataFrame,
    dit: pd.DataFrame,
    draws_file: Optional[str] = None,
    draws=None,
    categorical_threshold: int = 3,
    endpoint_type: Literal["items", "item_groups"] = "items",
    param_name: str = "ordered_prob_by_cat_qu_fit",
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Per-draw directional ``diff`` and ``ratio`` per item or item-group.

    Thin wrapper around :func:`_get_endpoints_per_draw`: resolves the posterior
    from ``draws`` or ``draws_file``, builds the long-form prob frame via
    :func:`_make_po`, then delegates aggregation.

    Parameters
    ----------
    dcati, dit : pd.DataFrame
        Pre-processed data and item metadata, same as :func:`get_endpoints`.
    draws_file : str, optional
        Path to the zarr holding the posterior. Provide this or ``draws``.
    draws : arviz.InferenceData, optional
        In-memory posterior (e.g. ``fit['draws']``). Provide this or ``draws_file``.
    categorical_threshold : int, default 3
    endpoint_type : {"items", "item_groups"}, default "items"
    param_name : {"ordered_prob_by_cat_qu_fit", "ordered_prob_by_cat_qu_pr"},
        default "ordered_prob_by_cat_qu_fit"
        Which ordered-prob parameter to read.

    Returns
    -------
    pd.DataFrame
        One row per (``draw``, item-or-group, time-pair) with columns
        ``draw``, the ``id_vars`` for the requested ``endpoint_type``,
        ``item_high_label``, ``Baseline``, ``Endline``, ``diff``, ``ratio``.
    """
    if endpoint_type not in ("items", "item_groups"):
        raise ValueError("endpoint_type must be either 'items' or 'item_groups'")
    if param_name not in ("ordered_prob_by_cat_qu_pr", "ordered_prob_by_cat_qu_fit"):
        raise ValueError(
            "param_name must be either 'ordered_prob_by_cat_qu_pr' or 'ordered_prob_by_cat_qu_fit'"
        )

    po = _make_po(_resolve_draws(draws, draws_file, param_name, verbose=verbose), dcati, dit)
    po = _get_endpoints_per_draw(
        po, dcati, dit,
        categorical_threshold=categorical_threshold,
        endpoint_type=endpoint_type,
    )
    po.rename(columns={'.draw': 'draw'}, inplace=True)
    return po
