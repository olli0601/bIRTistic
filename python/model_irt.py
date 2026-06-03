"""
:class:`IRTModel` -- mixin for the IRT family (partial credit, credit,
ordered logit). Owns the endpoint logic that was previously the free
functions in :mod:`get_endpoints` (now deleted from that module's
public surface): per-draw endpoint extraction and quantile summarisation
of posterior ordered-prob arrays.

Concrete IRT subclasses (:class:`model_pcm.PartialCreditModelNCats`,
:class:`model_credit.CreditModelNCats`,
:class:`model_ordered_logit.OrderedLogitNCats`) inline their own
``eval_loglik`` / fit drivers / ``make_stan_data`` etc. The two methods
shared across all three -- :meth:`get_endpoints_per_draw` and
:meth:`get_endpoints` -- live here because every IRT model records
``ordered_prob_by_cat_qu_fit`` (and ``..._pr``) in its posterior.
"""

from typing import Literal, Optional
import os

import arviz as az
import numpy as np
import pandas as pd

from model import Model
from utils import _map_cq_id_to_item_structure


class IRTModel(Model):
    """Marker mixin for the IRT family. Provides per-draw + summary
    endpoint extraction from posterior ordered-prob arrays."""

    # ------------------------------------------------------------------
    # Private helpers -- bodies verbatim from the legacy get_endpoints.py
    # ------------------------------------------------------------------

    @staticmethod
    def _resolve_draws(draws, draws_file: Optional[str], param_name: str,
                       verbose: bool = True) -> np.ndarray:
        """Return the ``(chain, draw, cq_id)`` array for ``param_name`` from
        either an in-memory ArviZ InferenceData (``draws``) or a zarr path
        (``draws_file``). Exactly one of the two must be provided."""
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

    def _make_po(self, po_arr: np.ndarray) -> pd.DataFrame:
        """Reshape a posterior probability array of shape ``(chain, draw,
        cq_id)`` into a long DataFrame with columns ``['.draw', 'cq_id',
        'prob']`` and merge in the item structure (``item_type_id``,
        ``item_time_id``, ``y``) recovered from ``cq_id`` via
        :func:`utils._map_cq_id_to_item_structure`."""
        dp1 = self.dcati
        dit = self.dit
        po = po_arr.reshape(-1, po_arr.shape[-1])
        po = pd.DataFrame.from_records(po)
        po = (
            po.melt(var_name='cq_id', value_name='prob', ignore_index=False)
            .reset_index()
            .rename(columns={'index': '.draw'})
        )
        cq_map = _map_cq_id_to_item_structure(dp1, dit)
        return po.merge(cq_map, on='cq_id')

    def _endpoints_per_draw_impl(
        self,
        po: pd.DataFrame,
        categorical_threshold: int = 3,
        endpoint_type: Literal["items", "item_groups"] = "items",
    ) -> pd.DataFrame:
        """Compute per-draw endpoint scalars + directional ``diff`` and
        ``ratio`` for both categorical and out-of-7 items. Body verbatim
        from legacy ``_get_endpoints_per_draw``; reads ``self.dcati``
        (``dp1``) and ``self.dit``."""
        dp1 = self.dcati
        dit = self.dit
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

    # ------------------------------------------------------------------
    # Public endpoint methods
    # ------------------------------------------------------------------

    def get_endpoints_per_draw(
        self,
        draws=None,
        draws_file: Optional[str] = None,
        categorical_threshold: int = 3,
        endpoint_type: Literal["items", "item_groups"] = "items",
        param_name: str = "ordered_prob_by_cat_qu_pr",
        verbose: bool = False,
    ) -> pd.DataFrame:
        """Per-draw directional ``diff`` and ``ratio`` per item or
        item-group, computed from ``self.dcati`` + ``self.dit`` and the
        supplied posterior."""
        if endpoint_type not in ("items", "item_groups"):
            raise ValueError("endpoint_type must be either 'items' or 'item_groups'")
        if param_name not in ("ordered_prob_by_cat_qu_pr", "ordered_prob_by_cat_qu_fit"):
            raise ValueError(
                "param_name must be either 'ordered_prob_by_cat_qu_pr' or 'ordered_prob_by_cat_qu_fit'"
            )
        po_arr = self._resolve_draws(draws, draws_file, param_name, verbose=verbose)
        po = self._make_po(po_arr)
        po = self._endpoints_per_draw_impl(
            po,
            categorical_threshold=categorical_threshold,
            endpoint_type=endpoint_type,
        )
        po.rename(columns={'.draw': 'draw'}, inplace=True)
        return po

    def get_endpoints(
        self,
        draws=None,
        draws_file: Optional[str] = None,
        categorical_threshold: int = 3,
        endpoint_type: Literal["items", "item_groups"] = "items",
        param_name: str = "ordered_prob_by_cat_qu_pr",
        verbose: bool = True,
    ) -> pd.DataFrame:
        """Endpoint summaries with quantile aggregation across draws
        (2.5%, 25%, 50%, 75%, 97.5%). Reads ``self.dcati`` + ``self.dit``
        and the supplied posterior."""
        vprint = print if verbose else (lambda *args, **kwargs: None)
        if endpoint_type not in ("items", "item_groups"):
            raise ValueError("endpoint_type must be either 'items' or 'item_groups'")
        if param_name not in ("ordered_prob_by_cat_qu_pr", "ordered_prob_by_cat_qu_fit"):
            raise ValueError(
                "param_name must be either 'ordered_prob_by_cat_qu_pr' or 'ordered_prob_by_cat_qu_fit'"
            )
        dit = self.dit

        vprint("Computing per-draw endpoints...")
        po_arr = self._resolve_draws(draws, draws_file, param_name, verbose=verbose)
        po = self._make_po(po_arr)
        po = self._endpoints_per_draw_impl(
            po,
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
