"""
:class:`IRTModel` -- mixin for the IRT family (partial credit, credit,
ordered logit). Owns the endpoint logic that was previously the free
functions in :mod:`get_endpoints` (now deleted from that module's
public surface): per-draw endpoint extraction and quantile summarisation
of posterior ordered-prob arrays.

Concrete IRT subclasses (:class:`model_pcm.PartialCreditModel`,
:class:`model_credit.CreditModel`,
:class:`model_ordered_logit.OrderedLogit`) inline their own
``eval_loglik`` / fit drivers / ``make_stan_data`` etc. The two methods
shared across all three -- :meth:`get_endpoints_per_draw` and
:meth:`get_endpoints` -- live here because every IRT model records
``ordered_prob_by_cat_qu_fit`` (and ``..._pr``) in its posterior.
"""

from typing import Literal, Optional
import os

import arviz as az
import jax
import numpy as np
import pandas as pd

from model import Model
from utils import _map_cq_id_to_item_structure


class IRTModel(Model):
    """Marker mixin for the IRT family. Provides per-draw + summary
    endpoint extraction from posterior ordered-prob arrays."""

    def __init__(self, dit: pd.DataFrame, dcati: pd.DataFrame,
                 x_formula: str = "~ time - 1", *,
                 seed: int = 123, categorical_threshold: int = 3):
        self.categorical_threshold = int(categorical_threshold)
        super().__init__(dit=dit, dcati=dcati, x_formula=x_formula, seed=seed)

    # ------------------------------------------------------------------
    # IS / SMC scoring helpers (IRT-specific overrides of Model defaults)
    # ------------------------------------------------------------------

    def make_stan_data_from_xi(self) -> dict:
        """stan_data for the x cohort. Sorts by (item_type_id, pid, time,
        item_label) and re-indexes oid/oidt before delegating to
        :meth:`make_stan_data`."""
        x_dcati = self.dcati.copy()
        if 'y_stan' not in x_dcati:
            x_dcati['y_stan'] = x_dcati['y'] + 1
        x_dcati = x_dcati.sort_values(
            ['item_type_id', 'pid', 'time', 'item_label']
        ).reset_index(drop=True)
        x_dcati['oid'] = np.arange(1, len(x_dcati) + 1)
        x_dcati['oidt'] = x_dcati.groupby('item_type').cumcount() + 1
        return self.make_stan_data(x_dcati, self.x_formula)

    def make_stan_data_from_zi(self, zi: pd.DataFrame, s_idx: int) -> dict:
        """stan_data for the future-data sample ``s_idx`` drawn from
        ``zi``. Promotes ``src_pid -> pid`` and ``ypred_s -> y_stan``,
        rebuilds oid/oidt by the IRT (item_type_id, pid, time, item_label)
        order, then delegates to :meth:`make_stan_data`."""
        zcol = f'ypred_{s_idx}'
        z_dcati = zi.assign(pid=zi['src_pid'], y_stan=zi[zcol].astype(int))
        z_dcati['y'] = z_dcati['y_stan'] - 1
        z_dcati = z_dcati.sort_values(
            ['item_type_id', 'pid', 'time', 'item_label']
        ).reset_index(drop=True)
        z_dcati['oid'] = np.arange(1, len(z_dcati) + 1)
        z_dcati['oidt'] = z_dcati.groupby('item_type').cumcount() + 1
        return self.make_stan_data(z_dcati, self.x_formula)

    # ------------------------------------------------------------------
    # IRT data builder
    # ------------------------------------------------------------------

    @staticmethod
    def get_interim_data_x(df: pd.DataFrame,
                           interim_date: Optional[str] = None
                           ) -> pd.DataFrame:
        """Slice a long-form panel ``df`` to the cohort observed on/before
        ``interim_date``, keep only (participant, item) responses present
        at BOTH baseline and endline, then re-index pids / oids / oidt.

        Callers operating on the instance pass ``self.dcati``; the
        nested-MC driver passes the augmented ``pd.concat([xi, z_s])``.
        Participants missing a timepoint for some items are retained for
        their complete items; only the incomplete (pid, item) responses
        are dropped (not the whole participant)."""
        dcati = df

        if interim_date is not None:
            dcati = dcati[dcati['submission_date'] <= interim_date]

        if dcati.empty:
            return dcati

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

    # ------------------------------------------------------------------
    # training-set summary W(z^(s))_t
    # ------------------------------------------------------------------

    def get_w(self, zi: pd.DataFrame,
              categorical_threshold: Optional[int] = None) -> pd.DataFrame:
        """Build the Strong-Oakley per-(item, draw) summary frame from
        the future-data block ``zi`` (must carry
        ``ypred_0 .. ypred_{S-1}`` columns).

        For each draw column:
        - ``out-of-7`` items -> mean of ypred_s (1..7),
        - ``categorical`` items -> proportion of ypred_s
          ``>= categorical_threshold``.

        Pivots to ``w_baseline`` / ``w_endline`` per (item, draw) and
        derives direction-aware ``w_diff`` / ``w_ratio`` consistent with
        the endpoint-ratio convention."""
        dit = self.dit
        ypred_cols = [c for c in zi.columns if c.startswith('ypred_')]

        cat_thresh = (categorical_threshold
                      if categorical_threshold is not None
                      else self.categorical_threshold)
        zi_o7 = zi[zi['item_type'] == 'out-of-7']
        zi_cat = zi[zi['item_type'] == 'categorical']
        w_o7 = (
            zi_o7.groupby(['item_label', 'item_type', 'time'])[ypred_cols]
            .mean().reset_index()
        )
        zi_cat_bin = zi_cat[ypred_cols].ge(cat_thresh).astype(float)
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
        return wa

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

    # ------------------------------------------------------------------
    # Public endpoint methods
    # ------------------------------------------------------------------

    def get_endpoints_per_draw_from_theta_batch(
        self, theta_batch, x_stan, endpoint_type: str = 'items',
    ) -> pd.DataFrame:
        """vmap ``eval_outcome_for_endpoint`` over the parameter batch on
        the x-cohort design, wrap into an arviz idata under
        ``ordered_prob_by_cat_qu_pr`` and score via
        :meth:`get_endpoints_per_draw`."""
        ordprob = np.asarray(
            jax.vmap(lambda p: self.eval_outcome_for_endpoint(x_stan, p))(theta_batch)
        )
        idata = az.from_dict(
            posterior={'ordered_prob_by_cat_qu_pr': ordprob[None, ...]},
        )
        return self.get_endpoints_per_draw(draws=idata, endpoint_type=endpoint_type)


    def get_endpoints_per_draw(
        self,
        draws=None,
        draws_file: Optional[str] = None,
        categorical_threshold: Optional[int] = None,
        endpoint_type: Literal["items", "item_groups"] = "items",
        param_name: str = "ordered_prob_by_cat_qu_pr",
        verbose: bool = False,
    ) -> pd.DataFrame:
        """Per-draw directional ``diff`` and ``ratio`` per item or
        item-group, computed from ``self.dcati`` + ``self.dit`` and the
        supplied posterior. ``categorical_threshold`` defaults to
        ``self.categorical_threshold`` (set in the constructor)."""
        if endpoint_type not in ("items", "item_groups"):
            raise ValueError("endpoint_type must be either 'items' or 'item_groups'")
        if param_name not in ("ordered_prob_by_cat_qu_pr", "ordered_prob_by_cat_qu_fit"):
            raise ValueError(
                "param_name must be either 'ordered_prob_by_cat_qu_pr' or 'ordered_prob_by_cat_qu_fit'"
            )
        if categorical_threshold is None:
            categorical_threshold = self.categorical_threshold
        po_arr = self._resolve_draws(draws, draws_file, param_name, verbose=verbose)
        po = self._make_po(po_arr)

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

        po['diff'] = np.nan
        po['ratio'] = np.nan
        tmp = po['item_high_label'] == 'lower_is_better'
        po.loc[tmp, 'diff'] = po.loc[tmp, 'Baseline'] - po.loc[tmp, 'Endline']
        po.loc[tmp, 'ratio'] = 1 - po.loc[tmp, 'Endline'] / po.loc[tmp, 'Baseline']
        tmp = po['item_high_label'] == 'higher_is_better'
        po.loc[tmp, 'diff'] = po.loc[tmp, 'Endline'] - po.loc[tmp, 'Baseline']
        po.loc[tmp, 'ratio'] = po.loc[tmp, 'Endline'] / po.loc[tmp, 'Baseline'] - 1

        po.rename(columns={'.draw': 'draw'}, inplace=True)
        return po

    def get_endpoints(
        self,
        draws=None,
        draws_file: Optional[str] = None,
        categorical_threshold: Optional[int] = None,
        endpoint_type: Literal["items", "item_groups"] = "items",
        param_name: str = "ordered_prob_by_cat_qu_pr",
        verbose: bool = True,
    ) -> pd.DataFrame:
        """Endpoint summaries with quantile aggregation across draws
        (2.5%, 25%, 50%, 75%, 97.5%). Reads ``self.dcati`` + ``self.dit``
        and the supplied posterior."""
        vprint = print if verbose else (lambda *args, **kwargs: None)
        dit = self.dit

        vprint("Computing per-draw endpoints...")
        po = self.get_endpoints_per_draw(
            draws=draws, draws_file=draws_file,
            categorical_threshold=categorical_threshold,
            endpoint_type=endpoint_type,
            param_name=param_name, verbose=verbose,
        )

        if endpoint_type == "item_groups":
            id_vars = ['item_type', 'group_label']
        else:
            id_vars = ['item_type_id', 'item_type', 'item_label', 'group_label']

        po = po.melt(
            id_vars=['draw'] + id_vars,
            value_vars=['diff', 'ratio', 'Baseline', 'Endline'],
            var_name='variable',
            value_name='value',
        )

        quantiles = [0.025, 0.25, 0.5, 0.75, 0.975]
        quantile_names = ['q_lower', 'iqr_lower', 'median', 'iqr_upper', 'q_upper']
        pos = po.groupby(id_vars + ['variable'])['value'].quantile(quantiles).unstack()
        pos.columns = quantile_names
        pos = pos.reset_index()

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
