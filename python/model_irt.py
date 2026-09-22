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
                 x_formula: str = "~ group - 1", *,
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
            ['item_type_id', 'pid', 'group', 'item_label']
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
            ['item_type_id', 'pid', 'group', 'item_label']
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
            dcati.groupby(['pid', 'item_label'])['group']
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

        dcati = dcati.sort_values(['item_type_id', 'pid', 'group', 'item_label']).reset_index(drop=True)
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
            zi_o7.groupby(['item_label', 'item_type', 'group'])[ypred_cols]
            .mean().reset_index()
        )
        zi_cat_bin = zi_cat[ypred_cols].ge(cat_thresh).astype(float)
        zi_cat_bin = pd.concat(
            [zi_cat[['item_label', 'item_type', 'group']].reset_index(drop=True),
             zi_cat_bin.reset_index(drop=True)], axis=1,
        )
        w_cat = (
            zi_cat_bin.groupby(['item_label', 'item_type', 'group'])[ypred_cols]
            .mean().reset_index()
        )
        wa = pd.concat([w_o7, w_cat], ignore_index=True).melt(
            id_vars=['item_label', 'item_type', 'group'],
            value_vars=ypred_cols, var_name='s_col', value_name='w',
        )
        wa['draw'] = wa['s_col'].str.replace('ypred_', '').astype(int)
        wa = wa.drop(columns='s_col')

        t_min, t_max = wa['group'].min(), wa['group'].max()
        wa = (
            wa.pivot_table(index=['item_label', 'item_type', 'draw'],
                           columns='group', values='w')
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
        ``item_group_id``, ``y``) recovered from ``cq_id`` via
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
        rho_specs: Optional[list] = None,
        contrast_col: str = "group",
        verbose: bool = False,
    ) -> pd.DataFrame:
        """Per-draw directional ``diff`` and ``ratio`` per item or
        item-group, computed from ``self.dcati`` + ``self.dit`` and the
        supplied posterior. ``categorical_threshold`` defaults to
        ``self.categorical_threshold`` (set in the constructor).

        If ``rho_specs`` is given, returns instead a long, MULTI-endpoint frame
        (one row per draw x item x rho) tagged by ``rho_id`` / ``rho_label`` — see
        :meth:`_rho_endpoints_per_draw`. Each spec picks a per-time-point reduction
        of the category probabilities (``reduction``) and a way to combine the
        Baseline / Endline values into the endpoint (``compare``). This is the
        general path (e.g. influenza HAI needs both a seroprotection-rate endpoint
        and a GMT-fold-rise endpoint); ``rho_specs=None`` keeps the legacy single
        ``diff``/``ratio`` behaviour used by every existing caller."""
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

        if rho_specs is not None:
            return self._rho_endpoints_per_draw(po, rho_specs, endpoint_type=endpoint_type,
                                                contrast_col=contrast_col)

        dp1 = self.dcati
        dit = self.dit
        parts = []
        # The two conditions being contrasted are keyed on the NEUTRAL index ``time`` (0/1),
        # not on the display string ``group_label``: time 0 -> ``group1``, time 1 -> ``group2``.
        # ``group_label`` is a free display label (Baseline/Endline for a paired design, or
        # e.g. adult/pediatric, mild/severe, vaccine/placebo for a between-group design) used
        # only by the plots. The endpoint (diff/ratio) is group2-vs-group1 either way.
        for item_type in ('categorical', 'out-of-7'):
            tmp = dp1[dp1['item_type'] == item_type][
                ['item_type_id', 'item_label', 'item_group_id', 'item_type', 'group']
            ].drop_duplicates()
            sub = po.merge(tmp, on=['item_type_id', 'item_group_id'])
            if sub.empty:
                continue
            if item_type == 'categorical':
                sub = sub[sub['y'] >= categorical_threshold]
                sub = sub.assign(_w=sub['prob'])
            else:
                sub = sub.assign(_w=sub['y'] * sub['prob'])
            sub = sub.groupby(
                ['.draw', 'item_type_id', 'item_label', 'item_group_id', 'item_type', 'group']
            ).agg(value=('_w', 'sum')).reset_index()
            parts.append(sub)
        po = pd.concat(parts, ignore_index=True)

        po = po.merge(
            dit[['item_type', 'item_label', 'construct']],
            on=['item_type', 'item_label'],
        )

        if endpoint_type == 'item_groups':
            id_vars = ['item_type', 'construct']
            po = po.groupby(['.draw', 'item_type', 'group', 'construct']).agg(
                value=('value', 'mean')
            ).reset_index()
        else:
            id_vars = ['item_type_id', 'item_type', 'item_label', 'construct']

        po = po.pivot_table(
            index=['.draw'] + id_vars,
            columns='group',
            values='value',
        ).rename(columns={0: 'group1', 1: 'group2'}).reset_index()
        po = po.dropna(subset=['group1', 'group2'])
        po = po.merge(
            dit[['item_type', 'construct', 'item_high_label']].drop_duplicates(),
            on=['item_type', 'construct'],
        )

        po['diff'] = np.nan
        po['ratio'] = np.nan
        tmp = po['item_high_label'] == 'lower_is_better'
        po.loc[tmp, 'diff'] = po.loc[tmp, 'group1'] - po.loc[tmp, 'group2']
        po.loc[tmp, 'ratio'] = 1 - po.loc[tmp, 'group2'] / po.loc[tmp, 'group1']
        tmp = po['item_high_label'] == 'higher_is_better'
        po.loc[tmp, 'diff'] = po.loc[tmp, 'group2'] - po.loc[tmp, 'group1']
        po.loc[tmp, 'ratio'] = po.loc[tmp, 'group2'] / po.loc[tmp, 'group1'] - 1

        po.rename(columns={'.draw': 'draw'}, inplace=True)
        return po

    @staticmethod
    def joint_any_met(endpoints_long: pd.DataFrame, group_keys, met_col: str = 'pps_H1_x',
                      draw_col: str = 'draw'):
        """Joint (composite) decision over several endpoints that are indexed on the SAME Monte-Carlo
        sample. `endpoints_long` is the multi-endpoint per-draw frame from
        :meth:`get_endpoints_per_draw` (one row per draw x item x rho, each carrying a per-rho
        criterion indicator `met_col` = 1{rho > its threshold}). Because the rhos share the draw
        index, the composite "**at least one criterion met**" (e.g. CHMP influenza: SCR *or* SPR
        *or* GMFR met) is evaluated WITHIN a draw as the OR of the per-rho indicators, then the PPS
        is the fraction of draws meeting it — the correct JOINT posterior probability, not a naive
        combination of the marginal PPS. Returns (per_draw, pps): `per_draw` has one
        `any_met` per (draw, *group_keys*); `pps` has `pps_any_met` per *group_keys*."""
        gk = list(group_keys)
        per_draw = (endpoints_long.groupby([draw_col] + gk)[met_col].max()
                    .reset_index().rename(columns={met_col: 'any_met'}))
        pps = (per_draw.groupby(gk)['any_met'].mean()
               .reset_index().rename(columns={'any_met': 'pps_any_met'}))
        return per_draw, pps

    @staticmethod
    def _apply_rho_compare(piv: pd.DataFrame, compare: str) -> "pd.Series":
        """Combine the two per-condition summaries (``group1`` = time 0, ``group2`` = time 1)
        into a scalar endpoint, respecting ``item_high_label`` direction. ``compare``:
          'endline'/'group2' -> the group2 level only (group1 computed but unused, e.g. a
                                seroprotection RATE at endline);
          'baseline'/'group1' -> the group1 level only;
          'fold'      -> geometric/level fold group2/group1 (good direction);
          'fold_log2' -> 2**(group2-group1) when the summary is on a log2 scale
                         (e.g. mean log2-titre -> GMT fold-rise);
          'ratio'     -> direction-aware relative change vs the group1 baseline (legacy
                         `get_endpoints`): higher_is_better group2/group1-1, lower_is_better
                         1-group2/group1;
          'diff'      -> group2 - group1.
        For 'lower_is_better' items the group1/group2 roles are swapped for the change-type
        compares (fold/fold_log2/ratio/diff). ('endline'/'baseline' kept as aliases so a
        paired design reads naturally; a between-group design can use 'group2'/'group1'.)"""
        hi = piv['item_high_label'].to_numpy() if 'item_high_label' in piv else np.array(['higher_is_better'] * len(piv))
        lower = hi == 'lower_is_better'
        B, E = piv['group1'].to_numpy(float), piv['group2'].to_numpy(float)
        if compare in ('endline', 'group2'):
            return pd.Series(E, index=piv.index)
        if compare in ('baseline', 'group1'):
            return pd.Series(B, index=piv.index)
        # change-type: orient so "good" is the numerator/positive
        num = np.where(lower, B, E)
        den = np.where(lower, E, B)
        if compare == 'fold':
            with np.errstate(divide='ignore', invalid='ignore'):
                return pd.Series(num / den, index=piv.index)
        if compare == 'fold_log2':
            return pd.Series(2.0 ** (num - den), index=piv.index)
        if compare == 'ratio':                                  # improvement relative to group1 (legacy)
            with np.errstate(divide='ignore', invalid='ignore'):
                good = np.where(lower, B - E, E - B)
                return pd.Series(good / B, index=piv.index)
        if compare == 'diff':
            return pd.Series(num - den, index=piv.index)
        raise ValueError(f"unknown compare {compare!r}")

    def _rho_endpoints_per_draw(self, po: pd.DataFrame, rho_specs: list,
                                endpoint_type: str = "items",
                                contrast_col: str = "group") -> pd.DataFrame:
        """Multi-endpoint per-draw frame. ``po`` is ``_make_po`` output (``.draw``,
        ``cq_id``, ``prob``, ``item_type_id``, ``item_group_id``, ``y``). Each spec in
        ``rho_specs`` is a dict with:
          rho_id, rho_label       identifiers carried into the output;
          reduction               'threshold' (P(y >= threshold), e.g. seroprotection)
                                  or 'mean' (E[y]; on the log2-dilution scale this is
                                  mean log2 titre -> pair with compare='fold_log2' for GMT);
          threshold               the y cut for reduction='threshold' (y is 0-indexed:
                                  y>=3 == titre>=1:40 on the HAI ladder start=5);
          compare                 how group1/group2 combine (see _apply_rho_compare);
          across/across_*         optional cross-stratum difference (see :meth:`_rho_across`).
        ``contrast_col`` is the 0/1 axis the endpoint contrasts (0->group1, 1->group2); it is
        ``'group'`` for the usual paired/between design, but a head-to-head fit whose flexible
        ``group`` carries arm x time passes ``contrast_col='phase'`` (the within-arm baseline/
        endline axis) and keeps ``arm`` as a stratum, so SPR/GMFR come out per (item, arm) and the
        cross-arm difference is a further ``across='arm'`` step. Any per-item_group_id descriptor
        columns present in ``dcati`` beyond the item (e.g. ``arm``) are carried as strata.
        Returns columns: draw, item_type_id, item_type, item_label, construct, item_high_label,
        [strata], rho_id, rho_label, reduction, compare, group1, group2, rho."""
        dp1, dit = self.dcati, self.dit
        # strata = per-item_group_id descriptors that are neither the item nor the contrast axis
        strata = [c for c in ('arm',) if c in dp1.columns and c != contrast_col]
        struct = dp1[['item_type_id', 'item_group_id', 'item_label', 'item_type',
                      contrast_col] + strata].drop_duplicates()
        base = po.merge(struct, on=['item_type_id', 'item_group_id'])
        base = base.merge(
            dit[['item_type', 'item_label', 'construct', 'item_high_label']].drop_duplicates(),
            on=['item_type', 'item_label'])
        if endpoint_type == 'item_groups':
            idx = ['.draw', 'item_type', 'construct', 'item_high_label']
        else:
            idx = ['.draw', 'item_type_id', 'item_type', 'item_label',
                   'construct', 'item_high_label'] + strata
        gkeys = idx + [contrast_col]
        out = []
        for spec in rho_specs:
            red = spec.get('reduction', 'mean')
            if red == 'threshold':
                k = int(spec['threshold'])
                b = base.assign(_w=np.where(base['y'] >= k, base['prob'], 0.0))
            elif red == 'mean':
                b = base.assign(_w=base['y'] * base['prob'])
            else:
                raise ValueError(f"unknown reduction {red!r}")
            val = b.groupby(gkeys, as_index=False).agg(value=('_w', 'sum'))
            if endpoint_type == 'item_groups':     # average the per-item summaries within a group
                val = val.groupby(idx + [contrast_col], as_index=False).agg(value=('value', 'mean'))
            piv = (val.pivot_table(index=idx, columns=contrast_col, values='value')
                   .rename(columns={0: 'group1', 1: 'group2'}).reset_index())
            for tcol in ('group1', 'group2'):
                if tcol not in piv.columns:
                    piv[tcol] = np.nan
            piv['rho'] = self._apply_rho_compare(piv, spec.get('compare', 'ratio'))
            if spec.get('across'):
                piv = self._rho_across(piv, spec)      # cross-arm (or other) difference endpoint
            piv['rho_id'] = spec.get('rho_id')
            piv['rho_label'] = spec.get('rho_label', str(spec.get('rho_id')))
            piv['reduction'] = red
            piv['compare'] = spec.get('compare', 'ratio')
            out.append(piv)
        res = pd.concat(out, ignore_index=True)
        res.rename(columns={'.draw': 'draw'}, inplace=True)
        return res

    @staticmethod
    def _rho_across(piv: pd.DataFrame, spec: dict) -> pd.DataFrame:
        """Contrast the per-item base endpoint `rho` ACROSS a descriptor column (e.g. `arm`),
        WITHIN a matching unit (e.g. `strain`), per draw — for a head-to-head difference such as
        GMFR_diff / SPR_diff on the shared strain. Because both arms come from the SAME joint fit,
        this is an exact per-draw contrast (no random draw-pairing). Spec keys: `across` (the
        column, e.g. 'arm'), `across_values` = [ref, foc] (endpoint oriented as foc-vs-ref),
        `across_within` (matching unit, default 'strain'), `across_compare` ('diff'|'fold'|'ratio',
        default 'diff'). Output keeps the standard schema with group1 = ref value, group2 = foc
        value, rho = their contrast; `item_label`/`strain` collapse to the matching unit and `arm`
        becomes '<foc>-<ref>'. Only units present for BOTH values survive."""
        ac = spec['across']; within = spec.get('across_within', 'item_label')
        ref, foc = spec['across_values']; comp = spec.get('across_compare', 'diff')
        w = piv.pivot_table(index=['.draw', within], columns=ac, values='rho')
        w = w.dropna(subset=[ref, foc]).reset_index()
        g1, g2 = w[ref].to_numpy(float), w[foc].to_numpy(float)
        if comp == 'diff':
            r = g2 - g1
        elif comp == 'fold':
            with np.errstate(divide='ignore', invalid='ignore'):
                r = g2 / g1
        elif comp == 'ratio':
            with np.errstate(divide='ignore', invalid='ignore'):
                r = g2 / g1 - 1.0
        else:
            raise ValueError(f"unknown across_compare {comp!r}")
        const = {c: piv[c].iloc[0] for c in
                 ('item_type_id', 'item_type', 'construct', 'item_high_label') if c in piv}
        res = pd.DataFrame({'.draw': w['.draw'], **{c: const[c] for c in const},
                            'item_label': w[within], within: w[within],
                            ac: f'{foc}-{ref}', 'group1': g1, 'group2': g2, 'rho': r})
        return res

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
            id_vars = ['item_type', 'construct']
        else:
            id_vars = ['item_type_id', 'item_type', 'item_label', 'construct']

        po = po.melt(
            id_vars=['draw'] + id_vars,
            value_vars=['diff', 'ratio', 'group1', 'group2'],
            var_name='variable',
            value_name='value',
        )

        quantiles = [0.025, 0.25, 0.5, 0.75, 0.975]
        quantile_names = ['q_lower', 'iqr_lower', 'median', 'iqr_upper', 'q_upper']
        pos = po.groupby(id_vars + ['variable'])['value'].quantile(quantiles).unstack()
        pos.columns = quantile_names
        pos = pos.reset_index()

        if endpoint_type == "item_groups":
            tmp = dit[['item_type', 'construct', 'construct_long', 'item_high_label']].drop_duplicates()
            pos = pos.merge(tmp, on=['item_type', 'construct'])
        else:
            tmp = dit[
                ['item_type', 'item_label', 'item_label_short', 'construct',
                 'construct_long', 'item_high_label']
            ].drop_duplicates()
            pos = pos.merge(tmp, on=['item_type', 'item_label', 'construct'])

        vprint(f"Computed endpoints for {len(pos)} item-variable combinations")
        return pos
