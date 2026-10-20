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
from utils import _map_cq_id_to_item_structure, define_effect_sizes


class IRTModel(Model):
    """Marker mixin for the IRT family. Provides per-draw + summary
    endpoint extraction from posterior ordered-prob arrays."""

    def __init__(self, dit: pd.DataFrame, dcati: pd.DataFrame,
                 x_formula: str = "~ group - 1", *,
                 seed: int = 123, categorical_threshold: int = 3,
                 drho: Optional[pd.DataFrame] = None):
        self.categorical_threshold = int(categorical_threshold)
        # Mirror the R side: the reduced-effect-size manifest `drho` is defined ONCE, up front,
        # and then consumed by both make_stan_data (-> the Stan/pyro reduction data block) and the
        # get_endpoints rho reader -- never recomputed inside a fit. Built from (dcati, dit) here so
        # it exists before super().__init__ triggers make_stan_data; callers may pass their own.
        if drho is None:
            drho = define_effect_sizes(
                dcati, dit, categorical_threshold=self.categorical_threshold)
        self.drho = drho
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

    def _rho_by_cat_qu_endpoints_per_draw(self, draws, draws_file, param_name: str,
                                          endpoint_type: str, verbose: bool = False,
                                          across=None) -> pd.DataFrame:
        """Per-draw endpoints read from the in-model reductions ``rho_by_cat_qu_*`` (R x Qcmp) -- the
        Python twin of R ``get_endpoints``. ``self.drho`` is the single source of
        truth: each (reduction r, comparison jc) carries its label ``rho`` in {Baseline, Endline,
        ratio, diff}, its ``scope`` (item / item_group) and item metadata. We read ``rho[r, jc]``,
        join the labels, pivot the four per (draw, jc), then form the directional endpoints exactly
        as R does: ``ratio`` is emitted by the model as a PERCENT (converted to a fraction here);
        ``diff``/``ratio`` keep the Baseline-Endline orientation for categorical OR lower-is-better
        items and flip (Endline-Baseline) for out-of-7 higher-is-better; ``diff2``/``ratio2`` are the
        unconditional negatives. group1 = Baseline level, group2 = Endline level (same schema as the
        ordered-prob path, so :meth:`get_endpoints` consumes it unchanged)."""
        drho = self.drho
        labels = set(drho['rho'].unique())

        arr = np.asarray(self._resolve_draws(draws, draws_file, param_name, verbose=verbose))
        arr = arr.reshape(-1, arr.shape[-2], arr.shape[-1])          # (chain, draw, R, Qcmp) -> (ndraws, R, Qcmp)
        ndraws, R, Qcmp = arr.shape

        rr, jj = np.meshgrid(np.arange(1, R + 1), np.arange(1, Qcmp + 1), indexing='ij')
        long = pd.DataFrame({
            'draw': np.repeat(np.arange(ndraws), R * Qcmp),
            'r': np.tile(rr.ravel(), ndraws),
            'jc': np.tile(jj.ravel(), ndraws),
            'val': arr.ravel(),
        })
        long = long.merge(drho[['jc', 'r', 'rho']].drop_duplicates(), on=['jc', 'r'])
        # one column per reduction label present in the manifest (Baseline/Endline/ratio/diff and/or
        # app-specific SPR/GMFR/...); arbitrary reductions are carried through, not just the fixed four.
        w = long.pivot_table(index=['draw', 'jc'], columns='rho', values='val').reset_index()
        w.columns.name = None

        cmpmeta = drho[['jc', 'scope', 'item_type_id', 'item_type', 'item_label',
                        'construct', 'item_high_label', 'stratum']].drop_duplicates('jc')
        w = w.merge(cmpmeta, on='jc')
        scope = 'item' if endpoint_type == 'items' else 'item_group'
        w = w[w['scope'] == scope].copy()

        rho_cols = sorted(labels)
        label2redfn = drho.drop_duplicates('rho').set_index('rho')['red_fn'].astype(int).to_dict()
        # Direction rule keyed on the REDUCTION (red_fn), so any app label is oriented correctly:
        #  - rel% (3 per-participant, 7 group): percent -> fraction, then direction;
        #  - diff  (4 per-participant, 8 group): direction;
        #  - direction = keep Baseline-Endline orientation for categorical OR lower-is-better, flip for
        #    out-of-7 higher-is-better. Each oriented label gets a `<label>2` = unconditional negative.
        #  - levels (5/6), GMFR (1/9), SCR (2), SPR (= endline level 6 @ thr): passed through raw.
        dir_pos = ((w['item_type'] == 'categorical')
                   | (w['item_high_label'] == 'lower_is_better')).to_numpy()
        extra2 = []
        for lab in rho_cols:
            rf = label2redfn.get(lab)
            col = w[lab].to_numpy(float)
            if rf in (3, 7):
                raw = col / 100.0
                w[lab] = np.where(dir_pos, raw, -raw)
                w[lab + '2'] = -raw
                extra2.append(lab + '2')
            elif rf in (4, 8):
                w[lab] = np.where(dir_pos, col, -col)
                w[lab + '2'] = -col
                extra2.append(lab + '2')
        # group1/group2 = raw Baseline/Endline levels (legacy alias), when the manifest carries them.
        if {'Baseline', 'Endline'}.issubset(labels):
            w['group1'] = w['Baseline']
            w['group2'] = w['Endline']

        # Stratum (e.g. 'arm'): exposed under its own name when the manifest is stratified, dropped
        # otherwise (all-NA) so the un-stratified output is byte-identical to the pre-stratum reader.
        strata_col = drho.attrs.get('strata_col', None)
        if strata_col:
            w = w.rename(columns={'stratum': strata_col})
        else:
            w = w.drop(columns=['stratum'])

        # Cross-stratum "across" endpoints (e.g. SDY269 TIV-vs-LAIV diffs): EXACT per-draw contrast of the
        # same reduction between two strata, matched on `within` (shared item_label). Each spec adds one new
        # rho column `out_label` carried on rows tagged stratum='foc-ref'; the per-stratum rows keep NaN there.
        if across:
            new_rows = []
            idx = ['draw', 'item_label', 'item_type_id', 'item_type', 'construct', 'item_high_label']
            for spec in across:
                rl, out = spec['rho_label'], spec['out_label']
                ref, foc, comp = spec['ref'], spec['foc'], spec.get('compare', 'diff')
                within = spec.get('within', 'item_label')
                sub = w[w[strata_col].isin([ref, foc])]
                piv = sub.pivot_table(index=[c for c in idx if c in sub.columns], columns=strata_col, values=rl)
                piv = piv.dropna(subset=[ref, foc])
                if comp == 'diff':
                    dv = piv[foc] - piv[ref]
                elif comp in ('fold', 'ratio'):
                    dv = piv[foc] / piv[ref]
                else:
                    raise ValueError(f"unknown across compare {comp!r}")
                dd = dv.reset_index(name=out)
                dd[strata_col] = f'{foc}-{ref}'
                new_rows.append(dd)
            if new_rows:
                from functools import reduce as _reduce
                mkeys = [c for c in idx if c in new_rows[0].columns] + [strata_col]   # shared id cols + stratum tag
                diff_df = _reduce(lambda a, b: a.merge(b, on=mkeys, how='outer'), new_rows)
                w = pd.concat([w, diff_df], ignore_index=True)
                rho_cols = rho_cols + [s['out_label'] for s in across if s['out_label'] not in rho_cols]

        if endpoint_type == 'items':
            keep = ['draw', 'item_type_id', 'item_type', 'item_label', 'construct', 'item_high_label']
        else:
            keep = ['draw', 'item_type', 'construct', 'item_high_label']
        if strata_col:
            keep.append(strata_col)
        extra = extra2 + [c for c in ('group1', 'group2') if c in w.columns]
        keep += rho_cols + [c for c in extra if c not in rho_cols]
        return w[[c for c in keep if c in w.columns]].reset_index(drop=True)

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
        category_thresholds: Optional[dict] = None,
        endpoint_type: Literal["items", "item_groups"] = "items",
        param_name: str = "ordered_prob_by_cat_qu_pr",
        across: Optional[list] = None,
        verbose: bool = False,
    ) -> pd.DataFrame:
        """Per-draw directional ``diff`` and ``ratio`` per item or
        item-group, computed from ``self.dcati`` + ``self.dit`` and the
        supplied posterior. ``categorical_threshold`` defaults to
        ``self.categorical_threshold`` (set in the constructor)."""
        if endpoint_type not in ("items", "item_groups"):
            raise ValueError("endpoint_type must be either 'items' or 'item_groups'")
        if param_name not in ("ordered_prob_by_cat_qu_pr", "ordered_prob_by_cat_qu_fit",
                              "rho_by_cat_qu_pr", "rho_by_cat_qu_fit"):
            raise ValueError(
                "param_name must be one of 'ordered_prob_by_cat_qu_pr', "
                "'ordered_prob_by_cat_qu_fit', 'rho_by_cat_qu_pr', 'rho_by_cat_qu_fit'"
            )
        # rho_by_cat_qu_* are the in-model reductions (R x Qcmp): read + label-decode them directly
        # from self.drho (the R get_endpoints path), bypassing the ordered-prob reduction.
        if param_name in ("rho_by_cat_qu_pr", "rho_by_cat_qu_fit"):
            return self._rho_by_cat_qu_endpoints_per_draw(
                draws, draws_file, param_name, endpoint_type, verbose=verbose, across=across)
        if categorical_threshold is None:
            categorical_threshold = self.categorical_threshold
        po_arr = self._resolve_draws(draws, draws_file, param_name, verbose=verbose)
        po = self._make_po(po_arr)

        dp1 = self.dcati
        dit = self.dit
        parts = []
        # The two conditions being contrasted are keyed on the NEUTRAL index ``time`` (0/1),
        # not on the display string ``group_label``: time 0 -> ``group1``, time 1 -> ``group2``.
        # ``group_label`` is a free display label (Baseline/Endline for a paired design, or
        # e.g. adult/pediatric, mild/severe, vaccine/placebo for a between-group design) used
        # only by the plots. The endpoint (diff/ratio) is group2-vs-group1 either way.
        # EVERY category type carries a threshold (default 0 => expected score E[y]; >=1 => caseness
        # P(y>=thr)), applied UNIFORMLY -- not gated on item_type=='categorical'. `category_thresholds`
        # {item_type: thr} overrides the per-type default (categorical -> categorical_threshold, any
        # other type -> 0). Also iterates whatever types are present (not a hardcoded pair), so no type
        # is silently dropped. (The in-model rho_by_cat_qu reduction already honours thr per type.)
        thr_map = dict(category_thresholds or {})
        for item_type in dp1['item_type'].drop_duplicates():
            thr_t = int(thr_map.get(item_type,
                                    categorical_threshold if item_type == 'categorical' else 0))
            tmp = dp1[dp1['item_type'] == item_type][
                ['item_type_id', 'item_label', 'item_group_id', 'item_type', 'group']
            ].drop_duplicates()
            sub = po.merge(tmp, on=['item_type_id', 'item_group_id'])
            if sub.empty:
                continue
            if thr_t >= 1:                                   # caseness P(y >= thr)
                sub = sub[sub['y'] >= thr_t].assign(_w=lambda d: d['prob'])
            else:                                            # thr == 0 -> expected score E[y]
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
    def joint_all_met(endpoints_long: pd.DataFrame, group_keys, met_col: str = 'pps_H1_x',
                      draw_col: str = 'draw'):
        """Composite decision requiring **ALL** endpoints met simultaneously (the conjunction), the
        strict counterpart of :meth:`joint_any_met`. Same shared-draw frame; evaluated WITHIN a draw
        as the AND of the per-rho indicators (product == min over 0/1), then the PPS is the fraction
        of draws meeting all of them. A much stricter bar than "at least one met" — it is the correct
        joint posterior probability that *every* criterion (e.g. SCR *and* SPR *and* GMFR) holds on
        the same draw. Returns (per_draw, pps): `per_draw` has one `all_met` per (draw, *group_keys*);
        `pps` has `pps_all_met` per *group_keys*."""
        gk = list(group_keys)
        per_draw = (endpoints_long.groupby([draw_col] + gk)[met_col].min()
                    .reset_index().rename(columns={met_col: 'all_met'}))
        pps = (per_draw.groupby(gk)['all_met'].mean()
               .reset_index().rename(columns={'all_met': 'pps_all_met'}))
        return per_draw, pps


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
