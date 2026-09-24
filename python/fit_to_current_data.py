"""Generic interim-SVI producer, shared by every application's *_startme.py.

`fit_interim_grid` runs one partial-credit SVI fit per interim over an accrual order of participants,
saving the standard per-interim files (dp1 / draws.zarr / prob_by_question_fit / the per-draw endpoint
pkl) into one output dir. The app-specific parts — which loader, how participants accrue, the rho
declarations + H1 rule — are passed in (built inline in the startme), so a new application needs no
new script, just a config + a few lines in its startme.

Helpers: assign_item_group_id (the common item x group difficulty index), the accrual orders
(shuffled / ordered / interleaved), the interim grids (linspace / step), and endpoint_frame (the
common one-call rho extraction -> the pkl frame). SDY269's one-by-one cross-arm extraction is done
inline in its startme via a custom `on_fit`."""
import os
import time
import itertools
import numpy as np
import pandas as pd

FILE_PREFIX = "pcm_1_interim"
_PKL_COLS = ['draw', 'item_label', 'item_type', 'item_high_label', 'rho_id', 'rho_label',
             'rho_label_long', 'reduction', 'compare', 'group1', 'group2', 'pps_rho_x',
             'pps_ratio_x', 'pps_H1_x']


# ---- structure + accrual + grid helpers ---------------------------------------------------
def assign_item_group_id(dp1):
    """Add item_group_id = the (item_label x group) difficulty index (sorted group, then item)."""
    dp1 = dp1.copy()
    it = (dp1[['item_label', 'group']].drop_duplicates()
          .sort_values(['group', 'item_label']).reset_index(drop=True))
    it['item_group_id'] = np.arange(1, len(it) + 1)
    return dp1.merge(it, on=['item_label', 'group'], how='left')


def shuffled_accrual(dp1, seed=123):
    """Random accrual order over participants (both arms present at every interim)."""
    return np.random.default_rng(seed).permutation(np.sort(dp1.pid.unique()))


def ordered_accrual(dp1):
    """Accrual pseudo-ordered by the participant-id string (deterministic, no arm balancing)."""
    return dp1[['pid', 'pid_label']].drop_duplicates().sort_values('pid_label')['pid'].to_numpy()


def interleaved_accrual(dp1, arms):
    """Pooled accrual that interleaves the arms so both accrue together (head-to-head)."""
    per = [dp1[dp1.arm == a][['pid', 'pid_label']].drop_duplicates().sort_values('pid_label')['pid'].tolist()
           for a in arms]
    return np.array([p for tup in itertools.zip_longest(*per) for p in tup if p is not None])


def linspace_grid(n_full, nint=10, floor=40):
    return np.unique(np.round(np.linspace(max(floor, n_full // nint), n_full, nint)).astype(int)).tolist()


def step_grid(n_full, step=10):
    return sorted(set(list(range(step, n_full, step)) + [n_full]))


# ---- the common one-call endpoint extraction --------------------------------------------
def endpoint_frame(model, fit, rho_specs, h1=0.5, contrast_col=None):
    """One get_endpoints_per_draw call -> the per-draw endpoint pkl frame. `h1` is a scalar or a
    {rho_id: threshold} map for pps_H1_x. Adds pps_ratio_x (= pps_rho_x, legacy deploy alias)."""
    kw = dict(draws=fit['draws'], rho_specs=rho_specs, endpoint_type='items')
    if contrast_col:
        kw['contrast_col'] = contrast_col
    xr = model.get_endpoints_per_draw(**kw).rename(columns={'rho': 'pps_rho_x'})
    rll = {s['rho_id']: s['rho_label_long'] for s in rho_specs}
    xr['rho_label_long'] = xr['rho_id'].map(rll)
    xr['pps_ratio_x'] = xr['pps_rho_x']
    thr = xr['rho_id'].map(h1) if isinstance(h1, dict) else h1
    xr['pps_H1_x'] = (xr['pps_rho_x'] > thr).astype(int)
    return xr[[c for c in _PKL_COLS if c in xr.columns]]


# ---- the generic interim loop -----------------------------------------------------------
def fit_interim_grid(dir_out, dp1, dit, pids, grid, x_formula, on_fit, *, seed=123,
                     nsteps=4000, output_samples=2000, prob_width=None, verbose=False, label=''):
    """SVI at each interim over `grid` (accrue `pids[:n]`). `on_fit(model, fit, k, n, xi) -> DataFrame`
    returns the per-draw endpoint frame to pickle (built with endpoint_frame or a custom extraction)."""
    from model_pcm import PartialCreditModel
    os.makedirs(dir_out, exist_ok=True)
    if prob_width is not None:
        os.environ.setdefault('PROB_FIT_WIDTH_MULT', str(prob_width))
    dit.to_csv(f"{dir_out}/{FILE_PREFIX}_1_data_dit.csv", index=False)
    for k, n in enumerate(grid, 1):
        obs = set(pids[:n]); xi = dp1[dp1.pid.isin(obs)].copy()
        xi = xi.sort_values(['item_type_id', 'pid', 'group', 'item_label']).reset_index(drop=True)
        xi['oid'] = range(1, len(xi) + 1); xi['oidt'] = xi.groupby('item_type').cumcount() + 1
        xi.to_csv(f"{dir_out}/{FILE_PREFIX}_{k}_data_dp1.csv", index=False)
        pre = f"{dir_out}/{FILE_PREFIX}_{k}"
        print(f"\n=== {label} interim {k}: n={n} ===")
        t0 = time.time()
        model = PartialCreditModel(dit=dit, dcati=xi, x_formula=x_formula, seed=seed)
        fit = model.fit_pyro_svi(output_file_prefix=pre, algorithm='AutoDiagonalNormal', lr=0.01,
                                 num_steps=nsteps, output_samples=output_samples, resume=True,
                                 with_core_analyses=True, with_additional_analyses=False, verbose=verbose)
        on_fit(model, fit, k, n, xi).to_pickle(f"{dir_out}/{FILE_PREFIX}_i{k}_regression_training.pkl")
        print(f"  done ({(time.time() - t0) / 60:.1f} min)")
    print(f"{label} SVI grid complete -> {dir_out}")
