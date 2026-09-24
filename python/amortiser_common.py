"""
Shared amortiser utilities (§11 of ``dev/amortised_decision_making.md``).

Contract:

- All amortiser :class:`flax.linen.Module` classes accept a **batch dict**
  ``batch`` (e.g. ``{'features': ...}`` for fixed pooled features, or
  ``{'x': ..., 'z': ..., 'mask_x': ..., 'mask_z': ..., 'sizes': ...}``
  for the learnable-features / DeepSets variant) and return a
  ``(batch, K)`` array of predicted quantiles at
  ``num_quantiles`` levels.

- ``sample_fn(rng, S)`` returns ``(batch_dict, rho)`` where
  ``batch_dict`` matches the schema expected by the amortiser class and
  ``rho`` is a ``(S,)`` continuous target.

- The shared functions
  :func:`train`, :func:`predict_amortised_p_h1_for_one_xz`,
  :func:`save_trained_model`, :func:`load_fitted_model` are
  amortiser-class-agnostic: they use ``net_class`` at construction
  time and treat the batch as an opaque pytree at call time.
"""

import importlib
import pickle
import time

import jax
import jax.numpy as jnp
import numpy as np
import optax
from flax import linen as nn


# ---------------------------------------------------------------------------
# Loss
# ---------------------------------------------------------------------------


def _pinball_loss(preds, target, pps_ProbH1_lwr_quantiles_mesh_arr):
    """Multi-quantile pinball / check loss.

    ``preds``: ``(batch, K)`` predicted quantiles.
    ``target``: ``(batch,)`` continuous label.
    ``pps_ProbH1_lwr_quantiles_mesh_arr``: ``(K,)`` quantile levels.
    Loss is the mean over ``batch x K``.
    """
    diff = target[:, None] - preds
    return jnp.mean(jnp.maximum(pps_ProbH1_lwr_quantiles_mesh_arr * diff, (pps_ProbH1_lwr_quantiles_mesh_arr - 1.0) * diff))


# ---------------------------------------------------------------------------
# Training
# ---------------------------------------------------------------------------


def _dummy_batch_from(peek_batch):
    """Return a size-1 replica of ``peek_batch`` for network init."""
    def _slice(x):
        x = jnp.asarray(x)
        return x[:1]
    return jax.tree_util.tree_map(_slice, peek_batch)


def train(
    sample_fn,
    net_class,
    pps_ProbH1_lwr_quantiles_mesh,
    *,
    net_kwargs=None,
    num_steps: int = 4000,
    batch_size: int = 1024,
    lr: float = 1e-3,
    seed: int = 0,
    verbose: bool = True,
):
    """Train the amortiser ``net_class`` on batches from ``sample_fn``.

    ``sample_fn(rng, S)`` returns ``(batch, rho)`` where ``batch`` is a
    pytree whose leaves are ``(S, ...)`` arrays (batch dimension first)
    and ``rho`` is a ``(S,)`` continuous target.

    ``net_class`` is a Flax module class. ``net_kwargs`` (default empty
    dict) is forwarded to the constructor.

    Returns a ``fit`` dict with keys:
      ``params``, ``apply_fn``, ``pps_ProbH1_lwr_quantiles_mesh``, ``net_class_module``,
      ``net_class_name``, ``net_kwargs``, ``loss_history``.
    """
    net_kwargs = dict(net_kwargs or {})
    net_kwargs.setdefault('num_quantiles', len(pps_ProbH1_lwr_quantiles_mesh))
    pps_ProbH1_lwr_quantiles_mesh_arr = jnp.asarray(pps_ProbH1_lwr_quantiles_mesh, dtype=jnp.float32)
    rng = np.random.default_rng(seed)

    peek_batch, _ = sample_fn(rng, 4)
    peek_batch = _jax_pytree(peek_batch)

    net = net_class(**net_kwargs)
    key = jax.random.PRNGKey(seed)
    dummy = _dummy_batch_from(peek_batch)
    params = net.init(key, dummy)

    if num_steps < 100:
        optimizer = optax.adam(learning_rate=lr)
    else:
        warmup_steps = min(max(200, num_steps // 20), num_steps // 2)
        schedule = optax.warmup_cosine_decay_schedule(
            init_value=0.0,
            peak_value=lr,
            warmup_steps=warmup_steps,
            decay_steps=num_steps,
            end_value=lr * 0.05,
        )
        optimizer = optax.adam(learning_rate=schedule)
    opt_state = optimizer.init(params)

    @jax.jit
    def train_step(params, opt_state, batch, target):
        def loss_fn(p):
            preds = net.apply(p, batch)
            return _pinball_loss(preds, target, pps_ProbH1_lwr_quantiles_mesh_arr)
        loss, grads = jax.value_and_grad(loss_fn)(params)
        updates, new_opt_state = optimizer.update(grads, opt_state)
        return optax.apply_updates(params, updates), new_opt_state, loss

    loss_history = []
    log_every = max(1, num_steps // 20)
    t_train_start = time.time()
    for step in range(num_steps):
        batch, rho = sample_fn(rng, batch_size)
        batch = _jax_pytree(batch)
        r = jnp.asarray(rho, dtype=jnp.float32)
        params, opt_state, loss = train_step(params, opt_state, batch, r)
        loss_history.append(float(loss))
        if verbose and ((step + 1) % log_every == 0 or step == 0):
            print(f"  step {step+1}/{num_steps}  loss={float(loss):.5f}")
    training_mins = (time.time() - t_train_start) / 60.0

    return {
        'params': params,
        'apply_fn': net.apply,
        'pps_ProbH1_lwr_quantiles_mesh': np.asarray(pps_ProbH1_lwr_quantiles_mesh, dtype=np.float32),
        'net_class_module': net_class.__module__,
        'net_class_name': net_class.__name__,
        'net_kwargs': dict(net_kwargs),
        'loss_history': loss_history,
        'training_mins': float(training_mins),
    }


def _jax_pytree(batch):
    """Convert numpy-array leaves to jax arrays, leaving jax arrays
    unchanged. Dict-only in practice; kept generic via ``tree_map``."""
    return jax.tree_util.tree_map(
        lambda x: jnp.asarray(x, dtype=jnp.float32), batch,
    )


# ---------------------------------------------------------------------------
# Save / load
# ---------------------------------------------------------------------------


def save_trained_model(fit: dict, path: str) -> None:
    """Save ``fit`` (params + config) to disk. ``apply_fn`` is dropped
    (rebuilt from ``net_class_module`` / ``net_class_name`` in
    :func:`load_fitted_model`)."""
    payload = {
        'params': jax.tree_util.tree_map(np.asarray, fit['params']),
        'pps_ProbH1_lwr_quantiles_mesh': np.asarray(fit['pps_ProbH1_lwr_quantiles_mesh']),
        'net_class_module': str(fit['net_class_module']),
        'net_class_name':   str(fit['net_class_name']),
        'net_kwargs':       dict(fit['net_kwargs']),
        'loss_history':     list(fit.get('loss_history', [])),
        'training_mins':    float(fit.get('training_mins', float('nan'))),
    }
    with open(path, 'wb') as f:
        pickle.dump(payload, f)


def load_fitted_model(path: str) -> dict:
    with open(path, 'rb') as f:
        payload = pickle.load(f)
    module = importlib.import_module(payload['net_class_module'])
    net_class = getattr(module, payload['net_class_name'])
    net = net_class(**payload['net_kwargs'])
    return {
        'params': jax.tree_util.tree_map(jnp.asarray, payload['params']),
        'apply_fn': net.apply,
        'pps_ProbH1_lwr_quantiles_mesh': np.asarray(payload['pps_ProbH1_lwr_quantiles_mesh']),
        'net_class_module': payload['net_class_module'],
        'net_class_name':   payload['net_class_name'],
        'net_kwargs':       dict(payload['net_kwargs']),
        'loss_history':     list(payload.get('loss_history', [])),
        'training_mins':    float(payload.get('training_mins', float('nan'))),
    }


# ---------------------------------------------------------------------------
# Prediction
# ---------------------------------------------------------------------------


def predict_amortised_p_h1_for_one_xz(
    fit: dict, batch, pps_H1_min_effect_size_thresh: float,
):
    """Forward-pass ``batch`` (pytree) through ``fit['apply_fn']`` and
    return ``(p_h1_xz, q_hat, quantiles)``.

    ``batch`` may be a dict (schema depends on the amortiser class) or a
    single array (for backwards-compatible callers passing pooled
    features directly, which are wrapped as ``{'features': ...}``).
    """
    apply_fn = fit['apply_fn']
    params = fit['params']
    pps_ProbH1_lwr_quantiles_mesh = np.asarray(fit['pps_ProbH1_lwr_quantiles_mesh'])
    # Backwards-compat: raw ndarray -> {'features': ...} for the
    # features-fixed amortiser class.
    if not isinstance(batch, dict):
        batch = {'features': batch}
    batch_jax = _jax_pytree(batch)
    preds = np.asarray(apply_fn(params, batch_jax))
    preds = np.maximum.accumulate(preds, axis=1)
    F = np.array([
        float(np.interp(pps_H1_min_effect_size_thresh, row, pps_ProbH1_lwr_quantiles_mesh,
                        left=0.0, right=1.0))
        for row in preds
    ])
    p_h1_xz = 1.0 - F
    med_idx = int(np.argmin(np.abs(pps_ProbH1_lwr_quantiles_mesh - 0.5)))
    q_hat = preds[:, med_idx]
    return p_h1_xz, q_hat, preds


# =============================================================================
# FEDERATED amortiser DEPLOY (§17.5 / §3.18-3.19). One generic loop, driven by the per-app config in
# the per-project *_startme.py CFG: for each declared endpoint it builds the per-endpoint reference (a
# `build` strategy re-using the SVI fit, no re-fit), runs the ragged BvM driver, and finally chains
# the combined diagnostics. Replaces the per-app *_amortiser_deploy.py / *_federated_deploy.py scripts.
# =============================================================================
import os as _os
import re as _re
import glob as _glob
import shutil as _shutil
import subprocess as _subprocess
import sys as _sys
import pandas as pd

_REPO = _os.path.dirname(_os.path.dirname(__file__))
DEFAULT_DRIVER = "scripts-py/Ukraine_interim_analysis_amortise_endpt_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead_contraction_bvm.py"
DEPLOY_SCRATCH = "/private/tmp/claude-501/-Users-or105-git-bIRTistic/49372e22-d11b-4882-a442-d1c60bcbdfb0/scratchpad"
# registry net-instance key -> (trained net dir under SB, widetok, case_c). warp is per-endpoint (cfg).
_NET = {
    'widetok-spr': ("py-ukraine-interim-amortise-deepsetXcompAtt-itemamortise-J64-widetok-spr-260919", 1, 3),
    'scale-feat':  ("py-ukraine-interim-amortise-deepsetXcompAtt-itemamortise-J64-scale-feat-260831",  0, 2),
    'groupdiff':   ("py-ukraine-interim-amortise-deepsetXcompAtt-itemamortise-J64-widetok-groupdiff-260922", 1, 3),
}
_DP_COLS = ['pid', 'pid_label', 'group', 'group_label', 'item_label', 'y', 'y_stan',
            'item_type', 'item_type_id', 'item_group_id', 'oid', 'oidt']


def _dp_interims(src, pfx):
    return sorted(int(_re.search(r'_(\d+)_data_dp1', f).group(1)) for f in _glob.glob(f"{src}/{pfx}_*_data_dp1.csv"))


def _symlink(s, d):
    if _os.path.islink(d) or _os.path.exists(d):
        (_os.remove if _os.path.islink(d) or _os.path.isfile(d) else _shutil.rmtree)(d)
    _os.symlink(s, d)


def _subset_zarr(src, pfx, k, orig_oid, ref):
    import xarray as xr
    d = xr.open_zarr(f"{src}/{pfx}_{k}_draws.zarr", group='posterior')
    yps = d['ypred'].values[:, :, orig_oid - 1]
    zp = f"{ref}/{pfx}_{k}_draws.zarr"; _shutil.rmtree(zp, ignore_errors=True)
    xr.Dataset({'ypred': (('chain', 'draw', 'ypred_dim_0'), yps)}).to_zarr(zp, group='posterior', mode='w')


def _write_dit(src, pfx, ref, endpoint_measure, strains, it):
    ditj = pd.read_csv(f"{src}/{pfx}_1_data_dit.csv").set_index('item_label')
    rows = [dict(item_label=s, item_type=it, item_type_id=1, cat_length=int(ditj.loc[s].cat_length),
                 item_label_short=np.nan, construct='HAI titre', construct_long=s,
                 item_high_label='higher_is_better', endpoint_measure=endpoint_measure,
                 cat_labels=ditj.loc[s].cat_labels) for s in strains]
    pd.DataFrame(rows).to_csv(f"{ref}/{pfx}_1_data_dit.csv", index=False)


# ---- build strategies: fill <fed>/<rho>/ref from the SVI fit (no re-fit) --------------------
def _ref_svi(src, pfx, ref, ep, interims):
    """Whole SVI fit symlinked as-is (single rho, all items; endpoint pkl already carries pps_ratio_x)."""
    _symlink(f"{src}/{pfx}_1_data_dit.csv", f"{ref}/{pfx}_1_data_dit.csv")
    for k in interims:
        for suf in (f"{k}_data_dp1.csv", f"{k}_draws.zarr", f"i{k}_regression_training.pkl"):
            _symlink(f"{src}/{pfx}_{suf}", f"{ref}/{pfx}_{suf}")


def _ref_rho_id(src, pfx, ref, ep, interims):
    """Symlink dp1/draws/dit; rewrite the endpoint pkl = the rho_id's per-draw value -> pps_ratio_x."""
    _symlink(f"{src}/{pfx}_1_data_dit.csv", f"{ref}/{pfx}_1_data_dit.csv")
    for k in interims:
        _symlink(f"{src}/{pfx}_{k}_data_dp1.csv", f"{ref}/{pfx}_{k}_data_dp1.csv")
        _symlink(f"{src}/{pfx}_{k}_draws.zarr", f"{ref}/{pfx}_{k}_draws.zarr")
        x = pd.read_pickle(f"{src}/{pfx}_i{k}_regression_training.pkl")
        r = x[x.rho_id == ep['rho_id']][['draw', 'item_label', 'item_type', 'item_high_label',
                                         'pps_rho_x', 'pps_H1_x']].copy()
        r.rename(columns={'pps_rho_x': 'pps_ratio_x'}).to_pickle(f"{ref}/{pfx}_i{k}_regression_training.pkl")


def _ref_level(src, pfx, ref, ep, interims):
    """SDY269 within-arm level: subset the joint fit to one arm's strains, group = phase (0/1)."""
    arm = ep['arm']; it = 'categorical' if 'spr' in ep['rho'] else 'out-of-7'
    dpf = pd.read_csv(f"{src}/{pfx}_{interims[-1]}_data_dp1.csv")
    strains = [s for s in pd.read_csv(f"{src}/{pfx}_1_data_dit.csv").item_label
               if s in set(dpf[dpf.group_label.str.startswith(arm)].item_label)]
    _write_dit(src, pfx, ref, ep['rho'], strains, it)
    for k in interims:
        dp = pd.read_csv(f"{src}/{pfx}_{k}_data_dp1.csv")
        e = dp[dp.group_label.str.startswith(arm)].copy()
        e['group'] = (e.group_label.str.contains('endline')).astype(int)
        e['item_label'] = pd.Categorical(e.item_label, categories=strains, ordered=True)
        e = e.sort_values(['group', 'item_label', 'pid']).reset_index(drop=True)
        e['item_label'] = e.item_label.astype(str); orig_oid = e['oid'].to_numpy()
        e['group_label'] = np.where(e.group == 0, 'Baseline', 'Endline')
        e['item_group_id'] = pd.factorize(e.group.astype(str) + '|' + e.item_label)[0] + 1
        e['item_type'] = it; e['item_type_id'] = 1
        e['oid'] = np.arange(1, len(e) + 1); e['oidt'] = np.arange(1, len(e) + 1)
        e[_DP_COLS].to_csv(f"{ref}/{pfx}_{k}_data_dp1.csv", index=False)
        _subset_zarr(src, pfx, k, orig_oid, ref)
        x = pd.read_pickle(f"{src}/{pfx}_i{k}_regression_training.pkl")
        r = x[(x.rho_label == ep['joint']) & (x.item_label.isin(strains))][
            ['draw', 'item_label', 'pps_rho_x', 'pps_H1_x']].copy()
        r['item_type'] = it; r['item_high_label'] = 'higher_is_better'
        r.rename(columns={'pps_rho_x': 'pps_ratio_x'}).to_pickle(f"{ref}/{pfx}_i{k}_regression_training.pkl")


def _ref_diffstd(src, pfx, ref, ep, interims):
    """SDY269 between-arm standardised SPR-difference on the shared strain(s); group = arm (0/1)."""
    lo_lab, hi_lab = ep['joint']
    dpf = pd.read_csv(f"{src}/{pfx}_{interims[-1]}_data_dp1.csv"); dpf['arm'] = dpf.group_label.str.split('_').str[0]
    per_arm = {a: set(g.item_label) for a, g in dpf.groupby('arm')}
    shared = sorted(set.intersection(*per_arm.values()))
    _write_dit(src, pfx, ref, 'SPR', shared, 'categorical')
    for k in interims:
        dp = pd.read_csv(f"{src}/{pfx}_{k}_data_dp1.csv")
        e = dp[dp.item_label.isin(shared) & dp.group_label.str.contains('endline')].copy()
        e['arm'] = np.where(e.group_label.str.contains('LAIV'), 0, 1)
        e = e.sort_values(['arm', 'item_label', 'pid']).reset_index(drop=True)
        orig_oid = e['oid'].to_numpy()
        e['group'] = e['arm']; e['group_label'] = np.where(e.arm == 0, 'LAIV', 'TIV')
        e['item_group_id'] = pd.factorize(e.arm.astype(str) + '|' + e.item_label)[0] + 1
        e['item_type'] = 'categorical'; e['item_type_id'] = 1
        e['oid'] = np.arange(1, len(e) + 1); e['oidt'] = np.arange(1, len(e) + 1)
        e[_DP_COLS].to_csv(f"{ref}/{pfx}_{k}_data_dp1.csv", index=False)
        _subset_zarr(src, pfx, k, orig_oid, ref)
        x = pd.read_pickle(f"{src}/{pfx}_i{k}_regression_training.pkl")
        rows = []
        for strain in shared:
            L = x[(x.rho_label == lo_lab) & (x.item_label == strain)].set_index('draw')['pps_rho_x']
            T = x[(x.rho_label == hi_lab) & (x.item_label == strain)].set_index('draw')['pps_rho_x']
            j = L.index.intersection(T.index); Lv, Tv = L.loc[j].to_numpy(), T.loc[j].to_numpy()
            pbar = np.clip((Lv + Tv) / 2, 1e-3, 1 - 1e-3); dstd = (Tv - Lv) / np.sqrt(pbar * (1 - pbar))
            rows.append(pd.DataFrame(dict(draw=j, item_label=strain, item_type='categorical',
                                          item_high_label='higher_is_better', pps_ratio_x=dstd,
                                          pps_H1_x=(dstd > 0).astype(int))))
        pd.concat(rows, ignore_index=True).to_pickle(f"{ref}/{pfx}_i{k}_regression_training.pkl")


_BUILD = {'svi': _ref_svi, 'rho_id': _ref_rho_id, 'level': _ref_level, 'diffstd': _ref_diffstd}


def federated_deploy(sb, cfg, driver=DEFAULT_DRIVER, scratch=DEPLOY_SCRATCH,
                     run_diagnostics=True, verbose=True):
    """Deploy every endpoint in `cfg` (the per-project *_startme.py CFG dict) into its federated subdir,
    then chain the combined diagnostics. `cfg` fields used: svi, fed, endpoints[net/warp/eta0/grid/build
    + build-specific], headmode?, ctag_prefix?, nref?."""
    fed = f"{sb}/{cfg['fed']}"; src = f"{sb}/{cfg['svi']}"; pfx = cfg.get('pfx', 'pcm_1_interim')
    _os.makedirs(fed, exist_ok=True); interims = _dp_interims(src, pfx)
    nref = cfg.get('nref') or max(pd.read_csv(f"{src}/{pfx}_{interims[-1]}_data_dp1.csv").pid.nunique() + 8, 64)
    headmode = cfg.get('headmode', ''); tag = cfg.get('ctag_prefix', _os.path.basename(fed)[:24])
    for ep in cfg['endpoints']:
        netdir, widetok, case_c = _NET[ep['net']]
        ref = f"{fed}/{ep['rho']}/ref"; _os.makedirs(ref, exist_ok=True)
        _BUILD[ep.get('build', 'rho_id')](src, pfx, ref, ep, interims)
        out = f"{fed}/{ep['rho']}"; ctag = f"{tag}-{ep['rho']}"
        _shutil.rmtree(f"{scratch}/ragged_deploy_cells_{ctag}", ignore_errors=True)
        env = dict(_os.environ, RAGD_BASE=f"{sb}/{netdir}", RAGD_RGE=ref, RAGD_OUT=out,
                   RAGD_WARP=ep['warp'], RAGD_WIDETOK=str(widetok), RAGD_KMAX='10', RAGD_CASE_C=str(case_c),
                   RAGD_BVM='1', RAGD_HEADFT='expand', RAGD_NREF=str(nref), RAGD_CTAG=ctag,
                   RAG_S='200', RAGD_PDFS='1', RAGD_PERINTERIM='0', RAGD_PLOTDATA='1',
                   RAGD_ETA0=str(ep['eta0']), RAGD_ETAH='0.89', RAGD_ETA0GRID='1', RAGD_ETA0GRID_VALS=ep['grid'])
        if headmode:
            env['RAGD_HEADMODE'] = headmode
        if verbose:
            print(f"\n===== deploy {ep['rho']} (net={ep['net']}, warp={ep['warp']}, head={headmode or 'plain'}, "
                  f"eta0={ep['eta0']}, NREF={nref}) -> {out}")
        _subprocess.run([_sys.executable, driver], env=env, cwd=_REPO)
    if verbose:
        print(f"\nfederated deploy done -> {fed}")
    if run_diagnostics:
        from amortiser_diag_plots import FederatedDiagnostics
        FederatedDiagnostics(sb, cfg, verbose=verbose).run()
