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
