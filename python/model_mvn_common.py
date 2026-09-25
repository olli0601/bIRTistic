"""Shared MVN-benchmark deploy helpers (§13).

The MVN drivers (the plain deepsetXcompAtt deploy, the contraction+BvM driver, and the by-architecture
diagnostics) all build the SAME per-interim deploy inputs from the interim pkl: the wide x cohort +
the (S, m, J) posterior-predictive z draws, and the padded nested-DeepSet batch fed to the amortiser.
This module holds that construction once so the drivers don't each re-implement it.

These are pure array builders (no net, no randomness): given identical inputs they return byte-identical
arrays, so the amortiser forward pass is unchanged. The net forward + output handling (predict vs the
head_in-capturing forward, quantile accumulation, PPS emission) stays in each driver.

Design note: the padded batch is split into `deepset_static_batch` (everything independent of the z
draw — build once per interim) + `with_z` (inject one posterior draw's z), so the S-loop only rebuilds
the z tensor, matching the plain driver's memory optimisation (the x tensor is the large one)."""
import numpy as np


def xz_arrays(dpi, zi, J, S, m):
    """From one interim block: x_wide (n, J) observed responses + ypred (S, m, J) posterior-predictive
    z draws, aligned to the zi ypred columns. Identical extraction used by every MVN deploy path."""
    x_wide = (dpi.pivot_table(index='pid', columns='j', values='y').sort_index()
              .reindex(columns=range(J)).to_numpy(np.float64))
    cols = sorted([c for c in zi.columns if c.startswith('ypred_')],
                  key=lambda c: int(c.split('_')[1]))[:S]
    zis = zi.sort_values(['pid', 'j']).reset_index(drop=True)
    ypred = zis[cols].to_numpy().T.reshape(S, m, J)
    return x_wide, ypred


def deepset_static_batch(x_wide, K, K_diag, n, m, N_max, J):
    """The z-independent part of the padded nested-DeepSet batch (build once per interim): the padded +
    broadcast x cohort, the x/z masks, per-item metadata K, the query index, and the aux vector
    ([n/N_max, m/N_max, K_diag]). Combine with each posterior draw via `with_z`. Byte-identical to the
    inline construction in the MVN deploy drivers."""
    x_pad = np.zeros((N_max, J), np.float32); x_pad[:n] = x_wide
    mx = np.zeros(N_max, np.float32); mx[:n] = 1.0
    mz = np.zeros(N_max, np.float32); mz[:m] = 1.0
    sizes = np.array([n / N_max, m / N_max], np.float32)
    return dict(
        x_responses=np.broadcast_to(x_pad[None, :, :, None], (J, N_max, J, 1)).astype(np.float32),
        mask_x=np.broadcast_to(mx[None], (J, N_max)).astype(np.float32),
        mask_z=np.broadcast_to(mz[None], (J, N_max)).astype(np.float32),
        item_metadata=K.astype(np.float32)[..., None],
        query_idx=np.arange(J, dtype=np.int32),
        aux=np.concatenate([np.broadcast_to(sizes[None, :], (J, 2)), K_diag[:, None]], -1).astype(np.float32),
        _N_max=int(N_max), _J=int(J),
    )


def with_z(static, z_s):
    """Add one posterior draw's z responses z_s (m, J) to a `deepset_static_batch` -> full net batch.
    Returns a fresh dict (the shared static arrays are reused as views; only z_responses is per-draw)."""
    N_max, J = static['_N_max'], static['_J']
    z_pad = np.zeros((N_max, J), np.float32); z_pad[:z_s.shape[0]] = z_s
    b = {k: v for k, v in static.items() if not k.startswith('_')}
    b['z_responses'] = np.broadcast_to(z_pad[None, :, :, None], (J, N_max, J, 1)).astype(np.float32)
    return b
