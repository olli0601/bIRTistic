"""
Amortiser variant 2: end-to-end learnable DeepSets encoder ``q_tau``
(per-item MLP) + sum-pooling + MLP head ``q_psi`` + multi-quantile
pinball loss (§7-§11 of ``dev/amortised_decision_making.md``, extension
of §11.6).

Batch schema (dict):
    ``x``:       ``(B, N_max, item_dim)`` — padded current-data items.
    ``mask_x``:  ``(B, N_max)`` — 1 for real item, 0 for padded.
    ``z``:       ``(B, N_max, item_dim)`` — padded future-data items.
    ``mask_z``:  ``(B, N_max)`` — 1 for real item, 0 for padded.
    ``sizes``:   ``(B, 2)`` — normalised cohort sizes ``(n / N_max,
                             m / N_max)``.

The same ``q_tau`` MLP is applied to every item in both ``x`` and ``z``
(shared weights → permutation invariance). Masked positions are zeroed
before sum-pooling. Pooled ``x`` and ``z`` embeddings are concatenated
with ``sizes`` and passed to a second MLP ``q_psi`` producing the
predicted quantile grid.

Shared amortiser utilities (``train``,
``predict_amortised_p_h1_for_one_xz``, ``save_trained_model`` /
``load_fitted_model``) live in :mod:`amortiser_common`; the same
``train`` / ``predict_amortised_p_h1_for_one_xz`` calls work for this
class and for the features-fixed sibling.
"""

import jax.numpy as jnp
from flax import linen as nn


class _MLP(nn.Module):
    """SiLU-activated MLP with a linear read-out. Shared between
    ``q_tau`` (per-item encoder) and ``q_psi`` (head)."""
    dims: tuple

    @nn.compact
    def __call__(self, x):
        for d in self.dims[:-1]:
            x = nn.silu(nn.Dense(d)(x))
        return nn.Dense(self.dims[-1])(x)


class Amortiser_PPS_features_MLP_qpsi_MLP_loss_multiquantilehead(nn.Module):
    """Learnable per-item encoder + sum-pooling + MLP head + multi-
    quantile pinball loss."""
    # Inner q_tau encoder.
    q_tau_hidden_dims: tuple = (32, 32)
    embed_dim: int = 16
    # Outer q_psi head.
    hidden_dims: tuple = (128, 128, 64)
    num_quantiles: int = 11
    default_pps_ProbH1_lwr_quantiles_mesh: tuple = (
        0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95,
    )

    def setup(self):
        self.q_tau = _MLP(dims=(*self.q_tau_hidden_dims, self.embed_dim))
        self.q_psi = _MLP(dims=(*self.hidden_dims, self.num_quantiles))

    def __call__(self, batch):
        # Per-item embeddings (SAME q_tau for x and z -> shared weights).
        emb_x = self.q_tau(batch['x'])                     # (B, N, embed_dim)
        emb_z = self.q_tau(batch['z'])
        # Sum-pool with padding mask; expand mask to (B, N, 1).
        mask_x = batch['mask_x'][..., None]
        mask_z = batch['mask_z'][..., None]
        pooled_x = (emb_x * mask_x).sum(axis=1)            # (B, embed_dim)
        pooled_z = (emb_z * mask_z).sum(axis=1)
        features = jnp.concatenate(
            [pooled_x, pooled_z, batch['sizes']], axis=-1,
        )                                                  # (B, 2*embed_dim + 2)
        return self.q_psi(features)
