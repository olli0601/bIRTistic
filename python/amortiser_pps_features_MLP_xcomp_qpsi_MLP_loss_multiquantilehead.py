"""
Amortiser variant 3: nested DeepSets encoder over
(participant i, component j) tokens + K-row query + multi-quantile
pinball loss (§12.10 of ``dev/amortised_decision_making.md``).

**xcomp variant** -- exploits the known cross-component covariance
``K`` via a two-scalar per-token feature ``(K[j*, j], y_{i, j})``.

Nested DeepSets:

- Inner ``q_tau_inner``: per (participant i, component j) token
  ``(K[j*, j], y_{i, j}) in R^2`` -> embedding in R^embed_inner.
  Same per-token MLP applied everywhere (shared weights ->
  permutation equivariance across both axes).
- Sum-pool over j (all J present per participant) -> per-i embedding
  ``h_i in R^embed_inner``.
- Sum-pool over i with participant mask -> set-of-participants
  embedding ``pooled in R^embed_inner``. Done twice, once for ``x``
  once for ``z``, with shared inner encoder.

Head ``q_psi``: concatenates ``(pooled_x, pooled_z, sizes,
K[j*, j*])`` and predicts K quantiles of ``mu_{j*} - mu_0``.

J-invariance: the K-row query is the *only* thing the head learns
about the queried component -- there is no positional or lookup
embedding of j*, so a single trained net handles any J and any K
family the training distribution covers.

Batch schema (dict):
    ``x``:       ``(B, N_max, J)`` -- padded observed participant vectors.
    ``mask_x``:  ``(B, N_max)``    -- 1 for real participant, 0 padded.
    ``z``:       ``(B, N_max, J)`` -- padded future participant vectors.
    ``mask_z``:  ``(B, N_max)``    -- 1 real, 0 padded.
    ``sizes``:   ``(B, 2)``        -- ``(n / N_max, m / N_max)``.
    ``k_row``:   ``(B, J)``        -- ``K[j*, :]``, one row per batch elem.
    ``k_diag``:  ``(B, 1)``        -- ``K[j*, j*]`` scalar.

Shared amortiser utilities (``train``,
``predict_amortised_p_h1_for_one_xz``, ``save_trained_model`` /
``load_fitted_model``) from :mod:`amortiser_common` work unchanged.
"""

import jax.numpy as jnp
from flax import linen as nn


class _MLP(nn.Module):
    """SiLU MLP with a linear read-out. Shared between the per-token
    inner encoder ``q_tau_inner`` and the head ``q_psi``."""
    dims: tuple

    @nn.compact
    def __call__(self, x):
        for d in self.dims[:-1]:
            x = nn.silu(nn.Dense(d)(x))
        return nn.Dense(self.dims[-1])(x)


class Amortiser_PPS_features_MLP_xcomp_qpsi_MLP_loss_multiquantilehead(nn.Module):
    """Nested DeepSets over (i, j) + K-row query + MLP head."""

    q_tau_inner_hidden: tuple = (32, 32)
    embed_inner: int = 32
    hidden_dims: tuple = (128, 128, 64)
    num_quantiles: int = 5
    default_pps_ProbH1_lwr_quantiles_mesh: tuple = (
        0.05, 0.25, 0.5, 0.75, 0.95,
    )

    def setup(self):
        self.q_tau_inner = _MLP(
            dims=(*self.q_tau_inner_hidden, self.embed_inner),
        )
        self.q_psi = _MLP(dims=(*self.hidden_dims, self.num_quantiles))

    def _encode(self, y, mask_participants, k_row):
        """y: (B, N, J). mask_participants: (B, N). k_row: (B, J)."""
        B, N, J = y.shape
        k_row_bnj = jnp.broadcast_to(k_row[:, None, :], (B, N, J))
        tok = jnp.stack([k_row_bnj, y], axis=-1)          # (B, N, J, 2)
        h_ij = self.q_tau_inner(tok)                       # (B, N, J, E)
        h_i = h_ij.sum(axis=2)                             # (B, N, E)
        h_i = h_i * mask_participants[..., None]           # (B, N, E)
        return h_i.sum(axis=1)                             # (B, E)

    def __call__(self, batch):
        pooled_x = self._encode(
            batch['x'], batch['mask_x'], batch['k_row'],
        )
        pooled_z = self._encode(
            batch['z'], batch['mask_z'], batch['k_row'],
        )
        head_in = jnp.concatenate(
            [pooled_x, pooled_z, batch['sizes'], batch['k_diag']],
            axis=-1,
        )                                                   # (B, 2E + 3)
        return self.q_psi(head_in)
