"""
Ragged-participant-axis variant of the nested-DeepSets + cross-attention
amortiser (§14.4.11 of dev/amortised_decision_making.md, component D).

Unlike `..._deepsetXcompAtt_...` which pads every cohort to a fixed
N_max and carries a 0/1 participant mask, this variant consumes a
**ragged** participant axis: participants from all batch elements are
concatenated into one flat axis and pooled with a segmented mean
(`jax.ops.segment_sum` over per-participant segment ids). There is no
padding and no mask — the encoder sees exactly the real n participants,
and the 1/n pooling delivers the correct posterior precision (BvM rate)
to the head. This is the exchangeable summary-network construction used
by SBI to learn contraction across dataset size (Radev et al. 2020;
Chan et al. 2018).

Cross-attention over items and the multi-quantile head are unchanged
from the padded variant. Q queries per batch element are handled inside
the module (query_idx is (B, Q)).

Batch schema (B = num batch elements, T = total x-participants across
the batch, U = total z-participants, J items, R response features,
M item-metadata features, A aux scalars):

  x_flat        (T, J, R)     concatenated observed-cohort participants
  x_seg         (T,)          int32, batch-element id in [0, B) per x-participant
  z_flat        (U, J, R)     concatenated future-cohort participants
  z_seg         (U,)          int32 in [0, B)
  item_metadata (B, J, M)
  query_idx     (B, Q)        int32 queried item per (batch element, query)
  aux           (B, A)

Output: (B, Q, num_quantiles).
"""

import jax
import jax.numpy as jnp
from flax import linen as nn


class _MLP(nn.Module):
    """SiLU MLP with linear read-out."""
    dims: tuple

    @nn.compact
    def __call__(self, x):
        for d in self.dims[:-1]:
            x = nn.silu(nn.Dense(d)(x))
        return nn.Dense(self.dims[-1])(x)


class Amortiser_PPS_features_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead(nn.Module):
    """Nested DeepSets (ragged participant axis) + cross-attention over items."""

    q_tau_hidden: tuple = (32, 32)
    q_tok_hidden: tuple = (32, 32)
    q_query_hidden: tuple = (32, 32)
    embed_dim: int = 32
    hidden_dims: tuple = (64, 64)
    num_quantiles: int = 5
    default_pps_ProbH1_lwr_quantiles_mesh: tuple = (0.05, 0.25, 0.5, 0.75, 0.95)

    def setup(self):
        self.q_tau = _MLP(dims=(*self.q_tau_hidden, self.embed_dim))
        self.q_tok = _MLP(dims=(*self.q_tok_hidden, self.embed_dim))
        self.q_query = _MLP(dims=(*self.q_query_hidden, self.embed_dim))
        self.q_psi = _MLP(dims=(*self.hidden_dims, self.num_quantiles))

    def _segment_mean_pool(self, flat, seg, meta, B):
        """flat (T, J, R), seg (T,), meta (B, J, M) -> (B, J, E) ragged mean.

        No padding, no mask: T is the true total participant count. Each
        participant is augmented with its batch element's item metadata,
        embedded by the shared q_tau, then averaged within its segment.
        """
        meta_t = meta[seg]                                   # (T, J, M)
        feats = jnp.concatenate([flat, meta_t], axis=-1)     # (T, J, R+M)
        e = self.q_tau(feats)                                # (T, J, E)
        summ = jax.ops.segment_sum(e, seg, num_segments=B)   # (B, J, E)
        cnt = jax.ops.segment_sum(
            jnp.ones((flat.shape[0], 1, 1), e.dtype), seg, num_segments=B,
        )                                                    # (B, 1, 1)
        return summ / jnp.maximum(cnt, 1.0)

    def __call__(self, batch):
        meta = batch['item_metadata']                        # (B, J, M)
        B, J = meta.shape[0], meta.shape[1]
        pool_x = self._segment_mean_pool(batch['x_flat'], batch['x_seg'], meta, B)
        pool_z = self._segment_mean_pool(batch['z_flat'], batch['z_seg'], meta, B)
        tok = jnp.concatenate([pool_x, pool_z], axis=-1)     # (B, J, 2E)
        h = self.q_tok(tok)                                  # (B, J, E)

        qidx = batch['query_idx'].astype(jnp.int32)          # (B, Q)
        Q = qidx.shape[1]
        E = h.shape[-1]
        # gather the queried item's embedding per (b, q): (B, Q, E)
        gather_idx = jnp.broadcast_to(qidx[:, :, None], (B, Q, E))
        h_query = jnp.take_along_axis(h, gather_idx, axis=1)  # (B, Q, E)
        q = self.q_query(h_query)                            # (B, Q, E)

        # cross-attention: 1 query per (b, q) vs the J item embeddings of b
        scores = jnp.einsum('bqe,bje->bqj', q, h) / jnp.sqrt(float(self.embed_dim))
        w = jax.nn.softmax(scores, axis=-1)                  # (B, Q, J)
        attn_out = jnp.einsum('bqj,bje->bqe', w, h)          # (B, Q, E)

        parts = [attn_out, h_query]
        aux = batch.get('aux', None)
        if aux is not None:
            parts.append(jnp.broadcast_to(aux[:, None, :], (B, Q, aux.shape[-1])))
        head_in = jnp.concatenate(parts, axis=-1)            # (B, Q, E+E+A)
        # expose for head-only recalibration (§14.2.5)
        self.sow('intermediates', 'head_in', head_in)
        return self.q_psi(head_in)                           # (B, Q, num_quantiles)
