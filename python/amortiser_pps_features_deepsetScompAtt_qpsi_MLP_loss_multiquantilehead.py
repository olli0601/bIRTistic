"""
Amortiser: nested DeepSets over participants + attention over items
(§9.2 nested-DeepSets IRT pattern + §13.5 attention-over-items on top).
Fully data-agnostic and fully amortised over arbitrary participant
cohort sizes ``(n, m)``: caller passes raw participant tensors + masks,
class does the aggregation over participants internally.

**Architecture.**

  1. Per ``(participant i, item j)`` token = raw response(s) at that
     ``(i, j)`` position, optionally concatenated with the item
     metadata vector (item_type indicator, item_high indicator, ...).
  2. ``q_tau`` -- per-token MLP shared across ``(i, j)`` -- maps each
     token to an embedding in ``R^embed_dim``. Analog of ``q_tau`` in
     §7 DeepSets.
  3. Sum-pool over participants (mask-aware) PER ITEM -> per-item
     summary embedding of shape ``(B, J, embed_dim)``. This is the
     "learned per-item sufficient statistic".
  4. ``q_tok`` -- MLP mixing the (pool_x, pool_z) per-item summaries
     into a single per-item embedding.
  5. Self-attention across items -> per-item embedding conditioned on
     all other items.
  6. Gather at the queried item + optional aux (sizes, etc.) ->
     ``q_psi`` head -> multi-quantile prediction of ``rho`` for that
     item.

**Batch schema (dict).**

    ``x_responses``:   ``(B, N_max_x, J, R)`` -- observed cohort's
                        raw per-(participant, item) feature vector.
    ``mask_x``:        ``(B, N_max_x)``        -- 1 real, 0 padded.
    ``z_responses``:   ``(B, N_max_z, J, R)``  -- future cohort's raw
                        per-(participant, item) feature vector.
    ``mask_z``:        ``(B, N_max_z)``        -- 1 real, 0 padded.
    ``item_metadata``: ``(B, J, M)`` optional  -- fixed per-item scalar
                        features (item_type indicator, item_high
                        indicator, ...); if present, broadcast over the
                        participant axis and concatenated with the raw
                        response before ``q_tau``.
    ``item_mask``:     ``(B, J)`` optional     -- 1 for real item, 0 for
                        padded; defaults to all ones.
    ``query_idx``:     ``(B,)`` int            -- which item to predict.
    ``aux``:           ``(B, A)`` optional     -- scalar extras fed to
                        the head (e.g. ``(n / N_max, m / N_max)``).

The response feature dimension ``R`` is arbitrary. Examples:

  - Ukraine PCM: two responses per (participant, item) -- ``R = 2`` =
    (baseline_response, endline_response).
  - MVN (§3.3): one continuous scalar per (participant, item) -- ``R = 1``.
  - Categorical: one-hot response per (participant, item) -- ``R = K``.
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


class Amortiser_PPS_features_deepsetScompAtt_qpsi_MLP_loss_multiquantilehead(nn.Module):
    """Nested DeepSets over participants + self-attention over items."""

    q_tau_hidden: tuple = (32, 32)
    q_tok_hidden: tuple = (32, 32)
    embed_dim: int = 32
    hidden_dims: tuple = (64, 64)
    num_quantiles: int = 5
    default_pps_ProbH1_lwr_quantiles_mesh: tuple = (
        0.05, 0.25, 0.5, 0.75, 0.95,
    )

    def setup(self):
        self.q_tau = _MLP(dims=(*self.q_tau_hidden, self.embed_dim))
        self.q_tok = _MLP(dims=(*self.q_tok_hidden, self.embed_dim))
        self.q_psi = _MLP(dims=(*self.hidden_dims, self.num_quantiles))
        # Fix-B: learned scalar bias on attention score at queried key
        # (see itemScompAtt / §13.6).
        self.query_bias_alpha = self.param(
            'query_bias_alpha', nn.initializers.zeros, (),
        )

    def _augment_with_meta(self, resp, meta):
        """resp: (B, N, J, R). meta: (B, J, M). Return (B, N, J, R+M)."""
        if meta is None:
            return resp
        B, N, J, R = resp.shape
        meta_b = jnp.broadcast_to(meta[:, None, :, :], (B, N, J, meta.shape[-1]))
        return jnp.concatenate([resp, meta_b], axis=-1)

    def _sum_pool_participants(self, e_ij, mask):
        """e_ij: (B, N, J, E). mask: (B, N). Return (B, J, E).
        Mean-pool (mask-aware): sum / count. Keeps the pooled magnitude
        scale-invariant to N so q_tok doesn't saturate as the cohort
        size varies from tens to thousands."""
        numer = jnp.einsum('bnje,bn->bje', e_ij, mask)
        denom = jnp.maximum(mask.sum(axis=1, keepdims=True), 1.0)
        return numer / denom[:, :, None]

    def _self_attention(self, h, mask, query_idx):
        """h: (B, J, E). mask: (B, J). Self-attention with fix-B
        learned scalar bias on queried key."""
        scores = jnp.einsum('bqd,bkd->bqk', h, h)
        scores = scores / jnp.sqrt(float(self.embed_dim))
        is_query = jax.nn.one_hot(
            query_idx.astype(jnp.int32), h.shape[1], dtype=h.dtype,
        )
        scores = scores + self.query_bias_alpha * is_query[:, None, :]
        mask_k = mask[:, None, :]
        scores = jnp.where(mask_k > 0, scores, -1e9)
        w = jax.nn.softmax(scores, axis=-1)
        return jnp.einsum('bqk,bkd->bqd', w, h)

    def __call__(self, batch):
        x_resp = batch['x_responses']                           # (B, N_x, J, R)
        z_resp = batch['z_responses']                           # (B, N_z, J, R)
        mask_x = batch['mask_x']                                # (B, N_x)
        mask_z = batch['mask_z']                                # (B, N_z)
        meta = batch.get('item_metadata', None)                 # (B, J, M) or None
        x_feats = self._augment_with_meta(x_resp, meta)         # (B, N_x, J, R+M)
        z_feats = self._augment_with_meta(z_resp, meta)         # (B, N_z, J, R+M)
        # Per-(participant, item) embedding.
        e_x = self.q_tau(x_feats)                               # (B, N_x, J, E)
        e_z = self.q_tau(z_feats)                               # (B, N_z, J, E)
        # Sum-pool over participants per item.
        pool_x_j = self._sum_pool_participants(e_x, mask_x)     # (B, J, E)
        pool_z_j = self._sum_pool_participants(e_z, mask_z)     # (B, J, E)
        # Per-item token embedding.
        tok = jnp.concatenate([pool_x_j, pool_z_j], axis=-1)    # (B, J, 2E)
        h = self.q_tok(tok)                                     # (B, J, E)
        # Attention across items.
        item_mask = batch.get(
            'item_mask', jnp.ones(h.shape[:2], dtype=h.dtype),
        )
        idx = batch['query_idx'].astype(jnp.int32)              # (B,)
        h_attn = self._self_attention(h, item_mask, query_idx=idx)  # (B, J, E)
        # Gather at query.
        B = h_attn.shape[0]
        b_range = jnp.arange(B, dtype=jnp.int32)
        h_query = h_attn[b_range, idx]                          # (B, E)
        tok_query = h[b_range, idx]                             # (B, E) -- pre-attn per-item token
        parts = [h_query, tok_query]
        aux = batch.get('aux', None)
        if aux is not None:
            parts.append(aux)
        head_in = jnp.concatenate(parts, axis=-1)
        return self.q_psi(head_in)                              # (B, num_quantiles)
