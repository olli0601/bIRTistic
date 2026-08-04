"""
Amortiser: nested DeepSets over participants + cross-attention across
items (approach ii of §12.13). Fully data-agnostic; fully amortised
over arbitrary participant cohort sizes ``(n, m)``.

**Architecture.**

  1. Per (participant i, item j) token = raw response(s) at that
     (i, j) position, optionally concatenated with per-item metadata
     broadcast over the participant axis.
  2. ``q_tau`` -- shared per-token MLP -- maps to an embedding
     in ``R^embed_dim``.
  3. Sum-pool over participants (mask-aware) PER ITEM -> per-item
     summary embedding of shape ``(B, J, embed_dim)`` -- the learned
     per-item sufficient statistic aggregated across participants.
  4. ``q_tok`` -- MLP mixing ``(pool_x, pool_z)`` per item -> per-item
     embedding ``h_j = q_tok(concat(pool_x_j, pool_z_j)) in R^embed_dim``.
  5. **Cross-attention (approach ii of §12.13).** Query built from
     the queried item's own per-item embedding via ``q_query`` MLP;
     attention over the J per-item embeddings. Number of scores per
     batch elem: J.
  6. Gather-free head: attention output + queried item's per-item
     embedding + optional aux -> ``q_psi`` -> multi-quantile prediction.

**Batch schema (dict; same as deepsetScompAtt so callers can swap).**

    ``x_responses``:   ``(B, N_max_x, J, R)`` -- raw per-(i, j) features.
    ``mask_x``:        ``(B, N_max_x)``       -- 1 real, 0 padded.
    ``z_responses``:   ``(B, N_max_z, J, R)`` -- raw per-(i, j) features.
    ``mask_z``:        ``(B, N_max_z)``       -- 1 real, 0 padded.
    ``item_metadata``: ``(B, J, M)`` optional -- fixed per-item features.
    ``item_mask``:     ``(B, J)`` optional    -- 1 real item, 0 padded.
    ``query_idx``:     ``(B,)`` int           -- which item to predict.
    ``aux``:           ``(B, A)`` optional    -- scalar extras for head.
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


class Amortiser_PPS_features_deepsetXcompAtt_qpsi_MLP_loss_multiquantilehead(nn.Module):
    """Nested DeepSets over participants + cross-attention over items.

    Cross-attention (1 query per batch elem vs J keys/values) built from
    the queried item's per-item embedding. Number of attention scores
    per batch elem is J (vs J^2 for deepsetScompAtt), so more
    sample-efficient at large J.
    """

    q_tau_hidden: tuple = (32, 32)
    q_tok_hidden: tuple = (32, 32)
    q_query_hidden: tuple = (32, 32)
    embed_dim: int = 32
    hidden_dims: tuple = (64, 64)
    num_quantiles: int = 5
    default_pps_ProbH1_lwr_quantiles_mesh: tuple = (
        0.05, 0.25, 0.5, 0.75, 0.95,
    )

    def setup(self):
        self.q_tau = _MLP(dims=(*self.q_tau_hidden, self.embed_dim))
        self.q_tok = _MLP(dims=(*self.q_tok_hidden, self.embed_dim))
        self.q_query = _MLP(dims=(*self.q_query_hidden, self.embed_dim))
        self.q_psi = _MLP(dims=(*self.hidden_dims, self.num_quantiles))

    def _augment_with_meta(self, resp, meta):
        if meta is None:
            return resp
        B, N, J, R = resp.shape
        meta_b = jnp.broadcast_to(meta[:, None, :, :], (B, N, J, meta.shape[-1]))
        return jnp.concatenate([resp, meta_b], axis=-1)

    def _sum_pool_participants(self, e_ij, mask):
        numer = jnp.einsum('bnje,bn->bje', e_ij, mask)
        denom = jnp.maximum(mask.sum(axis=1, keepdims=True), 1.0)
        return numer / denom[:, :, None]

    def _cross_attention(self, q, h, mask):
        """q: (B, E). h: (B, J, E). mask: (B, J). 1 query vs J keys/values."""
        scores = jnp.einsum('bd,bjd->bj', q, h) / jnp.sqrt(float(self.embed_dim))
        scores = jnp.where(mask > 0, scores, -1e9)
        w = jax.nn.softmax(scores, axis=-1)
        return jnp.einsum('bj,bjd->bd', w, h)

    def __call__(self, batch):
        x_resp = batch['x_responses']                           # (B, N_x, J, R)
        z_resp = batch['z_responses']                           # (B, N_z, J, R)
        mask_x = batch['mask_x']                                # (B, N_x)
        mask_z = batch['mask_z']                                # (B, N_z)
        meta = batch.get('item_metadata', None)                 # (B, J, M) or None
        x_feats = self._augment_with_meta(x_resp, meta)         # (B, N_x, J, R+M)
        z_feats = self._augment_with_meta(z_resp, meta)         # (B, N_z, J, R+M)
        e_x = self.q_tau(x_feats)                               # (B, N_x, J, E)
        e_z = self.q_tau(z_feats)                               # (B, N_z, J, E)
        pool_x_j = self._sum_pool_participants(e_x, mask_x)     # (B, J, E)
        pool_z_j = self._sum_pool_participants(e_z, mask_z)     # (B, J, E)
        tok = jnp.concatenate([pool_x_j, pool_z_j], axis=-1)    # (B, J, 2E)
        h = self.q_tok(tok)                                     # (B, J, E)
        idx = batch['query_idx'].astype(jnp.int32)              # (B,)
        B = h.shape[0]
        b_range = jnp.arange(B, dtype=jnp.int32)
        # Cross-attn query: built from queried item's per-item embedding
        # (approach ii of §12.13, deepset version). The queried item's
        # embedding already encodes both cohorts' per-item pool via q_tok.
        h_query = h[b_range, idx]                               # (B, E)
        q = self.q_query(h_query)                               # (B, E)
        item_mask = batch.get(
            'item_mask', jnp.ones(h.shape[:2], dtype=h.dtype),
        )
        attn_out = self._cross_attention(q, h, item_mask)       # (B, E)
        parts = [attn_out, h_query]
        aux = batch.get('aux', None)
        if aux is not None:
            parts.append(aux)
        head_in = jnp.concatenate(parts, axis=-1)
        return self.q_psi(head_in)                              # (B, num_quantiles)
