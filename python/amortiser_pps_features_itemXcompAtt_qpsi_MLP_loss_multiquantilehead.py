"""
Amortiser: cross-attention over items on caller-supplied per-item
feature tokens (§12.13 of ``dev/amortised_decision_making.md``).
Data-agnostic sibling of the self-attention variant in
``amortiser_pps_features_itemScompAtt_qpsi_MLP_loss_multiquantilehead``.

**Design.** Single learned query per batch element, built from the
queried item's own raw token features via a small MLP ``q_query``.
Cross-attention scores that query against the J per-item embeddings
``q_tok(t_j)``. Attention output = weighted sum over items with weights
telling "how much information does item j contribute to the queried
item j*". The head then reads out the multi-quantile grid from
``(attn_output, tok_query, aux)``.

**Where the K-conditioning goes.** Per-token features ``t_j`` still
carry ``K_{j*, j}`` (baked in by the MVN caller). The queried token
``t_{j*}`` carries ``K_{j*, j*} = 1`` (self-correlation, maximum);
its projection through ``q_query`` becomes the attention query. Tokens
whose embeddings are similar to the queried one (large K-row) attract
high attention weight -- the "K tells us how much information can be
gleaned from item j to item j*" intuition is preserved and made
explicit.

Batch schema (dict; same as itemScompAtt so callers can swap classes
without changing batch construction):

    ``tokens``:    ``(B, J, F)``    -- per-item feature vector.
    ``mask``:      ``(B, J)``       -- 1 for real item, 0 for padded.
    ``query_idx``: ``(B,)`` int     -- which item to predict.
    ``aux``:       ``(B, A)`` float -- optional scalar extras fed to
                                       the head. If absent, no aux.
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


class Amortiser_PPS_features_itemXcompAtt_qpsi_MLP_loss_multiquantilehead(nn.Module):
    """Item-axis cross-attention amortiser (data-agnostic).

    Attention is 1 query (built from the queried item's own token
    features) vs J keys/values (per-item token embeddings). Number of
    attention scores per batch elem is ``J`` (vs ``J^2`` for
    itemScompAtt), so this variant is more sample-efficient at large J.
    """

    q_tok_hidden: tuple = (32, 32)
    q_query_hidden: tuple = (32, 32)
    embed_dim: int = 32
    hidden_dims: tuple = (64, 64)
    num_quantiles: int = 5
    default_pps_ProbH1_lwr_quantiles_mesh: tuple = (
        0.05, 0.25, 0.5, 0.75, 0.95,
    )

    def setup(self):
        self.q_tok = _MLP(dims=(*self.q_tok_hidden, self.embed_dim))
        self.q_query = _MLP(dims=(*self.q_query_hidden, self.embed_dim))
        self.q_psi = _MLP(dims=(*self.hidden_dims, self.num_quantiles))

    def _cross_attention(self, q, h, mask):
        """q: (B, E). h: (B, J, E). mask: (B, J). Single query vs J
        keys/values."""
        # scores[b, j] = <q_b, h_{b, j}> / sqrt(E).
        scores = jnp.einsum('bd,bjd->bj', q, h) / jnp.sqrt(float(self.embed_dim))
        scores = jnp.where(mask > 0, scores, -1e9)
        w = jax.nn.softmax(scores, axis=-1)                     # (B, J)
        return jnp.einsum('bj,bjd->bd', w, h)                    # (B, E)

    def __call__(self, batch):
        tokens = batch['tokens']                                # (B, J, F)
        mask = batch['mask']                                    # (B, J)
        idx = batch['query_idx'].astype(jnp.int32)              # (B,)
        B = tokens.shape[0]
        b_range = jnp.arange(B, dtype=jnp.int32)
        # Per-item token embeddings (keys / values).
        h = self.q_tok(tokens)                                  # (B, J, E)
        # Query built from the queried item's raw features (approach ii of
        # §12.13). J-invariant (no fixed lookup table).
        tok_query = tokens[b_range, idx]                        # (B, F)
        q = self.q_query(tok_query)                             # (B, E)
        # Cross-attention.
        attn_out = self._cross_attention(q, h, mask)            # (B, E)
        # Head.
        parts = [attn_out, tok_query]
        aux = batch.get('aux', None)
        if aux is not None:
            parts.append(aux)
        head_in = jnp.concatenate(parts, axis=-1)
        return self.q_psi(head_in)                              # (B, num_quantiles)
