"""
Amortiser: self-attention over items on caller-supplied per-item feature
tokens (§12.11 MVN + §12.12 IRT extension of
``dev/amortised_decision_making.md``). Data-agnostic: works for any
case study where the caller has pre-pooled the raw participant
observations into per-item summary features (e.g. Strong-Oakley
``w_baseline`` / ``w_endline`` for IRT, ``sum_i y_{i, j}`` for MVN
with the K-row query baked into the token).

**Where the attention lives.** Between the per-item encoder ``q_tok``
and the per-item head ``q_psi`` -- it replaces the sum-pool that would
sit in a plain DeepSets over items with an attention-weighted mix. The
head then gathers the queried item's post-attention embedding and
predicts its quantile grid. ``q_tau`` (per-observation encoder over
participants) is NOT part of this class -- the caller is assumed to
have run whatever per-item aggregation over participants is appropriate
for the case study upstream. For a fully-amortised per-participant
variant see ``amortiser_pps_features_deepsetXcompAtt_...``.

Batch schema (dict):
    ``tokens``:    ``(B, J, F)``    -- per-item feature vector, F
                                       configurable per caller.
    ``mask``:      ``(B, J)``       -- 1 for real item, 0 for padded.
    ``query_idx``: ``(B,)`` int     -- which item to predict.
    ``aux``:       ``(B, A)`` float -- optional scalar extras fed to
                                       the head (cohort sizes, K-diag,
                                       ...). If absent, no aux.

Loss target: per-batch-elem prediction of ``rho`` at the queried
item. Same multi-quantile pinball loss as the other amortiser variants.
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


class Amortiser_PPS_features_itemScompAtt_qpsi_MLP_loss_multiquantilehead(nn.Module):
    """Item-axis self-attention amortiser (data-agnostic)."""

    q_tok_hidden: tuple = (32, 32)
    embed_dim: int = 32
    hidden_dims: tuple = (64, 64)
    num_quantiles: int = 5
    default_pps_ProbH1_lwr_quantiles_mesh: tuple = (
        0.05, 0.25, 0.5, 0.75, 0.95,
    )

    def setup(self):
        self.q_tok = _MLP(dims=(*self.q_tok_hidden, self.embed_dim))
        self.q_psi = _MLP(dims=(*self.hidden_dims, self.num_quantiles))
        # Learned scalar bias added to the attention score at the queried
        # key position (query-conditional attention weights; §12.11 fix B).
        self.query_bias_alpha = self.param(
            'query_bias_alpha', nn.initializers.zeros, (),
        )

    def _self_attention(self, h, mask, query_idx=None):
        """h: (B, J, E). mask: (B, J). Single-head self-attention with
        key-side masking + learned scalar bias on the queried key so
        every query preferentially attends to the queried token."""
        scores = jnp.einsum('bqd,bkd->bqk', h, h)
        scores = scores / jnp.sqrt(float(self.embed_dim))
        # Query-bias: learn a scalar alpha, add alpha * one_hot(query_idx)
        # to every row of the scores. Effect: attention softmax at every
        # query row is nudged toward the queried key. Restores the
        # query-conditional attention weights that the old cross-attn
        # design got for free.
        if query_idx is not None:
            is_query = jax.nn.one_hot(
                query_idx.astype(jnp.int32), h.shape[1], dtype=h.dtype,
            )                                                    # (B, J)
            scores = scores + self.query_bias_alpha * is_query[:, None, :]
        mask_k = mask[:, None, :]
        scores = jnp.where(mask_k > 0, scores, -1e9)
        w = jax.nn.softmax(scores, axis=-1)
        return jnp.einsum('bqk,bkd->bqd', w, h)

    def __call__(self, batch):
        tokens = batch['tokens']                                # (B, J, F)
        mask = batch['mask']                                    # (B, J)
        h = self.q_tok(tokens)                                  # (B, J, E)
        idx = batch['query_idx'].astype(jnp.int32)              # (B,)
        h_attn = self._self_attention(h, mask, query_idx=idx)   # (B, J, E)
        B = h_attn.shape[0]
        b_range = jnp.arange(B, dtype=jnp.int32)
        h_query = h_attn[b_range, idx]                          # (B, E)
        tok_query = tokens[b_range, idx]                        # (B, F)
        parts = [h_query, tok_query]
        aux = batch.get('aux', None)
        if aux is not None:
            parts.append(aux)
        head_in = jnp.concatenate(parts, axis=-1)
        return self.q_psi(head_in)                              # (B, num_quantiles)
