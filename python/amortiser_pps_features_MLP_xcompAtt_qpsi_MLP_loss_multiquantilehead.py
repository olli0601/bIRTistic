"""
Amortiser variant 4: attention over the component axis + K-row query
(§12.11 of ``dev/amortised_decision_making.md``).

**xcompAtt** replaces the inner sum-pool of the nested-DeepSets `xcomp`
variant (§12.10) with a single-head cross-attention block over the
J components. The query is derived from the queried component's own
K-diag + cohort sizes; keys / values are per-component summary
embeddings. Motivation: when the intra-block correlation is strong,
sum-pooling gives every component equal weight; attention lets the
network re-weight adaptively conditioned on the K-row profile.

Batch schema (dict):
    ``x_sum``:   ``(B, J)``   -- ``sum_i y_{i, j}`` over observed cohort.
    ``z_sum``:   ``(B, J)``   -- ``sum_i z_{i, j}`` over future cohort.
    ``k_row``:   ``(B, J)``   -- ``K[j*, :]``, one row per batch element.
    ``k_diag``:  ``(B, 1)``   -- ``K[j*, j*]`` scalar.
    ``sizes``:   ``(B, 2)``   -- ``(n / N_max, m / N_max)``.

Sufficient stats for the MVN posterior of ``mu_j`` given ``(x, z)``
are the per-component sums ``(sum_i y_{i, j}, sum_i z_{i, j})`` plus
sizes, so passing pre-pooled component summaries costs no information
vs the full participant matrix (nested-DeepSets in §12.10 sums over
participants first, then over components). This drops the per-forward
compute from ``O(B * N_max * J * embed)`` to ``O(B * J * embed)`` --
roughly ``N_max = 1050`` times cheaper -- and lets us afford the
larger training budget the target scale demands.

Standard Flax + jax.nn.softmax. No exotic layers.
"""

import jax
import jax.numpy as jnp
from flax import linen as nn


class _MLP(nn.Module):
    """SiLU MLP with a linear read-out."""
    dims: tuple

    @nn.compact
    def __call__(self, x):
        for d in self.dims[:-1]:
            x = nn.silu(nn.Dense(d)(x))
        return nn.Dense(self.dims[-1])(x)


class Amortiser_PPS_features_MLP_xcompAtt_qpsi_MLP_loss_multiquantilehead(nn.Module):
    """Attention over the component axis + K-row query + MLP head."""

    q_tok_hidden: tuple = (32, 32)
    embed_dim: int = 32
    hidden_dims: tuple = (64, 64)
    num_quantiles: int = 5
    default_pps_ProbH1_lwr_quantiles_mesh: tuple = (
        0.05, 0.25, 0.5, 0.75, 0.95,
    )

    def setup(self):
        self.q_tok = _MLP(dims=(*self.q_tok_hidden, self.embed_dim))
        self.q_query = _MLP(dims=(*self.q_tok_hidden, self.embed_dim))
        self.q_psi = _MLP(dims=(*self.hidden_dims, self.num_quantiles))

    def __call__(self, batch):
        # Per-component token: three scalars per (b, j).
        tok = jnp.stack(
            [batch['x_sum'], batch['z_sum'], batch['k_row']], axis=-1,
        )                                                       # (B, J, 3)
        h = self.q_tok(tok)                                     # (B, J, E)
        # Query: derive from (sizes, K[j*, j*]).
        q_feat = jnp.concatenate(
            [batch['sizes'], batch['k_diag']], axis=-1,
        )                                                       # (B, 3)
        q = self.q_query(q_feat)[:, None, :]                    # (B, 1, E)
        # Scaled dot-product attention over j.
        scores = jnp.einsum('bqd,bjd->bqj', q, h)               # (B, 1, J)
        scores = scores / jnp.sqrt(float(self.embed_dim))
        weights = jax.nn.softmax(scores, axis=-1)               # (B, 1, J)
        attn_out = jnp.einsum('bqj,bjd->bqd', weights, h)       # (B, 1, E)
        attn_out = attn_out.squeeze(axis=1)                     # (B, E)
        head_in = jnp.concatenate(
            [attn_out, batch['sizes'], batch['k_diag']], axis=-1,
        )                                                       # (B, E + 3)
        return self.q_psi(head_in)
