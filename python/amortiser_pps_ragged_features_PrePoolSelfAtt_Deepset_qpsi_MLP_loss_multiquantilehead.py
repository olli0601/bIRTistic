"""
Pre-pool self-attention variant of the nested-DeepSets amortiser
(§14.6.1 of dev/amortised_decision_making.md).

This is the `deepsetXcompAtt` ragged encoder
(`amortiser_pps_features_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead`)
with ONE change: a per-participant single-head self-attention block
(SAB) over the J items is inserted BEFORE the participant mean-pool, so
the pooled item embedding carries cross-item (within-participant)
covariance — the channel the shuffle test (§14.5.1) proved the plain
deep-set throws away. By default the §14.4.4 POST-pool cross-attention
is dropped (`post_pool_attn=False`): the head reads the queried item's
deep-set embedding e^{j*,s} directly, so this is the "pre-pool-only"
variant of §14.6.2. Set `post_pool_attn=True` for the "both" config
(pre-pool SAB AND post-pool cross-attention).

The post-pool module is kept untouched; this file reuses only its
`_MLP` read-out. Batch schema, ragged pooling and output are identical.

Batch schema (B = num batch elements, T = total x-participants, U =
total z-participants, J items, R response features, M item-metadata
features, A aux scalars):

  x_flat        (T, J, R)     concatenated observed-cohort participants
  x_seg         (T,)          int32 batch-element id in [0, B) per x-participant
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

from amortiser_pps_features_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead import _MLP


class Amortiser_PPS_ragged_features_PrePoolSelfAtt_Deepset_qpsi_MLP_loss_multiquantilehead(nn.Module):
    """Nested DeepSets (ragged participant axis) with per-participant
    self-attention over items BEFORE the pool; post-pool cross-attention
    optional (§14.6.1)."""

    q_tau_hidden: tuple = (32, 32)
    q_tok_hidden: tuple = (32, 32)
    q_query_hidden: tuple = (32, 32)
    embed_dim: int = 32
    hidden_dims: tuple = (64, 64)
    num_quantiles: int = 5
    head_mode: str = 'plain'
    precision_pool: bool = False
    z_contrast: bool = False
    raw_pool: bool = False
    perpart_map: bool = False
    linear_tau: bool = False
    head_raw: bool = False
    head_scale: bool = False
    query_from_token: bool = True            # used only when post_pool_attn (§14.5.2)
    separate_values: bool = True             # used only when post_pool_attn (§14.5.3)
    # §14.6.1: single-head self-attention over the J items of each participant, before the pool.
    prepool_attn: bool = True
    # §14.6.2: pre-pool-only (False, default here) feeds e^{j*,s} straight to the head; True adds the
    #   §14.4.4 post-pool cross-attention on top ("both").
    post_pool_attn: bool = False
    default_pps_ProbH1_lwr_quantiles_mesh: tuple = (0.05, 0.25, 0.5, 0.75, 0.95)

    def setup(self):
        self.q_tau = (_MLP(dims=(self.embed_dim,)) if self.linear_tau
                      else _MLP(dims=(*self.q_tau_hidden, self.embed_dim)))
        if self.prepool_attn:                                 # §14.6.1: single-head SAB projections
            self.sab_q = nn.Dense(self.embed_dim)
            self.sab_k = nn.Dense(self.embed_dim)
            self.sab_v = nn.Dense(self.embed_dim)
            self.sab_o = nn.Dense(self.embed_dim)
        if self.post_pool_attn:                               # §14.4.4 post-pool cross-attention
            self.q_tok = _MLP(dims=(*self.q_tok_hidden, self.embed_dim))
            self.q_query = _MLP(dims=(*self.q_query_hidden, self.embed_dim))
            if self.separate_values:
                self.q_values = _MLP(dims=(*self.q_tok_hidden, self.embed_dim))
        self.q_psi = _MLP(dims=(*self.hidden_dims, self.num_quantiles))
        if self.head_mode == 'factored':
            self.g_head = _MLP(dims=(16, 1))

    def _segment_mean_pool(self, flat, seg, meta, B):
        """flat (T, J, R), seg (T,), meta (B, J, M) -> (B, J, E) [or (B,J,2E)
        with precision_pool]. Embed each (participant, item), optionally
        self-attend over the J items within each participant (§14.6.1),
        then ragged mean (+ variance) over participants."""
        meta_t = meta[seg]                                   # (T, J, M)
        parts = [flat, meta_t]
        if self.perpart_map:                                 # per-participant effect proxy
            yb = flat[..., 0:1]; ye = flat[..., 1:2]
            sgn = 2.0 * meta_t[..., 1:2] - 1.0
            d = sgn * (ye - yb)
            rr = sgn * (ye - yb) / (yb + 0.1)
            parts += [d, rr]
        feats = jnp.concatenate(parts, axis=-1)              # (T, J, R+M[+2])
        e = self.q_tau(feats)                                # (T, J, E)
        if self.prepool_attn:                                # §14.6.1: per-participant self-attn over J items
            qh = self.sab_q(e); kh = self.sab_k(e); vh = self.sab_v(e)    # (T,J,E)
            asc = jnp.einsum('tje,tke->tjk', qh, kh) / jnp.sqrt(float(self.embed_dim))
            aw = jax.nn.softmax(asc, axis=-1)                # (T,J,J) cross-item weights within participant
            e = e + self.sab_o(jnp.einsum('tjk,tke->tje', aw, vh))       # residual, single head
        cnt = jax.ops.segment_sum(
            jnp.ones((flat.shape[0], 1, 1), e.dtype), seg, num_segments=B)
        cnt = jnp.maximum(cnt, 1.0)                          # (B, 1, 1)
        mean = jax.ops.segment_sum(e, seg, num_segments=B) / cnt
        if self.precision_pool:
            sq = jax.ops.segment_sum(e ** 2, seg, num_segments=B) / cnt
            var = jnp.maximum(sq - mean ** 2, 0.0)           # (B, J, E)
            return jnp.concatenate([mean, var], axis=-1)     # (B, J, 2E)
        return mean

    def _raw_mean(self, flat, seg, B):
        """Segment mean of the RAW responses (no embedding): (B, J, R)."""
        cnt = jax.ops.segment_sum(
            jnp.ones((flat.shape[0], 1, 1), flat.dtype), seg, num_segments=B)
        cnt = jnp.maximum(cnt, 1.0)
        return jax.ops.segment_sum(flat, seg, num_segments=B) / cnt

    def _change_std(self, flat, seg, B):
        """Per-item SD of the per-participant endpoint change (endline-baseline): (B, J, 1)."""
        ch = (flat[..., 1] - flat[..., 0])[..., None]               # (T, J, 1)
        cnt = jnp.maximum(jax.ops.segment_sum(jnp.ones_like(ch), seg, num_segments=B), 1.0)
        m = jax.ops.segment_sum(ch, seg, num_segments=B) / cnt
        m2 = jax.ops.segment_sum(ch ** 2, seg, num_segments=B) / cnt
        return jnp.sqrt(jnp.maximum(m2 - m ** 2, 0.0))              # (B, J, 1)

    def __call__(self, batch):
        meta = batch['item_metadata']                        # (B, J, M)
        B, J = meta.shape[0], meta.shape[1]
        pool_x = self._segment_mean_pool(batch['x_flat'], batch['x_seg'], meta, B)
        pool_z = self._segment_mean_pool(batch['z_flat'], batch['z_seg'], meta, B)
        parts_tok = [pool_x, pool_z]
        if self.z_contrast:
            parts_tok += [pool_z - pool_x, pool_z * pool_x]
        if self.raw_pool:
            rx = self._raw_mean(batch['x_flat'], batch['x_seg'], B)   # (B,J,2)=(wb,we)
            rz = self._raw_mean(batch['z_flat'], batch['z_seg'], B)
            ratio_x = rx[..., 1:2] / (rx[..., 0:1] + 0.1)
            ratio_z = rz[..., 1:2] / (rz[..., 0:1] + 0.1)
            parts_tok += [rx, rz, ratio_x, ratio_z, ratio_z - ratio_x]
        tok = jnp.concatenate(parts_tok, axis=-1)            # (B, J, 2E) or wider: item embeddings e^{j,s}

        qidx = batch['query_idx'].astype(jnp.int32)          # (B, Q)
        Q = qidx.shape[1]
        Et = tok.shape[-1]
        # queried item's deep-set embedding e^{j*,s}: (B, Q, 2E+) — head skip, and (default) query input
        tok_query = jnp.take_along_axis(tok, jnp.broadcast_to(qidx[:, :, None], (B, Q, Et)), axis=1)

        if self.post_pool_attn:                              # §14.4.4 cross-attention over the J items ("both")
            h = self.q_tok(tok)                              # (B, J, E) keys
            E = h.shape[-1]
            h_query = jnp.take_along_axis(h, jnp.broadcast_to(qidx[:, :, None], (B, Q, E)), axis=1)  # (B,Q,E)
            q_src = tok_query if self.query_from_token else h_query   # §14.5.2
            q = self.q_query(q_src)                          # (B, Q, E)
            v = self.q_values(tok) if self.separate_values else h     # §14.5.3
            scores = jnp.einsum('bqe,bje->bqj', q, h) / jnp.sqrt(float(self.embed_dim))
            w = jax.nn.softmax(scores, axis=-1)              # (B, Q, J)
            attn_out = jnp.einsum('bqj,bje->bqe', w, v)      # (B, Q, E)
            parts = [attn_out, h_query]
        else:                                                # §14.6.2 pre-pool-only: head reads e^{j*,s}
            parts = [tok_query]

        aux = batch.get('aux', None)
        aux_b = (jnp.broadcast_to(aux[:, None, :], (B, Q, aux.shape[-1]))
                 if aux is not None else None)
        if aux_b is not None:
            parts.append(aux_b)
        if self.head_raw:                                    # inject raw (wb,we,ratio) AT THE HEAD
            rx = self._raw_mean(batch['x_flat'], batch['x_seg'], B)   # (B,J,2)
            rz = self._raw_mean(batch['z_flat'], batch['z_seg'], B)
            ratx = rx[..., 1:2] / (rx[..., 0:1] + 0.1)
            ratz = rz[..., 1:2] / (rz[..., 0:1] + 0.1)
            rawf = jnp.concatenate([rx, rz, ratx, ratz, ratz - ratx], axis=-1)  # (B,J,7)
            gidx = jnp.broadcast_to(qidx[:, :, None], (B, Q, rawf.shape[-1]))
            parts.append(jnp.take_along_axis(rawf, gidx, axis=1))              # (B,Q,7)
        if self.head_scale:                                  # data-driven per-item scale at the head
            sc = self._change_std(batch['x_flat'], batch['x_seg'], B)          # (B,J,1)
            parts.append(jnp.take_along_axis(sc, jnp.broadcast_to(qidx[:, :, None], (B, Q, 1)), axis=1))
        head_in = jnp.concatenate(parts, axis=-1)            # (B, Q, 2E+A[...]) or (E+E+A[...]) with post-pool
        self.sow('intermediates', 'head_in', head_in)        # head-only recalibration (§14.2.5)
        raw = self.q_psi(head_in)                            # (B, Q, num_quantiles)

        if self.head_mode == 'plain':
            return raw
        if self.head_mode == 'successprob':                  # P(rho > eta0 | x, z)
            return jax.nn.sigmoid(raw)
        prec = aux_b[..., 0:1]                               # (B, Q, 1) = n^{-1/2}
        zq = jnp.asarray([-1.6449, -0.6745, 0.0, 0.6745, 1.6449])[:self.num_quantiles]
        zq = zq[None, None, :]
        if self.head_mode == 'factored':
            med = raw[..., self.num_quantiles // 2:self.num_quantiles // 2 + 1]
            g = jax.nn.softplus(self.g_head(aux_b))          # (B, Q, 1) > 0
            return med + (raw - med) * g
        if self.head_mode == 'semiparam':
            m = raw[..., 0:1]; s = jax.nn.softplus(raw[..., 1:2])
            return m + zq * s * prec
        if self.head_mode == 'powerlaw':
            m = raw[..., 0:1]; C = jax.nn.softplus(raw[..., 1:2])
            p = jax.nn.sigmoid(raw[..., 2:3])
            self.sow('intermediates', 'p_exp', p)
            return m + zq * (C * prec ** (2.0 * p))
        if self.head_mode == 'floor':
            m = raw[..., 0:1]; a = jax.nn.softplus(raw[..., 1:2]); b = jax.nn.softplus(raw[..., 2:3])
            return m + zq * jnp.sqrt(a ** 2 + (b * prec) ** 2)
        return raw
