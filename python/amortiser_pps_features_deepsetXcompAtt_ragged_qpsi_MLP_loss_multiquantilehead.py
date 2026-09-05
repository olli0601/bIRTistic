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
    # §14.4.14 capacity fixes (defaults preserve the plain build):
    #   head_mode: 'plain' | 'factored' | 'semiparam'
    #   precision_pool: append per-item pooled variance to the token
    head_mode: str = 'plain'
    precision_pool: bool = False
    # §14.4.18 sharper-z options:
    #   z_contrast: append (pool_z - pool_x) and (pool_z * pool_x) in EMBED space
    #   raw_pool: append the per-item RAW pooled response means (wb,we) for x and
    #     z plus the plug-in ratio we/wb — the exact statistic the oracle uses.
    #     Learned embed-then-pool does not preserve this ratio; hand it directly.
    z_contrast: bool = False
    raw_pool: bool = False
    # §14.4.21: per-participant rho-mapping channels. Every participant has BOTH
    #   baseline and endline, so the effect proxy can be formed PER PARTICIPANT
    #   (a nonlinear map) BEFORE pooling. Mean-pooling these recovers an effect
    #   statistic the raw [y_b, y_e] mean cannot -- notably the caseness/ratio
    #   shape the plain mean throws away. Fixes the reason-1/2 laundering at the
    #   participant level rather than patching pooled ratios (cf. raw_pool).
    perpart_map: bool = False
    # §14.4.21 diagnostic: force q_tau to a single LINEAR layer (no SiLU), so the
    #   pooled embedding is an affine image of the exact per-item means
    #   (w_b, w_e) -- i.e. q_tau can be the identity. Tests whether the deepset's
    #   weak z is the nonlinear-pool (Jensen) gap or a downstream one.
    linear_tau: bool = False
    # §14.4.21 diagnostic: inject the raw per-item (wb, we, ratio) of the QUERIED
    #   item straight into the head input, bypassing q_tok + cross-attention.
    #   Tests whether the z-signal dies in the token->attention path or at the
    #   head itself (raw_pool routes it through the token; this routes it direct).
    head_raw: bool = False
    head_scale: bool = False   # §14.4.9: inject a data-driven per-item scale (endpoint-change SD) at the head
    default_pps_ProbH1_lwr_quantiles_mesh: tuple = (0.05, 0.25, 0.5, 0.75, 0.95)

    def setup(self):
        self.q_tau = (_MLP(dims=(self.embed_dim,)) if self.linear_tau
                      else _MLP(dims=(*self.q_tau_hidden, self.embed_dim)))
        self.q_tok = _MLP(dims=(*self.q_tok_hidden, self.embed_dim))
        self.q_query = _MLP(dims=(*self.q_query_hidden, self.embed_dim))
        self.q_psi = _MLP(dims=(*self.hidden_dims, self.num_quantiles))
        if self.head_mode == 'factored':          # width-scaling head g(1/sqrt n)
            self.g_head = _MLP(dims=(16, 1))

    def _segment_mean_pool(self, flat, seg, meta, B):
        """flat (T, J, R), seg (T,), meta (B, J, M) -> (B, J, E) [or (B,J,2E)
        with precision_pool] ragged mean (+ variance).

        No padding, no mask: T is the true total participant count. Each
        participant is augmented with its batch element's item metadata,
        embedded by the shared q_tau, then averaged within its segment.
        precision_pool also appends the per-item pooled VARIANCE, giving
        the head the pooled mean AND its uncertainty.
        """
        meta_t = meta[seg]                                   # (T, J, M)
        parts = [flat, meta_t]
        if self.perpart_map:                                 # per-participant effect proxy
            yb = flat[..., 0:1]; ye = flat[..., 1:2]         # baseline, endline in [0,1]
            sgn = 2.0 * meta_t[..., 1:2] - 1.0               # +1 higher-better, -1 lower
            d = sgn * (ye - yb)                              # dir-aware per-part change
            rr = sgn * (ye - yb) / (yb + 0.1)                # dir-aware per-part rel. effect
            parts += [d, rr]                                 # mean pools -> effect statistic
        feats = jnp.concatenate(parts, axis=-1)              # (T, J, R+M[+2])
        e = self.q_tau(feats)                                # (T, J, E)
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
        """Per-item SD of the per-participant endpoint change (endline-baseline): (B, J, 1).
        A data-driven per-item scale signal (proxy for the effect sigma)."""
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
        if self.z_contrast:                                  # explicit z-vs-x contrast
            parts_tok += [pool_z - pool_x, pool_z * pool_x]
        if self.raw_pool:                                    # explicit plug-in ratio channel
            rx = self._raw_mean(batch['x_flat'], batch['x_seg'], B)   # (B,J,2)=(wb,we)
            rz = self._raw_mean(batch['z_flat'], batch['z_seg'], B)
            ratio_x = rx[..., 1:2] / (rx[..., 0:1] + 0.1)    # we/wb, x-cohort
            ratio_z = rz[..., 1:2] / (rz[..., 0:1] + 0.1)
            parts_tok += [rx, rz, ratio_x, ratio_z, ratio_z - ratio_x]
        tok = jnp.concatenate(parts_tok, axis=-1)            # (B, J, 2E) or wider
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
        head_in = jnp.concatenate(parts, axis=-1)            # (B, Q, E+E+A[+7][+1])
        # expose for head-only recalibration (§14.2.5)
        self.sow('intermediates', 'head_in', head_in)
        raw = self.q_psi(head_in)                            # (B, Q, num_quantiles)

        if self.head_mode == 'plain':
            return raw
        if self.head_mode == 'successprob':                  # P(rho > eta0 | x, z); eta0 is the
            return jax.nn.sigmoid(raw)                        # last aux feature (num_quantiles=1)
        # precision feature = 1/sqrt(n) (aux[...,0]); requires aux present
        prec = aux_b[..., 0:1]                               # (B, Q, 1) = n^{-1/2}
        zq = jnp.asarray([-1.6449, -0.6745, 0.0, 0.6745, 1.6449])[:self.num_quantiles]
        zq = zq[None, None, :]
        if self.head_mode == 'factored':
            # median from q_psi; offsets scaled by a learned g(1/sqrt n).
            med = raw[..., self.num_quantiles // 2:self.num_quantiles // 2 + 1]
            g = jax.nn.softplus(self.g_head(aux_b))          # (B, Q, 1) > 0
            return med + (raw - med) * g
        if self.head_mode == 'semiparam':                    # width ~ 1/sqrt(n)
            m = raw[..., 0:1]; s = jax.nn.softplus(raw[..., 1:2])
            return m + zq * s * prec
        if self.head_mode == 'powerlaw':                     # width = C * n^{-p}, learn C,p per item
            m = raw[..., 0:1]; C = jax.nn.softplus(raw[..., 1:2])
            p = jax.nn.sigmoid(raw[..., 2:3])                # p in (0,1)
            self.sow('intermediates', 'p_exp', p)            # for the toward-1/2 regulariser
            return m + zq * (C * prec ** (2.0 * p))          # prec^{2p} = n^{-p}
        if self.head_mode == 'floor':                        # width = sqrt(a^2 + b^2/n), no exponent
            m = raw[..., 0:1]; a = jax.nn.softplus(raw[..., 1:2]); b = jax.nn.softplus(raw[..., 2:3])
            return m + zq * jnp.sqrt(a ** 2 + (b * prec) ** 2)
        return raw
