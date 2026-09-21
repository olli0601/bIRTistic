"""
Shared deployment-calibration + Bernstein-von Mises correction for the amortised
PPS estimators (§14.4.6, §14.4.7). Used identically by the Ukraine PCM deploy
driver and the MVN drivers, so both case studies share one implementation of the
calibration layer (the reference posterior differs -- SVI for Ukraine, closed-form
for MVN -- but the operations on the deployed quantiles are the same).

The calibration operates on three per-interim arrays, keyed by interim id k:
  HD[k]  : (S, J, H)  the frozen-encoder head input (sown intermediate) per draw
  TGT[k] : (S, J)     reference draws of the endpoint rho (joint-draw target)
  NOBS[k]: int        observed cohort size n
and a per-item standardiser ``sig`` (J,) (1 for MVN; per-item SD for Ukraine).

Public functions
----------------
head_transform(raw, prec, head_mode)                 head-mode quantile map
expanding_head_ft(...)   -> HK                        §14.4.6 head fine-tune
affine_median_shift(HK, TGT, interims) -> HK          §14.4.6 location shift
bvm_correction(HK, TGT, NOBS, interims, N_ref) -> HKV §14.4.7 self-consistency
calibration_summary(PPSSRC, TGT, NOBS, interims, labels, ...) -> DataFrame
"""

import numpy as np
import jax
import jax.numpy as jnp
import optax
from flax import linen as nn

TAUS_DEFAULT = (0.05, 0.25, 0.5, 0.75, 0.95)
ZQ_DEFAULT = np.array([-1.6449, -0.6745, 0.0, 0.6745, 1.6449])


class _MLP(nn.Module):
    """SiLU MLP with linear read-out (same shape as the amortiser heads)."""
    dims: tuple

    @nn.compact
    def __call__(self, x):
        for d in self.dims[:-1]:
            x = nn.silu(nn.Dense(d)(x))
        return nn.Dense(self.dims[-1])(x)


# ---------------------------------------------------------------------------
# Head-mode quantile map (deploy-side mirror of the net's head_mode branches).
# ---------------------------------------------------------------------------


def head_transform(raw, prec, head_mode, zq=None, hd=None, g_head=None,
                   g_params=None, a_aux=None):
    """Standardised quantiles from a raw q_psi output under ``head_mode``. Single
    deploy-side head map shared by MVN and Ukraine.

    raw:  (..., K) q_psi output; prec: scalar or (...,1) = n^{-1/2}.
    Modes (§14.4.7): 'plain' monotone quantiles; 'powerlaw' m + z*C*n^{-p};
    'floor' m + z*sqrt(a^2 + b^2/n); 'semiparam' m + z*softplus(s)*n^{-1/2};
    'factored' med + (raw-med)*g where g = softplus(g_head(head_in aux tail)) is a
    frozen width-scaler (needs ``hd`` head input + ``g_head``/``g_params``/``a_aux``).
    """
    zq = jnp.asarray(ZQ_DEFAULT if zq is None else zq)[None, :]
    if head_mode == 'powerlaw':
        m = raw[..., 0:1]; C = jax.nn.softplus(raw[..., 1:2]); p = jax.nn.sigmoid(raw[..., 2:3])
        return m + zq * (C * prec ** (2.0 * p))
    if head_mode == 'floor':
        m = raw[..., 0:1]; a = jax.nn.softplus(raw[..., 1:2]); b = jax.nn.softplus(raw[..., 2:3])
        return m + zq * jnp.sqrt(a ** 2 + (b * prec) ** 2)
    if head_mode == 'freeq':
        # §16.3 skew-capable head: median + independent lower/upper monotone gaps (softplus),
        # scaled by the BvM rate n^{-1/2}; BvM then rescales the width to the power law. The
        # lower/upper gaps are free, so the quantile set can be asymmetric (any skew).
        m = raw[..., 0:1]
        d = jax.nn.softplus(raw[..., 1:5])                       # (...,4) positive gaps
        zeros = jnp.zeros_like(m)
        off = jnp.concatenate([-(d[..., 0:1] + d[..., 1:2]), -d[..., 1:2],
                               zeros, d[..., 2:3], d[..., 2:3] + d[..., 3:4]], axis=-1)
        return m + off * prec
    if head_mode == 'semiparam':
        return raw[..., 0:1] + zq * jax.nn.softplus(raw[..., 1:2]) * prec
    if head_mode == 'factored':
        K = raw.shape[-1]; med = raw[..., K // 2:K // 2 + 1]
        g = jax.nn.softplus(g_head.apply(g_params, hd[..., -a_aux:]))
        return med + (raw - med) * g
    return jnp.maximum.accumulate(raw, -1)


# ---------------------------------------------------------------------------
# §14.4.6 Expanding-window head-only fine-tune of q_psi (encoder frozen).
# ---------------------------------------------------------------------------


def expanding_head_ft(HD, TGT, NOBS, interims, q_psi_params, hidden_dims,
                      head_mode, sig, *, taus=TAUS_DEFAULT, steps=800, lr=3e-4,
                      seed=1, expanding=True, verbose=False):
    """Refit the quantile head on the reference draws of interims 1..k, reusing
    the cached head input HD (encoder never re-run). Returns HK: {k: (S,J,5)}
    real-scale calibrated quantiles. ``sig`` (J,) rescales to real units."""
    J = TGT[interims[0]].shape[1]
    hidden = tuple(hidden_dims)
    qpsi = _MLP(dims=(*hidden, len(taus)))
    qpsi0 = {'params': q_psi_params}
    taus_j = jnp.asarray(taus); opt = optax.adam(lr)
    sig = np.asarray(sig)

    def pinball(params, hd, y, prec):
        pr = head_transform(qpsi.apply(params, hd), prec, head_mode)
        e = y[:, None] - pr
        return jnp.mean(jnp.maximum(taus_j[None] * e, (taus_j[None] - 1) * e))

    @jax.jit
    def step(params, ost, h, y, pr):
        l, g = jax.value_and_grad(lambda p: pinball(p, h, y, pr))(params)
        up, ost2 = opt.update(g, ost, params)
        return optax.apply_updates(params, up), ost2, l

    def _pool(ks):
        HDp = np.concatenate([HD[i].reshape(-1, HD[i].shape[-1]) for i in ks])
        TGp = np.concatenate([(TGT[i] / sig[None]).reshape(-1) for i in ks])
        PRp = np.concatenate([np.full(HD[i].shape[0] * J, 1.0 / np.sqrt(NOBS[i]), np.float32)
                              for i in ks])
        ok = np.isfinite(TGp)
        return (HDp[ok].astype(np.float32), TGp[ok].astype(np.float32),
                PRp[ok, None].astype(np.float32))

    def _fit(ks):
        HDp, TGp, PRp = _pool(ks)
        params = qpsi0; ost = opt.init(params); r = np.random.default_rng(seed); l = 0.0
        for _ in range(steps):
            ix = r.integers(0, HDp.shape[0], min(2048, HDp.shape[0]))
            params, ost, l = step(params, ost, jnp.asarray(HDp[ix]), jnp.asarray(TGp[ix]),
                                   jnp.asarray(PRp[ix]))
        return params, float(l)

    def _predict(params, k):
        flat = jnp.asarray(HD[k].reshape(-1, HD[k].shape[-1]))
        pr = np.asarray(head_transform(qpsi.apply(params, flat), 1.0 / np.sqrt(NOBS[k]), head_mode))
        return pr.reshape(HD[k].shape[0], J, len(taus)) * sig[None, :, None]

    HK = {}
    if expanding:
        for k in interims:
            pk, l = _fit([i for i in interims if i <= k])
            HK[k] = _predict(pk, k)
            if verbose:
                print(f"    head-ft expand 1..{k} (loss {l:.4f})")
    else:
        pk, l = _fit([interims[0]])
        for k in interims:
            HK[k] = _predict(pk, k)
        if verbose:
            print(f"    head-ft i1-only (loss {l:.4f})")
    return HK


# ---------------------------------------------------------------------------
# §14.4.6 Per-item expanding affine median-shift (pure location move).
# ---------------------------------------------------------------------------


def _marg_median(qs, taus):
    lo, hi = float(qs.min()), float(qs.max())
    if hi <= lo:
        return float(np.median(qs))
    gr = np.linspace(lo, hi, 200); F = np.zeros_like(gr)
    for s in range(qs.shape[0]):
        F += np.interp(gr, qs[s], taus, 0., 1.)
    return float(np.interp(0.5, F / qs.shape[0], gr))


def affine_median_shift(HK, TGT, interims, *, taus=TAUS_DEFAULT):
    """Delta_j(k) = mean_{i<=k}[median(reference) - amortiser marginal median],
    applied as a location shift to all quantiles. Returns a new HK."""
    J = TGT[interims[0]].shape[1]
    resid = {}
    for k in interims:
        rj = np.full(J, np.nan)
        for j in range(J):
            y = TGT[k][:, j]; ok = np.isfinite(y); yv = y[ok]
            if yv.size >= 10:
                rj[j] = np.median(yv) - _marg_median(HK[k][ok, j, :], taus)
        resid[k] = rj
    HK0 = {k: HK[k].copy() for k in interims}
    out = {}
    for k in interims:
        dj = np.nan_to_num(np.nanmean(np.stack([resid[i] for i in interims if i <= k]), 0))
        out[k] = HK0[k] + dj[None, :, None]
    return out


# ---------------------------------------------------------------------------
# §14.4.7 Bernstein-von Mises between-first self-consistency correction.
# ---------------------------------------------------------------------------


def bvm_correction(HK, TGT, NOBS, interims, N_ref, *, verbose=False):
    """Per item fit the marginal power law SD=C n^-p to the reference SD, then
    reshape each interim's per-draw quantiles so the mixture variance equals the
    law read at n (T) with the conditional pinned at n+m=N_ref. Returns HKV."""
    J = TGT[interims[0]].shape[1]
    Cj = np.full(J, np.nan); pj = np.full(J, np.nan)
    for j in range(J):
        ns, sds = [], []
        for k in interims:
            yv = TGT[k][:, j]; yv = yv[np.isfinite(yv)]
            if yv.size >= 10 and np.std(yv) > 1e-9:
                ns.append(NOBS[k]); sds.append(np.std(yv))
        if len(ns) >= 3:
            b = np.polyfit(np.log(ns), np.log(sds), 1); pj[j] = -b[0]; Cj[j] = float(np.exp(b[1]))
    if verbose:
        print(f"    BvM law fit: median p={np.nanmedian(pj):.3f}")
    HKV = {}
    for k in interims:
        n = NOBS[k]; Hk = HK[k].copy()
        for j in range(J):
            if not np.isfinite(Cj[j]):
                continue
            C, p = Cj[j], pj[j]
            T = (C * n ** (-p)) ** 2
            Wt = (C * N_ref ** (-p)) ** 2
            Bs = max(T - Wt, 0.0)
            qs = HK[k][:, j, :]; mu = qs[:, 2]; mub = float(mu.mean())
            sg = (qs[:, 4] - qs[:, 0]) / 3.2897
            W = float(np.mean(sg ** 2)); B = float(np.var(mu))
            a = np.sqrt(Wt / W) if W > 1e-12 else 1.0
            be = np.sqrt(Bs / B) if B > 1e-12 else 0.0
            Hk[:, j, :] = (mub + be * (mu - mub))[:, None] + a * (qs - mu[:, None])
        HKV[k] = Hk
    return HKV, dict(Cj=Cj, pj=pj)


# ---------------------------------------------------------------------------
# Scalar diagnostics: PIT-KS (conditional) + marg-KS (marginal vs reference).
# ---------------------------------------------------------------------------


def _ks_uniform(u):
    u = np.sort(u[np.isfinite(u)])
    if u.size < 5:
        return np.nan
    return float(np.max(np.abs(np.arange(1, u.size + 1) / u.size - u)))


def _marg_cdf(qs, grid, taus):
    F = np.zeros_like(grid)
    for s in range(qs.shape[0]):
        F += np.interp(grid, qs[s], taus, 0., 1.)
    return F / qs.shape[0]


def _marg_ks_item(qs, yv, taus):
    lo = float(min(qs.min(), yv.min())); hi = float(max(qs.max(), yv.max()))
    if hi <= lo:
        return np.nan
    gr = np.linspace(lo, hi, 400)
    Fa = _marg_cdf(qs, gr, np.asarray(taus))
    ys = np.sort(yv); Fr = np.searchsorted(ys, gr, side='right') / ys.size
    return float(np.max(np.abs(Fa - Fr)))


def calibration_summary(PPSSRC, TGT, NOBS, interims, labels, *, taus=TAUS_DEFAULT,
                        extra=None):
    """Per (interim, item) PIT-KS + marg-KS from calibrated quantiles PPSSRC and
    reference draws TGT. Returns a long DataFrame. ``extra`` merges scalar columns."""
    import pandas as pd
    J = len(labels); rows = []
    for k in interims:
        for j in range(J):
            y = TGT[k][:, j]; ok = np.isfinite(y); yv = y[ok]
            if yv.size < 10:
                continue
            qs = PPSSRC[k][ok, j, :]
            u = np.array([np.interp(yv[s], qs[s], taus, 0., 1.) for s in range(len(yv))])
            qsa = PPSSRC[k][:, j, :]; qsa = qsa[np.isfinite(qsa).all(1)]
            r = dict(interim_id=k, n=NOBS[k], item_label=labels[j],
                     pit_ks=_ks_uniform(u), marg_ks=_marg_ks_item(qsa, yv, taus))
            if extra:
                r.update(extra)
            rows.append(r)
    return pd.DataFrame(rows)


def contraction_slope(TGT, NOBS, interims, labels):
    """Per-item power-law exponent p (log-log slope of reference SD vs n) +
    median. Returns (per_item_dict, median_p)."""
    J = len(labels); out = {}
    for j in range(J):
        ns, sds = [], []
        for k in interims:
            yv = TGT[k][:, j]; yv = yv[np.isfinite(yv)]
            if yv.size >= 10 and np.std(yv) > 1e-9:
                ns.append(NOBS[k]); sds.append(np.std(yv))
        if len(ns) >= 3:
            out[labels[j]] = float(-np.polyfit(np.log(ns), np.log(sds), 1)[0])
    med = float(np.nanmedian(list(out.values()))) if out else float('nan')
    return out, med
