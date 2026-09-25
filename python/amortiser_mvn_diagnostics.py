#!/usr/bin/env python3
"""
Unified MVN diagnostics harness (§13.5 results). Loads each trained amortiser
architecture, re-deploys it capturing the FULL per-draw quantile predictions (the
deploy scripts discard all but the median), and against the closed-form / HMC
reference draws (`mu_draws` in the interim-data pkls) computes, per (architecture,
J, interim): conditional calibration (PIT-KS), marginal calibration (marg-KS),
posterior-contraction exponent p (amortiser vs reference), and PPS-MSE vs the
closed-form Phi-tail PPS. Diagnostics are RAW (prior-trained net, no head fine-tune)
so they reflect the intrinsic encoder quality; the deepset family's CALIBRATED
diagnostics are in §13.8.

All six architectures share the batch-construction extracted from their deploy
scripts and the shared `amortiser_calibration` metrics (same code as Ukraine).

Output: py-mvn-interim-diagnostics-260916/mvn_arch_diagnostics.csv (+ .pdf).
Usage: pixi run python scripts-py/MVN_interim_diagnostics_by_architecture.py
Env: DIAG_S (200), DIAG_ARCHS (comma list of suffixes), DIAG_SMOKE (1).
"""
# ---- boilerplate ----

import os, sys, time, warnings
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'python'))
import numpy as np, pandas as pd
warnings.filterwarnings('ignore')
from amortiser_common import load_fitted_model, predict_amortised_p_h1_for_one_xz
import amortiser_calibration as cal

_sb = "/Users/or105/sandbox/bIRTistic"
DIR_SIM = f"{_sb}/py-mvn-interim-simulations-260609"
OUT = f"{_sb}/py-mvn-interim-diagnostics-260916"; os.makedirs(OUT, exist_ok=True)
SMOKE = os.environ.get('DIAG_SMOKE', '0') == '1'
# Native deploy budget per architecture kind (matches the deploy scripts + §13.3):
# idcomp S=4000, xcomp/hand-token S=500, deepset S=200. xcomp kept at 200 (it is
# the undertrained set-aside baseline of §13.4 and the compute bottleneck).
S_BY_KIND = {'idcomp': 4000, 'xcomp': 200, 'tokens': 500, 'deepset': 200}
_S_OVERRIDE = os.environ.get('DIAG_S')          # if set, forces a common S
ETA0, MU0 = 0.0, 1.0

sim = pd.read_pickle(f"{DIR_SIM}/mvn_sim_data.pkl"); sp = sim['simu_params']
N_MAX = int(sp['N_full']); J_GRID = list(sp['J_grid'])
ETAH = float(sp.get('pps_ProbH1_target_lwr_quantile', sp.get('pps_ProbH1_thresh', 0.89)))
pps_cf = pd.read_pickle(f"{DIR_SIM}/mvn_pps_closed_form.pkl")['pps_cf']

# arch registry: (label, suffix, dir, builder, J-list)
D = lambda s: f"{_sb}/py-mvn-interim-amortise-endptx-on-wz-with-features-{s}"
ARCHS = [
    ('idcomp (fixed)',   'RGEA', D("fixed-idcomp-qpsi-MLP-loss-multiquantilehead-260714"),   'idcomp',  [20, 60, 100]),
    ('xcomp (MLP)',      'RGEC', D("MLP-xcomp-qpsi-MLP-loss-multiquantilehead_15k_260716"),   'xcomp',   [20, 60, 100]),
    ('itemScompAtt',     'RGED', D("itemScompAtt-qpsi-MLP-loss-multiquantilehead_15k_260716"),'tokens',  [20, 60, 100]),
    ('itemXcompAtt',     'RGEF', D("itemXcompAtt-qpsi-MLP-loss-multiquantilehead_15k_260716"),'tokens',  [20, 60, 100]),
    ('deepsetScompAtt',  'RGDS', D("deepsetScompAtt-qpsi-MLP-loss-multiquantilehead_260803"), 'deepset', [20]),
    ('deepsetXcompAtt',  'RGDX', D("deepsetXcompAtt-qpsi-MLP-loss-multiquantilehead_260803"), 'deepset', [20]),
]
_only = os.environ.get('DIAG_ARCHS')
if _only:
    keep = set(_only.split(',')); ARCHS = [a for a in ARCHS if a[1] in keep]
if SMOKE:
    ARCHS = [a for a in ARCHS if a[1] in ('RGEA', 'RGEF', 'RGDX')]
    J_GRID = [20]


def _xz_arrays(dpi, zi, J, S, m):
    """x_wide (n,J); ypred_arr (S,m,J) aligned to zi.ypred draws."""
    x_wide = (dpi.pivot_table(index='pid', columns='j', values='y').sort_index()
              .reindex(columns=range(J)).to_numpy(np.float64))
    cols = sorted([c for c in zi.columns if c.startswith('ypred_')],
                  key=lambda c: int(c.split('_')[1]))[:S]
    zis = zi.sort_values(['pid', 'j']).reset_index(drop=True)
    ypred = zis[cols].to_numpy().T.reshape(S, m, J)
    return x_wide, ypred


def deploy(fit, kind, dpi, zi, K, Kd, n, m, J, S, taus):
    """Return qs (S,J,nq) capturing the full quantile matrix per draw."""
    nq = len(taus); x_wide, ypred = _xz_arrays(dpi, zi, J, S, m)
    qs = np.empty((S, J, nq), np.float64)
    if kind == 'idcomp':
        sum_x = dpi.groupby('j')['y'].sum().reindex(range(J), fill_value=0.0).to_numpy()
        sum_z = ypred.sum(1)                                              # (S,J)
        ntot = np.full(S, (n + m) / N_MAX, np.float32)
        for j in range(J):
            feats = np.stack([((sum_x[j] + sum_z[:, j]) / N_MAX).astype(np.float32), ntot], -1)
            _p, _q, preds = predict_amortised_p_h1_for_one_xz(fit, {'features': feats}, ETA0)
            qs[:, j, :] = preds
        return qs
    # per-s builders (batched over J queries)
    sizes = np.array([n / N_MAX, m / N_MAX], np.float32)
    sizes_b = np.broadcast_to(sizes[None, :], (J, 2)).astype(np.float32)
    Kd_b = Kd[:, None].astype(np.float32); qidx = np.arange(J, dtype=np.int32)
    for s in range(S):
        if kind == 'xcomp':
            x_pad = np.zeros((N_MAX, J), np.float32); x_pad[:n] = x_wide
            z_pad = np.zeros((N_MAX, J), np.float32); z_pad[:m] = ypred[s]
            mx = np.zeros(N_MAX, np.float32); mx[:n] = 1.0
            mz = np.zeros(N_MAX, np.float32); mz[:m] = 1.0
            batch = dict(x=np.broadcast_to(x_pad[None], (J, N_MAX, J)).astype(np.float32),
                         mask_x=np.broadcast_to(mx[None], (J, N_MAX)).astype(np.float32),
                         z=np.broadcast_to(z_pad[None], (J, N_MAX, J)).astype(np.float32),
                         mask_z=np.broadcast_to(mz[None], (J, N_MAX)).astype(np.float32),
                         sizes=sizes_b, k_row=K.astype(np.float32), k_diag=Kd_b)
        elif kind == 'tokens':
            xs = np.broadcast_to((x_wide.sum(0) / N_MAX).astype(np.float32)[None, :], (J, J))
            zs = np.broadcast_to((ypred[s].sum(0) / N_MAX).astype(np.float32)[None, :], (J, J))
            tok = np.stack([xs, zs, K.astype(np.float32)], -1).astype(np.float32)
            batch = dict(tokens=tok, mask=np.ones((J, J), np.float32), query_idx=qidx,
                         aux=np.concatenate([sizes_b, Kd_b], -1).astype(np.float32))
        else:   # deepset
            x_pad = np.zeros((N_MAX, J), np.float32); x_pad[:n] = x_wide
            z_pad = np.zeros((N_MAX, J), np.float32); z_pad[:m] = ypred[s]
            mx = np.zeros(N_MAX, np.float32); mx[:n] = 1.0
            mz = np.zeros(N_MAX, np.float32); mz[:m] = 1.0
            batch = dict(x_responses=np.broadcast_to(x_pad[None, :, :, None], (J, N_MAX, J, 1)).astype(np.float32),
                         mask_x=np.broadcast_to(mx[None], (J, N_MAX)).astype(np.float32),
                         z_responses=np.broadcast_to(z_pad[None, :, :, None], (J, N_MAX, J, 1)).astype(np.float32),
                         mask_z=np.broadcast_to(mz[None], (J, N_MAX)).astype(np.float32),
                         item_metadata=K.astype(np.float32)[..., None], query_idx=qidx,
                         aux=np.concatenate([sizes_b, Kd_b], -1).astype(np.float32))
        _p, _q, preds = predict_amortised_p_h1_for_one_xz(fit, batch, ETA0)
        qs[s] = preds
    return qs


def amo_contraction_p(QS, NOBS, interims, J):
    """Median amortiser posterior-contraction exponent (log-log slope of the
    marginal-predictive SD vs n)."""
    ps = []
    for j in range(J):
        ns, sds = [], []
        for k in interims:
            q = QS[k][:, j, :]; q = q[np.isfinite(q).all(1)]
            mu = q[:, q.shape[1] // 2]; sdw = (q[:, -1] - q[:, 0]) / 3.2897
            sd = np.sqrt(np.mean(sdw ** 2) + np.var(mu))
            if sd > 1e-9:
                ns.append(NOBS[k]); sds.append(sd)
        if len(ns) >= 3:
            ps.append(-np.polyfit(np.log(ns), np.log(sds), 1)[0])
    return float(np.nanmedian(ps)) if ps else np.nan


rows = []
for label, suf, ddir, kind, jlist in ARCHS:
    ckpt = f"{ddir}/mvn_interim_amortised_pps_net.pkl"
    if not os.path.exists(ckpt):
        print(f"[skip] {label}: no net at {ckpt}"); continue
    fit = load_fitted_model(ckpt); taus = np.asarray(fit['pps_ProbH1_lwr_quantiles_mesh'])
    for J in [j for j in jlist if j in J_GRID]:
        t0 = time.time()
        ipkl = pd.read_pickle(f"{DIR_SIM}/mvn_J{J}_interim_data.pkl")
        cell = sim['cells'][J]; K = (cell['R_chol'] @ cell['R_chol'].T); Kd = np.diag(K)
        labels = cell['dit']['item_label'].to_numpy()
        interims = sorted(ipkl.keys())
        if SMOKE:
            interims = interims[:3]
        _Sk = 20 if SMOKE else (int(_S_OVERRIDE) if _S_OVERRIDE else S_BY_KIND[kind])
        QS, TGT, NOBS = {}, {}, {}
        for k in interims:
            blk = ipkl[k]; n, m = int(blk['n_obs']), int(blk['interim_m'])
            ncols = sum(c.startswith('ypred_') for c in blk['zi'].columns)
            S = min(_Sk, blk['mu_draws'].shape[0], ncols)
            QS[k] = deploy(fit, kind, blk['dpi'], blk['zi'], K, Kd, n, m, J, S, taus)
            TGT[k] = (np.asarray(blk['mu_draws'][:S]) - MU0)
            NOBS[k] = n
        summ = cal.calibration_summary(QS, TGT, NOBS, interims, labels, taus=taus)
        _, ref_p = cal.contraction_slope(TGT, NOBS, interims, labels)
        amo_p = amo_contraction_p(QS, NOBS, interims, J)
        # PPS-MSE vs closed form (recomputed from captured quantiles, this S)
        cfJ = pps_cf[pps_cf.J == J].set_index(['interim_id', 'j'])['pps']
        se = []
        for k in interims:
            for j in range(J):
                q = QS[k][:, j, :]
                ph1 = 1.0 - np.array([np.interp(ETA0, q[s], taus, 0., 1.) for s in range(q.shape[0])])
                if (k, j) in cfJ.index:
                    se.append((float(np.mean(ph1 > ETAH)) - float(cfJ.loc[(k, j)])) ** 2)
        rows.append(dict(arch=label, suf=suf, J=J, S=S,
                         pps_mse=float(np.nanmean(se)),
                         pit_ks=float(summ.pit_ks.mean()), marg_ks=float(summ.marg_ks.mean()),
                         amo_contraction_p=amo_p, ref_contraction_p=ref_p,
                         mins=round((time.time() - t0) / 60, 2)))
        print(f"  {label:16s} J={J:3d}: PPS-MSE={rows[-1]['pps_mse']:.5f} PIT-KS={rows[-1]['pit_ks']:.3f}"
              f" marg-KS={rows[-1]['marg_ks']:.3f} amo-p={amo_p:.2f} ref-p={ref_p:.2f} ({rows[-1]['mins']}m)")

res = pd.DataFrame(rows)
res.to_csv(f"{OUT}/mvn_arch_diagnostics.csv", index=False)
print("\n===== architecture diagnostics (RAW, prior-trained) =====")
print(res.to_string(index=False))
print(f"\nSaved -> {OUT}/mvn_arch_diagnostics.csv")
