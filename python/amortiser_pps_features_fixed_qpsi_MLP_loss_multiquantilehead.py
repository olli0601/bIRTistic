"""
Amortiser variant 1: DeepSets-pooled features fixed by the data model,
MLP head ``q_psi``, multi-quantile pinball loss (§7-§11 of
``dev/amortised_decision_making.md``, prototype §11.6).

The pooled sufficient statistics ``Sum T(x_i)``, ``Sum T(z_i^*)`` +
cohort sizes ``(n, m)`` are stored on ``batch['features']`` (typically
already normalised to ``[0, 1]`` by the caller). The head ``q_psi`` is
a plain multi-layer perceptron (MLP) with SiLU non-linearities.

For default configuration ``hidden_dims = (128, 128, 64)`` and
``feature_dim = 2`` (Binomial), the MLP head has
``2*128 + 128 + 128*128 + 128 + 128*64 + 64 + 64*K + K`` free
parameters (``K = num_quantiles``). At ``K = 11`` that is
``25 867`` weights.

**Multi-quantile head + pinball loss.** With target grid
``0 < tau_1 < ... < tau_K < 1`` the network outputs
``q_hat(features) in R^K`` interpreted as the conditional quantiles of
``rho | features``. Training minimises the summed pinball / check loss

    L(q_hat, rho) = 1/K sum_{k=1..K} L_{tau_k}(rho - q_hat_k),
    L_tau(u)      = u * (tau - 1{u < 0}),

averaged over the mini-batch. See
:func:`amortiser_common._pinball_loss` for the JAX implementation.

Shared amortiser utilities (``train``,
``predict_amortised_p_h1_for_one_xz``, ``save_trained_model`` /
``load_fitted_model``) live in :mod:`amortiser_common`. This module
only defines the network class; ``predict_amortised_p_h1_for_many_xz``
below is a thin deployment wrapper for the ``wa`` DataFrame convention
used by the Binomial deployment scripts.
"""

import numpy as np
import pandas as pd
from flax import linen as nn

from amortiser_common import predict_amortised_p_h1_for_one_xz


class Amortiser_PPS_features_fixed_qpsi_MLP_loss_multiquantilehead(nn.Module):
    """DeepSets pooled features + MLP head + multi-quantile pinball loss.

    Expects ``batch = {'features': (B, feature_dim)}``. The network has
    no internal normalisation layer so absolute-scale information is
    preserved; callers are responsible for scaling features to a
    well-behaved range (e.g. dividing counts by ``N_max``).
    """
    hidden_dims: tuple = (128, 128, 64)
    num_quantiles: int = 11
    default_pps_ProbH1_lwr_quantiles_mesh: tuple = (
        0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95,
    )

    @nn.compact
    def __call__(self, batch):
        h = batch['features']
        for w in self.hidden_dims:
            h = nn.Dense(w)(h)
            h = nn.silu(h)
        return nn.Dense(self.num_quantiles)(h)


# ---------------------------------------------------------------------------
# Deployment wrapper (mquantile-schema ``wa`` interface)
# ---------------------------------------------------------------------------


def predict_amortised_p_h1_for_many_xz(
    wa: pd.DataFrame,
    fit: dict,
    *,
    pps_H1_min_effect_size_thresh: float = 0.5,
    pps_ProbH1_target_lwr_quantile: float = 0.89,
    feature_cols=('kn', 'km', 'n', 'm'),
):
    """
    Deployment wrapper matching the signature convention of
    :func:`fit_interim.fit_interim_regress_endptx_on_wz_with_mquantile_regr`.

    Expects ``wa`` to carry per-draw feature columns (``feature_cols``)
    already computed by the calling interim script. Returns
    ``(p_h1_xz_df, perf_df)`` with the same schema as the mquantile
    function.
    """
    p_rows = []
    perf_rows = []
    for item in wa['item_label'].unique():
        g = wa[wa['item_label'] == item].copy()
        features = g[list(feature_cols)].to_numpy(dtype=np.float32)
        p_h1_xz, q_hat, _ = predict_amortised_p_h1_for_one_xz(
            fit, {'features': features},
            pps_H1_min_effect_size_thresh=pps_H1_min_effect_size_thresh,
        )
        gout = g[['item_label', 'item_type', 'item_high_label', 'draw']].copy()
        gout['p_h1_xz'] = p_h1_xz.astype(np.float64)
        gout['q_hat'] = q_hat.astype(np.float64)
        gout['s'] = gout['draw'].astype(int) + 1
        p_rows.append(gout.drop(columns='draw'))

        row = {
            'item_label': item,
            'item_type': g['item_type'].iloc[0],
            'rho': float('nan'),
            'r2': float('nan'),
        }
        if 'w_ratio' in g.columns and 'pps_ratio_x' in g.columns:
            from fit_interim import (
                fit_interim_regress_H1x_on_wz_per_item_summary,
            )
            s = fit_interim_regress_H1x_on_wz_per_item_summary(g)
            row['rho'] = float(s['rho'])
            row['r2'] = float(s['r2'])
        perf_rows.append(row)

    p_h1_xz_df = pd.concat(p_rows, ignore_index=True)
    perf_df = pd.DataFrame(perf_rows)
    return p_h1_xz_df, perf_df
