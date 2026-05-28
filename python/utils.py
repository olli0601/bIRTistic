"""Utility functions for Bayesian IRT analysis."""

from __future__ import annotations

from typing import Optional

#%%
def _map_cq_id_to_item_structure(dp1: pd.DataFrame, dit: pd.DataFrame) -> pd.DataFrame:
    """
    Map probability question IDs (cq_id) back to item structure.
    
    Takes item metadata and generates a lookup table that maps each cq_id 
    to its corresponding item_type_id, item_time_id, and category (y).
    
    Parameters
    ----------
    dp1 : pd.DataFrame
        Pre-processed data containing item_type_id, item_label, item_time_id columns.
    dit : pd.DataFrame
        Item metadata containing item_type_id, item_label, cat_length columns.
    
    Returns
    -------
    pd.DataFrame
        Lookup table with columns: item_type_id, item_time_id, y, cq_id.
        One row per (item_type_id, item_time_id, category) combination.
    """
    import numpy as np

    tmp = dp1[['item_type_id', 'item_label', 'item_time_id']].drop_duplicates()
    tmp = tmp.merge(dit[['item_type_id', 'item_label', 'cat_length']], 
                   on=['item_type_id', 'item_label'])
    tmp = tmp.sort_values(['item_type_id', 'item_time_id']).reset_index(drop=True)
    tmp['cq_id'] = tmp['cat_length'].cumsum() - tmp['cat_length']
    tmp = tmp.assign(y=tmp['cat_length'].map(lambda n: np.arange(n))).explode('y', ignore_index=True)
    tmp['y'] = tmp['y'].astype(int)
    tmp['cq_id'] = tmp['cq_id'] + tmp['y']
    tmp.drop(columns=['cat_length', 'item_label'], inplace=True)
    
    return tmp

#%%
def _make_idata_from_advi_fit(fit) -> object:
    """
    Build an ArviZ InferenceData object from a CmdStan ADVI fit.

    Parameters
    ----------
    fit : CmdStanVB
        CmdStan variational fit object.

    Returns
    -------
    arviz.InferenceData
        InferenceData containing posterior draws reconstructed from
        flattened CmdStan ADVI columns.
    """
    import re
    import numpy as np
    import pandas as pd
    import arviz as az

    # CmdStanVB API differs across cmdstanpy versions.
    if hasattr(fit, "draws_pd"):
        draws_df = fit.draws_pd()
    elif hasattr(fit, "variational_sample_pd"):
        draws_df = fit.variational_sample_pd
        if callable(draws_df):
            draws_df = draws_df()
    elif hasattr(fit, "variational_sample") and hasattr(fit, "column_names"):
        draws_df = pd.DataFrame(fit.variational_sample, columns=list(fit.column_names))
    else:
        raise AttributeError(
            "Cannot extract ADVI draws from CmdStanVB fit: expected one of "
            "draws_pd, variational_sample_pd, or variational_sample+column_names."
        )

    if not isinstance(draws_df, pd.DataFrame):
        draws_df = pd.DataFrame(draws_df)

    draws_df = draws_df.drop(columns=["chain__", "iter__", "draw__"], errors="ignore")

    n_draws = len(draws_df)
    col_pattern = re.compile(r"^(?P<name>[^\[]+)(?:\[(?P<idx>[0-9,]+)\])?$")

    grouped_cols = {}
    for col in draws_df.columns:
        match = col_pattern.match(col)
        if match is None:
            raise ValueError(f"Unrecognized draw column format: {col}")

        base_name = match.group("name")
        idx_str = match.group("idx")
        idx_tuple = None if idx_str is None else tuple(int(i) for i in idx_str.split(","))
        grouped_cols.setdefault(base_name, []).append((idx_tuple, col))

    posterior_dict = {}
    dims = {}

    for base_name, entries in grouped_cols.items():
        first_idx = entries[0][0]

        if first_idx is None:
            posterior_dict[base_name] = draws_df[entries[0][1]].to_numpy().reshape(1, n_draws)
            continue

        ndim = len(first_idx)
        if any(idx is None for idx, _ in entries):
            raise ValueError(f"Mixed scalar/indexed columns for parameter '{base_name}'")
        if any(len(idx) != ndim for idx, _ in entries):
            raise ValueError(f"Inconsistent index dimensionality for parameter '{base_name}'")

        shape = tuple(max(idx[d] for idx, _ in entries) for d in range(ndim))
        arr = np.full((n_draws, *shape), np.nan, dtype=float)

        for idx, col in entries:
            arr[(slice(None),) + tuple(i - 1 for i in idx)] = draws_df[col].to_numpy()

        if np.isnan(arr).any():
            raise ValueError(f"Incomplete indexed draws for parameter '{base_name}'")

        posterior_dict[base_name] = arr[np.newaxis, ...]
        dims[base_name] = [f"{base_name}_dim_{i}" for i in range(ndim)]

    return az.from_dict(posterior=posterior_dict, dims=dims)

#%%
def _pdf_or_parts_exist(output_file_stem: str) -> bool:
    """
    Return True if a combined ``{stem}.pdf`` exists, or if at least the first
    split panel ``{stem}_part1.pdf`` exists. Used by the ``_plot_prob_barplots``
    callers to skip re-rendering when prior output is already on disk under
    either naming scheme (cowpatch combined vs cowpatch-missing per-panel).
    """
    import os
    return (
        os.path.exists(f"{output_file_stem}.pdf")
        or os.path.exists(f"{output_file_stem}_part1.pdf")
    )

#%%
def _get_autoguide_factory(algorithm: str):
    """
    Get the NumPyro autoguide factory class for the specified algorithm.

    Shared across the SVI fit functions for credit, ordered logit, and partial
    credit models.

    Parameters
    ----------
    algorithm : str
        One of 'AutoLaplaceApproximation', 'AutoMultivariateNormal',
        'AutoLowRankMultivariateNormal', 'AutoDiagonalNormal', 'AutoIAFNormal'.

    Returns
    -------
    Callable
        Autoguide factory class.
    """
    from numpyro.infer.autoguide import (
        AutoLaplaceApproximation,
        AutoMultivariateNormal,
        AutoLowRankMultivariateNormal,
        AutoDiagonalNormal,
        AutoIAFNormal,
    )

    factories = {
        'AutoLaplaceApproximation': AutoLaplaceApproximation,
        'AutoMultivariateNormal': AutoMultivariateNormal,
        'AutoLowRankMultivariateNormal': AutoLowRankMultivariateNormal,
        'AutoDiagonalNormal': AutoDiagonalNormal,
        'AutoIAFNormal': AutoIAFNormal,
    }
    if algorithm not in factories:
        raise ValueError(
            f"Unknown algorithm '{algorithm}'. Must be one of: {', '.join(factories.keys())}"
        )
    return factories[algorithm]

#%%
def _make_idata_from_svi_posterior(stan_data: dict, posterior_samples: dict, predictions: dict):
    """
    Build an ArviZ InferenceData object from NumPyro SVI posterior samples.

    Reshapes (n_draws, *shape) sample arrays into (1, n_draws, *shape) so the
    resulting idata matches the chain/draw layout used by the Stan back-ends.
    Generated quantities in ``predictions`` are placed into the same
    ``posterior`` group as the latent parameters to keep the on-disk zarr
    layout consistent across SVI/ADVI/HMC.
    """
    import numpy as np
    import arviz as az

    _ = stan_data
    posterior_dict = {}
    dims = {}

    for var_name, samples in posterior_samples.items():
        if var_name.startswith('_'):
            continue
        samples_arr = np.asarray(samples)
        if samples_arr.ndim == 1:
            posterior_dict[var_name] = samples_arr.reshape(1, -1)
        else:
            posterior_dict[var_name] = samples_arr[np.newaxis, ...]
            ndim = samples_arr.ndim - 1
            dims[var_name] = [f"{var_name}_dim_{i}" for i in range(ndim)]

    for var_name, samples in predictions.items():
        samples_arr = np.asarray(samples)
        if samples_arr.ndim == 1:
            posterior_dict[var_name] = samples_arr.reshape(1, -1)
        else:
            posterior_dict[var_name] = samples_arr[np.newaxis, ...]
            ndim = samples_arr.ndim - 1
            dims[var_name] = [f"{var_name}_dim_{i}" for i in range(ndim)]

    return az.from_dict(posterior=posterior_dict, dims=dims)

#%%
def _summarize_ordered_prob_quantiles(
    ordered_prob_values: np.ndarray,
    dcati: pd.DataFrame,
    dit: pd.DataFrame,
) -> pd.DataFrame:
    """
    Summarize ordered probability quantiles.

    Parameters
    ----------
    ordered_prob_values : np.ndarray
        Array of shape (chain, draw, cq_id) from posterior ordered probability samples.
    dcati : pd.DataFrame
        Pre-processed data containing item_type_id, item_time_id, item_label,
        time_label, y, and y_label columns.
    dit : pd.DataFrame
        Item metadata containing item_type_id, item_label, group_label_long,
        item_label_short, and endpoint_measure columns.
    
    Returns
    -------
    pd.DataFrame
        Summary DataFrame with columns: cq_id, q_lower, iqr_lower, median, iqr_upper,
        q_upper, item_type_id, item_time_id, y, item_label, time_label,
        group_label_long, item_label_short, endpoint_measure.
    """
    import numpy as np
    import pandas as pd

    po = ordered_prob_values  # shape: (chain, draw, cq_id)
    n_chains, n_draws, n_cq = po.shape
    po = po.reshape(-1, n_cq)  # (chain*draw, cq_id)

    pos = pd.DataFrame({
        'cq_id': range(1, n_cq + 1),
        'q_lower': np.percentile(po, 2.5, axis=0),
        'iqr_lower': np.percentile(po, 25, axis=0),
        'median': np.percentile(po, 50, axis=0),
        'iqr_upper': np.percentile(po, 75, axis=0),
        'q_upper': np.percentile(po, 97.5, axis=0),
    })

    tmp = _map_cq_id_to_item_structure(dcati, dit)
    tmp['cq_id'] = tmp['cq_id'] + 1
    tmp2 = dcati[['item_type_id', 'item_time_id', 'item_label', 'time_label']].drop_duplicates()
    tmp = tmp.merge(tmp2, on=['item_type_id', 'item_time_id'], how='left')
    tmp2 = dcati[['item_type_id', 'y', 'y_label']].drop_duplicates()
    tmp = tmp.merge(tmp2, on=['item_type_id', 'y'], how='left')
    pos = pos.merge(tmp, on='cq_id', how='left')
    pos = pos.merge(dit[['item_type_id', 'item_label', 'group_label_long', 
                         'item_label_short', 'endpoint_measure']],
                    on=['item_type_id', 'item_label'], how='left')

    return pos

#%%
def _plot_ppcheck(
    ypred_values: np.ndarray,
    dcati: pd.DataFrame,
    output_file_stem: str,
    verbose: bool = True,
) -> None:
    """
    Generate and save a posterior predictive check plot.

    Parameters
    ----------
    ypred_values : np.ndarray
        Array of shape (chain, draw, oid) from ``idata.posterior['ypred'].values``.
    dcati : pd.DataFrame
        Pre-processed data containing oid, item_type, item_type_id, oidt,
        y_stan, item_label, and time_label columns.
    output_file_stem : str
        Output file path without extension (e.g. ``".../prefix_ppcheck"``).
    """
    vprint = print if verbose else (lambda *args, **kwargs: None)
    import numpy as np
    import pandas as pd

    from plotnine import (
        ggplot, aes, geom_linerange, geom_point, facet_grid,
        scale_y_continuous, labs, theme_bw, theme,
        element_blank,
    )
    from ggsci import scale_color_npg 

    n_chains, n_draws, n_obs = ypred_values.shape
    ypred_flat = ypred_values.reshape(-1, n_obs)  # (chain*draw, oid)

    pos = pd.DataFrame({
        'oid': range(1, n_obs + 1),
        'q_lower': np.percentile(ypred_flat, 2.5, axis=0),
        'iqr_lower': np.percentile(ypred_flat, 25, axis=0),
        'median': np.percentile(ypred_flat, 50, axis=0),
        'iqr_upper': np.percentile(ypred_flat, 75, axis=0),
        'q_upper': np.percentile(ypred_flat, 97.5, axis=0),
    })

    tmp = dcati[['oid', 'item_type', 'item_type_id', 'oidt', 'y_stan',
                 'item_label', 'time_label']].copy()
    pos = pos.merge(tmp, on='oid')
    pos['in_ppi'] = (pos['y_stan'] >= pos['q_lower']) & (pos['y_stan'] <= pos['q_upper'])
    pos['in_ppi_label'] = np.where(pos['in_ppi'], 'TRUE', 'FALSE')

    vprint(f"Proportion in 95% PPI: {pos['in_ppi'].mean():.4f}")

    p = (
        ggplot(pos, aes(x='oid', group='oid')) +
        geom_linerange(aes(ymin='q_lower', ymax='q_upper'), color='#BDBDBD', size=0.35) +
        geom_linerange(aes(ymin='iqr_lower', ymax='iqr_upper'), color='#4D4D4D', size=0.7) +
        geom_point(aes(y='median'), color='#1A1A1A', size=0.7) +
        geom_point(aes(y='y_stan', color='in_ppi_label'), size=0.9) +
        facet_grid('item_label ~ time_label', scales='free') +
        scale_y_continuous() +
        scale_color_npg() +
        labs(
            x='',
            y='outcome',
            color='within\n95% posterior\nprediction\ninterval',
        ) +
        theme_bw() +
        theme(
            axis_text_x=element_blank(),
            axis_ticks_major_x=element_blank(),
            figure_size=(20, 20),
        )
    )
    p.save(
        f"{output_file_stem}.pdf",
        width=20,
        height=30,
        units='in',
        dpi=150,
        verbose=False,
        limitsize=False,
    )

#%%
def _compute_ordinal_brier_scores(
    ordinal_brier_score_values: np.ndarray,
    dcati: pd.DataFrame,
    output_file_stem: str,
    verbose: bool = True,
) -> None:
    """
    Compute ordinal Brier scores by question and save to CSV.

    Parameters
    ----------
    ordinal_brier_score_values : np.ndarray
        Array of shape (chain, draw, c, q) from
        ``idata.posterior['ordinal_brier_score'].values``.
    dcati : pd.DataFrame
        Pre-processed data containing item_type_id, item_time_id, item_type,
        item_label, and time_label columns.
    output_file_stem : str
        Output file path without extension (e.g. ``".../prefix_ordered_brierscore"``).
    """
    vprint = print if verbose else (lambda *args, **kwargs: None)
    import numpy as np
    import pandas as pd

    po = ordinal_brier_score_values  # shape: (chain, draw, c, q)
    n_chains, n_draws, n_cat, n_q = po.shape
    po = po.reshape(-1, n_cat, n_q)  # (chain*draw, c, q)

    rows = []
    for c in range(n_cat):
        for q in range(n_q):
            brier_vals = po[:, c, q]
            rows.append({
                'item_type_id': c + 1,
                'item_time_id': q + 1,
                'median_brier_score': np.median(brier_vals),
                'mean_brier_score': np.mean(brier_vals),
                'q025_brier_score': np.percentile(brier_vals, 2.5),
                'q975_brier_score': np.percentile(brier_vals, 97.5),
            })

    pos = pd.DataFrame(rows)
    tmp = dcati[['item_type_id', 'item_time_id', 'item_type',
                 'item_label', 'time_label']].drop_duplicates()
    pos = pos.merge(tmp, on=['item_type_id', 'item_time_id'])
    pos = pos.sort_values(['item_type_id', 'item_time_id'])

    brier_file = f"{output_file_stem}.csv"
    pos.to_csv(brier_file, index=False)
    vprint(f"Ordinal Brier scores saved to {brier_file}")

#%%
def _plot_prob_barplots(
    ordered_prob_values: np.ndarray,
    dcati: pd.DataFrame,
    dit: pd.DataFrame,
    output_file_stem: str,
    verbose: bool = True,
) -> None:
    """
    Generate and save posterior probability bar plots per item type.

    Parameters
    ----------
    ordered_prob_values : np.ndarray
        Array of shape (chain, draw, cq_id) from
        ``idata.posterior['ordered_prob_by_cat_qu_fit'].values``.
    dcati : pd.DataFrame
        Pre-processed data containing item_type_id, item_time_id, item_label,
        time_label, y, and y_label columns.
    dit : pd.DataFrame
        Item metadata containing item_type_id, item_label, group_label_long,
        item_label_short, and endpoint_measure columns.
    output_file_stem : str
        Output file path without extension (e.g. ``".../prefix_prob_by_question"``).
    """
    vprint = print if verbose else (lambda *args, **kwargs: None)
    import numpy as np
    import pandas as pd
    import ggsci
    try:
        import cowpatch as cow
    except ImportError:
        cow = None
    import matplotlib.colors as mcolors

    from plotnine import (
        ggplot, aes, geom_col, geom_boxplot, facet_grid,
        scale_y_continuous, scale_fill_manual, theme_bw, theme, labs, guides,
        element_text, position_dodge, guide_legend,
    )
    from plotnine import ggplot as _p9ggplot

    po = ordered_prob_values  # shape: (chain, draw, cq_id)
    n_chains, n_draws, n_cq = po.shape
    po = po.reshape(-1, n_cq)  # (chain*draw, cq_id)

    pos = pd.DataFrame({
        'cq_id': range(1, n_cq + 1),
        'q_lower': np.percentile(po, 2.5, axis=0),
        'iqr_lower': np.percentile(po, 25, axis=0),
        'median': np.percentile(po, 50, axis=0),
        'iqr_upper': np.percentile(po, 75, axis=0),
        'q_upper': np.percentile(po, 97.5, axis=0),
    })

    tmp = _map_cq_id_to_item_structure(dcati, dit)
    tmp['cq_id'] = tmp['cq_id'] + 1
    tmp2 = dcati[['item_type_id', 'item_time_id', 'item_label', 'time_label']].drop_duplicates()
    tmp = tmp.merge(tmp2, on=['item_type_id', 'item_time_id'], how='left')
    tmp2 = dcati[['item_type_id', 'y', 'y_label']].drop_duplicates()
    tmp = tmp.merge(tmp2, on=['item_type_id', 'y'], how='left')
    pos = pos.merge(tmp, on='cq_id', how='left')

    tmp = dcati.groupby(['time_label', 'item_label', 'y_label']).size().reset_index(name='n')
    tmp2 = dcati.groupby(['time_label', 'item_label']).size().reset_index(name='total')
    tmp = tmp.merge(tmp2, on=['time_label', 'item_label'])
    tmp['p_emp'] = tmp['n'] / tmp['total']
    pos = pos.merge(
        tmp[['time_label', 'item_label', 'y_label', 'p_emp']],
        on=['time_label', 'item_label', 'y_label'],
        how='left',
    )
    pos['p_emp'] = pos['p_emp'].fillna(0)

    pos = pos.merge(dit, on=['item_type_id', 'item_label'])
    pos['item_label_long'] = pos['group_label_long'] + np.where(
        pos['item_label_short'].notna(), '\n' + pos['item_label_short'], ''
    )

    csv_out = f"{output_file_stem}.csv"
    pos.to_csv(csv_out, index=False)
    vprint(f"Saving posterior probabilities to {csv_out}")

    all_items = list(pd.unique(pos['item_label_long']))
    base_pal = ggsci.pal_futurama("planetexpress")(12)
    pal = [
        mcolors.to_hex(c)
        for c in mcolors.LinearSegmentedColormap.from_list("futurama", base_pal)(
            np.linspace(0, 1, len(all_items))
        )
    ]
    color_dict = dict(zip(all_items, pal))

    subplot_plots = []
    subplot_widths = []
    subplot_heights = []
    for item_type_id in sorted(pos['item_type_id'].dropna().unique()):
        pos_i = pos[pos['item_type_id'] == item_type_id].copy()
        if pos_i.empty:
            continue
        pos_i['plot_group'] = pos_i['item_label_long'] + ':' + pos_i['y_label'].astype(str)
        items_i = list(pd.unique(pos_i['item_label_long']))
        color_dict_i = {k: color_dict[k] for k in items_i}
        n_time = max(1, pos_i['time_label'].nunique())
        subplot_width = max(8, 4 * n_time)
        subplot_height = 15
        y_label_i = str(pos_i['endpoint_measure'].dropna().iloc[0])
        p_i = (
            ggplot(pos_i, aes(x='y_label', group='plot_group')) +
            geom_col(
                aes(fill='item_label_long', y='p_emp'),
                position=position_dodge(width=0.9, preserve='single'),
                alpha=0.8,
                width=0.8,
            ) +
            geom_boxplot(
                aes(
                    ymin='q_lower',
                    lower='iqr_lower',
                    middle='median',
                    upper='iqr_upper',
                    ymax='q_upper',
                ),
                position=position_dodge(width=0.9, preserve='single'),
                stat='identity',
                alpha=0,
                width=0.3,
            ) +
            scale_y_continuous(labels=lambda l: [f'{v:.0%}' for v in l]) +
            scale_fill_manual(values=color_dict_i, limits=items_i, drop=False) +
            facet_grid('group_label_long ~ time_label') +
            theme_bw() +
            theme(
                axis_text_x=element_text(angle=45, vjust=1, hjust=1),
                legend_position='bottom',
                plot_title=element_text(size=10),
                figure_size=(subplot_width, subplot_height),
            ) +
            guides(fill=guide_legend(ncol=3)) +
            labs(x='', y=y_label_i, fill='')
        )
        subplot_plots.append(p_i)
        subplot_widths.append(subplot_width)
        subplot_heights.append(subplot_height)

    if not subplot_plots:
        return

    if cow is None:
        vprint(
            "cowpatch is not installed; saving separate plot files with _partN suffix instead of a combined panel."
        )
        for idx, p_i in enumerate(subplot_plots, start=1):
            part_height = subplot_heights[idx - 1]
            if idx == 1:
                part_height *= 0.5
            elif idx == 2:
                part_height *= 1.4
            p_i.save(
                f"{output_file_stem}_part{idx}.png",
                width=subplot_widths[idx - 1],
                height=part_height,
                units='in',
                limitsize=False,
            )
            p_i.save(
                f"{output_file_stem}_part{idx}.pdf",
                width=subplot_widths[idx - 1],
                height=part_height,
                units='in',
                limitsize=False,
            )
        return

    n_plots = len(subplot_plots)
    rel_h = [1] + [3.5] * (n_plots - 1)
    combined_patch = cow.patch(*subplot_plots) + cow.layout(
        ncol=1, nrow=n_plots, rel_heights=rel_h
    )

    # cowpatch calls ggplot.save(..., limitsize=True) internally, so
    # we monkey-patch save_helper to force limitsize=False for the duration.
    _orig_save_helper = _p9ggplot.save_helper

    def _unlimited_save_helper(self, *args, **kwargs):
        kwargs['limitsize'] = False
        return _orig_save_helper(self, *args, **kwargs)

    _p9ggplot.save_helper = _unlimited_save_helper
    try:
        raw_width = max(subplot_widths)
        raw_height = sum(subplot_heights)
        combined_patch.save(
            f"{output_file_stem}.png",
            width=raw_width,
            height=raw_height,
            dpi=150,
            verbose=False,
        )
        combined_patch.save(
            f"{output_file_stem}.pdf",
            width=raw_width,
            height=raw_height,
            verbose=False,
        )
    finally:
        _p9ggplot.save_helper = _orig_save_helper

#%%
def _plot_worst_chain_traces(
    po: pd.DataFrame,
    iter_warmup: int,
    output_file_stem: str,
) -> None:
    """
    Plot lp__ and worst-variable traces for the worst HMC chains.

    Parameters
    ----------
    po : pd.DataFrame
        Draws table containing lp__, chain/draw identifiers, and parameter columns.
    iter_warmup : int
        Number of warmup iterations used in sampling.
    output_file_stem : str
        Output file path without extension.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    meta_cols = {'chain__', '.chain', 'chain', 'iter__', 'draw__', 'iteration', '__iter__'}
    worst_vars = [c for c in po.columns if c not in meta_cols and c != 'lp__']
    if not worst_vars:
        return

    po = po.copy()
    chain_col = next((col for col in ['chain__', '.chain', 'chain'] if col in po.columns), None)
    iter_col = next((col for col in ['iter__', 'draw__', 'iteration'] if col in po.columns), None)

    if chain_col is not None and 'lp__' in po.columns:
        chain_summary = po.groupby(chain_col)['lp__'].mean().sort_values()
        worst_chains = chain_summary.head(min(4, len(chain_summary))).index.tolist()
    else:
        worst_chains = []

    if chain_col is not None and not worst_chains:
        worst_chains = list(pd.unique(po[chain_col]))[:4]

    if chain_col is not None and worst_chains:
        po = po[po[chain_col].isin(worst_chains)].copy()

    if iter_col is None:
        if chain_col is not None:
            po['__iter__'] = po.groupby(chain_col).cumcount() + 1
        else:
            po['__iter__'] = np.arange(len(po)) + 1
        iter_col = '__iter__'

    n_var_rows = int(np.ceil(len(worst_vars) / 2))
    fig = plt.figure(figsize=(20, 5 + 3.8 * n_var_rows), constrained_layout=True)
    outer_gs = fig.add_gridspec(2, 1, height_ratios=[1, 4], hspace=0.45)

    ax_lp = fig.add_subplot(outer_gs[0, 0])
    if 'lp__' in po.columns:
        if chain_col is not None and worst_chains:
            for chain_id in worst_chains:
                chain_data = po[po[chain_col] == chain_id]
                ax_lp.plot(chain_data[iter_col], chain_data['lp__'], linewidth=0.7, alpha=0.85, label=f"chain {chain_id}")
        else:
            ax_lp.plot(po[iter_col], po['lp__'], linewidth=0.8, alpha=0.9, color='C0')
    ax_lp.axvline(iter_warmup, color='0.4', linestyle='--', linewidth=1)
    ax_lp.set_title('lp__ trace for worst chains')
    ax_lp.set_xlabel('iteration')
    ax_lp.set_ylabel('lp__')
    if chain_col is not None and worst_chains:
        ax_lp.legend(loc='upper right', fontsize=8, frameon=False)
    ax_lp.grid(alpha=0.2)

    inner_gs = outer_gs[1, 0].subgridspec(n_var_rows, 2, hspace=0.95, wspace=0.25)
    for idx, var_name in enumerate(worst_vars):
        ax = fig.add_subplot(inner_gs[idx])
        if chain_col is not None and worst_chains:
            for chain_id in worst_chains:
                chain_data = po[po[chain_col] == chain_id]
                ax.plot(chain_data[iter_col], chain_data[var_name], linewidth=0.6, alpha=0.8)
        else:
            ax.plot(po[iter_col], po[var_name], linewidth=0.6, alpha=0.8)
        ax.axvline(iter_warmup, color='0.4', linestyle='--', linewidth=0.8)
        ax.set_title(var_name, fontsize=9, pad=10)
        if idx >= (n_var_rows - 1) * 2:
            ax.set_xlabel('iteration')
            ax.tick_params(axis='x', labelbottom=True)
        else:
            ax.set_xlabel('')
            ax.tick_params(axis='x', labelbottom=False)
        ax.grid(alpha=0.2)

    for idx in range(len(worst_vars), n_var_rows * 2):
        fig.add_subplot(inner_gs[idx]).axis('off')

    fig.savefig(f"{output_file_stem}.pdf", bbox_inches='tight')
    plt.close(fig)

#%%
def _futurama_palette(n: int) -> list[str]:
    """Return n hex colours interpolating the 12-colour ggsci futurama palette."""
    from ggsci import pal_futurama
    from matplotlib.colors import LinearSegmentedColormap, to_hex

    base = list(pal_futurama()(12))
    if n <= 12:
        return base[:n]
    cmap = LinearSegmentedColormap.from_list('futurama_x', base, N=n)
    return [to_hex(cmap(i / max(n - 1, 1))) for i in range(n)]
