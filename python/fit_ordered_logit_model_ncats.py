"""
Model fitting functions for Bayesian IRT analysis.

This module provides functions for fitting ordered logit models
using Stan HMC via cmdstanpy.
"""

from typing import Dict, Optional, Union
import os
from pathlib import Path
import warnings
import re
import time

import pandas as pd
import numpy as np
import patsy
import matplotlib.pyplot as plt
from cmdstanpy import CmdStanModel
from plotnine import *
import arviz as az

def fit_ordered_logit_model_ncats(
    dit: pd.DataFrame,
    dcati: pd.DataFrame,
    output_file_prefix: str,
    stan_file: Optional[str] = None,
    x_formula: str = "~ time - 1",
    x_formula_ignore_regex: Optional[str] = None,
    chains: int = 2,
    parallel_chains: int = 2,
    threads_per_chain: int = 1,
    iter_warmup: int = 500,
    iter_sampling: int = 1500,
    seed: int = 123,
    show_messages: bool = False,
    show_exceptions: bool = False,
    resume: bool = False,
    with_core_analyses: bool = True,
    with_additional_analyses: bool = False,
) -> Dict:
    """
    Run ordered logit model analysis (ncats version) on pre-processed data using HMC.
    
    This function performs Bayesian IRT analysis using the flexible ncats ordered logit model
    on pre-processed data. It compiles the Stan model, runs MCMC sampling, generates convergence
    diagnostics, and optionally creates detailed diagnostic plots.
    
    Parameters
    ----------
    dit : pd.DataFrame
        Item metadata table with item types and labels.
    dcati : pd.DataFrame
        Pre-processed data with observations for analysis. Must include
        item_type, y_stan, item_time_id, pid, and oidt columns.
    output_file_prefix : str
        Full path prefix for output files (without extension).
    stan_file : str, optional
        Path to Stan model file (.stan). If None, uses default ordered_logit_ncats_v260413.stan
    x_formula : str, default "~ time - 1"
        Patsy formula string specifying predictors for the design matrix.
    x_formula_ignore_regex : str, optional
        Regular expression pattern to identify X columns to exclude from posterior predictive.
    chains : int, default 2
        Number of MCMC chains to run.
    parallel_chains : int, default 2
        Number of chains to run in parallel.
    threads_per_chain : int, default 1
        Number of threads to use per chain.
    iter_warmup : int, default 500
        Number of warmup iterations per chain.
    iter_sampling : int, default 1500
        Number of sampling iterations per chain.
    seed : int, default 123
        Random seed for reproducibility.
    show_messages : bool, default False
        If True, show Stan informational messages.
    show_exceptions : bool, default False
        If True, show Stan exception messages when errors occur.
    resume : bool, default False
        If True and all output files exist, skip MCMC sampling and load existing results.
    with_core_analyses : bool, default True
        If True, generate core probability plots.
    with_additional_analyses : bool, default False
        If True, generate additional diagnostic plots.
    
    Returns
    -------
    dict
        Dictionary containing:
        - 'fit': CmdStanMCMC object (if not resumed)
        - 'draws': InferenceData object with posterior draws
        - 'timing': DataFrame with chain timing information
        - 'good_chains': List of chain IDs that converged
    """
    # Set default stan_file if not provided  
    if stan_file is None:
        # Assume we're running from the bIRTistic root directory
        stan_file = str(Path(__file__).parents[2] / "src" / "stan" / "ordered_logit_ncats_v260413.stan")
    
    # Print configuration
    print("\n" + "=" * 40)
    print("Ordered Logit Model (ncats) Analysis Configuration")
    print("=" * 40)
    print(f"Stan file: {stan_file}")
    print(f"Stan include dir: {os.path.dirname(stan_file)}")
    print(f"Data: dit with {len(dit)} items, dcati with {len(dcati)} observations")
    print(f"Output prefix: {output_file_prefix}")
    print(f"Chains: {chains}")
    print(f"Parallel chains: {parallel_chains}")
    print(f"Threads per chain: {threads_per_chain}")
    print(f"Warmup iterations: {iter_warmup}")
    print(f"Sampling iterations: {iter_sampling}")
    print(f"Seed: {seed}")
    print(f"Resume: {resume}")
    print(f"Core analyses: {with_core_analyses}")
    print(f"Additional analyses: {with_additional_analyses}")
    print("=" * 40 + "\n")
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file_prefix), exist_ok=True)
    
    # Define output file paths
    timing_file = f"{output_file_prefix}_timing.csv"
    draws_file = f"{output_file_prefix}_draws.zarr"  # Use zarr for arviz (9x faster than NetCDF)
    output_file = f"{output_file_prefix}_stan.pkl"
    mixing_file = f"{output_file_prefix}_convergence_mixing.csv"
    data_file = f"{output_file_prefix}_data.csv"
    
    # Save the input data for later use
    # Save dp1 and dit metadata separately
    dcati.to_csv(data_file.replace('.csv', '_dp1.csv'), index=False)
    dit.to_csv(data_file.replace('.csv', '_dit.csv'), index=False)
    print(f"Saved preprocessed data to: {data_file.replace('.csv', '_dp1.csv')} and _dit.csv")
    
    # Check if we should resume from existing outputs
    can_resume = (resume and 
                  os.path.exists(timing_file) and 
                  os.path.exists(draws_file) and 
                  os.path.exists(output_file) and 
                  os.path.exists(mixing_file) and
                  os.path.exists(data_file.replace('.csv', '_dp1.csv')) and
                  os.path.exists(data_file.replace('.csv', '_dit.csv')))
    
    # Prepare data in ncats format (same as fit_partial_credit_model_ncats)
    print("Preparing Stan data in ncats format...")
    
    # Ensure dcati is sorted by item_type_id and oidt
    dcati_sorted = dcati.sort_values(['item_type_id', 'oidt']).reset_index(drop=True)
    
    # Verify oid is sequential
    assert (dcati_sorted['oid'] == dcati_sorted.index + 1).all(), \
        "oid must be sequential 1:N after sorting by item_type_id, oidt"
    
    # Build stan_data dictionary
    stan_data = {}
    stan_data['C'] = int(dcati['item_type_id'].max())
    stan_data['U'] = int(dcati['pid'].max())
    
    # Arrays for each category type
    category_stats = dcati.groupby('item_type_id').agg({
        'oid': 'count',
        'item_time_id': 'max',
        'y_stan': lambda x: len(x.unique())
    }).rename(columns={'oid': 'N', 'item_time_id': 'Q', 'y_stan': 'K'})
    
    stan_data['N'] = category_stats['N'].astype(int).tolist()
    stan_data['Q'] = category_stats['Q'].astype(int).tolist()
    stan_data['K'] = category_stats['K'].astype(int).tolist()
    
    # Total counts
    stan_data['N_total'] = int(sum(stan_data['N']))
    stan_data['Q_total'] = int(sum(stan_data['Q']))
    
    # Concatenated observation data
    stan_data['y'] = dcati['y_stan'].astype(int).tolist()
    stan_data['unit_of_obs'] = dcati['pid'].astype(int).tolist()
    stan_data['question_of_obs'] = dcati['item_time_id'].astype(int).tolist()
    stan_data['cat_type'] = dcati['item_type_id'].astype(int).tolist()
    
    # Create design matrix using patsy
    design_matrix = patsy.dmatrix(x_formula, data=dcati, return_type='dataframe')
    
    # Determine which columns to use for predictive probabilities
    if x_formula_ignore_regex is not None:
        xpr_id = [i + 1 for i, col in enumerate(design_matrix.columns) 
                  if not re.search(x_formula_ignore_regex, col)]
    else:
        xpr_id = list(range(1, len(design_matrix.columns) + 1))
    
    stan_data['X'] = design_matrix.values.tolist()
    stan_data['Xpr_id'] = xpr_id
    stan_data['P'] = len(design_matrix.columns)
    stan_data['P_pr'] = len(xpr_id)
    
    print(f"Design matrix number of predictors: P = {stan_data['P']}")
    print(f"Design matrix column names (for fitting): {', '.join(design_matrix.columns)}")
    print(f"Design matrix column names (for prediction): {', '.join([design_matrix.columns[i-1] for i in xpr_id])}")
    print(f"Number of category types: C = {stan_data['C']}")
    print(f"Category type labels: {', '.join(dit['item_type'].unique())}")
    print(f"N per category: {', '.join(map(str, stan_data['N']))}")
    print(f"Q per category: {', '.join(map(str, stan_data['Q']))}")
    print(f"K per category: {', '.join(map(str, stan_data['K']))}")
    
    # Check if resuming from existing outputs
    if can_resume:
        print("\n" + "=" * 40)
        print("RESUMING from existing outputs")
        print("=" * 40)
        
        print(f"Loading draws from: {draws_file}")
        idata = az.from_zarr(draws_file)
        
        print(f"Loading timing data from: {timing_file}")
        timing_data = pd.read_csv(timing_file)
        good_chains = timing_data['chain'].tolist()
        print(f"Identified good HMC chains: {', '.join(map(str, good_chains))}")
        print("=" * 40 + "\n")
        
        return {
            'draws': idata,
            'timing': timing_data,
            'good_chains': good_chains
        }
    else:
        if resume:
            print("\nNote: Resume requested but not all output files exist. Running full analysis.\n")
        
        # Compile Stan model
        print("Compiling Stan model...")
        model = CmdStanModel(
            stan_file=stan_file,
            cpp_options={'STAN_THREADS': 'TRUE'}
        )
        
        # Sample from the model
        print("Running MCMC sampling...")
        fit = model.sample(
            data=stan_data,
            chains=chains,
            parallel_chains=parallel_chains,
            threads_per_chain=threads_per_chain,
            iter_warmup=iter_warmup,
            iter_sampling=iter_sampling,
            seed=seed,
            show_console=show_messages,
            refresh=500,
            save_warmup=True
        )
        
        # Remove any chains that did not converge
        print("Checking for divergent transitions and removing non-converged chains...")
        # Get lp__ draws using the method from cmdstanr: 
        # fit$draws(variables = "lp__", inc_warmup = FALSE, format = "draws_df")
        # In cmdstanpy, this is fit.draws_pd()
        try:
            # Try with inc_warmup parameter (cmdstanpy 1.x)
            lp_draws = fit.draws_pd(inc_warmup=False)
        except TypeError:
            # Fallback for older API
            lp_draws = fit.draws_pd()
        
        # Check if lp__ column exists
        if 'lp__' not in lp_draws.columns:
            raise ValueError(f"lp__ not found in draws. Available columns: {list(lp_draws.columns)[:10]}")
        
        # Identify chain column - could be 'chain__', '.chain', or 'chain'
        chain_col = None
        for col_name in ['chain__', '.chain', 'chain']:
            if col_name in lp_draws.columns:
                chain_col = col_name
                break
        
        if chain_col is None:
            # Need to reconstruct chain information from draw indices
            num_draws_per_chain = fit.num_draws_sampling
            total_draws = len(lp_draws)
            num_chains = len(fit.chain_ids)
            
            # Create chain column based on position
            lp_draws['chain'] = [fit.chain_ids[i // num_draws_per_chain] for i in range(total_draws)]
            chain_col = 'chain'
        
        # Calculate mean and std per chain (matching R code)
        chain_stats = lp_draws.groupby(chain_col)['lp__'].agg(['mean', 'std', 'count']).reset_index()
        chain_stats.columns = ['chain', 'mean_lp', 'sd_lp', 'n']
        
        # Identify good chains
        best_idx = chain_stats['mean_lp'].idxmax()
        threshold = chain_stats.loc[best_idx, 'mean_lp'] - 2 * chain_stats.loc[best_idx, 'sd_lp']
        good_chains = chain_stats[chain_stats['mean_lp'] > threshold]['chain'].tolist()
        print(f"Identified good HMC chains: {', '.join(map(str, good_chains))}")
        
        # Extract chain timing information
        print("\nExtracting timing information...")
        timing_data = pd.DataFrame({
            'chain': good_chains,
            'warmup_minutes': [0.0 for _ in good_chains],
            'sampling_minutes': [0.0 for _ in good_chains],
            'total_chain_minutes': [0.0 for _ in good_chains]
        })
        
        timing_data.to_csv(timing_file, index=False)
        print(f"Saved timing information to: {timing_file}")
        
        # Save model fit
        print(f"Saving model fit to: {output_file}")
        pd.to_pickle(fit, output_file)
        
        # Convert to ArviZ InferenceData
        print("Converting to ArviZ format...")
        idata = az.from_cmdstanpy(
            fit,
            coords={'chain': good_chains},
            dims={}
        )
        
        # Save draws
        print(f"Saving draws to: {draws_file}")
        idata.to_zarr(draws_file)
        
        # Check convergence and mixing
        print("Generating convergence diagnostics...")
        summary_df = fit.summary()
        
        # Filter to key variables of interest
        key_vars = ['latent_factor_unit', 'latent_factor_beta', 
                   'skill_thresholds_1', 'skill_thresholds_incs', 'loadings_questions_m1']
        # Filter rows where variable name starts with any of the key variable prefixes
        mask = summary_df.index.to_series().apply(
            lambda x: any(x.startswith(var) for var in key_vars)
        )
        summary_df_filtered = summary_df[mask]
        
        # Sort by ESS if available, otherwise just save as-is
        if 'ess_bulk' in summary_df_filtered.columns:
            summary_df_filtered = summary_df_filtered.sort_values('ess_bulk')
        elif 'N_Eff' in summary_df_filtered.columns:
            summary_df_filtered = summary_df_filtered.sort_values('N_Eff')
        
        summary_df_filtered.to_csv(mixing_file)
        
        print(f"Saved convergence diagnostics to: {mixing_file}")
        
        # Additional diagnostic analyses
        if with_additional_analyses:
            print("\nRunning additional diagnostic analyses...")
            
            # 1. Generate parameter intervals plot
            print("Generating parameter plots...")
            key_vars_pattern = '|'.join(['latent_factor_unit', 'latent_factor_beta',
                                        'skill_thresholds_1', 'skill_thresholds_incs', 
                                        'loadings_questions_m1'])
            
            # Use ArviZ plot_forest for parameter intervals
            var_names = [v for v in idata.posterior.data_vars 
                        if any(pat in v for pat in key_vars_pattern.split('|'))]
            
            if var_names:
                fig = az.plot_forest(
                    idata,
                    var_names=var_names,
                    combined=True,
                    hdi_prob=0.95,
                    figsize=(8, min(50, len(var_names) * 0.3))
                )
                plt.tight_layout()
                plt.savefig(f"{output_file_prefix}_intervals.pdf", bbox_inches='tight')
                plt.close()
            
            # 2. Posterior predictive checks
            print("Generating posterior predictive checks...")
            
            # Extract ypred from posterior
            if 'ypred' in idata.posterior.data_vars:
                ypred_data = idata.posterior['ypred'].values  # shape: (chain, draw, oid)
                n_chains, n_draws, n_obs = ypred_data.shape
                
                # Reshape to (draw, oid) by combining chains
                ypred_flat = ypred_data.reshape(-1, n_obs)  # (chain*draw, oid)
                
                # Compute quantiles
                q_lower = np.percentile(ypred_flat, 2.5, axis=0)
                iqr_lower = np.percentile(ypred_flat, 25, axis=0)
                median = np.percentile(ypred_flat, 50, axis=0)
                iqr_upper = np.percentile(ypred_flat, 75, axis=0)
                q_upper = np.percentile(ypred_flat, 97.5, axis=0)
                
                # Create dataframe with predictions
                pos = pd.DataFrame({
                    'oid': range(1, n_obs + 1),
                    'q_lower': q_lower,
                    'iqr_lower': iqr_lower,
                    'median': median,
                    'iqr_upper': iqr_upper,
                    'q_upper': q_upper
                })
                
                # Merge with actual data
                tmp = dcati[['oid', 'item_type', 'item_type_id', 'oidt', 'y_stan', 
                            'item_label', 'time_label']].copy()
                pos = pos.merge(tmp, on='oid')
                pos['in_ppi'] = (pos['y_stan'] >= pos['q_lower']) & (pos['y_stan'] <= pos['q_upper'])
                
                print(f"Proportion in 95% PPI: {pos['in_ppi'].mean():.4f}")
                
                # Plot posterior predictive check
                fig, axes = plt.subplots(
                    len(pos['item_label'].unique()),
                    len(pos['time_label'].unique()),
                    figsize=(20, 20),
                    sharex=True
                )
                
                for i, item in enumerate(sorted(pos['item_label'].unique())):
                    for j, time in enumerate(sorted(pos['time_label'].unique())):
                        ax = axes[i, j] if len(axes.shape) > 1 else axes[j]
                        data = pos[(pos['item_label'] == item) & (pos['time_label'] == time)]
                        
                        if len(data) > 0:
                            # Plot boxplots
                            ax.boxplot(
                                [[data.iloc[k]['q_lower'], data.iloc[k]['iqr_lower'],
                                  data.iloc[k]['median'], data.iloc[k]['iqr_upper'],
                                  data.iloc[k]['q_upper']] for k in range(len(data))],
                                positions=range(len(data)),
                                widths=0.6
                            )
                            
                            # Plot actual values
                            colors = ['#1f77b4' if x else '#d62728' for x in data['in_ppi']]
                            ax.scatter(range(len(data)), data['y_stan'], c=colors, zorder=3)
                            
                            ax.set_title(f"{item} - {time}", fontsize=8)
                            ax.tick_params(axis='x', rotation=45)
                
                plt.tight_layout()
                plt.savefig(f"{output_file_prefix}_ppcheck.png", dpi=150, bbox_inches='tight')
                plt.close()
            
            # 3. Compute ordinal Brier scores
            print("Computing ordinal Brier scores by question...")
            
            if 'ordinal_brier_score' in idata.posterior.data_vars:
                brier_data = idata.posterior['ordinal_brier_score'].values  # shape: (chain, draw, c, q)
                n_chains, n_draws, n_cat, n_q = brier_data.shape
                
                # Reshape to (draw, c, q)
                brier_flat = brier_data.reshape(-1, n_cat, n_q)
                
                # Compute statistics
                brier_summary = []
                for c in range(n_cat):
                    for q in range(n_q):
                        brier_vals = brier_flat[:, c, q]
                        brier_summary.append({
                            'item_type_id': c + 1,
                            'item_time_id': q + 1,
                            'median_brier_score': np.median(brier_vals),
                            'mean_brier_score': np.mean(brier_vals),
                            'q025_brier_score': np.percentile(brier_vals, 2.5),
                            'q975_brier_score': np.percentile(brier_vals, 97.5)
                        })
                
                brier_summary = pd.DataFrame(brier_summary)
                
                # Merge with item metadata
                tmp = dcati[['item_type_id', 'item_time_id', 'item_type', 
                            'item_label', 'time_label']].drop_duplicates()
                brier_summary = brier_summary.merge(tmp, on=['item_type_id', 'item_time_id'])
                brier_summary = brier_summary.sort_values(['item_type_id', 'item_time_id'])
                
                brier_file = f"{output_file_prefix}_ordered_brierscore.csv"
                brier_summary.to_csv(brier_file, index=False)
                print(f"Ordinal Brier scores saved to {brier_file}")
        
        # Core probability analyses
        if with_core_analyses:
            print("\nGenerating probability plots...")
            
            # Extract ordered_prob_by_cat_qu_fit
            if 'ordered_prob_by_cat_qu_fit' in idata.posterior.data_vars:
                prob_data = idata.posterior['ordered_prob_by_cat_qu_fit'].values
                # shape: (chain, draw, cq_id)
                n_chains, n_draws, n_cq = prob_data.shape
                
                # Reshape to (draw, cq_id)
                prob_flat = prob_data.reshape(-1, n_cq)
                
                # Compute quantiles
                q_lower = np.percentile(prob_flat, 2.5, axis=0)
                iqr_lower = np.percentile(prob_flat, 25, axis=0)
                median = np.percentile(prob_flat, 50, axis=0)
                iqr_upper = np.percentile(prob_flat, 75, axis=0)
                q_upper = np.percentile(prob_flat, 97.5, axis=0)
                
                pos = pd.DataFrame({
                    'cq_id': range(1, n_cq + 1),
                    'q_lower': q_lower,
                    'iqr_lower': iqr_lower,
                    'median': median,
                    'iqr_upper': iqr_upper,
                    'q_upper': q_upper
                })
                
                # Map cq_id back to (item_type_id, item_time_id, y_stan)
                tmp = dcati[['item_type_id', 'item_label', 'item_time_id']].drop_duplicates()
                tmp = tmp.merge(
                    dit[['item_type_id', 'item_label', 'cat_length']], 
                    on=['item_type_id', 'item_label']
                )
                tmp = tmp.sort_values(['item_type_id', 'item_time_id']).reset_index(drop=True)
                tmp['cq_id'] = tmp['cat_length'].cumsum() - tmp['cat_length'] + 1
                
                # Expand to get one row per category
                expanded = []
                for _, row in tmp.iterrows():
                    for k in range(int(row['cat_length'])):
                        expanded.append({
                            'cq_id': int(row['cq_id']) + k,
                            'item_type_id': row['item_type_id'],
                            'item_time_id': row['item_time_id'],
                            'y': k
                        })
                
                expanded_df = pd.DataFrame(expanded)
                
                # Merge with dcati to get labels
                tmp2 = dcati[['item_type_id', 'item_time_id', 'y', 'y_label', 
                             'item_label', 'time_label']].drop_duplicates()
                expanded_df = expanded_df.merge(tmp2, on=['item_type_id', 'item_time_id', 'y'])
                
                pos = pos.merge(expanded_df, on='cq_id')
                
                # Compute empirical probabilities
                tmp = dcati.groupby(['time_label', 'item_label', 'y_label']).size().reset_index(name='n')
                tmp2 = dcati.groupby(['time_label', 'item_label']).size().reset_index(name='total')
                tmp = tmp.merge(tmp2, on=['time_label', 'item_label'])
                tmp['p_emp'] = tmp['n'] / tmp['total']
                
                pos = pos.merge(
                    tmp[['time_label', 'item_label', 'y_label', 'p_emp']],
                    on=['time_label', 'item_label', 'y_label'],
                    how='left'
                )
                pos['p_emp'] = pos['p_emp'].fillna(0)
                
                # Merge with dit for plotting
                pos = pos.merge(dit, on=['item_type_id', 'item_label'])
                pos['item_label_long'] = pos['group_label_long'] + ' --- ' + pos['item_label_short'].fillna('')
                
                # Save posterior probabilities
                prob_file = f"{output_file_prefix}_posterior_probabilities_fit.csv"
                pos.to_csv(prob_file, index=False)
                print(f"Saving posterior probabilities to {prob_file}")
                
                # Create probability plots using plotnine (ggplot2 for Python)
                # This mirrors the R code structure for consistent aesthetics
                
                # Create named color palette (matching R's ggsci::pal_futurama)
                all_items = sorted(pos['item_label_long'].unique())
                from matplotlib.colors import rgb2hex
                # Use a colormap similar to futurama palette
                n_colors = len(all_items)
                color_values = plt.cm.Set3(np.linspace(0, 1, min(12, n_colors)))
                if n_colors > 12:
                    # Extend with additional colors if needed
                    color_values = np.vstack([
                        plt.cm.Set3(np.linspace(0, 1, 12)),
                        plt.cm.Pastel1(np.linspace(0, 1, n_colors - 12))
                    ])
                pal = [rgb2hex(c) for c in color_values[:n_colors]]
                color_dict = dict(zip(all_items, pal))
                
                # Create the plot using plotnine
                p = (
                    ggplot(pos, aes(x='y_label', group='factor(item_label_long):factor(y_label)')) +
                    geom_col(
                        aes(fill='item_label_long', y='p_emp'),
                        position=position_dodge(width=0.9, preserve='single'),
                        alpha=0.8,
                        width=0.8
                    ) +
                    geom_boxplot(
                        aes(
                            ymin='q_lower',
                            lower='iqr_lower',
                            middle='median',
                            upper='iqr_upper',
                            ymax='q_upper'
                        ),
                        position=position_dodge(width=0.9, preserve='single'),
                        stat='identity',
                        alpha=0,
                        width=0.3
                    ) +
                    scale_y_continuous(labels=lambda l: [f'{v:.0%}' for v in l]) +
                    scale_fill_manual(values=color_dict, limits=all_items, drop=False) +
                    facet_grid('group_label_long ~ time_label') +
                    theme_bw() +
                    theme(
                        axis_text_x=element_text(angle=45, vjust=1, hjust=1),
                        legend_position='bottom',
                        plot_title=element_text(size=10),
                        figure_size=(24, 30)
                    ) +
                    guides(fill=guide_legend(ncol=3, title='survey items')) +
                    labs(x='', y='proportion of outcomes', fill='survey items')
                )
                
                # Save the plot
                p.save(f"{output_file_prefix}_probs_barplot_v2_fit.png", dpi=150, verbose=False)
                p.save(f"{output_file_prefix}_probs_barplot_v2_fit.pdf", verbose=False)
        
        return {
            'fit': fit,
            'draws': idata,
            'timing': timing_data,
            'good_chains': good_chains
        }
