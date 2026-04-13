#' Run partial credit model analysis (ncats version) on pre-processed data
#'
#' This function performs Bayesian IRT analysis using the flexible ncats partial credit model
#' on pre-processed data. It compiles the Stan model, runs MCMC sampling, generates convergence
#' diagnostics, and optionally creates detailed diagnostic plots and posterior predictive checks.
#'
#' @param dit data.table. Item metadata table with item types and labels
#' @param dcati data.table. Pre-processed data with observations for analysis. Must include
#'   item_type, y_stan, item_time_id, pid, and oidt columns.
#' @param output_file_prefix Character. Full path prefix for output files (without extension)
#' @param stan_file Character. Path to Stan model file (.stan). Default: partial_credit_model_ncats_v260413.stan
#' @param x_formula Formula. Formula specifying predictors for the design matrix
#' @param chains Integer. Number of MCMC chains to run (default: 2)
#' @param parallel_chains Integer. Number of chains to run in parallel (default: 2)
#' @param threads_per_chain Integer. Number of threads to use per chain (default: 1)
#' @param iter_warmup Integer. Number of warmup iterations per chain (default: 500)
#' @param iter_sampling Integer. Number of sampling iterations per chain (default: 1500)
#' @param seed Integer. Random seed for reproducibility (default: 123)
#' @param show_messages Logical. If TRUE, show all Stan informational messages during sampling (default: FALSE)
#' @param show_exceptions Logical. If TRUE, show detailed Stan exception messages when errors occur (default: FALSE)
#' @param with_core_analyses Logical. If TRUE, generate core probability plots (default: TRUE)
#' @param with_additional_analyses Logical. If TRUE, generate additional diagnostic plots (default: FALSE)
#'
#' @return Invisibly returns the fitted cmdstanr model object
#'
#' @import data.table
#' @import ggplot2
#' @import cmdstanr
#'
#' @export
fit_partial_credit_model_ncats <- function(
  dit,
  dcati,
  output_file_prefix,
  stan_file = here::here("src", "stan", "partial_credit_model_ncats_v260413.stan"),
  x_formula = ~ time - 1,
  chains = 2L,
  parallel_chains = 2L,
  threads_per_chain = 1L,
  iter_warmup = 500L,
  iter_sampling = 1500L,
  seed = 123L,
  show_messages = FALSE,
  show_exceptions = FALSE,
  with_core_analyses = TRUE,
  with_additional_analyses = FALSE
) {
    # Suppress data.table NSE warnings
    item_type <- y_stan <- item_time_id <- pid <- time <- oidt <- variable <-
        ypred <- oid <- in_ppi <- item_label <- time_label <- y_label <- y <-
        n <- total <- p_emp <- group_label_long <- item_label_short <-
        item_label_long <- ess_bulk <- q_lower <- q_upper <- iqr_lower <-
        iqr_upper <- prob <- cat_type_id <- NULL

    require(data.table)
    require(ggplot2)
    require(ggsci)
    require(bayesplot)
    require(cmdstanr)

    # Print configuration
    cat("\n========================================\n")
    cat("Partial Credit Model (ncats) Analysis Configuration\n")
    cat("========================================\n")
    cat("Stan file:", stan_file, "\n")
    cat("Stan include dir:", dirname(stan_file), "\n")
    cat("Data: dit with", nrow(dit), "items, dcati with", nrow(dcati), "observations\n")
    cat("Output prefix:", output_file_prefix, "\n")
    cat("Chains:", chains, "\n")
    cat("Parallel chains:", parallel_chains, "\n")
    cat("Threads per chain:", threads_per_chain, "\n")
    cat("Warmup iterations:", iter_warmup, "\n")
    cat("Sampling iterations:", iter_sampling, "\n")
    cat("Seed:", seed, "\n")
    cat("Core analyses:", with_core_analyses, "\n")
    cat("Additional analyses:", with_additional_analyses, "\n")
    cat("========================================\n\n")

    # Create output directory if it doesn't exist
    dir.create(dirname(output_file_prefix), showWarnings = FALSE, recursive = TRUE)

    # Compile Stan model
    cat("Compiling Stan model...\n")
    pcm_compiled <- cmdstanr::cmdstan_model(
        stan_file,
        include_paths = dirname(stan_file),
        cpp_options = list(stan_threads = TRUE)
    )

    # Prepare data in ncats format
    cat("Preparing Stan data in ncats format...\n")
    
    # Get unique item types (sorted)
    item_types <- sort(unique(dcati$item_type))
    C <- length(item_types)
    
    # Assign cat_type_id based on sorted item types
    dcati_copy <- copy(dcati)
    dcati_copy[, cat_type_id := match(item_type, item_types)]
    
    # Sort data by cat_type_id to ensure proper ordering
    setkey(dcati_copy, cat_type_id, oidt)
    
    # Build stan_data
    stan_data <- list()
    stan_data$C <- C
    stan_data$U <- max(dcati_copy$pid)
    
    # Arrays for each category type
    stan_data$N <- integer(C)
    stan_data$Q <- integer(C)
    stan_data$K <- integer(C)
    
    for (i in 1:C) {
        item_type_i <- item_types[i]
        stan_data$N[i] <- nrow(dcati_copy[item_type == item_type_i])
        stan_data$Q[i] <- max(dcati_copy[item_type == item_type_i, item_time_id])
        stan_data$K[i] <- length(unique(dcati_copy[item_type == item_type_i, y_stan]))
    }
    
    # Total counts
    stan_data$N_total <- sum(stan_data$N)
    stan_data$Q_total <- sum(stan_data$Q)
    
    # Concatenated observation data (data is already sorted by cat_type_id)
    stan_data$y <- dcati_copy$y_stan
    stan_data$unit_of_obs <- dcati_copy$pid
    stan_data$question_of_obs <- dcati_copy$item_time_id
    stan_data$cat_type <- dcati_copy$cat_type_id
    
    # Create design matrix for all observations
    stan_data$X <- model.matrix(
        x_formula,
        data = as.data.frame(dcati_copy)
    )
    stan_data$P <- ncol(stan_data$X)
    
    cat("Design matrix number of predictors: P =", stan_data$P, "\n")
    cat("Design matrix column names:", paste(colnames(stan_data$X), collapse = ", "), "\n")
    cat("Number of category types: C =", C, "\n")
    cat("Category type labels:", paste(item_types, collapse = ", "), "\n")
    cat("N per category:", paste(stan_data$N, collapse = ", "), "\n")
    cat("Q per category:", paste(stan_data$Q, collapse = ", "), "\n")
    cat("K per category:", paste(stan_data$K, collapse = ", "), "\n")

    # Sample from the model
    cat("Running MCMC sampling...\n")
    flush.console()
    pcm_fit <- pcm_compiled$sample(
        data = stan_data,
        seed = seed,
        chains = chains,
        parallel_chains = parallel_chains,
        threads_per_chain = threads_per_chain,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        refresh = 500,
        save_warmup = TRUE,
        show_messages = show_messages,
        show_exceptions = show_exceptions
    )
    
    # Remove any trajectories that did not converge
    cat("Checking for divergent transitions and removing non-converged chains...\n")
    tmp <- as.data.table(pcm_fit$draws(variables = "lp__", inc_warmup = FALSE, format = "draws_df"))
    tmp <- tmp[, .(
        mean_lp = mean(lp__),
        sd_lp = sd(lp__),
        n = .N
    ), by = .chain]
    
    threshold <- tmp[which.max(mean_lp), mean_lp - 2 * sd_lp]
    pcm_fit_good_chains <- tmp[mean_lp > threshold, .chain]
    cat("Identified good HMC chains:", paste(pcm_fit_good_chains, collapse = ", "), "\n")

    # Extract chain timing information
    cat("\nExtracting timing information...\n")
    chain_times <- pcm_fit$time()
    timing_data <- data.table(
        chain = pcm_fit_good_chains,
        warmup_minutes = chain_times$chains$warmup[pcm_fit_good_chains] / 60,
        sampling_minutes = chain_times$chains$sampling[pcm_fit_good_chains] / 60,
        total_chain_minutes = chain_times$chains$total[pcm_fit_good_chains] / 60
    )
    
    timing_file <- paste0(output_file_prefix, "_timing.csv")
    write.csv(timing_data, file = timing_file, row.names = FALSE)
    cat("Saved timing information to:", timing_file, "\n")

    # Extract and save draws
    cat("Extracting draws...\n")
    draws_file <- paste0(output_file_prefix, "_draws.rds")
    draws <- pcm_fit$draws(format = "draws_array")
    draws <- draws[, pcm_fit_good_chains, , drop = FALSE]
    saveRDS(draws, file = draws_file)
    cat("Saved draws to:", draws_file, "\n")
    draws <- NULL
    gc()

    # Save output to RDS
    output_file <- paste0(output_file_prefix, "_stan.rds")
    cat("Saving model fit to:", output_file, "\n")
    pcm_fit$save_object(file = output_file)

    # Check convergence and mixing
    cat("Generating convergence diagnostics...\n")
    tmp <- pcm_fit$summary(
        variables = c(
            "latent_factor_unit", "latent_factor_beta",
            "skill_thresholds", "loadings_questions_m1"
        )
    )
    tmp <- as.data.table(tmp)
    tmp <- tmp[order(ess_bulk), ]
    write.csv(
        tmp,
        file = paste0(output_file_prefix, "_convergence_mixing.csv"),
        row.names = FALSE
    )

    # Additional diagnostic analyses (optional)
    if (with_additional_analyses) {
        cat("\nRunning additional diagnostic analyses...\n")

        # Worst parameters with lowest ess_bulk
        worst_var <- tmp$variable[1:9]

        # Make worst trace plot
        cat("Generating trace plots...\n")
        po <- pcm_fit$draws(
            variables = c("lp__", worst_var),
            inc_warmup = TRUE,
            format = "draws_array"
        )
        po <- po[, pcm_fit_good_chains, , drop = FALSE]
        
        lp_range <- range(po[(iter_warmup + 1):(iter_warmup + iter_sampling), , "lp__"])
        lp_ylim <- c(lp_range[1] - 0.05 * diff(lp_range), 
                     lp_range[2] + 0.05 * diff(lp_range))
        
        p_lp <- bayesplot:::mcmc_trace(po[,, "lp__", drop = FALSE],
            pars = "lp__",
            n_warmup = iter_warmup
        ) + 
        coord_cartesian(ylim = lp_ylim) +
        theme_bw()
        
        p_other <- bayesplot:::mcmc_trace(po[,, worst_var, drop = FALSE],
            pars = worst_var,
            n_warmup = iter_warmup,
            facet_args = list(ncol = 2)
        ) + 
        theme_bw()
        
        p <- patchwork::wrap_plots(p_lp, p_other, ncol = 1, heights = c(1, 4))
        
        ggsave(
            file = paste0(output_file_prefix, "_worsttrace.png"),
            plot = p,
            h = 10,
            w = 20
        )

        # Make intervals plot
        cat("Generating parameter plots...\n")
        po <- pcm_fit$draws(
            variables = c(
                "latent_factor_unit", "latent_factor_beta",
                "skill_thresholds", "loadings_questions_m1"
            ),
            inc_warmup = FALSE,
            format = "draws_array"
        )
        po <- po[, pcm_fit_good_chains, , drop = FALSE]

        color_scheme_set("teal")
        p <- bayesplot::mcmc_intervals(
            po,
            prob = 0.5, prob_outer = 0.95, outer_size = 1, point_size = 2
        ) + theme_bw()
        ggsave(
            file = paste0(output_file_prefix, "_intervals.pdf"),
            plot = p,
            h = 50,
            w = 8,
            limitsize = FALSE
        )

        # Make posterior predictive check
        cat("Generating posterior predictive checks...\n")
        po <- pcm_fit$draws(
            variables = c("ypred"),
            inc_warmup = FALSE,
            format = "draws_array"
        )
        po <- po[, pcm_fit_good_chains, , drop = FALSE]
        po <- posterior::as_draws_df(po)
        po <- as.data.table(po)
        po <- data.table::melt(
            po,
            id.vars = c(".draw", ".chain", ".iteration"),
            value.name = "ypred"
        )
        po[, obs_idx := as.integer(gsub("ypred\\[([0-9]+)\\]", "\\1", variable))]
        set(po, NULL, "variable", NULL)

        # Map back to original observation structure using sorted dcati_copy
        dcati_copy[, obs_idx := .I]
        po <- merge(po, dcati_copy[, .(obs_idx, item_type, oidt, y_stan, oid, item_label, time_label)], 
                    by = "obs_idx")

        pos <- po[,
            list(
                summary_value = quantile(
                    ypred,
                    prob = c(0.025, 0.25, 0.5, 0.75, 0.975)
                ),
                summary_name = c("q_lower", "iqr_lower", "median", "iqr_upper", "q_upper")
            ),
            by = c("item_type", "oidt")
        ]
        pos <- data.table::dcast(pos,
            item_type + oidt ~ summary_name,
            value.var = "summary_value"
        )

        pos <- merge(pos, unique(dcati_copy[, .(item_type, oidt, y_stan, oid, item_label, time_label)]), 
                     by = c("item_type", "oidt"))
        pos[, in_ppi := y_stan >= q_lower & y_stan <= q_upper]
        cat("Proportion in 95% PPI:", pos[, mean(in_ppi)], "\n")

        # Plot posterior predictive check
        p <- ggplot(pos, aes(x = oid, group = oid)) +
            geom_boxplot(
                aes(
                    ymin = q_lower,
                    lower = iqr_lower,
                    middle = median,
                    upper = iqr_upper,
                    ymax = q_upper
                ),
                stat = "identity"
            ) +
            geom_point(aes(y = y_stan, colour = in_ppi)) +
            facet_grid(item_label ~ time_label, scales = "free") +
            scale_x_discrete() +
            scale_y_continuous() +
            ggsci::scale_color_npg() +
            labs(
                x = "",
                y = "outcome",
                colour = "within\n95% posterior\nprediction\ninterval"
            ) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        ggsave(
            file = paste0(output_file_prefix, "_ppcheck.png"),
            p,
            w = 20,
            h = 20
        )

        # Compute and save ordinal Brier scores by question
        cat("Computing ordinal Brier scores by question...\n")
        
        # Extract Brier scores (now in format ordinal_brier_score[c,q])
        brier <- pcm_fit$draws(
            variables = c("ordinal_brier_score"),
            inc_warmup = FALSE,
            format = "draws_array"
        )
        brier <- brier[, pcm_fit_good_chains, , drop = FALSE]
        brier <- posterior::as_draws_df(brier)
        brier <- as.data.table(brier)
        brier <- data.table::melt(brier,
            id.vars = c(".draw", ".chain", ".iteration"),
            value.name = "brier_score"
        )
        
        # Parse indices: ordinal_brier_score[c,q]
        brier[, cat_id := as.integer(gsub(
            "ordinal_brier_score\\[([0-9]+),([0-9]+)\\]", "\\1",
            variable
        ))]
        brier[, question_id := as.integer(gsub(
            "ordinal_brier_score\\[([0-9]+),([0-9]+)\\]", "\\2",
            variable
        ))]
        set(brier, NULL, "variable", NULL)
        
        # Map cat_id back to item_type
        cat_type_map <- data.table(cat_id = 1:C, item_type = item_types)
        brier <- merge(brier, cat_type_map, by = "cat_id")
        
        # Compute median
        brier_summary <- brier[,
            list(
                median_brier_score = median(brier_score),
                mean_brier_score = mean(brier_score),
                q025_brier_score = quantile(brier_score, 0.025),
                q975_brier_score = quantile(brier_score, 0.975)
            ),
            by = c("item_type", "question_id")
        ]
        
        data.table::fwrite(
            brier_summary,
            file = paste0(output_file_prefix, "_ordered_brierscore.csv")
        )
        cat("Ordinal Brier scores saved to", paste0(output_file_prefix, "_ordered_brierscore.csv"), "\n")
    }

    # Generate probability plots
    if (with_core_analyses) {
        cat("Generating probability plots...\n")
        
        # Extract ordered_prob_by_obs (long vector format)
        # NOTE: In ncats model, ordered_prob_by_obs is a single long vector [eta_length]
        # where probabilities for each observation are stored contiguously.
        # This differs from the 2cats model which has separate matrices
        # cat1_ordered_prob_by_obs[Ncat1, Kcat1] and cat2_ordered_prob_by_obs[Ncat2, Kcat2]
        po <- pcm_fit$draws(
            variables = c("ordered_prob_by_obs"),
            inc_warmup = FALSE,
            format = "draws_array"
        )
        po <- po[, pcm_fit_good_chains, , drop = FALSE]
        po <- posterior::as_draws_df(po)
        po <- as.data.table(po)
        po <- data.table::melt(po,
            id.vars = c(".draw", ".chain", ".iteration"),
            value.name = "prob"
        )
        
        # Parse index from ordered_prob_by_obs[i]
        po[, prob_idx := as.integer(gsub(
            "ordered_prob_by_obs\\[([0-9]+)\\]", "\\1",
            variable
        ))]
        set(po, NULL, "variable", NULL)
        
        # Map prob_idx back to observation and response category
        # Need to reconstruct eta_obs_offset from stan_data
        K_vec <- stan_data$K[stan_data$cat_type]  # K for each observation
        eta_obs_offset <- cumsum(c(1, K_vec[-length(K_vec)]))  # Starting position for each obs
        
        # Create mapping table
        mapping <- data.table()
        for (obs_idx in 1:stan_data$N_total) {
            K_obs <- K_vec[obs_idx]
            offset <- eta_obs_offset[obs_idx]
            obs_mapping <- data.table(
                prob_idx = offset:(offset + K_obs - 1),
                obs_idx = obs_idx,
                y_stan = 1:K_obs
            )
            mapping <- rbind(mapping, obs_mapping)
        }
        
        # Merge with probability data
        po <- merge(po, mapping, by = "prob_idx")
        
        # Merge with observation metadata
        dcati_copy[, obs_idx := .I]
        po <- merge(po, dcati_copy[, .(obs_idx, item_type, oidt, pid, item_time_id, time, 
                                       time_label, item_label)], 
                    by = "obs_idx")
        
        # Average over observations at same time/item/response
        po <- po[,
            list(prob = mean(prob)),
            by = c(".draw", "time", "item_time_id", "y_stan", "item_type")
        ]
        
        # Compute quantiles
        pos <- po[,
            list(
                summary_value = quantile(
                    prob,
                    prob = c(0.025, 0.25, 0.5, 0.75, 0.975)
                ),
                summary_name = c("q_lower", "iqr_lower", "median", "iqr_upper", "q_upper")
            ),
            by = c("time", "item_time_id", "y_stan", "item_type")
        ]
        pos <- data.table::dcast(pos,
            time + item_time_id + y_stan + item_type ~ summary_name,
            value.var = "summary_value"
        )
        
        # Add labels
        tmp <- unique(dcati_copy[, .(time, time_label)])
        pos <- merge(pos, tmp, by = "time")
        tmp <- unique(dcati_copy[, .(item_time_id, item_label, item_type)])
        pos <- merge(pos, tmp, by = c("item_time_id", "item_type"))
        tmp <- unique(dcati_copy[, .(y, y_stan, y_label)])
        pos <- merge(pos, tmp, by = "y_stan")
        
        # Compute empirical probabilities
        tmp <- dcati_copy[,
            list(n = .N),
            by = c("time_label", "item_label", "y_label")
        ]
        tmp2 <- tmp[, list(total = sum(n)), by = c("time_label", "item_label")]
        tmp <- merge(tmp, tmp2, by = c("time_label", "item_label"))
        tmp[, p_emp := n / total]
        pos <- merge(pos,
            subset(tmp, select = -c(n, total)),
            by = c("time_label", "item_label", "y_label"),
            all.x = TRUE
        )
        set(pos, pos[, which(is.na(p_emp))], "p_emp", 0.)
        pos[, y_stan := NULL]
        
        # Merge with dit for plotting
        pos <- merge(pos, dit, by = "item_label")
        pos[, item_label_long := paste0(group_label_long, " --- ", item_label_short)]

        cat(
            "\nSaving posterior probabilities to RDS with name",
            paste0(output_file_prefix, "__posterior_probabilities.rds"),
            "...\n"
        )
        saveRDS(
            pos,
            file = paste0(output_file_prefix, "__posterior_probabilities.rds")
        )

        # Make plots
        pal <- colorRampPalette(ggsci::pal_futurama("planetexpress")(12))(
            pos[, length(unique(item_label))]
        )

        p <- ggplot(
            pos,
            aes(x = y_label, group = interaction(item_label_long, y_label))
        ) +
            geom_col(aes(fill = item_label_long, y = p_emp),
                position = position_dodge(width = 0.9, preserve = "single"),
                alpha = 0.8,
                width = 0.8
            ) +
            geom_boxplot(
                aes(
                    ymin = q_lower, lower = iqr_lower,
                    middle = median, upper = iqr_upper,
                    ymax = q_upper
                ),
                position = position_dodge(0.9, preserve = "single"),
                stat = "identity",
                alpha = 0,
                width = 0.3
            ) +
            scale_y_continuous(labels = scales::percent) +
            scale_fill_manual(values = pal) +
            facet_wrap(group_label_long ~ time_label, scales = "free_x", ncol = 2) +
            theme_bw() +
            theme(
                axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                legend.position = "bottom"
            ) +
            labs(x = "", y = "proportion of outcomes", fill = "survey items") +
            guides(fill = guide_legend(ncol = 3))
        ggsave(
            file = paste0(output_file_prefix, "_probs_barplot_v2.png"),
            plot = p,
            h = 40,
            w = 12
        )
        ggsave(
            file = paste0(output_file_prefix, "_probs_barplot_v2.pdf"),
            plot = p,
            h = 40,
            w = 12
        )
    }

    invisible(pcm_fit)
}
