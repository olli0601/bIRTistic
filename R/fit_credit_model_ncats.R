#' Run credit model analysis (ncats version) on pre-processed data
#'
#' This function performs Bayesian IRT analysis using the flexible ncats credit model
#' on pre-processed data. It compiles the Stan model, runs MCMC sampling, generates convergence
#' diagnostics, and optionally creates detailed diagnostic plots and posterior predictive checks.
#'
#' @param dit data.table. Item metadata table with item types and labels
#' @param dcati data.table. Pre-processed data with observations for analysis. Must include
#'   item_type, y_stan, item_time_id, pid, and oidt columns.
#' @param output_file_prefix Character. Full path prefix for output files (without extension)
#' @param stan_file Character. Path to Stan model file (.stan). Default: credit_model_ncats_v260413.stan
#' @param x_formula Formula. Formula specifying predictors for the design matrix
#' @param x_formula_ignore_regex Character. Regular expression pattern to identify X columns to exclude 
#'   from posterior predictive probability calculations (ordered_prob_by_cat_qu_pr). Columns matching this 
#'   pattern will be excluded. If NA (default), all X columns are used for both fitted and predictive probabilities.
#' @param chains Integer. Number of MCMC chains to run (default: 2)
#' @param parallel_chains Integer. Number of chains to run in parallel (default: 2)
#' @param threads_per_chain Integer. Number of threads to use per chain (default: 1)
#' @param iter_warmup Integer. Number of warmup iterations per chain (default: 500)
#' @param iter_sampling Integer. Number of sampling iterations per chain (default: 1500)
#' @param seed Integer. Random seed for reproducibility (default: 123)
#' @param show_messages Logical. If TRUE, show all Stan informational messages during sampling (default: FALSE)
#' @param show_exceptions Logical. If TRUE, show detailed Stan exception messages when errors occur (default: FALSE)
#' @param resume Logical. If TRUE and all output files exist, skip MCMC sampling and load existing results (default: FALSE)
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
fit_credit_model_ncats <- function(
  dit,
  dcati,
  output_file_prefix,
  stan_file = here::here("src", "stan", "credit_model_ncats_v260413.stan"),
  x_formula = ~ time - 1,
  x_formula_ignore_regex = NA_character_,
  chains = 2L,
  parallel_chains = 2L,
  threads_per_chain = 1L,
  iter_warmup = 500L,
  iter_sampling = 1500L,
  seed = 123L,
  show_messages = FALSE,
  show_exceptions = FALSE,
  resume = FALSE,
  with_core_analyses = TRUE,
  with_additional_analyses = FALSE
) {
    # Suppress data.table NSE warnings
    item_type <- y_stan <- item_time_id <- pid <- time <- oidt <- variable <-
        ypred <- oid <- in_ppi <- item_label <- time_label <- y_label <- y <-
        n <- total <- p_emp <- group_label_long <- item_label_short <-
        item_label_long <- ess_bulk <- q_lower <- q_upper <- iqr_lower <-
        iqr_upper <- prob <- item_type_id <- stat <- NULL

    require(data.table)
    require(ggplot2)
    require(ggsci)
    require(bayesplot)
    require(cmdstanr)
    require(cowplot)

    # Print configuration
    cat("\n========================================\n")
    cat("Credit Model (ncats) Analysis Configuration\n")
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
    cat("Resume:", resume, "\n")
    cat("Core analyses:", with_core_analyses, "\n")
    cat("Additional analyses:", with_additional_analyses, "\n")
    cat("========================================\n\n")

    # Create output directory if it doesn't exist
    dir.create(dirname(output_file_prefix), showWarnings = FALSE, recursive = TRUE)
    
    # Define output file paths
    timing_file <- paste0(output_file_prefix, "_timing.csv")
    draws_file <- paste0(output_file_prefix, "_draws.rds")
    output_file <- paste0(output_file_prefix, "_stan.rds")
    mixing_file <- paste0(output_file_prefix, "_convergence_mixing.csv")
    
    # Check if we should resume from existing outputs
    can_resume <- resume && 
                  file.exists(timing_file) && 
                  file.exists(draws_file) && 
                  file.exists(output_file) && 
                  file.exists(mixing_file)

    # Compile Stan model
    cat("Compiling Stan model...\n")
    cm_compiled <- cmdstanr::cmdstan_model(
        stan_file,
        include_paths = dirname(stan_file),
        cpp_options = list(stan_threads = TRUE),
        force_recompile = TRUE
    )

    # Prepare data in ncats format
    cat("Preparing Stan data in ncats format...\n")
    
    # Ensure oid is sequential and ordered by item_type_id and oidt
    setkey(dcati, item_type_id, oidt)
    stopifnot(all( 1:nrow(dcati) == dcati[,oid] ))  

    # Build stan_data
    stan_data <- list()
    stan_data$C <- dcati[, max(item_type_id)]
    stan_data$U <- dcati[, max(pid)] 
    
    # Arrays for each category type
    tmp <- dcati[, .(
        N = .N,
        Q = max(item_time_id),
        K = length(unique(y_stan))
    ), by = item_type_id]
    stan_data$N <- tmp$N
    stan_data$Q <- tmp$Q
    stan_data$K <- tmp$K
    
    # Total counts
    stan_data$N_total <- sum(stan_data$N)
    stan_data$Q_total <- sum(stan_data$Q)
    
    # Concatenated observation data (data is already sorted by item_type_id)
    stan_data$y <- dcati$y_stan
    stan_data$unit_of_obs <- dcati$pid
    stan_data$question_of_obs <- dcati$item_time_id
    stan_data$cat_type <- dcati$item_type_id
    
    # Create design matrix for all observations
    stan_data$X <- model.matrix(
        x_formula,
        data = as.data.frame(dcati)
    )
    if(!is.na(x_formula_ignore_regex)) {
        stan_data$Xpr_id <- which(!grepl(x_formula_ignore_regex, colnames(stan_data$X)))
    }
    else {
        stan_data$Xpr_id <- seq_len(ncol(stan_data$X))
    }
    stan_data$P <- ncol(stan_data$X)
    stan_data$P_pr <- length(stan_data$Xpr_id)
    
    cat("Design matrix number of predictors: P =", stan_data$P, "\n")
    cat("Design matrix column names (for fitting):", paste(colnames(stan_data$X), collapse = ", "), "\n")
    cat("Design matrix column names (for prediction):", paste(colnames(stan_data$X)[stan_data$Xpr_id], collapse = ", "), "\n")
    cat("Number of category types: C =", stan_data$C, "\n")
    cat("Category type labels:", paste(unique(dit$item_type), collapse = ", "), "\n")
    cat("N per category:", paste(stan_data$N, collapse = ", "), "\n")
    cat("Q per category:", paste(stan_data$Q, collapse = ", "), "\n")
    cat("K per category:", paste(stan_data$K, collapse = ", "), "\n")

    # Check if resuming from existing outputs
    if (can_resume) {
        cat("\n========================================\n")
        cat("RESUMING from existing outputs\n")
        cat("========================================\n")
        
        cat("Loading draws into poa from:", draws_file, "\n")
        poa <- readRDS(draws_file)
        
        cat("Loading timing data from:", timing_file, "\n")
        timing_data <- read.csv(timing_file)
        cm_fit_good_chains <- timing_data$chain
        cat("Identified good HMC chains:", paste(cm_fit_good_chains, collapse = ", "), "\n")
        cat("========================================\n\n")
    } else {
        if (resume) {
            cat("\nNote: Resume requested but not all output files exist. Running full analysis.\n\n")
        }
        
        # Sample from the model
        cat("Running MCMC sampling...\n")
        flush.console()
        cm_fit <- cm_compiled$sample(
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
        tmp <- as.data.table(cm_fit$draws(variables = "lp__", inc_warmup = FALSE, format = "draws_df"))
        tmp <- tmp[, .(
            mean_lp = mean(lp__),
            sd_lp = sd(lp__),
            n = .N
        ), by = .chain]
        
        threshold <- tmp[which.max(mean_lp), mean_lp - 2 * sd_lp]
        cm_fit_good_chains <- tmp[mean_lp > threshold, .chain]
        cat("Identified good HMC chains:", paste(cm_fit_good_chains, collapse = ", "), "\n")

        # Extract chain timing information
        cat("\nExtracting timing information...\n")
        chain_times <- cm_fit$time()
        timing_data <- data.table(
            chain = cm_fit_good_chains,
            warmup_minutes = chain_times$chains$warmup[cm_fit_good_chains] / 60,
            sampling_minutes = chain_times$chains$sampling[cm_fit_good_chains] / 60,
            total_chain_minutes = chain_times$chains$total[cm_fit_good_chains] / 60
        )
        
        write.csv(timing_data, file = timing_file, row.names = FALSE)
        cat("Saved timing information to:", timing_file, "\n")

        # Save output to RDS
        cat("Saving model fit to:", output_file, "\n")
        cm_fit$save_object(file = output_file)

        # Check convergence and mixing
        cat("Generating convergence diagnostics...\n")
        tmp <- cm_fit$summary(
            variables = c(
                "latent_factor_unit", "latent_factor_beta",
                "skill_thresholds", "loadings_questions_m1"
            )
        )
        tmp <- as.data.table(tmp)
        tmp <- tmp[order(ess_bulk), ]
        write.csv(
            tmp,
            file = mixing_file,
            row.names = FALSE
        )

        # Worst parameters with lowest ess_bulk
        worst_var <- tmp$variable[1:9]

        # Make worst trace plot
        cat("Generating trace plots...\n")        
        po <- cm_fit$draws(
            variables = c("lp__", worst_var),
            inc_warmup = TRUE,
            format = "draws_array"
        )
        po <- po[, cm_fit_good_chains, , drop = FALSE]
        
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

        # Extract and save draws
        cat("Extracting draws...\n")
        poa <- cm_fit$draws(format = "draws_array")
        poa <- poa[, cm_fit_good_chains, , drop = FALSE]
        saveRDS(poa, file = draws_file)
        cat("Saved draws to:", draws_file, "\n")

        cm_fit <- NULL
        gc()

    }

    # Additional diagnostic analyses (optional)
    if (with_additional_analyses) {
        cat("\nRunning additional diagnostic analyses...\n")

        # Make intervals plot
        cat("Generating parameter plots...\n")
        tmp <- c("latent_factor_unit", "latent_factor_beta","skill_thresholds", "loadings_questions_m1")                
        po <- poa[,, grepl(paste(tmp, collapse = "|"), dimnames(poa)[[3]]), drop = FALSE]
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
        po <- poa[,, grepl("ypred", dimnames(poa)[[3]]), drop = FALSE]
        po <- posterior::as_draws_df(po)
        po <- as.data.table(po)
        po <- data.table::melt(
            po,
            id.vars = c(".draw", ".chain", ".iteration"),
            value.name = "ypred"
        )
        po[, oid := as.integer(gsub("ypred\\[([0-9]+)\\]", "\\1", variable))]
        set(po, NULL, "variable", NULL)
        pos <- po[,
            list(
                summary_value = quantile(
                    ypred,
                    prob = c(0.025, 0.25, 0.5, 0.75, 0.975)
                ),
                summary_name = c("q_lower", "iqr_lower", "median", "iqr_upper", "q_upper")
            ),
            by = c("oid")
        ]
        pos <- data.table::dcast(pos,
            oid ~ summary_name,
            value.var = "summary_value"
        )
        tmp <- dcati[, .(oid, item_type, item_type_id, oidt, y_stan, item_label, time_label)]
        pos <- merge(pos, tmp, by = "oid")
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
        brier <- poa[,, grepl("ordinal_brier_score", dimnames(poa)[[3]]), drop = FALSE]
        brier <- posterior::as_draws_df(brier)
        brier <- as.data.table(brier)
        brier <- data.table::melt(brier,
            id.vars = c(".draw", ".chain", ".iteration"),
            value.name = "brier_score"
        )
        
        # Parse indices: ordinal_brier_score[c,q]
        brier[, item_type_id := as.integer(gsub(
            "ordinal_brier_score\\[([0-9]+),([0-9]+)\\]", "\\1",
            variable
        ))]
        brier[, item_time_id := as.integer(gsub(
            "ordinal_brier_score\\[([0-9]+),([0-9]+)\\]", "\\2",
            variable
        ))]
        set(brier, NULL, "variable", NULL)
        
        # Compute median
        brier_summary <- brier[,
            list(
                median_brier_score = median(brier_score),
                mean_brier_score = mean(brier_score),
                q025_brier_score = quantile(brier_score, 0.025),
                q975_brier_score = quantile(brier_score, 0.975)
            ),
            by = c("item_type_id", "item_time_id")
        ]
        tmp <- unique(subset(dcati, select = c(item_type_id, item_time_id, item_type, item_label, time_label)))
        brier_summary <- merge(brier_summary, tmp, by = c("item_type_id", "item_time_id"))
        brier_summary <- brier_summary[order(item_type_id, item_time_id), ]

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
        # NOTE: In ncats model, ordered_prob_by_obs stores averaged probabilities
        # per (category type, question) combination.
        # Format: for each category c, for each question q in that category,
        # store K[c] probabilities (one per response category)        
        po <- poa[,, grepl("ordered_prob_by_cat_qu_fit", dimnames(poa)[[3]]), drop = FALSE]
        po <- posterior::as_draws_df(po)
        po <- as.data.table(po)
        po <- data.table::melt(po,
            id.vars = c(".draw", ".chain", ".iteration"),
            value.name = "prob"
        )
        
        # Parse index from ordered_prob_by_cat_qu_fit[i]
        po[, cq_id := as.integer(gsub(
            "ordered_prob_by_cat_qu_fit\\[([0-9]+)\\]", "\\1",
            variable
        ))]
        set(po, NULL, "variable", NULL)

        # Compute quantiles
        pos <- po[,
            list(
                summary_value = quantile(
                    prob,
                    prob = c(0.025, 0.25, 0.5, 0.75, 0.975)
                ),
                summary_name = c("q_lower", "iqr_lower", "median", "iqr_upper", "q_upper")
            ),
            by = c("cq_id")
        ]
        po <- NULL
        gc()

        pos <- data.table::dcast(pos,
            cq_id ~ summary_name,
            value.var = "summary_value"
        )

        # Map cq_id back to (item_type_id, item_time_id, y_stan)
        tmp <- unique(subset(dcati, select = c(item_type_id, item_label, item_time_id)))
        tmp <- merge(tmp, subset(dit, select = c(item_type_id, item_label, cat_length)), by = c("item_type_id", "item_label"))
        setkey(tmp, item_type_id, item_time_id)
        tmp[, cq_id := cumsum(cat_length) - cat_length + 1L]
        tmp <- tmp[, list(cq_id = cq_id - 1L + 1:cat_length, y = 1:cat_length - 1L), by = c("item_type_id", "item_time_id")]
        tmp <- merge(tmp, 
                     unique(subset(dcati, 
                                    select = c("item_type_id", "item_time_id", "y", "y_label", "item_label", "time_label")
                                    )), 
                     by = c("item_type_id", "item_time_id", "y")
                     )
        pos <- merge(tmp, pos, by = "cq_id")              
        
        # Compute empirical probabilities
        tmp <- dcati[,
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
        
        # Merge with dit for plotting
        pos <- merge(pos, dit, by = c("item_type_id","item_label"))
        pos[, item_label_long := paste0(group_label_long, " --- ", item_label_short)]

        cat(
            "\nSaving posterior probabilities to RDS with name",
            paste0(output_file_prefix, "_posterior_probabilities_fit.rds"),
            "...\n"
        )
        saveRDS(
            pos,
            file = paste0(output_file_prefix, "_posterior_probabilities_fit.rds")
        )

        # Make plots
        # Create named color palette to ensure consistent colors across all subplots
        all_items <- unique(pos$item_label_long)
        pal <- colorRampPalette(ggsci::pal_futurama("planetexpress")(12))(length(all_items))
        names(pal) <- all_items

        # Get unique groups and times
        groups <- unique(pos$group_label_long)
        times <- unique(pos$time_label)
        
        # Create individual plots for each group-time combination
        plot_list <- list()
        for (g in groups) {
            for (t in times) {
                plot_data <- pos[group_label_long == g & time_label == t]
                
                p_sub <- ggplot(
                    plot_data,
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
                    scale_fill_manual(values = pal, limits = all_items, drop = FALSE) +
                    ggtitle(paste0(g, " - ", t)) +
                    theme_bw() +
                    theme(
                        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                        legend.position = "none",
                        plot.title = element_text(size = 10)
                    ) +
                    labs(x = "", y = "proportion of outcomes", fill = "survey items")
                
                plot_list[[paste(g, t, sep = "_")]] <- p_sub
            }
        }
        
        # Extract legend from one plot
        p_legend <- ggplot(
            pos,
            aes(x = y_label, group = interaction(item_label_long, y_label))
        ) +
            geom_col(aes(fill = item_label_long, y = p_emp)) +
            scale_fill_manual(values = pal, limits = all_items, drop = FALSE) +
            theme_bw() +
            theme(legend.position = "bottom") +
            labs(fill = "survey items") +
            guides(fill = guide_legend(ncol = 3))
        
        legend <- cowplot::get_legend(p_legend)
        
        # Arrange plots in grid: 4 columns (2 groups × 2 times each)
        p <- cowplot::plot_grid(
            plotlist = plot_list,
            ncol = 4,
            align = "hv",
            axis = "tb"
        )
        
        # Add legend below
        p <- cowplot::plot_grid(p, legend, ncol = 1, rel_heights = c(1, 0.1))
        ggsave(
            file = paste0(output_file_prefix, "_probs_barplot_v2_fit.png"),
            plot = p,
            h = 30,
            w = 24
        )
        ggsave(
            file = paste0(output_file_prefix, "_probs_barplot_v2_fit.pdf"),
            plot = p,
            h = 30,
            w = 24
        )
        
        # Generate plots for posterior predictive probabilities if x_formula_ignore_regex is specified
        if (!is.na(x_formula_ignore_regex)) {
            cat("Generating posterior predictive probability plots (ordered_prob_by_cat_qu_pr)...\n")
            
            # Extract ordered_prob_by_cat_qu_pr
            po2 <- poa[,, grepl("ordered_prob_by_cat_qu_pr", dimnames(poa)[[3]]), drop = FALSE]            
            po2 <- posterior::as_draws_df(po2)
            po2 <- as.data.table(po2)
            po2 <- data.table::melt(po2,
                id.vars = c(".draw", ".chain", ".iteration"),
                value.name = "prob"
            )
            
            # Parse index from ordered_prob_by_cat_qu_pr[i]
            po2[, cq_id := as.integer(gsub(
                "ordered_prob_by_cat_qu_pr\\[([0-9]+)\\]", "\\1",
                variable
            ))]
            set(po2, NULL, "variable", NULL)

            # Compute quantiles
            pos2 <- po2[,
                list(
                    summary_value = quantile(
                        prob,
                        prob = c(0.025, 0.25, 0.5, 0.75, 0.975)
                    ),
                    summary_name = c("q_lower", "iqr_lower", "median", "iqr_upper", "q_upper")
                ),
                by = c("cq_id")
            ]
            po2 <- NULL
            gc()

            pos2 <- data.table::dcast(pos2,
                cq_id ~ summary_name,
                value.var = "summary_value"
            )

            # Map cq_id back to (item_type_id, item_time_id, y_stan)
            tmp <- unique(subset(dcati, select = c(item_type_id, item_label, item_time_id)))
            tmp <- merge(tmp, subset(dit, select = c(item_type_id, item_label, cat_length)), by = c("item_type_id", "item_label"))
            setkey(tmp, item_type_id, item_time_id)
            tmp[, cq_id := cumsum(cat_length) - cat_length + 1L]
            tmp <- tmp[, list(cq_id = cq_id - 1L + 1:cat_length, y = 1:cat_length - 1L), by = c("item_type_id", "item_time_id")]
            tmp <- merge(tmp, 
                         unique(subset(dcati, 
                                        select = c("item_type_id", "item_time_id", "y", "y_label", "item_label", "time_label")
                                        )), 
                         by = c("item_type_id", "item_time_id", "y")
                         )
            pos2 <- merge(tmp, pos2, by = "cq_id")                            
            
            # Merge with dit for plotting
            pos2 <- merge(pos2, dit, by = c("item_type_id","item_label"))
            pos2[, item_label_long := paste0(group_label_long, " --- ", item_label_short)]

            cat(
                "\nSaving posterior predictive probabilities to RDS with name",
                paste0(output_file_prefix, "_posterior_predictive_probabilities_pr.rds"),
                "...\n"
            )
            saveRDS(
                pos2,
                file = paste0(output_file_prefix, "_posterior_predictive_probabilities_pr.rds")
            )

            pos[, stat:= 'fitted']
            pos2[, stat:= 'predictive']
            pos <- rbind(pos, pos2, use.names = TRUE, fill = TRUE)
            
            # Create comparison plots (fitted vs predictive)
            all_items <- unique(pos$item_label_long)
            pal <- colorRampPalette(ggsci::pal_futurama("planetexpress")(12))(length(all_items))
            names(pal) <- all_items

            # Get unique groups and times
            groups <- unique(pos$group_label_long)
            times <- unique(pos$time_label)      
            
            # Create individual comparison plots for each group-time combination
            ps <- list()
            for (g in groups) {
                for (t in times) {
                    tmp <- pos[group_label_long == g & time_label == t]                    
                    p <- ggplot(
                        tmp,
                        aes(x = y_label, group = interaction(stat, item_label_long, y_label))
                    ) +
                    geom_boxplot(
                        aes(
                                ymin = q_lower, lower = iqr_lower,
                                middle = median, upper = iqr_upper,
                                ymax = q_upper,
                                fill = item_label_long,
                                linetype = stat,
                                alpha = stat
                            ),
                            position = position_dodge2(width = 0.9, preserve = "single"),
                            stat = "identity",
                            outlier.shape = NA
                        ) +
                        scale_alpha_manual(values = c("fitted" = 0.7, "predictive" = 0.3)) +
                        scale_linetype_manual(values = c("fitted" = "solid", "predictive" = "dashed")) +
                        scale_y_continuous(labels = scales::percent) +
                        scale_fill_manual(values = pal, limits = all_items, drop = FALSE) +
                        ggtitle(paste0(g, " - ", t)) +
                        theme_bw() +
                        theme(
                            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                            legend.position = "none",
                            plot.title = element_text(size = 10)
                        ) +
                        labs(x = "", y = "proportion of outcomes", fill = "survey items", linetype = "Model", alpha = "Model")
                    
                    ps[[paste(g, t, sep = "_")]] <- p
                }
            }
            
            # Extract combined legend showing both fill and linetype/alpha
            # Create minimal dataset for legend (one row per item_label_long x stat x y_label combination)
            pos_legend <- unique(pos[, .(item_label_long, stat, y_label, q_lower, iqr_lower, median, iqr_upper, q_upper)])            
            p_legend_combined <- ggplot(
                pos_legend,
                aes(x = y_label, group = interaction(stat, item_label_long, y_label))
            ) +
            geom_col(
                    aes(
                        y = median, 
                        fill = item_label_long,
                        linetype = stat,
                        alpha = stat
                    )                    
            ) +
            scale_alpha_manual(values = c("fitted" = 0.7, "predictive" = 0.3)) +
            scale_linetype_manual(values = c("fitted" = "solid", "predictive" = "dashed")) +
            scale_fill_manual(values = pal, limits = all_items, drop = FALSE) +
            theme_bw() +
            theme(legend.position = "bottom") +
            labs(fill = "survey items", linetype = "Model", alpha = "Model") +
            guides(
                fill = guide_legend(ncol = 3, order = 1),
                linetype = guide_legend(order = 2),
                alpha = guide_legend(order = 2)
            )            
            legend_combined <- cowplot::get_legend(p_legend_combined)
            
            # Arrange comparison plots in grid
            p <- cowplot::plot_grid(
                plotlist = ps,
                ncol = 4,
                align = "hv",
                axis = "tb"
            )
            
            # Add legend below
            p <- cowplot::plot_grid(p, legend_combined, ncol = 1, rel_heights = c(1, 0.15))
            ggsave(
                file = paste0(output_file_prefix, "_probs_comparison_fit_vs_pr.png"),
                plot = p,
                h = 30,
                w = 24
            )
            ggsave(
                file = paste0(output_file_prefix, "_probs_comparison_fit_vs_pr.pdf"),
                plot = p,
                h = 30,
                w = 24
            )
        }
    }
}
