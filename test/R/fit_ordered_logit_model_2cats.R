#' Run ordered logit model analysis on pre-processed data
#'
#' This function performs Bayesian IRT analysis using an ordered logit model on pre-processed
#' data. It compiles the Stan model, runs MCMC sampling, generates convergence
#' diagnostics, and optionally creates detailed diagnostic plots and posterior
#' predictive checks.
#'
#' @param dit data.table. Item metadata table with item types and labels
#' @param dcati data.table. Pre-processed data with observations for analysis
#' @param output_file_prefix Character. Full path prefix for output files (without extension)
#' @param stan_file Character. Path to Stan model file (.stan)
#' @param x_formula Formula. Formula specifying predictors for the design matrix. The formula 
#'   is used with \code{model.matrix()} to create predictor matrices. Typically should include 
#'   \code{-1} to remove the intercept term, as the model includes baseline parameters separately. 
#'   For example, \code{~ time - 1} for a single time predictor, or 
#'   \code{~ time + factor(group) - 1} for time and group effects. When factors are included, 
#'   one level is automatically dropped for identifiability. Default: \code{~ time - 1}
#' @param chains Integer. Number of MCMC chains to run (default: 2)
#' @param parallel_chains Integer. Number of chains to run in parallel (default: 2)
#' @param threads_per_chain Integer. Number of threads to use per chain (default: 1)
#' @param iter_warmup Integer. Number of warmup iterations per chain (default: 500)
#' @param iter_sampling Integer. Number of sampling iterations per chain (default: 1500)
#' @param seed Integer. Random seed for reproducibility (default: 123)
#' @param show_messages Logical. If TRUE, show all Stan informational messages during sampling (default: FALSE)
#' @param show_exceptions Logical. If TRUE, show detailed Stan exception messages when errors occur (default: FALSE)
#' @param with_core_analyses Logical. If TRUE, generate core probability plots (default: TRUE)
#' @param with_additional_analyses Logical. If TRUE, generate additional diagnostic
#'   plots including trace plots, parameter intervals/areas, and posterior predictive
#'   checks (default: FALSE)
#'
#' @return Invisibly returns the fitted cmdstanr model object. Saves the following
#'   files to output_dir:
#'   \item{*_stan.rds}{Fitted Stan model object}
#'   \item{*_convergence_mixing.csv}{Convergence diagnostics and parameter summaries}
#'   \item{*_probs_barplot_v2.png/pdf}{Probability distribution plots}
#'   \item{*_worsttrace.png}{Trace plots for worst-mixing parameters (if with_additional_analyses=TRUE)}
#'   \item{*_intervals.png}{Parameter credible intervals (if with_additional_analyses=TRUE)}
#'   \item{*_areas.png}{Parameter posterior densities (if with_additional_analyses=TRUE)}
#'   \item{*_ppcheck.png}{Posterior predictive check plots (if with_additional_analyses=TRUE)}
#'
#' @import data.table
#' @import ggplot2
#' @import cmdstanr
#' @importFrom ggsci pal_futurama scale_color_npg
#' @importFrom bayesplot mcmc_trace mcmc_intervals mcmc_areas color_scheme_set
#' @importFrom posterior default_summary_measures default_convergence_measures quantile2
#'
#' @examples
#' \dontrun{
#' # Load data
#' dit <- readRDS("data/dit.rds")
#' dcati <- readRDS("data/dcati_processed.rds")
#'
#' # Run analysis with default settings
#' fit_ordered_logit_model(
#'     dit = dit,
#'     dcati = dcati,
#'     output_file_prefix = "output/analysis_job1",
#'     stan_file = "src/stan/ordered_logit_2cats_v251120.stan"
#' )
#'
#' # Run with additional diagnostics
#' fit_ordered_logit_model_2cats(
#'     dit = dit,
#'     dcati = dcati,
#'     output_file_prefix = "output/analysis_job1",
#'     stan_file = "test/stan/ordered_logit_2cats_v251222.stan",
#'     with_additional_analyses = TRUE
#' )
#' }
#'
#' @export
fit_ordered_logit_model_2cats <- function(
  dit,
  dcati,
  output_file_prefix,
  stan_file = here::here("src", "stan", "ordered_logit_2cats_v251216.stan"),
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
        iqr_upper <- prob <- NULL

    require(data.table)
    require(ggplot2)
    require(ggsci)
    require(bayesplot)
    require(cmdstanr)

    # Detect if using threading based on threads_per_chain
    use_threading <- threads_per_chain > 1L
    if (use_threading && !grepl("threading", stan_file)) {
        warning(
            "threads_per_chain > 1 but stan_file does not appear to be a threading version. ",
            "Consider using a Stan model with threading support."
        )
    }
    # Print configuration
    cat("\n========================================\n")
    cat("Ordered Logit Model Analysis Configuration\n")
    cat("========================================\n")
    cat("Stan file:", stan_file, "\n")
    cat("Stan include dir:", dirname(stan_file), "\n")
    cat("Data: dit with", nrow(dit), "items, dcati with", nrow(dcati), "observations\n")
    cat("Output prefix:", output_file_prefix, "\n")
    cat("Chains:", chains, "\n")
    cat("Parallel chains:", parallel_chains, "\n")
    cat("Threads per chain:", threads_per_chain, "\n")
    if (use_threading && grepl("threading", stan_file)) {
        cat("Threading enabled: using threaded Stan model\n")
    }
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
    ol_compiled <- cmdstanr::cmdstan_model(
        stan_file,
        include_paths = dirname(stan_file),
        cpp_options = list(stan_threads = TRUE)
    )

    # Define data in format needed for model specification
    cat("Preparing Stan data...\n")
    stan_data <- list()
    stan_data$U <- max(dcati$pid)
    stan_data$Ncat1 <- nrow(dcati[item_type == "categorical", ])
    stan_data$Qcat1 <- max(dcati[item_type == "categorical", item_time_id])
    stan_data$Kcat1 <- length(unique(dcati[item_type == "categorical", y_stan]))
    stan_data$cat1_y <- dcati[item_type == "categorical", y_stan]
    stan_data$cat1_question_of_obs <- dcati[
        item_type == "categorical", item_time_id
    ]
    stan_data$cat1_unit_of_obs <- dcati[item_type == "categorical", pid]
    stan_data$cat1_X <- model.matrix(
        x_formula,
        data = as.data.frame(dcati[item_type == "categorical", ])
    )
    stan_data$Ncat2 <- nrow(dcati[item_type == "out-of-7", ])
    stan_data$Qcat2 <- max(dcati[item_type == "out-of-7", item_time_id])
    stan_data$Kcat2 <- length(unique(dcati[item_type == "out-of-7", y_stan]))
    stan_data$cat2_y <- dcati[item_type == "out-of-7", y_stan]
    stan_data$cat2_question_of_obs <- dcati[
        item_type == "out-of-7", item_time_id
    ]
    stan_data$cat2_unit_of_obs <- dcati[item_type == "out-of-7", pid]
    stan_data$cat2_X <- model.matrix(
        x_formula,
        data = as.data.frame(dcati[item_type == "out-of-7", ])
    )
    
    # Set P based on actual number of columns in design matrix
    # This accounts for factors which drop one level for identifiability
    stan_data$P <- ncol(stan_data$cat1_X)
    cat("Design matrix number of predictors: P =", stan_data$P, "\n")
    cat("Design matrix column names:", paste(colnames(stan_data$cat1_X), collapse = ", "), "\n")
    
    # Add grainsize if using threading version
    if (use_threading && grepl("threading", stan_file)) {
        stan_data$grainsize <- 1
    }

    # Sample from the model
    cat("Running MCMC sampling...\n")
    flush.console()  
    ol_fit <- ol_compiled$sample(
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
    tmp <- as.data.table(ol_fit$draws(variables = "lp__", inc_warmup = FALSE, format = "draws_df"))
    tmp <- tmp[, .(
        mean_lp = mean(lp__),
        sd_lp = sd(lp__),
        n = .N
    ), by = .chain]
    
    # Test if any other chains have mean lp__ < avg - 2*SE of best chain    
    threshold <- tmp[which.max(mean_lp), mean_lp - 2 * sd_lp]
    ol_fit_good_chains <- tmp[mean_lp > threshold, .chain]
    cat("Identified good HMC chains:", paste(ol_fit_good_chains, collapse = ", "), "\n")

    # Extract chain timing information
    cat("\nExtracting timing information...\n")
    chain_times <- ol_fit$time()
    timing_data <- data.table(
        chain = ol_fit_good_chains,
        warmup_minutes = chain_times$chains$warmup[ol_fit_good_chains] / 60,
        sampling_minutes = chain_times$chains$sampling[ol_fit_good_chains] / 60,
        total_chain_minutes = chain_times$chains$total[ol_fit_good_chains] / 60
    )
    
    # Save timing information
    timing_file <- paste0(output_file_prefix, "_timing.csv")
    write.csv(timing_data, file = timing_file, row.names = FALSE)
    cat("Saved timing information to:", timing_file, "\n")

    # Extract and save draws immediately (while CSV files still exist)
    cat("Extracting draws...\n")
    draws_file <- paste0(output_file_prefix, "_draws.rds")
    draws <- ol_fit$draws(format = "draws_array")
    draws <- draws[, ol_fit_good_chains, , drop = FALSE]
    saveRDS(draws, file = draws_file)
    cat("Saved draws to:", draws_file, "\n")
    draws <- NULL
    gc()

    # Save output to RDS
    output_file <- paste0(output_file_prefix, "_stan.rds")
    cat("Saving model fit to:", output_file, "\n")
    ol_fit$save_object(file = output_file)

    # Check convergence and mixing - compute only rhat and ess for efficiency
    cat("Generating convergence diagnostics...\n")
    tmp <- ol_fit$summary(
        variables = c(
            "latent_factor_unit", "latent_factor_beta",
            "cat1_skill_thresholds_1", "cat1_skill_thresholds_incs",
            "cat1_loadings_questions_m1",
            "cat2_skill_thresholds_1", "cat2_skill_thresholds_incs",
            "cat2_loadings_questions_m1"
        )
    )
    tmp <- as.data.table(tmp)
    setorder(tmp, ess_bulk)
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
        po <- ol_fit$draws(
            variables = c("lp__", worst_var),
            inc_warmup = TRUE,
            format = "draws_array"
        )
        po <- po[, ol_fit_good_chains, , drop = FALSE]
        
        # Calculate lp__ range from post-warmup samples for y-axis scaling
        lp_range <- range(po[(iter_warmup + 1):(iter_warmup + iter_sampling), , "lp__"])
        lp_ylim <- c(lp_range[1] - 0.05 * diff(lp_range), 
                     lp_range[2] + 0.05 * diff(lp_range))
        
        # Create lp__ trace plot with constrained y-axis
        p_lp <- bayesplot:::mcmc_trace(po[,, "lp__", drop = FALSE],
            pars = "lp__",
            n_warmup = iter_warmup
        ) + 
            coord_cartesian(ylim = lp_ylim) +
            theme_bw()
        
        # Create trace plot for other parameters with free y-axis
        p_other <- bayesplot:::mcmc_trace(po[,, worst_var, drop = FALSE],
            pars = worst_var,
            n_warmup = iter_warmup,
            facet_args = list(ncol = 2)
        ) + 
            theme_bw()
        
        # Combine plots
        p <- patchwork::wrap_plots(p_lp, p_other, ncol = 1, heights = c(1, 4))
        
        ggsave(
            file = paste0(output_file_prefix, "_worsttrace.png"),
            plot = p,
            h = 10,
            w = 20
        )

        # Make intervals/areas plot
        cat("Generating parameter interval plots...\n")
        po <- ol_fit$draws(
            variables = c(
                "latent_factor_unit", "latent_factor_beta",
                "cat1_skill_thresholds_1", "cat1_skill_thresholds_incs",
                "cat1_loadings_questions_m1",
                "cat2_skill_thresholds_1", "cat2_skill_thresholds_incs",
                "cat2_loadings_questions_m1"
            ),
            inc_warmup = FALSE,
            format = "draws_array"
        )
        po <- po[, ol_fit_good_chains, , drop = FALSE]

        color_scheme_set("teal")
        p <- bayesplot::mcmc_intervals(
            po,
            prob = 0.5, prob_outer = 0.95, outer_size = 1, point_size = 2
        ) + theme_bw()
        ggsave(
            file = paste0(output_file_prefix, "_intervals.png"),
            plot = p,
            h = 70,
            w = 8,
            limitsize = FALSE
        )

        cat("Generating parameter density overlay plots...\n")
        po <- ol_fit$draws(
            variables = c(
                "latent_factor_unit", "latent_factor_beta"                
            ),
            inc_warmup = FALSE,
            format = "draws_array"
        )
        po <- po[, ol_fit_good_chains, , drop = FALSE]
        color_scheme_set("mix-teal-pink")
        p <- bayesplot::mcmc_dens_overlay(po,
                                          facet_args = list(ncol = 5)) + 
            theme_bw() +
            labs(color = "Chain")
        ggsave(
            file = paste0(output_file_prefix, "_areas_latent_factor_and_baselines.pdf"),
            plot = p,
            h = 100,
            w = 20,
            limitsize = FALSE
        )

        po <- ol_fit$draws(
            variables = c(
                "cat1_skill_thresholds_1", "cat1_skill_thresholds_incs",
                "cat1_loadings_questions_m1"
            ),
            inc_warmup = FALSE,
            format = "draws_array"
        )
        po <- po[, ol_fit_good_chains, , drop = FALSE]
        p <- bayesplot::mcmc_dens_overlay(po,
                                          facet_args = list(ncol = 5)) + 
            theme_bw() +
            labs(color = "Chain")
        ggsave(
            file = paste0(output_file_prefix, "_areas_cat1_parameters.pdf"),
            plot = p,
            h = 40,
            w = 20,
            limitsize = FALSE
        )

        po <- ol_fit$draws(
            variables = c(
                "cat2_skill_thresholds_1", "cat2_skill_thresholds_incs",
                "cat2_loadings_questions_m1"
            ),
            inc_warmup = FALSE,
            format = "draws_array"
        )
        po <- po[, ol_fit_good_chains, , drop = FALSE]
        p <- bayesplot::mcmc_dens_overlay(po,
                                          facet_args = list(ncol = 5)) + 
            theme_bw() +
            labs(color = "Chain")
        ggsave(
            file = paste0(output_file_prefix, "_areas_cat2_parameters.pdf"),
            plot = p,
            h = 60,
            w = 20,
            limitsize = FALSE
        )

        # Make pairs plots with 2D density contours colored by chain
        cat("Generating pairs plots...\n")    
        require(GGally)
        po <- ol_fit$draws(
            variables = c(
                "latent_factor_unit[10]", "latent_factor_unit[20]", "latent_factor_unit[30]", 
                "latent_factor_unit[40]", "latent_factor_unit[50]", "latent_factor_unit[60]", 
                "cat1_loadings_questions_m1[3]", "cat1_loadings_questions_m1[4]", 
                "cat1_skill_thresholds_1[3]",
                "cat1_skill_thresholds_incs[3,1]", "cat1_skill_thresholds_incs[3,2]"
            ),
            format = "draws_array"
            )
        po <- po[, ol_fit_good_chains, , drop = FALSE]
        # Convert to draws_df format for ggpairs
        po <- posterior::as_draws_df(po)
        
        # Convert to data.table and prepare for ggplot
        po <- as.data.table(po)
        po[, chain := factor(.chain)]
        
        # Select only parameter columns and chain
        param_names <- setdiff(names(po), c(".draw", ".chain", ".iteration", "chain"))
        po <- po[, c(param_names, "chain"), with = FALSE]
        
        # Create pairs plot with GGally
        p <- GGally::ggpairs(
            po,
            mapping = aes(color = chain, fill = chain),
            columns = param_names,
            lower = list(continuous = GGally::wrap("density", alpha = 0.5, contour = TRUE)),
            diag = list(continuous = GGally::wrap("densityDiag", alpha = 0.3)),
            upper = list(continuous = "blank")
        ) +
        scale_color_manual(values = c("1" = "#3B9AB2", "2" = "#F21A00")) +
        scale_fill_manual(values = c("1" = "#3B9AB2", "2" = "#F21A00")) +
        theme_bw()
        
        ggsave(file = paste0(output_file_prefix, "_pairs_plot_sampled.pdf"),
               plot = p,
               h = 20,
               w = 20,
               limitsize = FALSE
            )
        

        # Make posterior predictive check
        cat("Generating posterior predictive checks...\n")
        po <- ol_fit$draws(
            variables = c("cat1_ypred", "cat2_ypred"),
            inc_warmup = FALSE,
            format = "draws_array"
        )
        po <- po[, ol_fit_good_chains, , drop = FALSE]
        # Convert to draws_df format
        po <- posterior::as_draws_df(po)
        po <- as.data.table(po)
        po <- data.table::melt(
            po,
            id.vars = c(".draw", ".chain", ".iteration"),
            value.name = "ypred"
        )
        po[, item_type := gsub("([a-z0-9_]+)\\[([0-9]+)\\]", "\\1", variable)]
        po[, oidt := as.integer(gsub("([a-z0-9_]+)\\[([0-9]+)\\]", "\\2", variable))]
        tmp <- po[, which(grepl("cat2", item_type))]
        set(po, tmp, "item_type", "out-of-7")
        set(po, po[, which(grepl("ypred", item_type))], "item_type", "categorical")

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

        pos <- merge(pos, dcati, by = c("item_type", "oidt"))
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
        
        # Extract Brier scores for cat1
        brier_cat1 <- ol_fit$draws(
            variables = c("cat1_ordinal_brier_score"),
            inc_warmup = FALSE,
            format = "draws_array"
        )
        brier_cat1 <- brier_cat1[, ol_fit_good_chains, , drop = FALSE]
        brier_cat1 <- posterior::as_draws_df(brier_cat1)
        brier_cat1 <- as.data.table(brier_cat1)
        brier_cat1 <- data.table::melt(brier_cat1,
            id.vars = c(".draw", ".chain", ".iteration"),
            value.name = "brier_score"
        )
        brier_cat1[, question_id := as.integer(gsub(
            "cat1_ordinal_brier_score\\[([0-9]+)\\]", "\\1",
            variable
        ))]
        brier_cat1[, item_type := "categorical"]
        set(brier_cat1, NULL, "variable", NULL)
        
        # Extract Brier scores for cat2
        brier_cat2 <- ol_fit$draws(
            variables = c("cat2_ordinal_brier_score"),
            inc_warmup = FALSE,
            format = "draws_array"
        )
        brier_cat2 <- brier_cat2[, ol_fit_good_chains, , drop = FALSE]
        brier_cat2 <- posterior::as_draws_df(brier_cat2)
        brier_cat2 <- as.data.table(brier_cat2)
        brier_cat2 <- data.table::melt(brier_cat2,
            id.vars = c(".draw", ".chain", ".iteration"),
            value.name = "brier_score"
        )
        brier_cat2[, question_id := as.integer(gsub(
            "cat2_ordinal_brier_score\\[([0-9]+)\\]", "\\1",
            variable
        ))]
        brier_cat2[, item_type := "out-of-7"]
        set(brier_cat2, NULL, "variable", NULL)
        
        # Combine and compute median
        brier_all <- rbind(brier_cat1, brier_cat2)
        brier_summary <- brier_all[,
            list(
                median_brier_score = median(brier_score),
                mean_brier_score = mean(brier_score),
                q025_brier_score = quantile(brier_score, 0.025),
                q975_brier_score = quantile(brier_score, 0.975)
            ),
            by = c("item_type", "question_id")
        ]
        
        # Save to CSV
        data.table::fwrite(
            brier_summary,
            file = paste0(output_file_prefix, "_ordered_brierscore.csv")
        )
        cat("Ordinal Brier scores saved to", paste0(output_file_prefix, "_ordered_brierscore.csv"), "\n")
    } else {
        cat("\nSkipping additional diagnostic analyses (set with_additional_analyses=TRUE to enable)\n")
    }

    # Generate probability plots for cat1
    if (with_core_analyses) {
                cat("Generating probability plots for categorical outcomes...\n")
        po <- ol_fit$draws(
            variables = c("cat1_ordered_prob_by_obs"),
            inc_warmup = FALSE,
            format = "draws_array"
        )
        po <- po[, ol_fit_good_chains, , drop = FALSE]
        # Convert to draws_df format
        po <- posterior::as_draws_df(po)
        po <- as.data.table(po)
        po <- data.table::melt(po,
            id.vars = c(".draw", ".chain", ".iteration"),
            value.name = "prob"
        )
        po[, y_stan := as.integer(gsub(
            "([a-z0-9_]+)\\[([0-9]+),([0-9]+)\\]", "\\3",
            variable
        ))]
        po[, oidt := as.integer(gsub(
            "([a-z0-9_]+)\\[([0-9]+),([0-9]+)\\]", "\\2",
            variable
        ))]
        set(po, NULL, "variable", NULL)
        tmp <- unique(subset(dcati[item_type == "categorical", ],
            select = c(oidt, pid, item_time_id, time)
        ))
        po <- merge(po, tmp, by = c("oidt"))
        po <- po[,
            list(prob = mean(prob)),
            by = c(".draw", "time", "item_time_id", "y_stan")
        ]
        pos <- po[,
            list(
                summary_value = quantile(
                    prob,
                    prob = c(0.025, 0.25, 0.5, 0.75, 0.975)
                ),
                summary_name = c("q_lower", "iqr_lower", "median", "iqr_upper", "q_upper")
            ),
            by = c("time", "item_time_id", "y_stan")
        ]
        pos <- data.table::dcast(pos,
            time + item_time_id + y_stan ~ summary_name,
            value.var = "summary_value"
        )
        tmp <- unique(subset(dcati[item_type == "categorical", ],
            select = c(time, time_label)
        ))
        pos <- merge(pos, tmp, by = c("time"))
        tmp <- unique(subset(dcati[item_type == "categorical", ],
            select = c(item_time_id, item_label)
        ))
        pos <- merge(pos, tmp, by = c("item_time_id"))
        tmp <- unique(subset(dcati[item_type == "categorical", ],
            select = c(y, y_stan, y_label)
        ))
        pos <- merge(pos, tmp, by = c("y_stan"))
        tmp <- dcati[item_type == "categorical",
            list(n = length(pid)),
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
        pos_cat <- copy(pos)

        # Generate probability plots for cat2
        cat("Generating probability plots for out-of-7 outcomes...\n")
        po <- ol_fit$draws(
            variables = c("cat2_ordered_prob_by_obs"),
            inc_warmup = FALSE,
            format = "draws_array"
        )
        po <- po[, ol_fit_good_chains, , drop = FALSE]
        # Convert to draws_df format
        po <- posterior::as_draws_df(po)
        po <- as.data.table(po)
        po <- data.table::melt(po,
            id.vars = c(".draw", ".chain", ".iteration"),
            value.name = "prob"
        )
        po[, y_stan := as.integer(gsub(
            "([a-z0-9_]+)\\[([0-9]+),([0-9]+)\\]", "\\3",
            variable
        ))]
        po[, oidt := as.integer(gsub(
            "([a-z0-9_]+)\\[([0-9]+),([0-9]+)\\]", "\\2",
            variable
        ))]
        set(po, NULL, "variable", NULL)
        tmp <- unique(subset(dcati[item_type == "out-of-7", ],
            select = c(oidt, pid, item_time_id, time)
        ))
        po <- merge(po, tmp, by = c("oidt"))
        po <- po[,
            list(prob = mean(prob)),
            by = c(".draw", "time", "item_time_id", "y_stan")
        ]
        pos <- po[,
            list(
                summary_value = quantile(
                    prob,
                    prob = c(0.025, 0.25, 0.5, 0.75, 0.975)
                ),
                summary_name = c("q_lower", "iqr_lower", "median", "iqr_upper", "q_upper")
            ),
            by = c("time", "item_time_id", "y_stan")
        ]
        pos <- data.table::dcast(pos,
            time + item_time_id + y_stan ~ summary_name,
            value.var = "summary_value"
        )
        tmp <- unique(subset(dcati[item_type == "out-of-7", ],
            select = c(time, time_label)
        ))
        pos <- merge(pos, tmp, by = c("time"))
        tmp <- unique(subset(dcati[item_type == "out-of-7", ],
            select = c(item_time_id, item_label)
        ))
        pos <- merge(pos, tmp, by = c("item_time_id"))
        tmp <- unique(subset(dcati[item_type == "out-of-7", ],
            select = c(y, y_stan, y_label)
        ))
        pos <- merge(pos, tmp, by = c("y_stan"))
        tmp <- dcati[item_type == "out-of-7",
            list(n = length(pid)),
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
        pos_cat2 <- copy(pos)

        # Combine cat1 and cat2 probability data
        pos <- rbind(pos_cat, pos_cat2)
        pos <- merge(pos, dit, by = c("item_label"))
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
    } else {
        cat("\nSkipping core analyses (set with_core_analyses=TRUE to enable)\n")
    }

    invisible(ol_fit)
}
