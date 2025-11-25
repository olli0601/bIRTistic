#' Run credit model analysis on pre-processed data
#'
#' This function performs Bayesian IRT analysis using a credit model on pre-processed
#' data. It compiles the Stan model, runs MCMC sampling, generates convergence
#' diagnostics, and optionally creates detailed diagnostic plots and posterior
#' predictive checks.
#'
#' @param dit data.table. Item metadata table with item types and labels
#' @param dcati data.table. Pre-processed data with observations for analysis
#' @param output_file_prefix Character. Full path prefix for output files (without extension)
#' @param stan_file Character. Path to Stan model file (.stan)
#' @param chains Integer. Number of MCMC chains to run (default: 2)
#' @param parallel_chains Integer. Number of chains to run in parallel (default: 2)
#' @param iter_warmup Integer. Number of warmup iterations per chain (default: 500)
#' @param iter_sampling Integer. Number of sampling iterations per chain (default: 1500)
#' @param seed Integer. Random seed for reproducibility (default: 123)
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
#' credit_model_run_analysis(
#'     dit = dit,
#'     dcati = dcati,
#'     output_file_prefix = "output/analysis_job1",
#'     stan_file = "src/stan/credit_model_2cats_v251120.stan"
#' )
#'
#' # Run with additional diagnostics
#' credit_model_run_analysis(
#'     dit = dit,
#'     dcati = dcati,
#'     output_file_prefix = "output/analysis_job1",
#'     stan_file = "src/stan/credit_model_2cats_v251120.stan",
#'     with_additional_analyses = TRUE
#' )
#' }
#'
#' @export
credit_model_run_analysis <- function(
  dit,
  dcati,
  output_file_prefix,
  stan_file = here::here("src", "stan", "credit_model_2cats_v251120.stan"),
  chains = 2L,
  parallel_chains = 2L,
  iter_warmup = 500L,
  iter_sampling = 1500L,
  seed = 123L,
  with_core_analyses = TRUE,
  with_additional_analyses = FALSE
) {
    # Suppress data.table NSE warnings
    item_type <- y_stan <- item_time_id <- pid <- time <- oidt <- variable <-
        ypred <- oid <- IN_PPI <- item_label <- time_label <- y_label <- y <-
        n <- total <- p_emp <- group_label_long <- item_label_short <-
        item_label_long <- ess_bulk <- q_lower <- q_upper <- iqr_lower <-
        iqr_upper <- prob <- NULL

    require(data.table)
    require(ggplot2)
    require(ggsci)
    require(bayesplot)
    require(cmdstanr)

    # Print configuration
    cat("\n========================================\n")
    cat("Credit Model Analysis Configuration\n")
    cat("========================================\n")
    cat("Stan file:", stan_file, "\n")
    cat("Stan include dir:", dirname(stan_file), "\n")
    cat("Data: dit with", nrow(dit), "items, dcati with", nrow(dcati), "observations\n")
    cat("Output prefix:", output_file_prefix, "\n")
    cat("Chains:", chains, "\n")
    cat("Parallel chains:", parallel_chains, "\n")
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
    cm_compiled <- cmdstanr::cmdstan_model(
        stan_file,
        include_paths = dirname(stan_file)
    )

    # Define data in format needed for model specification
    cat("Preparing Stan data...\n")
    stan_data <- list()
    stan_data$P <- 1L
    stan_data$U <- max(dcati$pid)
    stan_data$Ncat1 <- nrow(dcati[item_type == "categorical", ])
    stan_data$Qcat1 <- max(dcati[item_type == "categorical", item_time_id])
    stan_data$Kcat1 <- length(unique(dcati[item_type == "categorical", y_stan]))
    stan_data$cat1_y <- dcati[item_type == "categorical", y_stan]
    stan_data$cat1_question_of_obs <- dcati[
        item_type == "categorical", item_time_id
    ]
    stan_data$cat1_unit_of_obs <- dcati[item_type == "categorical", pid]
    stan_data$cat1_X <- as.matrix(
        dcati[item_type == "categorical", time - 1L],
        ncol = 1
    )
    stan_data$Ncat2 <- nrow(dcati[item_type == "out-of-7", ])
    stan_data$Qcat2 <- max(dcati[item_type == "out-of-7", item_time_id])
    stan_data$Kcat2 <- length(unique(dcati[item_type == "out-of-7", y_stan]))
    stan_data$cat2_y <- dcati[item_type == "out-of-7", y_stan]
    stan_data$cat2_question_of_obs <- dcati[
        item_type == "out-of-7", item_time_id
    ]
    stan_data$cat2_unit_of_obs <- dcati[item_type == "out-of-7", pid]
    stan_data$cat2_X <- as.matrix(
        dcati[item_type == "out-of-7", time - 1L],
        ncol = 1
    )

    # Sample from the model
    cat("Running MCMC sampling...\n")
    cm_fit <- cm_compiled$sample(
        data = stan_data,
        seed = seed,
        chains = chains,
        parallel_chains = parallel_chains,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        refresh = 500,
        save_warmup = TRUE
    )

    # Save output to RDS

    output_file <- paste0(output_file_prefix, "_stan.rds")
    cat("Saving model fit to:", output_file, "\n")
    cm_fit$save_object(file = output_file)

    # Check convergence and mixing
    cat("Generating convergence diagnostics...\n")
    tmp <- cm_fit$summary(
        variables = c(
            "latent_factor_unit", "latent_factor_beta",
            "cat1_skill_thresholds_1", "cat1_skill_thresholds_incs",
            "cat1_loadings_questions",
            "cat2_skill_thresholds_1", "cat2_skill_thresholds_incs",
            "cat2_loadings_questions"
        ),
        posterior::default_summary_measures(),
        posterior::default_convergence_measures(),
        extra_quantiles = ~ posterior::quantile2(., probs = c(.0275, .975))
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
        po <- cm_fit$draws(
            variables = c("lp__", worst_var),
            inc_warmup = TRUE,
            format = "draws_array"
        )
        p <- bayesplot:::mcmc_trace(po,
            pars = c("lp__", worst_var),
            n_warmup = iter_warmup,
            facet_args = list(nrow = 2)
        )
        p <- p + theme_bw()
        ggsave(
            file = paste0(output_file_prefix, "_worsttrace.png"),
            plot = p,
            h = 10,
            w = 20
        )

        # Make intervals/areas plot
        cat("Generating parameter plots...\n")
        po <- cm_fit$draws(
            variables = c(
                "latent_factor_unit", "latent_factor_beta",
                "cat1_skill_thresholds_1", "cat1_skill_thresholds_incs",
                "cat1_loadings_questions",
                "cat2_skill_thresholds_1", "cat2_skill_thresholds_incs",
                "cat2_loadings_questions"
            ),
            inc_warmup = FALSE,
            format = "draws_array"
        )

        color_scheme_set("teal")
        p <- bayesplot::mcmc_intervals(
            po,
            prob = 0.5, prob_outer = 0.95, outer_size = 1, point_size = 2
        ) + theme_bw()
        ggsave(
            file = paste0(output_file_prefix, "_intervals.png"),
            plot = p,
            h = 50,
            w = 8,
            limitsize = FALSE
        )

        p <- bayesplot::mcmc_areas(
            po,
            prob = 0.5, prob_outer = 0.95, point_est = "median"
        ) + theme_bw()
        ggsave(
            file = paste0(output_file_prefix, "_areas.png"),
            plot = p,
            h = 150,
            w = 8,
            limitsize = FALSE
        )

        # Make posterior predictive check
        cat("Generating posterior predictive checks...\n")
        po <- cm_fit$draws(
            variables = c("cat1_ypred", "cat2_ypred"),
            inc_warmup = FALSE,
            format = "draws_df"
        )
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
        pos[, IN_PPI := y_stan >= q_lower & y_stan <= q_upper]
        cat("Proportion in 95% PPI:", pos[, mean(IN_PPI)], "\n")

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
            geom_point(aes(y = y_stan, colour = IN_PPI)) +
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
    } else {
        cat("\nSkipping additional diagnostic analyses (set with_additional_analyses=TRUE to enable)\n")
    }

    # Generate probability plots for cat1
    if (with_core_analyses) {
        cat("Generating probability plots for categorical outcomes...\n")
        po <- cm_fit$draws(
            variables = c("cat1_ordered_prob_by_obs"),
            inc_warmup = FALSE,
            format = "draws_df"
        )
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
        po <- cm_fit$draws(
            variables = c("cat2_ordered_prob_by_obs"),
            inc_warmup = FALSE,
            format = "draws_df"
        )
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

    invisible(cm_fit)
}
