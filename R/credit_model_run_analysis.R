#!/usr/bin/env Rscript

# Self-contained script to run analysis on arbitrary datasets
# Designed for execution on HPC/cloud systems
#
# Usage:
#   Rscript run_interim_analysis.R \
#     --stan_file <path_to_stan_model> \
#     --stan_include_dir <path_to_stan_includes> \
#     --dit_file <path_to_dit_rds> \
#     --dcati_file <path_to_dcati_rds> \
#     --job_id <integer> \
#     --output_dir <path_to_output> \
#     --output_prefix <prefix_string> \
#     --chains <integer> \
#     --parallel_chains <integer> \
#     --iter_warmup <integer> \
#     --iter_sampling <integer>

# Load required packages
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggsci)
    library(bayesplot)
    library(cmdstanr)
    library(argparse)
})

# Parse command line arguments using argparse
parser <- ArgumentParser(description = "Run analysis on arbitrary datasets")

# Required arguments
parser$add_argument("--stan_file",
    required = TRUE,
    help = "Path to Stan model file"
)
parser$add_argument("--stan_include_dir",
    required = TRUE,
    help = "Directory containing Stan include files"
)
parser$add_argument("--dit_file",
    required = TRUE,
    help = "Path to dit RDS file"
)
parser$add_argument("--dcati_file",
    required = TRUE,
    help = "Path to dcati RDS file (pre-processed data)"
)
parser$add_argument("--job_id",
    type = "integer", required = TRUE,
    help = "Integer ID for this analysis job"
)
parser$add_argument("--output_dir",
    required = TRUE,
    help = "Directory for output files"
)
parser$add_argument("--output_prefix",
    required = TRUE,
    help = "Prefix for output filenames"
)

# Optional arguments with defaults
parser$add_argument("--chains",
    type = "integer", default = 2L,
    help = "Number of MCMC chains (default: 2)"
)
parser$add_argument("--parallel_chains",
    type = "integer", default = 2L,
    help = "Number of parallel chains (default: 2)"
)
parser$add_argument("--iter_warmup",
    type = "integer", default = 500L,
    help = "Warmup iterations (default: 500)"
)
parser$add_argument("--iter_sampling",
    type = "integer", default = 1500L,
    help = "Sampling iterations (default: 1500)"
)
parser$add_argument("--seed",
    type = "integer", default = 123L,
    help = "Random seed for MCMC (default: 123)"
)
parser$add_argument("--with_additional_analyses",
    action = "store_true",
    default = FALSE,
    help = "Include additional diagnostic plots (trace, intervals, areas, ppcheck)"
)

args <- parser$parse_args()

# Extract arguments
stan_file <- args$stan_file
stan_include_dir <- args$stan_include_dir
dit_file <- args$dit_file
dcati_file <- args$dcati_file
job_id <- args$job_id
output_dir <- args$output_dir
output_prefix <- args$output_prefix
chains <- args$chains
parallel_chains <- args$parallel_chains
iter_warmup <- args$iter_warmup
iter_sampling <- args$iter_sampling
seed <- args$seed
with_additional_analyses <- args$with_additional_analyses

# Print configuration
cat("\n========================================\n")
cat("Analysis Configuration\n")
cat("========================================\n")
cat("Job ID:", job_id, "\n")
cat("Stan file:", stan_file, "\n")
cat("Stan include dir:", stan_include_dir, "\n")
cat("DIT file:", dit_file, "\n")
cat("Dcati file:", dcati_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Output prefix:", output_prefix, "\n")
cat("Chains:", chains, "\n")
cat("Parallel chains:", parallel_chains, "\n")
cat("Warmup iterations:", iter_warmup, "\n")
cat("Sampling iterations:", iter_sampling, "\n")
cat("Seed:", seed, "\n")
cat("Additional analyses:", with_additional_analyses, "\n")
cat("========================================\n\n")

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
cat("Loading data...\n")
dit <- readRDS(dit_file)
dcati <- readRDS(dcati_file)
cat("Loaded dcati with", nrow(dcati), "observations\n")

# Compile Stan model
cat("Compiling Stan model...\n")
ugpcm_m1_compiled <- cmdstanr::cmdstan_model(
    stan_file,
    include_paths = stan_include_dir
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
upgcm_m1_all2_fit <- ugpcm_m1_compiled$sample(
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
output_file <- file.path(
    output_dir,
    paste0(output_prefix, "_", interim_id, "_stan.rds")
)
cat("Saving model fit to:", output_file, "\n")
upgcm_m1_all2_fit$save_object(file = output_file)

# Check convergence and mixing
cat("Generating convergence diagnostics...\n")
tmp <- upgcm_m1_all2_fit$summary(
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
    file = file.path(
        output_dir,
        paste0(output_prefix, "_", job_id, "_convergence_mixing.csv")
    ),
    row.names = FALSE
)

# Additional diagnostic analyses (optional)
if (with_additional_analyses) {
    cat("\nRunning additional diagnostic analyses...\n")

    # Worst parameters with lowest ess_bulk
    worst_var <- tmp$variable[1:9]

    # Make worst trace plot
    cat("Generating trace plots...\n")
    po <- upgcm_m1_all2_fit$draws(
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
        file = file.path(
            output_dir,
            paste0(output_prefix, "_", job_id, "_worsttrace.png")
        ),
        plot = p,
        h = 10,
        w = 20
    )

    # Make intervals/areas plot
    cat("Generating parameter plots...\n")
    po <- upgcm_m1_all2_fit$draws(
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
        file = file.path(
            output_dir,
            paste0(output_prefix, "_", job_id, "_intervals.png")
        ),
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
        file = file.path(
            output_dir,
            paste0(output_prefix, "_", job_id, "_areas.png")
        ),
        plot = p,
        h = 150,
        w = 8,
        limitsize = FALSE
    )

    # Make posterior predictive check
    cat("Generating posterior predictive checks...\n")
    po <- upgcm_m1_all2_fit$draws(
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
        file = file.path(
            output_dir,
            paste0(output_prefix, "_", job_id, "_ppcheck.png")
        ),
        p,
        w = 20,
        h = 20
    )
} else {
    cat("\nSkipping additional diagnostic analyses (set --with_additional_analyses to enable)\n")
}

# Generate probability plots for cat1
cat("Generating probability plots for categorical outcomes...\n")
po <- upgcm_m1_all2_fit$draws(
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
po <- upgcm_m1_all2_fit$draws(
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

# Combine and plot
pos <- rbind(pos_cat, pos_cat2)
pos <- merge(pos, dit, by = c("item_label"))
pos[, item_label_long := paste0(group_label_long, " --- ", item_label_short)]

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
    file = file.path(
        output_dir,
        paste0(output_prefix, "_", job_id, "_probs_barplot_v2.png")
    ),
    plot = p,
    h = 40,
    w = 12
)
ggsave(
    file = file.path(
        output_dir,
        paste0(output_prefix, "_", job_id, "_probs_barplot_v2.pdf")
    ),
    plot = p,
    h = 40,
    w = 12
)

cat("\n========================================\n")
cat("Analysis complete for job", job_id, "!\n")
cat("Results saved to:", output_dir, "\n")
cat("========================================\n")
