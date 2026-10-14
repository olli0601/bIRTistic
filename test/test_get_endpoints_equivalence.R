#!/usr/bin/env Rscript
# Equivalence test: the new get_endpoints (reads the Stan reductions rho_by_cat_qu_*, group-mean ids 5/6/7/8)
# must return EXACTLY the same items pos as the legacy old/get_endpoints_from_probs_by_cat_qu (which computes
# Baseline/Endline levels and ratio-/diff-of-group-means in R from ordered_prob_by_cat_qu_*). Both read the
# SAME draws of a single PCM v261012 fit on the Colombia data, so the match is per-draw and does not depend on
# MCMC convergence (few iterations suffice).

suppressMessages({
    library(data.table)
    library(cmdstanr)
    library(posterior)
    library(here)
})

root <- here::here()
source(file.path(root, "R", "fit_partial_credit_model_ncats_stanhmc.R"))
source(file.path(root, "R", "get_endpoints.R"))
source(file.path(root, "old", "get_endpoints_from_probs_by_cat_qu.R"))

CAT_THR <- 3L

# ---- load Colombia fixtures (= read_data_colombia output) ----
fx <- file.path(root, "test", "python", "fixtures", "r_baselines", "phase2")
dp_col  <- fread(file.path(fx, "colombia_dp.csv"))
dit_col <- fread(file.path(fx, "colombia_dit.csv"))

# ---- preprocess (mirror scripts-R/Colombia_analysis_for_HGpaper_pcm_v251224.Rmd) ----
dp1 <- subset(dp_col, !grepl("agg", item_label))
dp1[, y_stan := y + 1L]
tmp <- CJ(item_label = sort(unique(dp1$item_label)), time = unique(dp1$time))
tmp[, item_type := ifelse(grepl("CG-MH", item_label), "categorical", "out-of-7")]
tmp <- tmp[, list(item_label = item_label, time = time, item_time_id = seq_along(item_label)), by = "item_type"]
dp1 <- merge(dp1, tmp, by = c("item_label", "time"))
dp1 <- merge(dp1, unique(subset(dit_col, select = c(item_type, item_type_id))), by = "item_type")
setkey(dp1, item_type_id, pid, time, item_label)
dp1[, oid := .I]
setkey(dp1, oid)
tmp <- dp1[, list(oid = oid, oidt = seq_along(y_stan)), by = "item_type"]
dp1 <- merge(dp1, tmp, by = c("item_type", "oid"))
setkey(dp1, pid, time, item_label)

# ---- fit v261012 (few iters; equivalence is per-draw so convergence is irrelevant) ----
outdir <- file.path(tempdir(), "colombia_pcm_equivtest")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
prefix <- file.path(outdir, "pcm_test")
draws_file <- paste0(prefix, "_draws.rds")

if (!file.exists(draws_file)) {
    fit_partial_credit_model_ncats_stanhmc(
        dit = dit_col, dcati = dp1, output_file_prefix = prefix,
        stan_file = file.path(root, "src", "stan", "partial_credit_model_ncats_v261012.stan"),
        x_formula = ~ time - 1, categorical_threshold = CAT_THR,
        chains = 2L, parallel_chains = 2L, iter_warmup = 150L, iter_sampling = 150L, seed = 1L,
        with_core_analyses = FALSE, with_additional_analyses = FALSE,
        show_messages = FALSE, show_exceptions = TRUE
    )
}
stopifnot(file.exists(draws_file))

# ---- run both endpoint functions (endpoint_type = "items") ----
pos_new <- get_endpoints(dp1 = dp1, dit = dit_col, draws_file = draws_file,
                         categorical_threshold = CAT_THR, endpoint_type = "items")
pos_old <- get_endpoints_from_probs_by_cat_qu(dp1 = dp1, dit = dit_col, draws_file = draws_file,
                         categorical_threshold = CAT_THR, endpoint_type = "items")

# ---- compare the quantile summaries ----
keys    <- c("item_label", "item_type", "variable")
valcols <- c("q_lower", "iqr_lower", "median", "iqr_upper", "q_upper")
m <- merge(pos_new[, c(keys, valcols), with = FALSE],
           pos_old[, c(keys, valcols), with = FALSE],
           by = keys, suffixes = c(".new", ".old"))

cat(sprintf("\nrows: new=%d old=%d matched=%d\n", nrow(pos_new), nrow(pos_old), nrow(m)))
cat("variables:", paste(sort(unique(pos_new$variable)), collapse = ", "), "\n\n")

maxdiff <- 0
for (v in valcols) {
    d <- max(abs(m[[paste0(v, ".new")]] - m[[paste0(v, ".old")]]), na.rm = TRUE)
    cat(sprintf("  %-10s max|new-old| = %.3e\n", v, d))
    maxdiff <- max(maxdiff, d)
}

tol <- 1e-6
ok  <- nrow(m) == nrow(pos_new) && nrow(m) == nrow(pos_old) && is.finite(maxdiff) && maxdiff < tol
if (ok) {
    cat(sprintf("\nPASS: new get_endpoints == old within tol %.1e (max diff %.3e)\n", tol, maxdiff))
} else {
    cat(sprintf("\nFAIL: maxdiff=%.3e tol=%.1e; unmatched new=%d old=%d\n",
                maxdiff, tol, nrow(pos_new) - nrow(m), nrow(pos_old) - nrow(m)))
    m[, dmax := pmax(abs(median.new - median.old), abs(q_lower.new - q_lower.old),
                     abs(q_upper.new - q_upper.old))]
    print(head(m[order(-dmax)], 12))
    quit(status = 1L)
}
