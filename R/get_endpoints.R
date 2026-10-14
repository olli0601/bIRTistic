#' Get endpoints from PCM model draws (reads the in-model reductions rho_by_cat_qu_*)
#'
#' This function reads the per-comparison reduced effect sizes computed inside the PCM generated quantities
#' (rho_by_cat_qu_*) and summarises them with credible intervals. Baseline/Endline group levels and the
#' differences/ratios are all read from that one matrix -- there is no ordered-probability reduction in R.
#' For categorical items the group level is P(y >= threshold); for out-of-7 items it is the expected mean.
#'
#' @param dp1 A data.table containing the input data with columns including:
#'   item_type_id, item_label, item_time_id, item_type, time_label, y, y_label.
#' @param dit A data.table containing item definitions with columns including:
#'   item_type_id, item_label, cat_length, item_type, group_label, item_label_short,
#'   group_label_long, item_high_label.
#' @param draws_file Character string specifying the path to the RDS file containing
#'   posterior draws (array format with rho_by_cat_qu_pr / rho_by_cat_qu_fit parameters from the PCM fit).
#' @param categorical_threshold Integer specifying the threshold for aggregating
#'   categorical responses. Categories >= this value are summed. Default is 3,
#'   which aggregates "most" and "all of the time" categories for 5-category items.
#' @param endpoint_type Character string specifying the aggregation level:
#'   - "items": Keeps individual items (item_type_id, item_time_id preserved)
#'   - "item_groups": Aggregates across items within groups (group_label level)
#'   Default is "items".
#' @param param_name Character string specifying which rho variant to read:
#'   - "rho_by_cat_qu_pr": posterior-predictive (confounding-adjusted) reductions
#'   - "rho_by_cat_qu_fit": fitted reductions
#'   Default is "rho_by_cat_qu_pr".
#' @param effect_sizes data.table. The reduced-effect-size manifest (drho) from define_effect_sizes_*; one row
#'   per (comparison jc, reduction r), tagging each with its reduction label (rho in Baseline/Endline/ratio/diff),
#'   scope and metadata. MUST be the SAME manifest passed to the fit. If NULL, built via define_effect_sizes_colombia.
#'
#' @return A data.table containing endpoint summaries with columns:
#'   - item_type, group_label, group_label_long: Item type and groupings
#'   - item_label, item_label_short (if endpoint_type="items"): Item identifiers
#'   - item_time_id (if endpoint_type="items"): Time point identifier (categorical items only)
#'   - item_high_label: Direction indicator (higher_is_better/lower_is_better)
#'   - variable: Endpoint type (Baseline, Endline, diff, ratio, diff2, ratio2)
#'   - q_lower, iqr_lower, median, iqr_upper, q_upper: Quantile summaries (2.5%, 25%, 50%, 75%, 97.5%)
#'
#' @details
#' The function processes two types of items:
#'
#' **Categorical items**: Aggregates probabilities for categories >= categorical_threshold.
#' - endpoint_type="items": Computes for each individual item
#' - endpoint_type="item_groups": Averages probabilities across items within each group
#'
#' **Out-of-7 items**: Computes expected means (weighted sum of category probabilities).
#' - endpoint_type="items": Computes for each individual item
#' - endpoint_type="item_groups": Averages means across items within each group
#'
#' Computed measures:
#' - diff: Baseline - Endline (reduction, for categorical or lower_is_better items)
#' - ratio: 1 - Endline/Baseline (proportional reduction)
#' - diff2: Endline - Baseline (increase)
#' - ratio2: Endline/Baseline - 1 (proportional increase)
#'
#' For out-of-7 items, direction of diff/ratio respects item_high_label when endpoint_type="items".
#'
#' Source of the reductions (both endpoint types read ONLY rho_by_cat_qu_*, never ordered_prob):
#' - endpoint_type="items": Baseline/Endline and diff/ratio/diff2/ratio2 come from rho_by_cat_qu_* [r, jc],
#'   where red_fn[r] names the reduction. The legacy group-mean pos uses ids 5 (baseline level), 6 (endline
#'   level), 7 (rel %, ratio-of-means), 8 (diff-of-means); the per-participant pre/post estimand (ids 3 rel %,
#'   4 diff) is in the same matrix. param_name selects the _pr or _fit variant; red_fn must match the fit.
#' - endpoint_type="item_groups": read from the per-group comparisons the fit schedules for each (item_type,
#'   group_label). The Stan reduction pools each participant's probs across the construct's items they answered
#'   (participant-weighted), so under item-level missingness this differs from the legacy item-weighted
#'   mean-of-per-item-group-means; it is identical when every participant answered every item in the group.
#'
#' @examples
#' \dontrun{
#' # Load data
#' tmp <- readRDS("path/to/data.rds")
#' dp1 <- copy(tmp$dp1)
#' dit <- copy(tmp$dit)
#'
#' # Get detailed endpoints by individual items using posterior predictions
#' pos <- get_endpoints(
#'   dp1 = dp1,
#'   dit = dit,
#'   draws_file = "path/to/draws.rds",
#'   categorical_threshold = 3,
#'   endpoint_type = "items",
#'   param_name = "rho_by_cat_qu_pr"
#' )
#'
#' # Get grouped endpoints aggregated across items using fitted values
#' pos_grouped <- get_endpoints(
#'   dp1 = dp1,
#'   dit = dit,
#'   draws_file = "path/to/draws.rds",
#'   categorical_threshold = 3,
#'   endpoint_type = "item_groups",
#'   param_name = "rho_by_cat_qu_fit"
#' )
#' }
#'
#' @importFrom data.table copy merge setkey dcast melt set rbind data.table as.data.table setnames setorder
#' @importFrom posterior as_draws_df
#' @export
get_endpoints <- function(dp1, dit, draws_file, categorical_threshold = 3, endpoint_type = "items", param_name = "rho_by_cat_qu_pr", effect_sizes = NULL) {
    
    # Ensure data.table is available
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Package 'data.table' is required but not installed.")
    }
    
    # Validate inputs
    if (!file.exists(draws_file)) {
        stop("Draws file not found: ", draws_file)
    }
    
    if (!endpoint_type %in% c("items", "item_groups")) {
        stop("endpoint_type must be either 'items' or 'item_groups'")
    }
    
    if (!param_name %in% c("rho_by_cat_qu_pr", "rho_by_cat_qu_fit")) {
        stop("param_name must be either 'rho_by_cat_qu_pr' or 'rho_by_cat_qu_fit'")
    }

    # =========================================================================
    # Read the reductions from the Stan model (rho_by_cat_qu_*), driven by `drho`
    # =========================================================================
    # drho (define_effect_sizes_*; one row per comparison jc x reduction r) is the single source of truth: it
    # tags each (r, jc) with its reduction label (rho in {Baseline, Endline, ratio, diff}), scope (item /
    # item_group) and output metadata. We read rho_by_cat_qu_*[r, jc], join drho by (r, jc), pivot the four
    # labels into Baseline/Endline/ratio/diff, form diff2 = -diff and ratio2 = -ratio (categorical not
    # direction-aware, out-of-7 respect item_high_label), and return the rows for the requested endpoint_type.
    # Label-driven, so it is agnostic to whether ratio/diff come from the per-participant ids (3,4) or the
    # group-mean ids (7,8). ratio is emitted by Stan as a percent and converted to a fraction here.
    if (is.null(effect_sizes)) {
        effect_sizes <- define_effect_sizes_colombia(dp1, dit, categorical_threshold = categorical_threshold)
    }
    drho <- data.table::as.data.table(effect_sizes)
    need_lbl <- c("Baseline", "Endline", "ratio", "diff")
    if (!all(need_lbl %in% drho$rho)) {
        stop("effect_sizes must contain reductions labelled ", paste(need_lbl, collapse = ", "),
             "; got: ", paste(sort(unique(drho$rho)), collapse = ", "))
    }

    # --- read rho_by_cat_qu_*[r, jc] ---
    draws <- readRDS(draws_file)
    rho <- draws[,, grepl(param_name, dimnames(draws)[[3]]), drop = FALSE]
    rho <- data.table::as.data.table(posterior::as_draws_df(rho))
    rho <- data.table::melt(rho, id.vars = c(".draw", ".chain", ".iteration"), value.name = "val")
    rho[, r  := as.integer(gsub(paste0(param_name, "\\[([0-9]+),([0-9]+)\\]"), "\\1", variable))]
    rho[, jc := as.integer(gsub(paste0(param_name, "\\[([0-9]+),([0-9]+)\\]"), "\\2", variable))]
    data.table::set(rho, NULL, "variable", NULL)

    # --- join reduction labels, pivot to Baseline/Endline/ratio/diff per (draw, jc) ---
    rho <- merge(rho, unique(drho[, .(jc, r, rho)]), by = c("jc", "r"))
    w <- data.table::dcast(rho[rho %in% need_lbl], .draw + jc ~ rho, value.var = "val")

    # --- per-draw Baseline/Endline + diff/ratio/diff2/ratio2 (direction-aware) ---
    cmpmeta <- unique(drho[, .(jc, scope, item_type_id, item_type, item_label, group_label,
                               group_label_long, item_high_label, item_label_short)])
    w <- merge(w, cmpmeta[, .(jc, item_type, item_high_label)], by = "jc")
    w[, ratio_lb := ratio / 100]                                        # Stan emits percent; (B-E)/B as fraction
    w[, is_cat := item_type == "categorical"]
    w[, lb := item_high_label == "lower_is_better"]
    w[, diff_v   := data.table::fifelse(is_cat | lb,  diff, -diff)]
    w[, ratio_v  := data.table::fifelse(is_cat | lb, ratio_lb, -ratio_lb)]
    w[, diff2_v  := -diff]
    w[, ratio2_v := -ratio_lb]
    allv <- data.table::melt(
        w[, .(.draw, jc, Baseline, Endline, diff = diff_v, ratio = ratio_v, diff2 = diff2_v, ratio2 = ratio2_v)],
        id.vars = c(".draw", "jc"), variable.name = "variable", value.name = "value")
    allv[, variable := as.character(variable)]

    # --- select rows for the requested endpoint, attach metadata, summarise quantiles ---
    if (endpoint_type == "items") {
        allv <- allv[jc %in% cmpmeta[scope == "item", jc]]
        keys <- c("item_type_id", "item_label", "item_type", "group_label",
                  "item_high_label", "item_label_short", "group_label_long", "variable")
    } else {
        allv <- allv[jc %in% cmpmeta[scope == "item_group", jc]]
        keys <- c("item_type", "group_label", "group_label_long", "item_high_label", "variable")
    }
    allv <- merge(allv, cmpmeta, by = "jc")
    pos <- allv[,
        list(
            summary_value = quantile(value, prob = c(0.025, 0.25, 0.5, 0.75, 0.975)),
            summary_name  = c("q_lower", "iqr_lower", "median", "iqr_upper", "q_upper")
        ),
        by = keys]
    pos <- data.table::dcast(pos, as.formula(paste(paste(keys, collapse = " + "), "~ summary_name")),
        value.var = "summary_value")
    if (endpoint_type == "items") {
        pos[item_type != "categorical", item_type_id := NA_integer_]
    }
    return(pos)
}
