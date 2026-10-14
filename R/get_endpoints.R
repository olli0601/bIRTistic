#' Get endpoints from ordered logit model draws
#'
#' This function processes posterior draws from an ordered logit model to compute
#' endpoint statistics for both categorical and out-of-7 response items.
#' For categorical items, it aggregates probabilities above a threshold. For out-of-7
#' items, it computes means. Results include baseline, endline, differences, and ratios
#' with credible intervals.
#'
#' @param dp1 A data.table containing the input data with columns including:
#'   item_type_id, item_label, item_time_id, item_type, time_label, y, y_label.
#' @param dit A data.table containing item definitions with columns including:
#'   item_type_id, item_label, cat_length, item_type, group_label, item_label_short,
#'   group_label_long, item_high_label.
#' @param draws_file Character string specifying the path to the RDS file containing
#'   posterior draws (array format with ordered_prob_by_cat_qu_pr or ordered_prob_by_cat_qu_fit parameters).
#' @param categorical_threshold Integer specifying the threshold for aggregating
#'   categorical responses. Categories >= this value are summed. Default is 3,
#'   which aggregates "most" and "all of the time" categories for 5-category items.
#' @param endpoint_type Character string specifying the aggregation level:
#'   - "items": Keeps individual items (item_type_id, item_time_id preserved)
#'   - "item_groups": Aggregates across items within groups (group_label level)
#'   Default is "items".
#' @param param_name Character string specifying which ordered probability parameters to use:
#'   - "ordered_prob_by_cat_qu_pr": Posterior predictions
#'   - "ordered_prob_by_cat_qu_fit": Fitted values
#'   Default is "ordered_prob_by_cat_qu_pr".
#' @param red_fn Integer vector. The reduction-id vector used in the fit (maps rho_by_cat_qu_* slot r ->
#'   reduction id; see reduction_functions.stan). Only used for endpoint_type="items"; must contain the
#'   group-mean ids 5,6,7,8. Default c(5,6,7,8) matches fit_partial_credit_model_ncats_stanhmc.
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
#' Source of the reductions:
#' - endpoint_type="items": Baseline/Endline and diff/ratio/diff2/ratio2 are read from the single PCM Stan
#'   output rho_by_cat_qu_* [r, jc], where red_fn[r] names the reduction. The legacy group-mean pos uses ids
#'   5 (baseline level), 6 (endline level), 7 (rel %, ratio-of-means), 8 (diff-of-means); the per-participant
#'   pre/post estimand (ids 3 rel %, 4 diff) is in the same matrix. param_name selects the _pr or _fit variant;
#'   red_fn must match the vector used in the fit.
#' - endpoint_type="item_groups": still computed in R from ordered_prob_by_cat_qu_* (cross-item aggregation is
#'   not done inside the model).
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
#'   param_name = "ordered_prob_by_cat_qu_pr"
#' )
#'
#' # Get grouped endpoints aggregated across items using fitted values
#' pos_grouped <- get_endpoints(
#'   dp1 = dp1,
#'   dit = dit,
#'   draws_file = "path/to/draws.rds",
#'   categorical_threshold = 3,
#'   endpoint_type = "item_groups",
#'   param_name = "ordered_prob_by_cat_qu_fit"
#' )
#' }
#'
#' @importFrom data.table copy merge setkey dcast melt set rbind data.table as.data.table setnames setorder
#' @importFrom posterior as_draws_df
#' @export
get_endpoints <- function(dp1, dit, draws_file, categorical_threshold = 3, endpoint_type = "items", param_name = "ordered_prob_by_cat_qu_pr", red_fn = c(5L, 6L, 7L, 8L)) {
    
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
    
    if (!param_name %in% c("ordered_prob_by_cat_qu_pr", "ordered_prob_by_cat_qu_fit")) {
        stop("param_name must be either 'ordered_prob_by_cat_qu_pr' or 'ordered_prob_by_cat_qu_fit'")
    }

    # =========================================================================
    # endpoint_type = "items": read ALL reductions straight from the Stan model
    # =========================================================================
    # The PCM generated quantities compute every reduction inside the model and emit one matrix
    # rho_by_cat_qu_* [r, jc] (r = reduction slot, jc = item paired Baseline->Endline). red_fn[r] names the
    # reduction (see reduction_functions.stan). To reproduce the legacy ratio-/diff-of-group-means pos we read
    # the GROUP-MEAN ids 5 (baseline level), 6 (endline level), 7 (rel %), 8 (diff). diff2 = -diff, ratio2 =
    # -ratio; categorical items are not direction-aware, out-of-7 respect item_high_label. The per-participant
    # ids (3 rel %, 4 diff) are in the same matrix for the correct pre/post estimand -- select them via red_fn.
    # (endpoint_type = "item_groups" keeps the legacy ordered_prob path below -- no cross-item aggregation in Stan.)
    if (endpoint_type == "items") {
        need <- c(5L, 6L, 7L, 8L)
        if (!all(need %in% red_fn)) {
            stop("red_fn must contain the group-mean ids 5,6,7,8 to build the items pos; got: ",
                 paste(red_fn, collapse = ","))
        }
        rho_param <- sub("ordered_prob_by_cat_qu", "rho_by_cat_qu", param_name)
        draws <- readRDS(draws_file)

        # comparison index jc in the SAME order fit_partial_credit_model_ncats_stanhmc built question_cmp:
        # items present at BOTH Baseline and Endline, keyed by (item_type_id, item_label).
        cmp_base <- unique(dp1[time_label == "Baseline", .(item_type_id, item_label, item_type, item_time_id)])
        cmp_end  <- unique(dp1[time_label == "Endline",  .(item_type_id, item_label)])
        cmp <- merge(cmp_base, cmp_end, by = c("item_type_id", "item_label"))
        data.table::setkey(cmp, item_type_id, item_label)
        cmp[, jc := seq_len(.N)]
        cmp <- merge(cmp,
            unique(subset(dit, select = c(item_type_id, item_label, item_high_label,
                                          item_label_short, group_label, group_label_long))),
            by = c("item_type_id", "item_label"))

        # --- read rho_by_cat_qu_*[r, jc], map slot r -> reduction id via red_fn, pivot id -> column ---
        rho <- draws[,, grepl(rho_param, dimnames(draws)[[3]]), drop = FALSE]
        rho <- data.table::as.data.table(posterior::as_draws_df(rho))
        rho <- data.table::melt(rho, id.vars = c(".draw", ".chain", ".iteration"), value.name = "val")
        rho[, r  := as.integer(gsub(paste0(rho_param, "\\[([0-9]+),([0-9]+)\\]"), "\\1", variable))]
        rho[, jc := as.integer(gsub(paste0(rho_param, "\\[([0-9]+),([0-9]+)\\]"), "\\2", variable))]
        data.table::set(rho, NULL, "variable", NULL)
        rho[, red := red_fn[r]]
        w <- data.table::dcast(rho[red %in% need], .draw + jc ~ red, value.var = "val")
        data.table::setnames(w, c("5", "6", "7", "8"), c("base", "end", "grp_relpct", "grp_diff"))

        # --- Baseline / Endline group levels (ids 5, 6) ---
        gm <- data.table::melt(w[, .(.draw, jc, Baseline = base, Endline = end)],
            id.vars = c(".draw", "jc"), variable.name = "variable", value.name = "value")
        gm[, variable := as.character(variable)]

        # --- diff / ratio / diff2 / ratio2 from the group-mean reductions (ids 7, 8) ---
        gr <- merge(w[, .(.draw, jc, grp_relpct, grp_diff)], cmp[, .(jc, item_type, item_high_label)], by = "jc")
        gr[, ratio_lb := grp_relpct / 100]                             # (Baseline-Endline)/Baseline as a fraction
        gr[, is_cat := item_type == "categorical"]
        gr[, lb := item_high_label == "lower_is_better"]
        gr[, diff   := data.table::fifelse(is_cat | lb,  grp_diff, -grp_diff)]
        gr[, ratio  := data.table::fifelse(is_cat | lb, ratio_lb, -ratio_lb)]
        gr[, diff2  := -grp_diff]
        gr[, ratio2 := -ratio_lb]
        gr <- data.table::melt(gr, id.vars = c(".draw", "jc"),
            measure.vars = c("diff", "ratio", "diff2", "ratio2"))
        gr[, variable := as.character(variable)]

        # --- combine, attach item metadata, summarise quantiles ---
        allv <- rbind(gm, gr, use.names = TRUE)
        allv <- merge(allv,
            cmp[, .(jc, item_type_id, item_time_id, item_label, item_type, group_label,
                    item_high_label, item_label_short, group_label_long)],
            by = "jc")
        pos <- allv[,
            list(
                summary_value = quantile(value, prob = c(0.025, 0.25, 0.5, 0.75, 0.975)),
                summary_name  = c("q_lower", "iqr_lower", "median", "iqr_upper", "q_upper")
            ),
            by = c("item_type_id", "item_time_id", "item_label", "item_type", "group_label",
                   "item_high_label", "item_label_short", "group_label_long", "variable")
        ]
        pos <- data.table::dcast(pos,
            item_type_id + item_time_id + item_label + item_type + group_label +
                item_high_label + item_label_short + group_label_long + variable ~ summary_name,
            value.var = "summary_value")
        # mirror the legacy schema: out-of-7 rows carry no item_type_id / item_time_id
        pos[item_type != "categorical", c("item_type_id", "item_time_id") := NA_integer_]
        return(pos)
    }

    # =========================================================================
    # Part 1: Categorical items - aggregate probabilities for high categories
    # =========================================================================
    
    # Extract ordered probability parameters
    po <- readRDS(draws_file)
    po <- po[,, grepl(param_name, dimnames(po)[[3]]), drop = FALSE]
    po <- posterior::as_draws_df(po)
    po <- data.table::as.data.table(po)
    po <- data.table::melt(po,
        id.vars = c(".draw", ".chain", ".iteration"),
        value.name = "prob"
    )
    
    po[, cq_id := as.integer(gsub(
        paste0(param_name, "\\[([0-9]+)\\]"), "\\1",
        variable
    ))]
    data.table::set(po, NULL, "variable", NULL)
    
    # Map cq_id back to item structure
    tmp <- unique(subset(dp1, select = c(item_type_id, item_label, item_time_id)))
    tmp <- merge(tmp, subset(dit, select = c(item_type_id, item_label, cat_length)), 
                 by = c("item_type_id", "item_label"))
    data.table::setkey(tmp, item_type_id, item_time_id)
    tmp[, cq_id := cumsum(cat_length) - cat_length + 1L]
    tmp <- tmp[, list(cq_id = cq_id - 1L + 1:cat_length, y = 1:cat_length - 1L), 
               by = c("item_type_id", "item_time_id")]
    po <- merge(tmp, po, by = "cq_id")
    
    # Filter for categorical items and aggregate high categories
    tmp <- unique(subset(dp1[item_type == "categorical", ],
        select = c(item_type_id, item_label, item_time_id, item_type, time_label)
    ))
    
    po <- merge(po, tmp, by = c("item_type_id", "item_time_id"))
    
    po <- subset(po, y >= categorical_threshold)
    
    po <- po[,
        list(y_stan = 45, prob = sum(prob)),
        by = c(".draw", "item_type_id", "item_time_id", "item_label", "item_type", "time_label")
    ]
    
    po <- merge(po, subset(dit, select = c(item_type, item_label, group_label)), 
                by = c("item_type", "item_label"))
    
    # Aggregate across items if endpoint_type="item_groups"
    if (endpoint_type == "item_groups") {
        po <- po[,
            list(prob = mean(prob)),
            by = c(".draw", "item_type", "time_label", "group_label")
        ]
        # Reshape to wide format (Baseline vs Endline)
        po <- data.table::dcast(po,
            .draw + item_type + group_label ~ time_label,
            value.var = "prob"
        )
    } else {
        # Reshape to wide format (Baseline vs Endline)
        # Note: item_time_id is removed from formula since Baseline and Endline have different item_time_ids
        # But we save it first to add back later
        item_time_mapping <- unique(po[time_label == "Baseline", .(item_type_id, item_label, item_time_id)])
        po <- data.table::dcast(po,
            .draw + item_type_id + item_label + item_type + group_label ~ time_label,
            value.var = "prob"
        )
        # Add back item_time_id from Baseline measurements
        po <- merge(po, item_time_mapping, by = c("item_type_id", "item_label"), all.x = TRUE)
    }
    
    po <- subset(po, !is.na(Baseline) & !is.na(Endline))
    
    # Compute differences and ratios
    po[, diff := Baseline - Endline]
    po[, ratio := 1 - Endline / Baseline]
    po[, diff2 := Endline - Baseline]
    po[, ratio2 := Endline / Baseline - 1]
    
    # Reshape to long format for summarization
    if (endpoint_type == "item_groups") {
        po <- data.table::melt(po,
            id.vars = c(".draw", "item_type", "group_label"),
            measure.vars = c("diff", "ratio", "diff2", "ratio2", "Baseline", "Endline")
        )
        # Compute quantile summaries
        pos <- po[,
            list(
                summary_value = quantile(value, prob = c(0.025, 0.25, 0.5, 0.75, 0.975)),
                summary_name = c("q_lower", "iqr_lower", "median", "iqr_upper", "q_upper")
            ),
            by = c("item_type", "group_label", "variable")
        ]
        pos <- data.table::dcast(pos,
            item_type + group_label + variable ~ summary_name,
            value.var = "summary_value"
        )
        # Merge with item metadata
        tmp <- unique(subset(dit, select = c(item_type, group_label, group_label_long, item_high_label)))
        pos <- merge(pos, tmp, by = c("item_type", "group_label"))
    } else {
        po <- data.table::melt(po,
            id.vars = c(".draw", "item_type_id", "item_time_id", "item_label", "item_type", "group_label"),
            measure.vars = c("diff", "ratio", "diff2", "ratio2", "Baseline", "Endline")
        )
        # Compute quantile summaries
        pos <- po[,
            list(
                summary_value = quantile(value, prob = c(0.025, 0.25, 0.5, 0.75, 0.975)),
                summary_name = c("q_lower", "iqr_lower", "median", "iqr_upper", "q_upper")
            ),
            by = c("item_type_id", "item_time_id", "item_label", "item_type", "group_label", "variable")
        ]
        pos <- data.table::dcast(pos,
            item_type_id + item_time_id + item_label + item_type + group_label + variable ~ summary_name,
            value.var = "summary_value"
        )
        # Merge with item metadata
        tmp <- unique(subset(dit, select = c(item_type, item_label, item_label_short, 
                                              group_label, group_label_long, item_high_label)))
        pos <- merge(pos, tmp, by = c("item_type", "item_label"))
    }
    pos_cat1 <- data.table::copy(pos)
    
    # Clean up
    po <- NULL
    gc()
    
    # =========================================================================
    # Part 2: Out-of-7 items - compute means
    # =========================================================================
    
    # Extract ordered probability parameters again
    po <- readRDS(draws_file)
    po <- po[,, grepl(param_name, dimnames(po)[[3]]), drop = FALSE]
    po <- posterior::as_draws_df(po)
    po <- data.table::as.data.table(po)
    po <- data.table::melt(po,
        id.vars = c(".draw", ".chain", ".iteration"),
        value.name = "prob"
    )
    po[, cq_id := as.integer(gsub(
        paste0(param_name, "\\[([0-9]+)\\]"), "\\1",
        variable
    ))]
    data.table::set(po, NULL, "variable", NULL)
    
    # Map cq_id back to item structure
    tmp <- unique(subset(dp1, select = c(item_type_id, item_label, item_time_id)))
    tmp <- merge(tmp, subset(dit, select = c(item_type_id, item_label, cat_length)), 
                 by = c("item_type_id", "item_label"))
    data.table::setkey(tmp, item_type_id, item_time_id)
    tmp[, cq_id := cumsum(cat_length) - cat_length + 1L]
    tmp <- tmp[, list(cq_id = cq_id - 1L + 1:cat_length, y = 1:cat_length - 1L), 
               by = c("item_type_id", "item_time_id")]
    po <- merge(tmp, po, by = "cq_id")
    
    # Filter for out-of-7 items and compute means
    tmp <- unique(subset(dp1[item_type == "out-of-7", ],
        select = c(item_type_id, item_label, item_time_id, item_type, time_label)
    ))
    po <- merge(po, tmp, by = c("item_type_id", "item_time_id"))
    po <- po[,
        list(mean = sum(y * prob)),
        by = c(".draw", "item_label", "item_time_id", "item_type", "time_label")
    ]
    po <- merge(po, subset(dit, select = c(item_type, item_label, group_label)), 
                by = c("item_type", "item_label"))
    
    # Aggregate across items if endpoint_type="item_groups"
    if (endpoint_type == "item_groups") {
        po <- po[,
            list(mean = mean(mean)),
            by = c(".draw", "item_type", "time_label", "group_label")
        ]
        # Reshape to wide format (Baseline vs Endline)
        po <- data.table::dcast(po,
            .draw + item_type + group_label ~ time_label,
            value.var = "mean"
        )
        po <- subset(po, !is.na(Baseline) & !is.na(Endline))
        po <- merge(po,
            unique(subset(dit, select = c(item_type, group_label, item_high_label))),
            by = c("item_type", "group_label")
        )
        # Compute differences and ratios (direction depends on item_high_label)
        po[, diff := NA_real_]
        po[, ratio := NA_real_]
        tmp <- po[, which(item_high_label == "lower_is_better")]
        data.table::set(po, tmp, "diff", po[tmp, Baseline - Endline])
        data.table::set(po, tmp, "ratio", po[tmp, 1 - Endline / Baseline])
        tmp <- po[, which(item_high_label == "higher_is_better")]
        data.table::set(po, tmp, "diff", po[tmp, Endline - Baseline])
        data.table::set(po, tmp, "ratio", po[tmp, Endline / Baseline - 1])
        po[, diff2 := Endline - Baseline]
        po[, ratio2 := Endline / Baseline - 1]
        # Reshape to long format for summarization
        po <- data.table::melt(po,
            id.vars = c(".draw", "item_type", "group_label", "item_high_label"),
            measure.vars = c("diff", "ratio", "diff2", "ratio2", "Baseline", "Endline")
        )
        # Compute quantile summaries
        pos <- po[,
            list(
                summary_value = quantile(value, prob = c(0.025, 0.25, 0.5, 0.75, 0.975)),
                summary_name = c("q_lower", "iqr_lower", "median", "iqr_upper", "q_upper")
            ),
            by = c("item_type", "group_label", "item_high_label", "variable")
        ]
        pos <- data.table::dcast(pos,
            item_type + group_label + item_high_label + variable ~ summary_name,
            value.var = "summary_value"
        )
        # Merge with item metadata
        tmp <- unique(subset(dit, select = c(item_type, group_label, group_label_long, item_high_label)))
        pos <- merge(pos, tmp, by = c("item_type", "group_label", "item_high_label"))
    } else {
        # Reshape to wide format (Baseline vs Endline)
        # Note: item_time_id is removed from formula since Baseline and Endline have different item_time_ids
        # But we save it first to add back later
        item_time_mapping <- unique(po[time_label == "Baseline", .(item_label, item_time_id)])
        po <- data.table::dcast(po,
            .draw + item_label + item_type + group_label ~ time_label,
            value.var = "mean"
        )
        # Add back item_time_id from Baseline measurements
        po <- merge(po, item_time_mapping, by = "item_label", all.x = TRUE)
        po <- subset(po, !is.na(Baseline) & !is.na(Endline))
        po <- merge(po,
            unique(subset(dit, select = c(item_type, group_label, item_high_label))),
            by = c("item_type", "group_label")
        )
        # Compute differences and ratios (direction depends on item_high_label)
        po[, diff := NA_real_]
        po[, ratio := NA_real_]
        tmp <- po[, which(item_high_label == "lower_is_better")]
        data.table::set(po, tmp, "diff", po[tmp, Baseline - Endline])
        data.table::set(po, tmp, "ratio", po[tmp, 1 - Endline / Baseline])
        tmp <- po[, which(item_high_label == "higher_is_better")]
        data.table::set(po, tmp, "diff", po[tmp, Endline - Baseline])
        data.table::set(po, tmp, "ratio", po[tmp, Endline / Baseline - 1])
        po[, diff2 := Endline - Baseline]
        po[, ratio2 := Endline / Baseline - 1]
        # Reshape to long format for summarization
        po <- data.table::melt(po,
            id.vars = c(".draw", "item_label", "item_type", "group_label", "item_high_label"),
            measure.vars = c("diff", "ratio", "diff2", "ratio2", "Baseline", "Endline")
        )
        # Compute quantile summaries
        pos <- po[,
            list(
                summary_value = quantile(value, prob = c(0.025, 0.25, 0.5, 0.75, 0.975)),
                summary_name = c("q_lower", "iqr_lower", "median", "iqr_upper", "q_upper")
            ),
            by = c("item_label", "item_type", "group_label", "item_high_label", "variable")
        ]
        pos <- data.table::dcast(pos,
            item_label + item_type + group_label + item_high_label + variable ~ summary_name,
            value.var = "summary_value"
        )
        # Merge with item metadata
        tmp <- unique(subset(dit, select = c(item_type, item_label, item_label_short, 
                                              group_label, group_label_long, item_high_label)))
        pos <- merge(pos, tmp, by = c("item_type", "item_label", "item_high_label"))
    }
    pos_cat2 <- data.table::copy(pos)
    
    # Clean up
    po <- NULL
    gc()
    
    # =========================================================================
    # Combine both parts
    # =========================================================================
    
    pos <- rbind(pos_cat1, pos_cat2, fill = TRUE)

    return(pos)
}
