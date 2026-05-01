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
#'   posterior draws (array format with ordered_prob_by_cat_qu_pr parameters).
#' @param categorical_threshold Integer specifying the threshold for aggregating
#'   categorical responses. Categories >= this value are summed. Default is 3,
#'   which aggregates "most" and "all of the time" categories for 5-category items.
#' @param endpoint_type Character string specifying the aggregation level:
#'   - "items": Keeps individual items (item_type_id, item_time_id preserved)
#'   - "item_groups": Aggregates across items within groups (group_label level)
#'   Default is "items".
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
#' @examples
#' \dontrun{
#' # Load data
#' tmp <- readRDS("path/to/data.rds")
#' dp1 <- copy(tmp$dp1)
#' dit <- copy(tmp$dit)
#'
#' # Get detailed endpoints by individual items
#' pos <- get_endpoints(
#'   dp1 = dp1,
#'   dit = dit,
#'   draws_file = "path/to/draws.rds",
#'   categorical_threshold = 3,
#'   endpoint_type = "items"
#' )
#'
#' # Get grouped endpoints aggregated across items
#' pos_grouped <- get_endpoints(
#'   dp1 = dp1,
#'   dit = dit,
#'   draws_file = "path/to/draws.rds",
#'   categorical_threshold = 3,
#'   endpoint_type = "item_groups"
#' )
#' }
#'
#' @importFrom data.table copy merge setkey dcast melt set rbind
#' @importFrom posterior as_draws_df
#' @export
get_endpoints <- function(dp1, dit, draws_file, categorical_threshold = 3, endpoint_type = "items") {
    
    # Validate inputs
    if (!file.exists(draws_file)) {
        stop("Draws file not found: ", draws_file)
    }
    
    if (!endpoint_type %in% c("items", "item_groups")) {
        stop("endpoint_type must be either 'items' or 'item_groups'")
    }
    
    # Load draws
    poa <- readRDS(draws_file)
    
    # =========================================================================
    # Part 1: Categorical items - aggregate probabilities for high categories
    # =========================================================================
    
    # Extract ordered probability parameters
    po <- poa[,, grepl("ordered_prob_by_cat_qu_pr", dimnames(poa)[[3]]), drop = FALSE]
    po <- posterior::as_draws_df(po)
    po <- data.table::as.data.table(po)
    po <- data.table::melt(po,
        id.vars = c(".draw", ".chain", ".iteration"),
        value.name = "prob"
    )
    po[, cq_id := as.integer(gsub(
        "ordered_prob_by_cat_qu_pr\\[([0-9]+)\\]", "\\1",
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
        po <- data.table::dcast(po,
            .draw + item_type_id + item_time_id + item_label + item_type + group_label ~ time_label,
            value.var = "prob"
        )
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
    po <- poa[,, grepl("ordered_prob_by_cat_qu_pr", dimnames(poa)[[3]]), drop = FALSE]
    po <- posterior::as_draws_df(po)
    po <- data.table::as.data.table(po)
    po <- data.table::melt(po,
        id.vars = c(".draw", ".chain", ".iteration"),
        value.name = "prob"
    )
    po[, cq_id := as.integer(gsub(
        "ordered_prob_by_cat_qu_pr\\[([0-9]+)\\]", "\\1",
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
        po <- data.table::dcast(po,
            .draw + item_label + item_time_id + item_type + group_label ~ time_label,
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
            id.vars = c(".draw", "item_label", "item_time_id", "item_type", "group_label", "item_high_label"),
            measure.vars = c("diff", "ratio", "diff2", "ratio2", "Baseline", "Endline")
        )
        # Compute quantile summaries
        pos <- po[,
            list(
                summary_value = quantile(value, prob = c(0.025, 0.25, 0.5, 0.75, 0.975)),
                summary_name = c("q_lower", "iqr_lower", "median", "iqr_upper", "q_upper")
            ),
            by = c("item_label", "item_time_id", "item_type", "group_label", "item_high_label", "variable")
        ]
        pos <- data.table::dcast(pos,
            item_label + item_time_id + item_type + group_label + item_high_label + variable ~ summary_name,
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
    
    pos <- data.table::rbind(pos_cat1, pos_cat2, fill = TRUE)
    
    return(pos)
}
