#' Define the effect-size manifest (drho) that schedules the in-model reductions
#'
#' `drho` is a LONG data.table with one row per emitted rho = (comparison jc, reduction r). It is the single
#' source of truth shared by fit_partial_credit_model_ncats_stanhmc (builds the Stan data entries R, Qcmp,
#' Scmp, cmp_cat_type, question_cmp, cmp_reduce_thr) and get_endpoints (reads rho_by_cat_qu_*[r, jc] and joins
#' drho for labels / metadata / direction). It is the R analogue of the Python `rho_specs` (model_irt.py), but
#' makes the COMPARISONS explicit (each item, then each item-group) because the Stan model reduces per
#' comparison. No list columns: a group's member items are recovered from the item rows sharing its construct.
#'
#' Columns:
#'   jc               int    comparison index (shared by its reduction rows); items 1..n_item, then item-groups
#'   r                int    reduction slot within a comparison (1..R); red_fn is constant per r across jc (cartesian)
#'   scope            chr    "item" | "item_group"
#'   item_type_id     int    Stan cmp_cat_type (category type)
#'   item_type        chr    e.g. "out-of-7" / "categorical"
#'   cmp_key          chr    identifies the comparison: item_label (items) or construct (item-groups)
#'   red_fn           int    reduction id (see reduction_functions.stan): 1 GMFR, 2 SCR, 3 rel% (per-participant),
#'                           4 diff (per-participant), 5 baseline level, 6 endline level, 7 rel% (group-mean),
#'                           8 diff (group-mean), 9 GMFR (group-mean)
#'   rho              chr    short reduction label (Baseline/Endline/ratio/diff/...)
#'   thr              int    threshold for the level/threshold family: 0 = E[y], >=1 = P(y>=thr)
#'   item_label, construct, construct_long, item_high_label, item_label_short  chr  output metadata
#'   (item_label NA / item_label_short NA for item-group rows)
#'
#' @param dp1 preprocessed data.table (item_type_id, item_label, item_type, item_time_id, group_label, ...)
#' @param dit item metadata (item_type_id, item_label, construct, construct_long, item_high_label, item_label_short)
#' @param categorical_threshold integer y-cut for categorical items (0-indexed); out-of-7 use expected level (thr 0)
#' @param reductions data.table(r, red_fn, rho) listing the reductions to compute for EVERY comparison. Default is
#'   the Colombia production set: Baseline/Endline levels (5,6) + the CORRECT per-participant pre/post change
#'   rel% (3) and diff (4). (For the legacy group-mean reproduction pass c(5,6,7,8).)
#' @return drho data.table, ordered by (jc, r)
#' @importFrom data.table data.table as.data.table setkey setorder setcolorder fifelse rbindlist
#' @export
define_effect_sizes_colombia <- function(dp1, dit, categorical_threshold = 3L,
                                         reductions = data.table::data.table(
                                             r = 1:4, red_fn = c(5L, 6L, 3L, 4L),
                                             rho = c("Baseline", "Endline", "ratio", "diff"))) {
    item_type_id <- item_label <- item_type <- group_label <- construct <-
        construct_long <- item_high_label <- item_label_short <- jc <- NULL
    stopifnot(all(c("r", "red_fn", "rho") %in% names(reductions)))
    reductions <- data.table::as.data.table(reductions)

    # items measured at BOTH timepoints
    both <- intersect(dp1[group_label == "Baseline", unique(item_label)],
                      dp1[group_label == "Endline",  unique(item_label)])
    it <- unique(dp1[item_label %in% both, .(item_type_id, item_label, item_type)])
    it <- merge(it, unique(dit[, .(item_type_id, item_label, construct,
                                   construct_long, item_high_label, item_label_short)]),
                by = c("item_type_id", "item_label"))
    data.table::setkey(it, item_type_id, item_label)
    it[, jc := seq_len(.N)]
    n_item <- nrow(it)

    # item-groups: one per (item_type_id, construct); members are the item rows sharing construct
    gr <- unique(it[, .(item_type_id, item_type, construct, construct_long, item_high_label)])
    data.table::setkey(gr, item_type_id, construct)
    gr[, jc := n_item + seq_len(.N)]

    cmp <- rbind(
        it[, .(jc, scope = "item", item_type_id, item_type, cmp_key = item_label,
               item_label, construct, construct_long, item_high_label, item_label_short)],
        gr[, .(jc, scope = "item_group", item_type_id, item_type, cmp_key = construct,
               item_label = NA_character_, construct, construct_long, item_high_label,
               item_label_short = NA_character_)],
        use.names = TRUE)
    cmp[, thr := data.table::fifelse(item_type == "categorical", as.integer(categorical_threshold), 0L)]
    data.table::setorder(cmp, jc)

    # cartesian comparison x reduction (red_fn constant per r across comparisons)
    nred <- nrow(reductions)
    drho <- cmp[rep(seq_len(.N), each = nred)]
    drho <- cbind(drho, reductions[rep(seq_len(nred), times = nrow(cmp)), .(r, red_fn, rho)])
    data.table::setcolorder(drho, c("jc", "r", "scope", "item_type_id", "item_type", "cmp_key",
                                     "red_fn", "rho", "thr", "item_label", "construct",
                                     "construct_long", "item_high_label", "item_label_short"))
    data.table::setorder(drho, jc, r)
    drho[]
}

#' Define the effect-size manifest for the Ukraine application.
#'
#' NOTE: the current UkraineRCT estimand is a BETWEEN-ARM (disjoint-participant) contrast, so only the
#' group-mean reductions (5-9) are meaningful (per-participant pairing has no both-sides units). This stub
#' reuses the Colombia builder with the group-mean set; the baseline/endline "questions" must index the ARM
#' contrast (not pre/post) once the arm-keyed dp1 encoding is wired -- revisit when porting Ukraine to v261012.
#' @export
define_effect_sizes_ukraine <- function(dp1, dit, categorical_threshold = 2L,
                                        reductions = data.table::data.table(
                                            r = 1:4, red_fn = c(5L, 6L, 7L, 8L),
                                            rho = c("Baseline", "Endline", "ratio", "diff"))) {
    define_effect_sizes_colombia(dp1, dit, categorical_threshold = categorical_threshold,
                                 reductions = reductions)
}
