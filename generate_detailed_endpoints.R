# Generate missing detailed endpoint files
library(data.table)
library(posterior)

dir_out_base <- "/Users/or105/sandbox"

# Process all three datasets in sequence
datasets_to_process <- list(
    list(
        dir = "colombia-plots-tables-for-paper-ordered-logit-260413-ps-tsd",
        data_file = "ol_1_ps_tsd_data.rds",
        prefix = "ol_1_fit_ps_tsd_",
        categorical_threshold = 3
    ),
    list(
        dir = "colombia-plots-tables-for-paper-ordered-logit-260413-matched-ps-tsd",
        data_file = "ol_1_matched_ps_tsd_data.rds",
        prefix = "ol_1_fit_matched_ps_tsd_",
        categorical_threshold = 3
    ),
    list(
        dir = "colombia-plots-tables-for-paper-ordered-logit-260413-ukraine-comparison",
        data_file = "ol_1_ukraine_vanilla_data.rds",
        prefix = "ol_1_ukraine_fit_vanilla_",
        categorical_threshold = 2
    )
)

for (ds in datasets_to_process) {
    cat("\nProcessing dataset (detailed):", ds$dir, "\n")
    
    tryCatch({
        # Load data
        dir_out_ol <- file.path(dir_out_base, ds$dir)
        
        # Check if detailed file already exists
        detailed_file <- file.path(dir_out_ol, paste0(ds$prefix, "primary_endpoints_detailed.rds"))
        if (file.exists(detailed_file)) {
            cat("File already exists, skipping:", detailed_file, "\n")
            next
        }
        
        # Check if draws file exists
        draws_file <- file.path(dir_out_ol, paste0(ds$prefix, "draws.rds"))
        if (!file.exists(draws_file)) {
            cat("Draws file not found, skipping:", draws_file, "\n")
            next
        }
        
        tmp <- readRDS(file.path(dir_out_ol, ds$data_file))
        dp1 <- copy(tmp$dp1)
        dit <- copy(tmp$dit)        
        output_file_prefix <- file.path(dir_out_ol, ds$prefix)    
        poa <- readRDS(paste0(output_file_prefix, "draws.rds"))
        categorical_threshold <- ds$categorical_threshold

        # Aggregate categorical probabilities for most + all of the time
        po <- poa[,, grepl("ordered_prob_by_cat_qu_pr", dimnames(poa)[[3]]), drop = FALSE]            
        po <- posterior::as_draws_df(po)
        po <- as.data.table(po)
        po <- data.table::melt(po,
            id.vars = c(".draw", ".chain", ".iteration"),
            value.name = "prob"
            )
        po[, cq_id := as.integer(gsub(
            "ordered_prob_by_cat_qu_pr\\[([0-9]+)\\]", "\\1",
            variable
            ))]
        set(po, NULL, "variable", NULL)
        tmp <- unique(subset(dp1, select = c(item_type_id, item_label, item_time_id)))
        tmp <- merge(tmp, subset(dit, select = c(item_type_id, item_label, cat_length)), by = c("item_type_id", "item_label"))
        setkey(tmp, item_type_id, item_time_id)
        tmp[, cq_id := cumsum(cat_length) - cat_length + 1L]
        tmp <- tmp[, list(cq_id = cq_id - 1L + 1:cat_length, y = 1:cat_length - 1L), by = c("item_type_id", "item_time_id")]
        po <- merge(tmp, po, by = "cq_id")  
        
        tmp <- unique(subset(dp1[item_type == "categorical", ],
            select = c(item_type_id, item_label, item_time_id, item_type, time_label)
        ))
        po <- merge(po, tmp, by = c("item_type_id", "item_time_id"))
        po <- subset(po, y >= categorical_threshold)
        po <- po[,
            list(y_stan = 45, prob = sum(prob)),
            by = c(".draw", "item_type_id", "item_time_id", "item_label", "item_type", "time_label")
        ]
        po <- merge(po, subset(dit, select = c(item_type, item_label, group_label)), by = c("item_type", "item_label"))
        po <- data.table::dcast(po,
            .draw + item_type_id + item_time_id + item_label + item_type + group_label ~ time_label,
            value.var = "prob"
        )
        po <- subset(po, !is.na(Baseline) & !is.na(Endline))
        po[, diff := Baseline - Endline]
        po[, ratio := 1 - Endline / Baseline]
        po[, diff2 := Endline - Baseline]
        po[, ratio2 := Endline / Baseline - 1]
        po <- melt(po,
            id.vars = c(".draw", "item_type_id", "item_time_id", "item_label", "item_type", "group_label"),
            measure.vars = c("diff", "ratio", "diff2", "ratio2", "Baseline", "Endline")
        )
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
        tmp <- unique(subset(dit, select = c(item_type, item_label, item_label_short, group_label, group_label_long, item_high_label)))
        pos <- merge(pos, tmp, by = c("item_type", "item_label"))
        pos_cat1 <- copy(pos)
        po <- NULL
        gc()
        
        # Compute mean for out-of-7 responses
        po <- poa[,, grepl("ordered_prob_by_cat_qu_pr", dimnames(poa)[[3]]), drop = FALSE]            
        po <- posterior::as_draws_df(po)
        po <- as.data.table(po)
        po <- data.table::melt(po,
            id.vars = c(".draw", ".chain", ".iteration"),
            value.name = "prob"
            )
        po[, cq_id := as.integer(gsub(
            "ordered_prob_by_cat_qu_pr\\[([0-9]+)\\]", "\\1",
            variable
            ))]
        set(po, NULL, "variable", NULL)
        tmp <- unique(subset(dp1, select = c(item_type_id, item_label, item_time_id)))
        tmp <- merge(tmp, subset(dit, select = c(item_type_id, item_label, cat_length)), by = c("item_type_id", "item_label"))
        setkey(tmp, item_type_id, item_time_id)
        tmp[, cq_id := cumsum(cat_length) - cat_length + 1L]
        tmp <- tmp[, list(cq_id = cq_id - 1L + 1:cat_length, y = 1:cat_length - 1L), by = c("item_type_id", "item_time_id")]
        po <- merge(tmp, po, by = "cq_id")  
        
        tmp <- unique(subset(dp1[item_type == "out-of-7", ],
            select = c(item_type_id, item_label, item_time_id, item_type, time_label)
        ))
        po <- merge(po, tmp, by = c("item_type_id", "item_time_id"))
        po <- po[,
            list(mean = sum(y * prob)),
            by = c(".draw", "item_label", "item_time_id", "item_type", "time_label")
        ]
        po <- merge(po, subset(dit, select = c(item_type, item_label, group_label)), by = c("item_type", "item_label"))
        po <- data.table::dcast(po,
            .draw + item_label + item_time_id + item_type + group_label ~ time_label,
            value.var = "mean"
        )
        po <- subset(po, !is.na(Baseline) & !is.na(Endline))
        po <- merge(po,
            unique(subset(dit, select = c(item_type, group_label, item_high_label))),
            by = c("item_type", "group_label")
        )
        po[, diff := NA_real_]
        po[, ratio := NA_real_]
        tmp <- po[, which(item_high_label == "lower_is_better")]
        set(po, tmp, "diff", po[tmp, Baseline - Endline])
        set(po, tmp, "ratio", po[tmp, 1 - Endline / Baseline])
        tmp <- po[, which(item_high_label == "higher_is_better")]
        set(po, tmp, "diff", po[tmp, Endline - Baseline])
        set(po, tmp, "ratio", po[tmp, Endline / Baseline - 1])
        po[, diff2 := Endline - Baseline]
        po[, ratio2 := Endline / Baseline - 1]
        po <- melt(po,
            id.vars = c(".draw", "item_label", "item_time_id", "item_type", "group_label", "item_high_label"),
            measure.vars = c("diff", "ratio", "diff2", "ratio2", "Baseline", "Endline")
        )
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
        tmp <- unique(subset(dit, select = c(item_type, item_label, item_label_short, group_label, group_label_long, item_high_label)))
        pos <- merge(pos, tmp, by = c("item_type", "item_label", "item_high_label"))
        pos_cat2 <- copy(pos)
        po <- NULL
        gc()
        
        pos <- rbind(pos_cat1, pos_cat2, fill = TRUE)
        
        cat("\nwrite endpoints to file=", paste0(output_file_prefix, "primary_endpoints_detailed.rds"), "\n")
        saveRDS(pos,
            file = paste0(output_file_prefix, "primary_endpoints_detailed.rds"))
        
        cat("Successfully created:", paste0(output_file_prefix, "primary_endpoints_detailed.rds"), "\n")
        
    }, error = function(e) {
        cat("Error processing", ds$dir, ":", conditionMessage(e), "\n")
        traceback()
    })
}

cat("\nDone processing all datasets\n")
