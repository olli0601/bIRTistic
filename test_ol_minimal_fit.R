library(data.table)
library(here)
source(here("R", "read_data_colombia.R"))
source(here("R", "fit_ordered_logit_model_ncats.R"))

# Define directories
dir_data <- "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/2025_project_Hope_Groups/data"
dir_out_base <- "/Users/or105/sandbox"
file_data <- file.path(dir_data, "Colombia_data_baseline_endline_itemised_250927.csv")

set.seed(42L)
seed <- 123L

# Read data
tmp <- read_data_colombia(file_data)
dp_col <- copy(tmp$dp)
dit_col <- copy(tmp$dit)
dmeta_col <- copy(tmp$dmeta)

# Preprocess data
dp1_col <- subset(dp_col, !grepl("agg", item_label))
dp1_col[, y_stan := y + 1L]
tmp <- CJ(
    item_label = sort(unique(dp1_col$item_label)),
    time = unique(dp1_col$time)
)
tmp[, item_type := ifelse(
    grepl("CG-MH", item_label), "categorical", "out-of-7"
)]
tmp <- tmp[
    , list(
        item_label = item_label, time = time,
        item_time_id = seq_along(item_label)
    ),
    by = "item_type"
]
dp1_col <- merge(dp1_col, tmp, by = c("item_label", "time"))
dp1_col <- merge(dp1_col, unique(subset(dit_col, select = c(item_type, item_type_id))), by = "item_type")
setkey(dp1_col, item_type_id, pid, time, item_label)
dp1_col[, oid := .I]
setkey(dp1_col, oid)
tmp <- dp1_col[, list(oid = oid, oidt = seq_along(y_stan)), by = "item_type"]
dp1_col <- merge(dp1_col, tmp, by = c("item_type", "oid"))
setkey(dp1_col, pid, time, item_label)

# Setup directories
dir_out_ol <- file.path(dir_out_base, "colombia-plots-tables-for-paper-ordered-logit-260413-vanilla")
dir_logs_ol <- file.path(dir_out_ol, "logs")
if (!dir.exists(dir_out_ol)) {
    dir.create(dir_out_ol, recursive = TRUE)
}
if (!dir.exists(dir_logs_ol)) {
    dir.create(dir_logs_ol, recursive = TRUE)
}

output_file_prefix <- file.path(dir_out_ol, "ol_1")

# Try fitting with minimal iterations first to test
cat("\n=== Testing ordered logit fit with minimal settings ===\n")
fit_ordered_logit_model_ncats(
    dit_col,
    dp1_col,
    output_file_prefix = paste0(output_file_prefix, "_fit_vanilla_test"),
    stan_file = here::here("src", "stan", "ordered_logit_ncats_v260413.stan"),
    chains = 1L,
    parallel_chains = 1L,
    threads_per_chain = 1L,
    iter_warmup = 50L,
    iter_sampling = 50L,
    seed = seed,
    x_formula = ~ time - 1,
    resume = FALSE,
    with_core_analyses = FALSE,
    with_additional_analyses = FALSE,
    show_messages = TRUE,
    show_exceptions = TRUE
)

cat("\n=== Test completed ===\n")
