library(data.table)
library(here)
source(here("R", "read_data_colombia.R"))

# Define directories
dir_data <- "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/2025_project_Hope_Groups/data"
file_data <- file.path(dir_data, "Colombia_data_baseline_endline_itemised_250927.csv")

# Read data
tmp <- read_data_colombia(file_data)
dp_col <- copy(tmp$dp)
dit_col <- copy(tmp$dit)

# Preprocess
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

cat("Data prepared successfully\n")
cat("Rows:", nrow(dp1_col), "\n")
cat("Item types:", paste(unique(dp1_col$item_type), collapse=", "), "\n")
cat("y_stan range:", paste(range(dp1_col$y_stan, na.rm=TRUE), collapse=" to "), "\n")

# Check for categorical item type
cat_items <- unique(subset(dp1_col, item_type == "categorical", select=c(item_label, y_stan)))
cat("Categorical items:", nrow(cat_items), "\n")
if(nrow(cat_items) > 0) {
  cat("Categorical y_stan values:", paste(sort(unique(cat_items$y_stan)), collapse=", "), "\n")
}
