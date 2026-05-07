#!/usr/bin/env Rscript
# Generate Colombia baseline outputs from R function
# Author: testing-engineer
# Date: 2026-05-05

# Load function
source("R/read_data_colombia.R")

# Set data file path
file_data <- "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/2025_project_Hope_Groups/data/Colombia_data_baseline_endline_itemised_250927.csv"

# Run function
cat("Running read_data_colombia()...\n")
tmp <- read_data_colombia(file_data)

# Inspect outputs
cat("\nOutput structure:\n")
cat("dp dimensions:", dim(tmp$dp), "\n")
cat("dit dimensions:", dim(tmp$dit), "\n")
cat("dmeta dimensions:", dim(tmp$dmeta), "\n")

# Save to CSV (for Python to load)
cat("\nSaving baseline CSV files...\n")
write.csv(tmp$dp, "python/tests/fixtures/r_baselines/phase2/colombia_dp.csv", row.names = FALSE)
write.csv(tmp$dit, "python/tests/fixtures/r_baselines/phase2/colombia_dit.csv", row.names = FALSE)
write.csv(tmp$dmeta, "python/tests/fixtures/r_baselines/phase2/colombia_dmeta.csv", row.names = FALSE)

# Verify files created
cat("\nGenerated files:\n")
files <- list.files("python/tests/fixtures/r_baselines/phase2", pattern = "colombia.*\\.csv$", full.names = TRUE)
for(f in files) {
    cat("  ", f, " (", file.size(f), " bytes)\n", sep="")
}

cat("\n✓ Colombia baseline generation complete\n")
