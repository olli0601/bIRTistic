# Generate R baselines for Ukraine data
# Generated: 2026-05-05 for Task 2.2 validation

cat("============================================================\n")
cat("Generating R baselines for Ukraine data\n")
cat("============================================================\n\n")

# Source the R function
source("R/read_data_ukraine.R")

# Data file path
file_data <- "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/2025_project_Hope_Groups/data/Ukraine_Hope_Groups_Baseline_Endline_Wide_Aug6.csv"

cat("Reading Ukraine data...\n")
tmp <- read_data_ukraine(file_data)

cat("Data loaded:\n")
cat(sprintf("  dp: %d rows x %d columns\n", nrow(tmp$dp), ncol(tmp$dp)))
cat(sprintf("  dit: %d rows x %d columns\n", nrow(tmp$dit), ncol(tmp$dit)))
cat(sprintf("  dmeta: %d rows x %d columns\n", nrow(tmp$dmeta), ncol(tmp$dmeta)))

# Create output directory if needed
dir.create("python/tests/fixtures/r_baselines/phase2", recursive = TRUE, showWarnings = FALSE)

# Write baselines
cat("\nWriting baselines to python/tests/fixtures/r_baselines/phase2/...\n")

write.csv(tmp$dp, "python/tests/fixtures/r_baselines/phase2/ukraine_dp.csv", row.names = FALSE)
cat("  ✓ ukraine_dp.csv written\n")

write.csv(tmp$dit, "python/tests/fixtures/r_baselines/phase2/ukraine_dit.csv", row.names = FALSE)
cat("  ✓ ukraine_dit.csv written\n")

write.csv(tmp$dmeta, "python/tests/fixtures/r_baselines/phase2/ukraine_dmeta.csv", row.names = FALSE)
cat("  ✓ ukraine_dmeta.csv written\n")

cat("\n============================================================\n")
cat("✓ R baselines generated successfully\n")
cat("============================================================\n")
