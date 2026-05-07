#!/usr/bin/env Rscript

# Helper script to prepare data for HPC analysis
# This script exports the necessary data objects from the main analysis
# into separate RDS files that can be used by the HPC script

library(data.table)

# Source: Run the main Rmd up to the data preparation step to get dcat6, dit, and di
# Then run this script to export the data

cat("Preparing data files for HPC analysis...\n")

# Define output directory for prepared data
data_export_dir <- "data_for_hpc"
dir.create(data_export_dir, showWarnings = FALSE, recursive = TRUE)

# Export dit (item type definitions)
cat("Saving dit...\n")
saveRDS(dit, file = file.path(data_export_dir, "dit.rds"))

# Create job configuration tracking
job_config <- data.table(JOB_ID = integer(), DCATI_FILE = character())

# Export dcati for each analysis
cat("\\nProcessing and saving dcati for each analysis...\\n")
for (idx in seq_len(nrow(di))) {
    job_id <- di[idx, interim_id]
    interim_date <- di[idx, interim_date]

    cat("  Processing job", job_id, "(", as.character(interim_date), ")...\\n")

    # Process data for this analysis
    dcati <- subset(dcat6, submission_date <= interim_date)
    tmp <- dcati[, list(n = length(unique(time))), by = c("pid", "item_label")]
    tmp <- tmp[, list(complete = all(n == 2)), by = "pid"]
    tmp <- subset(tmp, complete)
    tmp[, pid_new := seq_len(nrow(tmp))]
    dcati <- merge(dcati, tmp, by = "pid")
    setnames(
        dcati,
        c("pid", "oid", "oidt", "pid_new"),
        c("pid_orig", "oid_orig", "oidt_orig", "pid")
    )
    dcati[, oid := seq_len(nrow(dcati))]
    tmp <- dcati[, list(oid = oid, oidt = seq_along(y_stan)), by = "item_type"]
    dcati <- merge(dcati, tmp, by = c("item_type", "oid"))

    # Save dcati for this job
    dcati_filename <- paste0("dcati_job_", job_id, ".rds")
    output_file <- file.path(data_export_dir, dcati_filename)
    saveRDS(dcati, file = output_file)
    cat("    Saved", nrow(dcati), "observations to", basename(output_file), "\\n")

    # Add to job configuration
    job_config <- rbind(job_config, data.table(JOB_ID = job_id, DCATI_FILE = dcati_filename))
}

# Save job configuration CSV
config_file <- "job_config.csv"
fwrite(job_config, file = config_file)
cat("\\nSaved job configuration to:", config_file, "\\n")

cat("\\nData files prepared successfully!\\n")
cat("Files saved to:", data_export_dir, "\\n")
cat("\\nFiles created:\\n")
cat("  - dit.rds (item definitions)\\n")
cat("  - dcati_job_*.rds (processed data for each job, n=", nrow(di), ")\\n", sep = "")
cat("  - ../", config_file, " (job configuration CSV)\\n", sep = "")
