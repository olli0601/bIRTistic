#' Create and submit a PBS job submission script
#'
#' Creates a PBS job submission script that runs an R script with a JSON input file,
#' and optionally submits it to the HPC scheduler.
#'
#' @param script_path Character. Path to the R script to execute (e.g., "R/train_test_item_response_models.Rscript")
#' @param json_path Character. Path to the JSON file with input arguments
#' @param job_name Character. Name for the PBS job
#' @param output_dir Character. Directory where PBS output and error logs will be written
#' @param walltime Character. Wall time limit in format "HH:MM:SS" (default: "24:00:00")
#' @param memory Character. Memory requirement (default: "8gb")
#' @param ncpus Integer. Number of CPUs to request (default: 1)
#' @param email Character or NULL. Email address for job notifications (default: NULL)
#' @param email_options Character. PBS email options: "a" (abort), "b" (begin), "e" (end) (default: "abe")
#' @param modules Character vector. Environment modules to load (default: c("R/4.5.2"))
#' @param additional_commands Character vector. Additional shell commands to run before executing the script (default: NULL)
#' @param submit Logical. Whether to submit the job immediately using qsub (default: FALSE)
#' @param pbs_script_path Character or NULL. Path where the PBS script should be saved. If NULL, uses output_dir/job_name.pbs
#' @return Character. Path to the created PBS submission script
#' @export
#' @examples
#' \dontrun{
#' # Create PBS script without submitting
#' pbs_file <- hpc_submit_generic_job(
#'     script_path = "R/train_test_item_response_models.Rscript",
#'     json_path = "inputs/config.json",
#'     job_name = "validate_model",
#'     output_dir = "logs",
#'     walltime = "12:00:00",
#'     memory = "16gb",
#'     ncpus = 4
#' )
#'
#' # Create and submit PBS script
#' pbs_file <- hpc_submit_generic_job(
#'     script_path = "R/train_test_item_response_models.Rscript",
#'     json_path = "inputs/config.json",
#'     job_name = "validate_model",
#'     output_dir = "logs",
#'     submit = TRUE
#' )
#' }
hpc_submit_generic_job <- function(
  script_path,
  json_path,
  job_name,
  output_dir,
  walltime = "24:00:00",
  memory = "8gb",
  ncpus = 1L,
  email = NULL,
  email_options = "abe",
  modules = c("R/4.5.2"),
  additional_commands = NULL,
  submit = FALSE,
  pbs_script_path = NULL
) {
    # Validate inputs
    if (!file.exists(script_path)) {
        stop("Script file not found: ", script_path)
    }
    if (!file.exists(json_path)) {
        stop("JSON file not found: ", json_path)
    }

    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        cat("Created output directory:", output_dir, "\n")
    }

    # Set PBS script path
    if (is.null(pbs_script_path)) {
        pbs_script_path <- file.path(output_dir, paste0(job_name, ".pbs"))
    }

    # Get absolute paths
    script_path_abs <- normalizePath(script_path, mustWork = TRUE)
    json_path_abs <- normalizePath(json_path, mustWork = TRUE)
    output_dir_abs <- normalizePath(output_dir, mustWork = TRUE)

    # Build PBS script content
    pbs_content <- c(
        "#!/bin/bash",
        "",
        "# PBS job submission script",
        paste0("#PBS -N ", job_name),
        paste0("#PBS -o ", file.path(output_dir_abs, paste0(job_name, ".out"))),
        paste0("#PBS -e ", file.path(output_dir_abs, paste0(job_name, ".err"))),
        paste0("#PBS -l walltime=", walltime),
        paste0("#PBS -l mem=", memory),
        paste0("#PBS -l ncpus=", ncpus)
    )

    # Add email notifications if specified
    if (!is.null(email)) {
        pbs_content <- c(
            pbs_content,
            paste0("#PBS -M ", email),
            paste0("#PBS -m ", email_options)
        )
    }

    pbs_content <- c(
        pbs_content,
        "",
        "# Change to submission directory",
        "cd $PBS_O_WORKDIR",
        "",
        "# Print job information",
        "echo \"Job started on $(hostname) at $(date)\"",
        "echo \"Job ID: $PBS_JOBID\"",
        "echo \"Job name: $PBS_JOBNAME\"",
        "echo \"Working directory: $PBS_O_WORKDIR\"",
        ""
    )

    # Add module loading
    if (length(modules) > 0) {
        pbs_content <- c(
            pbs_content,
            "# Load required modules",
            paste0("module load ", modules),
            ""
        )
    }

    # Add additional commands if specified
    if (!is.null(additional_commands) && length(additional_commands) > 0) {
        pbs_content <- c(
            pbs_content,
            "# Additional setup commands",
            additional_commands,
            ""
        )
    }

    # Add main execution command
    pbs_content <- c(
        pbs_content,
        "# Execute R script with JSON input",
        paste0("Rscript ", script_path_abs, " --json_file ", json_path_abs),
        "",
        "# Print job completion",
        "echo \"Job completed at $(date)\"",
        "exit 0"
    )

    # Write PBS script
    writeLines(pbs_content, pbs_script_path)
    cat("Created PBS submission script:", pbs_script_path, "\n")

    # Make script executable
    Sys.chmod(pbs_script_path, mode = "0755")

    # Submit job if requested
    if (submit) {
        cat("Submitting job to PBS scheduler...\n")
        submit_cmd <- paste("qsub", pbs_script_path)
        result <- system(submit_cmd, intern = TRUE)
        cat("Job submitted. Job ID:", result, "\n")
    } else {
        cat("PBS script created. To submit, run: qsub", pbs_script_path, "\n")
    }

    invisible(pbs_script_path)
}
