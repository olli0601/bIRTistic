#' Run ordered logit model analysis (ncats version) using ADVI
#'
#' This function performs Bayesian IRT analysis using the flexible ncats ordered logit model
#' with ADVI (Automatic Differentiation Variational Inference) instead of HMC.
#'
#' @param dit data.table. Item metadata table with item types and labels
#' @param dcati data.table. Pre-processed data with observations for analysis. Must include
#'   item_type, y_stan, item_time_id, pid, and oidt columns.
#' @param output_file_prefix Character. Full path prefix for output files (without extension)
#' @param stan_file Character. Path to Stan model file (.stan). Default: ordered_logit_ncats_v260413.stan
#' @param x_formula Formula. Formula specifying predictors for the design matrix
#' @param x_formula_ignore_regex Character. Regular expression pattern to identify X columns to exclude 
#'   from posterior predictive probability calculations (ordered_prob_by_cat_qu_pr). Columns matching this 
#'   pattern will be excluded. If NA (default), all X columns are used for both fitted and predictive probabilities.
#' @param iter Integer. Number of ADVI iterations (default: 10000)
#' @param grad_samples Integer. Number of samples for ELBO gradient estimation (default: 1)
#' @param elbo_samples Integer. Number of samples for ELBO evaluation (default: 100)
#' @param output_samples Integer. Number of samples to draw from the approximate posterior (default: 4000)
#' @param seed Integer. Random seed for reproducibility (default: 123)
#' @param show_messages Logical. If TRUE, show all Stan informational messages during sampling (default: FALSE)
#' @param show_exceptions Logical. If TRUE, show detailed Stan exception messages when errors occur (default: FALSE)
#' @param resume Logical. If TRUE and all output files exist, skip ADVI and load existing results (default: FALSE)
#'
#' @return Invisibly returns the fitted cmdstanr model object
#'
#' @import data.table
#' @import ggplot2
#' @import cmdstanr
#'
#' @export
fit_ordered_logit_model_ncats_advi <- function(
  dit,
  dcati,
  output_file_prefix,
  stan_file = here::here("src", "stan", "ordered_logit_ncats_v260413.stan"),
  x_formula = ~ time - 1,
  x_formula_ignore_regex = NA_character_,
  iter = 10000L,
  grad_samples = 1L,
  elbo_samples = 100L,
  output_samples = 4000L,
  seed = 123L,
  show_messages = FALSE,
  show_exceptions = FALSE,
  resume = FALSE
) {
    # Suppress data.table NSE warnings
    item_type <- y_stan <- item_time_id <- pid <- time <- oidt <- variable <-
        ypred <- oid <- in_ppi <- item_label <- time_label <- y_label <- y <-
        n <- total <- p_emp <- group_label_long <- item_label_short <-
        item_label_long <- ess_bulk <- q_lower <- q_upper <- iqr_lower <-
        iqr_upper <- prob <- item_type_id <- stat <- NULL

    require(data.table)
    require(ggplot2)
    require(cmdstanr)

    # Print configuration
    cat("\n========================================\n")
    cat("Ordered Logit Model (ncats) ADVI Analysis Configuration\n")
    cat("========================================\n")
    cat("Stan file:", stan_file, "\n")
    cat("Stan include dir:", dirname(stan_file), "\n")
    cat("Data: dit with", nrow(dit), "items, dcati with", nrow(dcati), "observations\n")
    cat("Output prefix:", output_file_prefix, "\n")
    cat("ADVI iterations:", iter, "\n")
    cat("Gradient samples:", grad_samples, "\n")
    cat("ELBO samples:", elbo_samples, "\n")
    cat("Output samples:", output_samples, "\n")
    cat("Seed:", seed, "\n")
    cat("Resume:", resume, "\n")
    cat("========================================\n\n")

    # Create output directory if it doesn't exist
    dir.create(dirname(output_file_prefix), showWarnings = FALSE, recursive = TRUE)
    
    # Define output file paths
    timing_file <- paste0(output_file_prefix, "_timing.csv")
    draws_file <- paste0(output_file_prefix, "_draws.rds")
    output_file <- paste0(output_file_prefix, "_stan.rds")
    
    # Check if we should resume from existing outputs
    can_resume <- resume && 
                  file.exists(timing_file) && 
                  file.exists(draws_file) && 
                  file.exists(output_file)

    # Compile Stan model
    cat("Compiling Stan model...\n")
    ol_compiled <- cmdstanr::cmdstan_model(
        stan_file,
        include_paths = dirname(stan_file),
        force_recompile = TRUE
    )

    # Prepare data in ncats format
    cat("Preparing Stan data in ncats format...\n")
    
    # Ensure oid is sequential and ordered by item_type_id and oidt
    setkey(dcati, item_type_id, oidt)
    stopifnot(all( 1:nrow(dcati) == dcati[,oid] ))  

    # Build stan_data
    stan_data <- list()
    stan_data$C <- dcati[, max(item_type_id)]
    stan_data$U <- dcati[, max(pid)] 
    
    # Arrays for each category type
    tmp <- dcati[, .(
        N = .N,
        Q = max(item_time_id),
        K = length(unique(y_stan))
    ), by = item_type_id]
    stan_data$N <- tmp$N
    stan_data$Q <- tmp$Q
    stan_data$K <- tmp$K
    
    # Total counts
    stan_data$N_total <- sum(stan_data$N)
    stan_data$Q_total <- sum(stan_data$Q)
    
    # Concatenated observation data (data is already sorted by item_type_id)
    stan_data$y <- dcati$y_stan
    stan_data$unit_of_obs <- dcati$pid
    stan_data$question_of_obs <- dcati$item_time_id
    stan_data$cat_type <- dcati$item_type_id
    
    # Create design matrix for all observations
    stan_data$X <- model.matrix(
        x_formula,
        data = as.data.frame(dcati)
    )
    if(!is.na(x_formula_ignore_regex)) {
        stan_data$Xpr_id <- which(!grepl(x_formula_ignore_regex, colnames(stan_data$X)))
    }
    else {
        stan_data$Xpr_id <- seq_len(ncol(stan_data$X))
    }
    stan_data$P <- ncol(stan_data$X)
    stan_data$P_pr <- length(stan_data$Xpr_id)
    
    cat("Design matrix number of predictors: P =", stan_data$P, "\n")
    cat("Design matrix column names (for fitting):", paste(colnames(stan_data$X), collapse = ", "), "\n")
    cat("Design matrix column names (for prediction):", paste(colnames(stan_data$X)[stan_data$Xpr_id], collapse = ", "), "\n")
    cat("Number of category types: C =", stan_data$C, "\n")
    cat("Category type labels:", paste(unique(dit$item_type), collapse = ", "), "\n")
    cat("N per category:", paste(stan_data$N, collapse = ", "), "\n")
    cat("Q per category:", paste(stan_data$Q, collapse = ", "), "\n")
    cat("K per category:", paste(stan_data$K, collapse = ", "), "\n")

    # Check if resuming from existing outputs
    if (can_resume) {
        cat("\n========================================\n")
        cat("RESUMING from existing outputs\n")
        cat("========================================\n")
        
        cat("Loading draws into poa from:", draws_file, "\n")
        poa <- readRDS(draws_file)
        
        cat("Loading timing data from:", timing_file, "\n")
        timing_data <- read.csv(timing_file)
        cat("========================================\n\n")
    } else {
        if (resume) {
            cat("\nNote: Resume requested but not all output files exist. Running full analysis.\n\n")
        }
        
        # Run ADVI
        cat("Running ADVI...\n")
        flush.console()
        start_time <- Sys.time()
        ol_fit <- ol_compiled$variational(
            data = stan_data,
            seed = seed,
            algorithm = "meanfield",
            iter = iter,
            grad_samples = grad_samples,
            elbo_samples = elbo_samples,
            output_samples = output_samples,
            show_messages = show_messages,
            show_exceptions = show_exceptions
        )
        end_time <- Sys.time()
        
        # Extract timing information
        cat("\nExtracting timing information...\n")
        total_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
        timing_data <- data.table(
            method = "ADVI",
            iterations = iter,
            total_minutes = total_minutes
        )
        
        write.csv(timing_data, file = timing_file, row.names = FALSE)
        cat("Saved timing information to:", timing_file, "\n")

        # Save output to RDS
        cat("Saving model fit to:", output_file, "\n")
        ol_fit$save_object(file = output_file)

        # Extract and save draws
        cat("Extracting draws...\n")
        poa <- ol_fit$draws(format = "draws_array")
        saveRDS(poa, file = draws_file)
        cat("Saved draws to:", draws_file, "\n")

        ol_fit <- NULL
        gc()
    }

    cat("\nADVI analysis complete!\n")
    cat("========================================\n\n")
    
    invisible(poa)
}
