#!/usr/bin/env Rscript
# Unit test: Verify ordered_logit_ncats_v260413.stan computational correctness
#
# This test verifies that the ncats model produces identical log-likelihood values
# to the reference 2cats model when run with identical fixed parameter values.
#
# Expected outcome: PASS if max absolute difference in log_lik < 1e-10

library(testthat)
library(cmdstanr)
library(data.table)

# Set working directory to project root
if (!grepl("bIRTistic$", getwd())) {
  if (file.exists("test/test_ordered_logit_ncats_correctness.R")) {
    # Already in project root
  } else if (file.exists("../R")) {
    setwd("..")
  } else {
    stop("Cannot determine project root directory")
  }
}

test_that("Ordered logit ncats model produces identical log_lik to 2cats model with fixed parameters", {
  
  # Compile models
  model_2cats <- cmdstan_model(
    "src/stan/ordered_logit_2cats_v251222.stan",
    include_paths = "src/stan",
    quiet = TRUE
  )
  
  model_ncats <- cmdstan_model(
    "src/stan/ordered_logit_ncats_v260413.stan",
    include_paths = "src/stan",
    quiet = TRUE
  )
  
  # Create simple test data
  # 2 units, 1 predictor
  # Category 1: 2 observations, 2 questions, 3 response options (K=3)
  # Category 2: 2 observations, 2 questions, 4 response options (K=4)
  
  stan_data_2cats <- list(
    U = 2, P = 1,
    Ncat1 = 2, Qcat1 = 2, Kcat1 = 3,
    cat1_y = c(1, 2),
    cat1_question_of_obs = c(1, 2),
    cat1_unit_of_obs = c(1, 2),
    cat1_X = matrix(c(0.5, -0.3), ncol = 1),
    Ncat2 = 2, Qcat2 = 2, Kcat2 = 4,
    cat2_y = c(2, 3),
    cat2_question_of_obs = c(1, 2),
    cat2_unit_of_obs = c(1, 2),
    cat2_X = matrix(c(0.1, -0.2), ncol = 1)
  )
  
  # Prepare ncats data (concatenated format)
  combined <- rbind(
    data.table(
      y = stan_data_2cats$cat1_y,
      unit = stan_data_2cats$cat1_unit_of_obs,
      question = stan_data_2cats$cat1_question_of_obs,
      cat_type = 1L,
      X = stan_data_2cats$cat1_X[, 1]
    ),
    data.table(
      y = stan_data_2cats$cat2_y,
      unit = stan_data_2cats$cat2_unit_of_obs,
      question = stan_data_2cats$cat2_question_of_obs,
      cat_type = 2L,
      X = stan_data_2cats$cat2_X[, 1]
    )
  )
  
  stan_data_ncats <- list(
    C = 2L, U = 2, P = 1,
    N = c(2L, 2L),
    Q = c(2L, 2L),
    K = c(3L, 4L),
    N_total = 4L,
    Q_total = 4L,
    y = combined$y,
    unit_of_obs = combined$unit,
    question_of_obs = combined$question,
    cat_type = combined$cat_type,
    X = matrix(combined$X, ncol = 1),
    P_pr = 1L,
    Xpr_id = c(1L)
  )
  
  # Define fixed parameter values for ordered logit
  # For ordered logit: skill_thresholds_1 + cumsum(skill_thresholds_incs)
  # Category 1: K=3, so 2 cutpoints per question (Km1=2)
  #   - For each question: 1 skill_thresholds_1 + 1 skill_thresholds_inc (Km2=1)
  # Category 2: K=4, so 3 cutpoints per question (Km1=3)  
  #   - For each question: 1 skill_thresholds_1 + 2 skill_thresholds_incs (Km2=2)
  
  init_2cats <- list(
    latent_factor_unit = c(-0.5, 0.5),
    latent_factor_beta = 1.0,
    cat1_skill_thresholds_1 = c(0.5, -0.5),
    cat1_skill_thresholds_incs = matrix(c(0.8, 0.6), nrow = 2, byrow = TRUE),
    cat1_loadings_questions_m1 = 1.5,
    cat2_skill_thresholds_1 = c(0.3, -0.3),
    cat2_skill_thresholds_incs = matrix(c(0.5, 0.7, 0.4, 0.6), nrow = 2, byrow = TRUE),
    cat2_loadings_questions_m1 = 0.8
  )
  
  # Convert to ncats format - flat vectors without padding
  # skill_thresholds_1: all Q questions (4 total)
  # skill_thresholds_incs: cat1 has 2*1=2 values, cat2 has 2*2=4 values, total 6
  
  skill_thresholds_1_vec <- c(
    init_2cats$cat1_skill_thresholds_1,  # 2 values
    init_2cats$cat2_skill_thresholds_1   # 2 values
  )
  
  skill_thresholds_incs_vec <- c(
    as.vector(t(init_2cats$cat1_skill_thresholds_incs)),  # cat1: 2 values (2 questions × 1 inc)
    as.vector(t(init_2cats$cat2_skill_thresholds_incs))   # cat2: 4 values (2 questions × 2 incs)
  )
  
  init_ncats <- list(
    latent_factor_unit = c(-0.5, 0.5),
    latent_factor_beta = 1.0,
    skill_thresholds_1 = skill_thresholds_1_vec,
    skill_thresholds_incs = skill_thresholds_incs_vec,
    loadings_questions_m1 = c(1.5, 0.8)
  )
  
  # Run models with fixed parameters
  fit_2cats <- model_2cats$sample(
    data = stan_data_2cats,
    init = list(init_2cats),
    chains = 1L,
    iter_warmup = 0L,
    iter_sampling = 1L,
    fixed_param = TRUE,
    refresh = 0,
    show_messages = FALSE
  )
  
  fit_ncats <- model_ncats$sample(
    data = stan_data_ncats,
    init = list(init_ncats),
    chains = 1L,
    iter_warmup = 0L,
    iter_sampling = 1L,
    fixed_param = TRUE,
    refresh = 0,
    show_messages = FALSE
  )
  
  # Extract log_lik values
  ll_2cats <- as.numeric(fit_2cats$draws("log_lik", format = "draws_matrix"))
  ll_ncats <- as.numeric(fit_ncats$draws("log_lik", format = "draws_matrix"))
  
  # Test expectations
  expect_equal(length(ll_2cats), 4, 
               label = "2cats model should return 4 log_lik values")
  expect_equal(length(ll_ncats), 4,
               label = "ncats model should return 4 log_lik values")
  
  # Main assertion: models produce identical log_lik values
  expect_equal(ll_ncats, ll_2cats, tolerance = 1e-10,
               label = "ncats and 2cats models should produce identical log_lik values")
  
  # Also check element-wise for better diagnostics
  for (i in seq_along(ll_2cats)) {
    expect_equal(ll_ncats[i], ll_2cats[i], tolerance = 1e-10,
                 info = sprintf("Observation %d: ncats=%.6f, 2cats=%.6f", 
                                i, ll_ncats[i], ll_2cats[i]))
  }
})
