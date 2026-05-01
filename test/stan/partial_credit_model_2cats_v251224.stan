functions {
  #include credit_model_functions.stan
}
data {
  int<lower=1> P; // number of predictors
  int<lower=1> U; // number of units
  
  int<lower=1> Ncat1; // number of ordered categorial observations
  int<lower=2> Qcat1; // number of ordered categorial questions
  int<lower=3> Kcat1; // number of categories
  array[Ncat1] int<lower=1, upper=Kcat1> cat1_y; // observations
  array[Ncat1] int<lower=1, upper=U> cat1_unit_of_obs;
  array[Ncat1] int<lower=1, upper=Qcat1> cat1_question_of_obs;
  matrix[Ncat1, P] cat1_X;
  
  int<lower=1> Ncat2;
  int<lower=2> Qcat2;
  int<lower=3> Kcat2;
  array[Ncat2] int<lower=1, upper=Kcat2> cat2_y; // observations
  array[Ncat2] int<lower=1, upper=U> cat2_unit_of_obs;
  array[Ncat2] int<lower=1, upper=Qcat2> cat2_question_of_obs;
  matrix[Ncat2, P] cat2_X;
}
transformed data {
  real s2z_sd_unit;
  matrix[Kcat1 - 1, Kcat1 - 1] cumsum_matrix_cat1;
  matrix[Kcat2 - 1, Kcat2 - 1] cumsum_matrix_cat2;
  matrix[Kcat1, Kcat1 - 1] cumsum_matrix_cat1_brier;
  matrix[Kcat2, Kcat2 - 1] cumsum_matrix_cat2_brier;
  matrix[Ncat1, Kcat1 - 1] cum_obs_cat1;
  matrix[Ncat2, Kcat2 - 1] cum_obs_cat2;
  array[Qcat1] int cat1_question_start;
  array[Qcat1] int cat1_question_end;
  array[Ncat1] int cat1_sorted_indices;
  array[Qcat2] int cat2_question_start;
  array[Qcat2] int cat2_question_end;
  array[Ncat2] int cat2_sorted_indices;
  
  s2z_sd_unit = inv(sqrt(1. - inv(U)));
  
  // Create upper triangular matrices for cumulative sums (PCM parameterization)
  cumsum_matrix_cat1 = rep_matrix(0., Kcat1 - 1, Kcat1 - 1);
  for (i in 1 : (Kcat1 - 1)) {
    for (j in i : (Kcat1 - 1)) {
      cumsum_matrix_cat1[i, j] = 1.;
    }
  }
  
  cumsum_matrix_cat2 = rep_matrix(0., Kcat2 - 1, Kcat2 - 1);
  for (i in 1 : (Kcat2 - 1)) {
    for (j in i : (Kcat2 - 1)) {
      cumsum_matrix_cat2[i, j] = 1.;
    }
  }
  
  // Create cumsum matrices for Brier score (K x K-1 format)
  for (i in 1 : Kcat1) {
    for (j in 1 : (Kcat1 - 1)) {
      cumsum_matrix_cat1_brier[i, j] = (i <= j) ? 1.0 : 0.0;
    }
  }
  
  for (i in 1 : Kcat2) {
    for (j in 1 : (Kcat2 - 1)) {
      cumsum_matrix_cat2_brier[i, j] = (i <= j) ? 1.0 : 0.0;
    }
  }
  
  // Precompute cumulative observation indicators
  for (n in 1 : Ncat1) {
    for (k in 1 : (Kcat1 - 1)) {
      cum_obs_cat1[n, k] = (cat1_y[n] <= k) ? 1.0 : 0.0;
    }
  }
  
  for (n in 1 : Ncat2) {
    for (k in 1 : (Kcat2 - 1)) {
      cum_obs_cat2[n, k] = (cat2_y[n] <= k) ? 1.0 : 0.0;
    }
  }
  
  // Build sorted indices and start/end positions for cat1
  {
    array[Qcat1] int count = rep_array(0, Qcat1);
    // Count observations per question
    for (n in 1 : Ncat1) {
      count[cat1_question_of_obs[n]] += 1;
    }
    // Compute start positions
    cat1_question_start[1] = 1;
    for (q in 2 : Qcat1) {
      cat1_question_start[q] = cat1_question_start[q - 1] + count[q - 1];
    }
    // Compute end positions
    for (q in 1 : Qcat1) {
      cat1_question_end[q] = cat1_question_start[q] + count[q] - 1;
    }
    // Fill sorted indices
    count = rep_array(0, Qcat1); // Reset to use as write position
    for (n in 1 : Ncat1) {
      int q = cat1_question_of_obs[n];
      int pos = cat1_question_start[q] + count[q];
      cat1_sorted_indices[pos] = n;
      count[q] += 1;
    }
  }
  
  // Build sorted indices and start/end positions for cat2
  {
    array[Qcat2] int count = rep_array(0, Qcat2);
    // Count observations per question
    for (n in 1 : Ncat2) {
      count[cat2_question_of_obs[n]] += 1;
    }
    // Compute start positions
    cat2_question_start[1] = 1;
    for (q in 2 : Qcat2) {
      cat2_question_start[q] = cat2_question_start[q - 1] + count[q - 1];
    }
    // Compute end positions
    for (q in 1 : Qcat2) {
      cat2_question_end[q] = cat2_question_start[q] + count[q] - 1;
    }
    // Fill sorted indices
    count = rep_array(0, Qcat2); // Reset to use as write position
    for (n in 1 : Ncat2) {
      int q = cat2_question_of_obs[n];
      int pos = cat2_question_start[q] + count[q];
      cat2_sorted_indices[pos] = n;
      count[q] += 1;
    }
  }
}
parameters {
  sum_to_zero_vector[U] latent_factor_unit;
  vector[P] latent_factor_beta;
  matrix[Qcat1, Kcat1 - 1] cat1_skill_thresholds;
  vector<lower=0>[Qcat1 - 1] cat1_loadings_questions_m1;
  matrix[Qcat2, Kcat2 - 1] cat2_skill_thresholds;
  vector<lower=0>[Qcat2 - 1] cat2_loadings_questions_m1;
}
transformed parameters {
  matrix[Ncat1, Kcat1] cat1_eta;
  matrix[Ncat2, Kcat2] cat2_eta;
  cat1_eta = pcm_get_etas(append_row(1.0, cat1_loadings_questions_m1),
                          cat1_skill_thresholds,
                          latent_factor_unit[cat1_unit_of_obs]
                          + cat1_X * latent_factor_beta,
                          cat1_question_of_obs, cat1_unit_of_obs,
                          cumsum_matrix_cat1);
  cat2_eta = pcm_get_etas(append_row(1.0, cat2_loadings_questions_m1),
                          cat2_skill_thresholds,
                          latent_factor_unit[cat2_unit_of_obs]
                          + cat2_X * latent_factor_beta,
                          cat2_question_of_obs, cat2_unit_of_obs,
                          cumsum_matrix_cat2);
}
model {
  // likelihood under partial credit model
  for (n in 1 : Ncat1) {
    target += categorical_logit_lupmf(cat1_y[n] | cat1_eta[n]');
  }
  for (n in 1 : Ncat2) {
    target += categorical_logit_lupmf(cat2_y[n] | cat2_eta[n]');
  }
  
  // priors for latent factors
  target += normal_lupdf(latent_factor_unit | 0, s2z_sd_unit);
  target += std_normal_lupdf(latent_factor_beta);
  
  // priors for skill thresholds
  target += normal_lupdf(to_vector(cat1_skill_thresholds) | 0, 3.5);
  target += normal_lupdf(to_vector(cat2_skill_thresholds) | 0, 3.5);
  
  // priors for loadings
  target += student_t_lupdf(cat1_loadings_questions_m1 | 3, 0, 1);
  target += student_t_lupdf(cat2_loadings_questions_m1 | 3, 0, 1);
}
generated quantities {
  matrix[Ncat1, Kcat1] cat1_ordered_prob_by_obs;
  matrix[Ncat2, Kcat2] cat2_ordered_prob_by_obs;
  array[Ncat1] int<lower=0> cat1_ypred;
  array[Ncat2] int<lower=0> cat2_ypred;
  array[Ncat1 + Ncat2] real log_lik;
  vector[Qcat1] cat1_ordinal_brier_score;
  vector[Qcat2] cat2_ordinal_brier_score;
  
  for (n in 1 : Ncat1) {
    cat1_ordered_prob_by_obs[n,  : ] = softmax(cat1_eta[n]')';
    cat1_ypred[n] = categorical_logit_rng(cat1_eta[n]');
  }
  for (n in 1 : Ncat2) {
    cat2_ordered_prob_by_obs[n,  : ] = softmax(cat2_eta[n]')';
    cat2_ypred[n] = categorical_logit_rng(cat2_eta[n]');
  }
  
  // Compute ordinal Brier scores by question
  cat1_ordinal_brier_score = pcm_get_ordered_brier_score(
                               cat1_ordered_prob_by_obs,
                               cat1_question_of_obs,
                               cumsum_matrix_cat1_brier, cum_obs_cat1,
                               cat1_question_start, cat1_question_end,
                               cat1_sorted_indices);
  cat2_ordinal_brier_score = pcm_get_ordered_brier_score(
                               cat2_ordered_prob_by_obs,
                               cat2_question_of_obs,
                               cumsum_matrix_cat2_brier, cum_obs_cat2,
                               cat2_question_start, cat2_question_end,
                               cat2_sorted_indices);
  
  for (n in 1 : Ncat1) {
    log_lik[n] = categorical_logit_lpmf(cat1_y[n] | cat1_eta[n]');
  }
  for (n in 1 : Ncat2) {
    log_lik[Ncat1 + n] = categorical_logit_lpmf(cat2_y[n] | cat2_eta[n]');
  }
}
