// Partial Credit Model with Arbitrary Number of Category Types
// Version: v250413
//
// This model extends partial_credit_model_2cats_v251224.stan to handle an arbitrary
// number of category types (C >= 1).
//
// Key differences from 2-category version:
// - Supports C category types instead of hardcoded cat1/cat2
// - Uses segment-based indexing for flexible category types
// - Each category type can have different N, Q, and K values
//
// Data structure:
// - C: number of category types  
// - N[c]: number of observations for category type c
// - Q[c]: number of questions for category type c  
// - K[c]: number of response categories for category type c
// - Data stored in concatenated arrays with obs_start[c]:obs_end[c] indexing
//
functions {
  #include credit_model_functions.stan
}
data {
  int<lower=1> C; // number of category types
  int<lower=1> P; // number of predictors
  int<lower=1> U; // number of units
  
  array[C] int<lower=1> N; // number of observations per category type
  array[C] int<lower=2> Q; // number of questions per category type
  array[C] int<lower=3> K; // number of categories per category type
  
  int<lower=1> N_total; // sum(N) - total observations across all category types
  int<lower=1> Q_total; // sum(Q) - total questions across all category types
  
  // Concatenated observation data for all category types
  array[N_total] int<lower=1> y;
  array[N_total] int<lower=1, upper=U> unit_of_obs;
  array[N_total] int<lower=1> question_of_obs; // local question id within category type
  array[N_total] int<lower=1, upper=C> cat_type; // which category type each observation belongs to
  matrix[N_total, P] X;
}
transformed data {
  real s2z_sd_unit;
  array[C] int cat_obs_start;
  array[C] int cat_obs_end;
  array[C] int cat_question_offset;
  int K_max = max(K);
  int Q_max = max(Q);
  int N_max = max(N);
  
  // Check that cat_type is ordered (all cat_type[i]==1, then all ==2, etc.)
  for (i in 2 : N_total) {
    if (cat_type[i] < cat_type[i - 1]) {
      reject("cat_type must be ordered: all observations for category type 1, then 2, etc.");
    }
  }
  
  // Derive observation ranges for each category type
  {
    int current_cat = 1;
    cat_obs_start[1] = 1;
    
    for (i in 2 : N_total) {
      if (cat_type[i] > cat_type[i - 1]) {
        cat_obs_end[current_cat] = i - 1;
        current_cat += 1;
        cat_obs_start[current_cat] = i;
      }
    }
    cat_obs_end[C] = N_total;
  }
  
  // Derive question offsets for each category type
  cat_question_offset[1] = 0;
  for (c in 2 : C) {
    cat_question_offset[c] = cat_question_offset[c - 1] + Q[c - 1];
  }
  
  // Compute total length for eta long vector and offsets for each observation
  int eta_length = 0; // total number of eta values across all observations
  array[N_total] int eta_obs_offset; // starting position in eta vector for each observation
  
  for (i in 1 : N_total) {
    eta_obs_offset[i] = eta_length + 1; // 1-indexed
    eta_length += K[cat_type[i]];
  }
  
  // Store processing structures for each category type
  array[C] matrix[K_max - 1, K_max - 1] cat_cumsum_matrix;
  array[C] matrix[K_max, K_max - 1] cat_cumsum_matrix_brier;
  array[C] matrix[N_max, K_max - 1] cat_cum_obs;
  array[C, Q_max] int cat_question_start;
  array[C, Q_max] int cat_question_end;
  array[C, N_max] int cat_sorted_indices;
  
  s2z_sd_unit = inv(sqrt(1. - inv(U)));
  
  // For each category type, create necessary matrices and indices
  for (c in 1 : C) {
    int Nc = N[c];
    int Qc = Q[c];
    int Kc = K[c];
    int Km1c = Kc - 1;
    int this_cat_obs_start = cat_obs_start[c];
    int this_cat_obs_end = cat_obs_end[c];
    
    // Create upper triangular matrix for cumulative sums (PCM parameterization)
    cat_cumsum_matrix[c] = rep_matrix(0., K_max - 1, K_max - 1);
    for (i in 1 : Km1c) {
      for (j in i : Km1c) {
        cat_cumsum_matrix[c][i, j] = 1.;
      }
    }
    
    // Create cumsum matrix for Brier score (K x K-1 format)
    cat_cumsum_matrix_brier[c] = rep_matrix(0., K_max, K_max - 1);
    for (i in 1 : Kc) {
      for (j in 1 : Km1c) {
        cat_cumsum_matrix_brier[c][i, j] = (i <= j) ? 1.0 : 0.0;
      }
    }
    
    // Precompute cumulative observation indicators
    cat_cum_obs[c] = rep_matrix(0., N_max, K_max - 1);
    for (n in 1 : Nc) {
      int obs_idx = this_cat_obs_start + n - 1;
      for (k in 1 : Km1c) {
        cat_cum_obs[c][n, k] = (y[obs_idx] <= k) ? 1.0 : 0.0;
      }
    }
    
    // Build sorted indices and start/end positions
    {
      array[Q_max] int count = rep_array(0, Q_max);
      
      // Count observations per question
      for (i in this_cat_obs_start : this_cat_obs_end) {
        count[question_of_obs[i]] += 1;
      }
      
      // Compute start positions
      for (q in 1 : Q_max) {
        cat_question_start[c, q] = 0;
        cat_question_end[c, q] = 0;
      }
      cat_question_start[c, 1] = 1;
      for (q in 2 : Qc) {
        cat_question_start[c, q] = cat_question_start[c, q - 1]
                                   + count[q - 1];
      }
      
      // Compute end positions
      for (q in 1 : Qc) {
        cat_question_end[c, q] = cat_question_start[c, q] + count[q] - 1;
      }
      
      // Fill sorted indices (relative to category type's observations, 1-indexed)
      for (n in 1 : N_max) {
        cat_sorted_indices[c, n] = 0;
      }
      count = rep_array(0, Q_max);
      for (n in 1 : Nc) {
        int obs_idx = this_cat_obs_start + n - 1;
        int pos = cat_question_start[c, question_of_obs[obs_idx]]
                  + count[question_of_obs[obs_idx]];
        cat_sorted_indices[c, pos] = n;
        count[question_of_obs[obs_idx]] += 1;
      }
    }
  }
}
parameters {
  sum_to_zero_vector[U] latent_factor_unit;
  vector[P] latent_factor_beta;
  
  // Parameters stored in concatenated format
  matrix[Q_total, K_max - 1] skill_thresholds;
  vector<lower=0>[Q_total - C] loadings_questions_m1; // all questions minus one per category type
}
transformed parameters {
  vector[eta_length] eta; // long vector storing all eta values without padding
  
  {
    int Qc;
    int Km1c;
    int q_off;
    int this_cat_obs_start;
    int this_cat_obs_end;
    for (c in 1 : C) {
      Qc = Q[c];
      Km1c = K[c] - 1;
      q_off = cat_question_offset[c];
      this_cat_obs_start = cat_obs_start[c];
      this_cat_obs_end = cat_obs_end[c];
      
      // Compute eta matrix for this category type (fully vectorized)
      matrix[N[c], K[c]] eta_matrix_c;
      eta_matrix_c = pcm_get_etas(
                                  append_row(1.0,
                                             segment(loadings_questions_m1,
                                                     q_off - (c - 1) + 1,
                                                     Qc - 1)),
                                  skill_thresholds[(q_off + 1) : (q_off + Qc), 1 : Km1c],
                                  latent_factor_unit[unit_of_obs[this_cat_obs_start : this_cat_obs_end]]
                                  + X[this_cat_obs_start : this_cat_obs_end,  : ]
                                    * latent_factor_beta,
                                  question_of_obs[this_cat_obs_start : this_cat_obs_end],
                                  unit_of_obs[this_cat_obs_start : this_cat_obs_end],
                                  cat_cumsum_matrix[c][1 : Km1c, 1 : Km1c]);
      
      // Flatten matrix to vector and insert into eta at correct positions
      eta[eta_obs_offset[this_cat_obs_start] : (eta_obs_offset[this_cat_obs_end]
                                                + Km1c)] = to_vector(
                                                                    eta_matrix_c');
    }
  }
}
model {
  // likelihood under partial credit model
  for (c in 1 : C) {
    for (obs_idx in cat_obs_start[c] : cat_obs_end[c]) {
      target += categorical_logit_lupmf(y[obs_idx] |
                  segment(eta, eta_obs_offset[obs_idx], K[c]));
    }
  }
  
  // priors for latent factors
  target += normal_lupdf(latent_factor_unit | 0, s2z_sd_unit);
  target += std_normal_lupdf(latent_factor_beta);
  
  // priors for skill thresholds
  target += normal_lupdf(to_vector(skill_thresholds) | 0, 3.5);
  
  // priors for loadings
  target += student_t_lupdf(loadings_questions_m1 | 3, 0, 1);
}
generated quantities {
  vector[eta_length] ordered_prob_by_obs; // long vector storing all probabilities without padding
  array[N_total] int<lower=0> ypred;
  array[N_total] real log_lik;
  array[C] vector[Q_max] ordinal_brier_score;
  
  for (c in 1 : C) {
    int Nc = N[c];
    int Qc = Q[c];
    int Kc = K[c];
    int Km1c = Kc - 1;
    int this_cat_obs_start = cat_obs_start[c];
    int this_cat_obs_end = cat_obs_end[c];
    matrix[Nc, Kc] this_cat_prob_obs;
    
    // Initialize category-specific structures
    ordinal_brier_score[c] = rep_vector(0., Q_max);
    
    // Compute probabilities and predictions
    for (n in 1 : Nc) {
      int obs_idx = this_cat_obs_start + n - 1;
      int eta_start = eta_obs_offset[obs_idx];
      vector[Kc] eta_obs = segment(eta, eta_start, Kc);
      vector[Kc] prob_obs = softmax(eta_obs);
      
      this_cat_prob_obs[n,  : ] = prob_obs';
      ypred[obs_idx] = categorical_logit_rng(eta_obs);
      log_lik[obs_idx] = categorical_logit_lpmf(y[obs_idx] | eta_obs);
      
      // Store probabilities in long vector
      ordered_prob_by_obs[(eta_start) : (eta_start + Kc - 1)] = prob_obs;
    }
    
    // Compute ordinal Brier scores by question
    {
      array[Qc] int q_start_c;
      array[Qc] int q_end_c;
      array[Nc] int sorted_idx_c;
      
      for (q in 1 : Qc) {
        q_start_c[q] = cat_question_start[c, q];
        q_end_c[q] = cat_question_end[c, q];
      }
      for (n in 1 : Nc) {
        sorted_idx_c[n] = cat_sorted_indices[c, n];
      }
      
      ordinal_brier_score[c][1 : Qc] = pcm_get_ordered_brier_score(
                                         this_cat_prob_obs,
                                         question_of_obs[this_cat_obs_start : (
                                         this_cat_obs_start + Nc - 1)],
                                         cat_cumsum_matrix_brier[c][1 : Kc, 1 : Km1c],
                                         cat_cum_obs[c][1 : Nc, 1 : Km1c],
                                         q_start_c, q_end_c, sorted_idx_c);
    }
  }
}
