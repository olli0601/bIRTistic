// Ordered Logit Model with Arbitrary Number of Category Types
// Version: v260413
//
// This model extends ordered_logit_2cats_v251222.stan to handle an arbitrary
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
  #include ordered_logit_functions.stan
}
data {
  int<lower=1> C; // number of category types
  int<lower=1> P; // number of predictors
  int<lower=1> U; // number of units
  int<lower=1, upper=P> P_pr; // number of predictors to include in the model for predicting the ordered probabilities
  
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
  array[P_pr] int<lower=1, upper=P> Xpr_id; // X indices to include in the model for predicting the ordered probabilities
}
transformed data {
  real s2z_sd_unit;
  real s2z_sd_questions;
  int skill_threshold_total_length;
  int cutpoint_total_length;
  int ordered_prob_length;
  int K_max = max(K);
  int Q_max = max(Q);
  int N_max = max(N);
  array[C] int skill_threshold_offset;
  array[C] int cutpoint_offset;
  array[C] int ordered_prob_offset;
  array[C] int cat_obs_start;
  array[C] int cat_obs_end;
  array[C] int cat_question_offset;
  array[C] matrix[K_max, K_max - 1] cat_cumsum_matrix;
  array[C] matrix[N_max, K_max - 1] cat_cum_obs;
  array[C, Q_max] int cat_question_start;
  array[C, Q_max] int cat_question_end;
  array[C, N_max] int cat_sorted_indices;
  
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
  
  // Compute total length for skill threshold increments flat vector and offsets
  skill_threshold_total_length = 0;
  for (c in 1 : C) {
    skill_threshold_offset[c] = skill_threshold_total_length;
    skill_threshold_total_length += Q[c] * (K[c] - 2);
  }
  
  // Compute total length for cutpoints flat vector
  cutpoint_total_length = 0;
  for (c in 1 : C) {
    cutpoint_offset[c] = cutpoint_total_length;
    cutpoint_total_length += Q[c] * (K[c] - 1);
  }
  
  // Compute total length for ordered probabilities (averaged by question)
  ordered_prob_length = 0;
  for (c in 1 : C) {
    ordered_prob_offset[c] = ordered_prob_length;
    ordered_prob_length += Q[c] * K[c];
  }
  
  s2z_sd_unit = inv(sqrt(1. - inv(U)));
  s2z_sd_questions = inv(sqrt(1. - inv(Q_total)));
  
  // For each category type, create necessary matrices and indices
  for (c in 1 : C) {
    int Nc = N[c];
    int Qc = Q[c];
    int Kc = K[c];
    int Km1c = Kc - 1;
    int this_cat_obs_start = cat_obs_start[c];
    int this_cat_obs_end = cat_obs_end[c];
    
    // Create cumsum matrix for Brier score (K x K-1 format)
    cat_cumsum_matrix[c] = rep_matrix(0., K_max, K_max - 1);
    for (i in 1 : Kc) {
      for (j in 1 : Km1c) {
        cat_cumsum_matrix[c][i, j] = (i <= j) ? 1.0 : 0.0;
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
  
  // Parameters stored in flat vector format (no padding)
  vector[Q_total] skill_thresholds_1; // first threshold for each question
  vector<lower=0>[skill_threshold_total_length] skill_thresholds_incs; // increments for subsequent thresholds
  vector<lower=0>[Q_total - C] loadings_questions_m1; // all questions minus one per category type
}
transformed parameters {
  vector[N_total] eta; // long vector storing all eta values
  vector[cutpoint_total_length] cutpoints; // flat vector of all cutpoints
  
  {
    int Qc;
    int Km1c;
    int Km2c;
    int q_off;
    int this_cat_obs_start;
    int this_cat_obs_end;
    for (c in 1 : C) {
      Qc = Q[c];
      Km1c = K[c] - 1;
      Km2c = K[c] - 2;
      q_off = cat_question_offset[c];
      this_cat_obs_start = cat_obs_start[c];
      this_cat_obs_end = cat_obs_end[c];
      
      // Extract and reshape skill thresholds for this category type
      vector[Qc] skill_thresholds_1_c = segment(skill_thresholds_1,
                                                q_off + 1, Qc);
      matrix[Qc, Km2c] skill_thresholds_incs_c;
      if (Km2c > 0) {
        skill_thresholds_incs_c = to_matrix(
                                            segment(skill_thresholds_incs,
                                                    skill_threshold_offset[c]
                                                    + 1, Qc * Km2c),
                                            Km2c, Qc)';
      }
      
      // Compute cutpoints for this category type
      vector[Qc] loadings_c = append_row(1.0,
                                         segment(loadings_questions_m1,
                                                 q_off - (c - 1) + 1, 
                                                 Qc - 1));
      matrix[Qc, Km1c] cutpoints_c;
      if (Km2c > 0) {
        cutpoints_c = cm_get_cutpoints_noncentered(skill_thresholds_1_c,
                        skill_thresholds_incs_c, loadings_c);
      } else {
        // K=3 case: only one cutpoint per question
        cutpoints_c = rep_matrix(skill_thresholds_1_c .* loadings_c, 1);
      }
      
      // Flatten and store cutpoints
      int cp_start = cutpoint_offset[c] + 1;
      int cp_end = cp_start + Qc * Km1c - 1;
      cutpoints[cp_start : cp_end] = to_vector(cutpoints_c');
      
      // Compute eta for this category type (fully vectorized)
      eta[this_cat_obs_start : this_cat_obs_end] = ol_get_etas(loadings_c,
                                                               latent_factor_unit[unit_of_obs[this_cat_obs_start : this_cat_obs_end]]
                                                               + X[this_cat_obs_start : this_cat_obs_end,  : ]
                                                                 * latent_factor_beta,
                                                               question_of_obs[this_cat_obs_start : this_cat_obs_end]);
    }
  }
}
model {
  // likelihood under ordered logit model
  for (c in 1 : C) {
    int Km1c = K[c] - 1;
    int cp_off = cutpoint_offset[c];
    for (obs_idx in cat_obs_start[c] : cat_obs_end[c]) {
      int q_local = question_of_obs[obs_idx];
      int cp_start = cp_off + (q_local - 1) * Km1c + 1;
      target += ordered_logistic_lupmf(y[obs_idx] | eta[obs_idx],
                  segment(cutpoints, cp_start, Km1c));
    }
  }
  
  // priors for latent factors
  target += normal_lupdf(latent_factor_unit | 0, s2z_sd_unit);
  target += std_normal_lupdf(latent_factor_beta);
  
  // priors for skill thresholds
  target += normal_lupdf(skill_thresholds_1 | 0, 3.5);
  target += normal_lupdf(skill_thresholds_incs | 0, 2.5);
  
  // priors for loadings
  target += student_t_lupdf(loadings_questions_m1 | 3, 0, 1);
}
generated quantities {
  vector[ordered_prob_length] ordered_prob_by_cat_qu_fit; // averaged probabilities per (cat_type, question), best fit
  vector[ordered_prob_length] ordered_prob_by_cat_qu_pr; // averaged probabilities per (cat_type, question), posterior predictive, potentially excluding X columns
  array[N_total] int<lower=0> ypred;
  array[N_total] real log_lik;
  array[C] vector[Q_max] ordinal_brier_score;
  
  for (c in 1 : C) {
    int Nc = N[c];
    int Qc = Q[c];
    int Kc = K[c];
    int Km1c = Kc - 1;
    int cp_off = cutpoint_offset[c];
    int this_cat_obs_start = cat_obs_start[c];
    int this_cat_obs_end = cat_obs_end[c];
    int obs_idx;
    int start;
    vector[Kc] tmp_Kc;
    vector[Kc] prob_obs;
    matrix[Nc, Kc] this_cat_prob_obs;
    matrix[Qc, Kc] thiscat_ordered_prob_by_question = rep_matrix(0., Qc, Kc);
    array[Qc] int count_by_question = rep_array(0, Qc);
    
    // Initialize category-specific structures
    ordinal_brier_score[c] = rep_vector(0., Q_max);
    
    // Compute probabilities and predictions
    for (n in 1 : Nc) {
      obs_idx = this_cat_obs_start + n - 1;
      int q_local = question_of_obs[obs_idx];
      int cp_start = cp_off + (q_local - 1) * Km1c + 1;
      vector[Km1c] cutpoints_q = segment(cutpoints, cp_start, Km1c);
      
      // Compute ordered probabilities
      tmp_Kc[1] = 1 - inv_logit(eta[obs_idx] - cutpoints_q[1]);
      for (k in 2 : Km1c) {
        tmp_Kc[k] = inv_logit(eta[obs_idx] - cutpoints_q[k - 1])
                    - inv_logit(eta[obs_idx] - cutpoints_q[k]);
      }
      tmp_Kc[Kc] = inv_logit(eta[obs_idx] - cutpoints_q[Km1c]);
      
      this_cat_prob_obs[n,  : ] = tmp_Kc';
      ypred[obs_idx] = ordered_logistic_rng(eta[obs_idx], cutpoints_q);
      log_lik[obs_idx] = ordered_logistic_lpmf(y[obs_idx] | eta[obs_idx],
                           cutpoints_q);
      thiscat_ordered_prob_by_question[q_local,  : ] += tmp_Kc';
      count_by_question[q_local] += 1;
    }
    
    // Compute averaged probabilities and store in long vector
    for (q in 1 : Qc) {
      start = ordered_prob_offset[c] + (q - 1) * Kc + 1;
      tmp_Kc = to_vector(thiscat_ordered_prob_by_question[q,  : ])
               / count_by_question[q];
      ordered_prob_by_cat_qu_fit[start : (start + Kc - 1)] = tmp_Kc;
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
      
      ordinal_brier_score[c][1 : Qc] = ol_get_ordered_brier_score(
                                         this_cat_prob_obs,
                                         segment(question_of_obs,
                                                 this_cat_obs_start, Nc),
                                         cat_cumsum_matrix[c][1 : Kc, 1 : Km1c],
                                         cat_cum_obs[c][1 : Nc, 1 : Km1c],
                                         q_start_c, q_end_c, sorted_idx_c);
    }
  }
  
  // Compute ordered_prob_by_cat_qu_pr
  if (P_pr == P) {
    // If using all predictors, just copy the fitted probabilities
    ordered_prob_by_cat_qu_pr = ordered_prob_by_cat_qu_fit;
  } else {
    // Recompute probabilities using only specified X columns
    for (c in 1 : C) {
      int Nc = N[c];
      int Qc = Q[c];
      int Kc = K[c];
      int Km1c = Kc - 1;
      int Km2c = Kc - 2;
      int q_off = cat_question_offset[c];
      int cp_off = cutpoint_offset[c];
      int this_cat_obs_start = cat_obs_start[c];
      int this_cat_obs_end = cat_obs_end[c];
      matrix[Qc, Kc] thiscat_ordered_prob_by_question_pr = rep_matrix(0., Qc,
                                                                    Kc);
      array[Qc] int count_by_question = rep_array(0, Qc);
      
      // Compute eta using only specified X columns
      vector[Qc] loadings_c = append_row(1.0,
                                         segment(loadings_questions_m1,
                                                 q_off - (c - 1) + 1, 
                                                 Qc - 1));
      vector[Nc] eta_pr = ol_get_etas(loadings_c,
                                      X[this_cat_obs_start : this_cat_obs_end, Xpr_id]
                                      * latent_factor_beta[Xpr_id],
                                      question_of_obs[this_cat_obs_start : this_cat_obs_end]);
      
      // Compute probabilities and average by question
      for (n in 1 : Nc) {
        int obs_idx = this_cat_obs_start + n - 1;
        int q_local = question_of_obs[obs_idx];
        int cp_start = cp_off + (q_local - 1) * Km1c + 1;
        vector[Km1c] cutpoints_q = segment(cutpoints, cp_start, Km1c);
        vector[Kc] prob_pr;
        
        prob_pr[1] = 1 - inv_logit(eta_pr[n] - cutpoints_q[1]);
        for (k in 2 : Km1c) {
          prob_pr[k] = inv_logit(eta_pr[n] - cutpoints_q[k - 1])
                       - inv_logit(eta_pr[n] - cutpoints_q[k]);
        }
        prob_pr[Kc] = inv_logit(eta_pr[n] - cutpoints_q[Km1c]);
        
        thiscat_ordered_prob_by_question_pr[q_local,  : ] += prob_pr';
        count_by_question[q_local] += 1;
      }
      
      // Store averaged probabilities in long vector
      for (q in 1 : Qc) {
        int start = ordered_prob_offset[c] + (q - 1) * Kc + 1;
        vector[Kc] tmp_Kc = to_vector(
                                      thiscat_ordered_prob_by_question_pr[q,  : ])
                            / count_by_question[q];
        ordered_prob_by_cat_qu_pr[start : (start + Kc - 1)] = tmp_Kc;
      }
    }
  }
}
