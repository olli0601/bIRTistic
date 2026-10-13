// Partial Credit Model with Arbitrary Number of Category Types
// Version: v261012  (adds per-participant PAIRED reduced effect sizes: participant_reduced_rho_qu_pr,
//   R x Qcmp, e.g. GMFR + SCR for influenza HAI; reduction done INSIDE the model so output stays tiny)
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
  #include reduction_functions.stan
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
  
  // --- per-participant PAIRED reduced effect sizes (R reductions x Qcmp pairwise question comparisons) ---
  // Each reduction is a per-participant functional of a (question_cmp[,1], question_cmp[,2]) pair, averaged
  // over units present at BOTH questions. Emitted both as _fit (full eta, unit RE retained -> exact per-
  // participant joint) and _pr (confounding-adjusted eta). Output R x Qcmp per MC draw. Reduction-specific
  // constants (e.g. SCR cut/floor) live INSIDE reduction_functions.stan, NOT in this data block.
  int<lower=0> R; // number of reduced individual effect sizes (0 = disabled)
  int<lower=0> Qcmp; // number of pairwise question comparisons
  array[Qcmp] int<lower=1, upper=C> cmp_cat_type; // category type of each comparison pair
  array[Qcmp, 2] int<lower=1> question_cmp; // [,1] = baseline question, [,2] = endline question (local ids)
  array[R] int<lower=1> red_fn; // reduction fn per r: 1 = GMFR (2^mean log2-fold), 2 = SCR (seroconversion), 3 = rel % change (1-indexed Likert)
}
transformed data {
  real s2z_sd_unit;
  int skill_threshold_total_length;
  int ordered_prob_length;
  int eta_length;
  int K_max = max(K);
  int Q_max = max(Q);
  int N_max = max(N);
  array[C] int skill_threshold_offset;
  array[C] int ordered_prob_offset;
  array[N_total] int eta_obs_offset;
  array[C] int cat_obs_start;
  array[C] int cat_obs_end;
  array[C] int cat_question_offset;
  array[C] matrix[K_max - 1, K_max - 1] cat_cumsum_matrix;
  array[C] matrix[K_max, K_max - 1] cat_cumsum_matrix_brier;
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
  
  // Compute total length for skill thresholds flat vector and offsets
  skill_threshold_total_length = 0;
  for (c in 1 : C) {
    skill_threshold_offset[c] = skill_threshold_total_length;
    skill_threshold_total_length += Q[c] * (K[c] - 1);
  }
  
  // Compute total length for ordered probabilities (averaged by question)
  ordered_prob_length = 0;
  for (c in 1 : C) {
    ordered_prob_offset[c] = ordered_prob_length;
    ordered_prob_length += Q[c] * K[c];
  }
  
  // Compute total length for eta long vector and offsets for each observation
  eta_length = 0; // total number of eta values across all observations
  for (i in 1 : N_total) {
    eta_obs_offset[i] = eta_length + 1; // 1-indexed
    eta_length += K[cat_type[i]];
  }
  
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
  
  // Parameters stored in flat vector format (no padding)
  vector[skill_threshold_total_length] skill_thresholds;
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
      
      // Extract and reshape skill thresholds for this category type
      // Reshape: Stan's to_matrix fills column-by-column, but we stored row-by-row
      // So we reshape to (Km1c x Qc) and transpose to get (Qc x Km1c)  
      matrix[Qc, Km1c] skill_thresholds_c;
      skill_thresholds_c = to_matrix(
                                     segment(skill_thresholds,
                                             skill_threshold_offset[c] + 1,
                                             Qc * Km1c),
                                     Km1c, Qc)';
      
      // Compute eta matrix for this category type (fully vectorized)
      matrix[N[c], K[c]] eta_matrix_c;
      eta_matrix_c = pcm_get_etas(
                                  append_row(1.0,
                                             segment(loadings_questions_m1,
                                                     q_off - (c - 1) + 1,
                                                     Qc - 1)),
                                  skill_thresholds_c,
                                  latent_factor_unit[unit_of_obs[this_cat_obs_start : this_cat_obs_end]]
                                  + X[this_cat_obs_start : this_cat_obs_end,  : ]
                                    * latent_factor_beta,
                                  question_of_obs[this_cat_obs_start : this_cat_obs_end],
                                  unit_of_obs[this_cat_obs_start : this_cat_obs_end],
                                  cat_cumsum_matrix[c][1 : Km1c, 1 : Km1c]);
      
      // Flatten matrix row-by-row and insert into eta at correct positions
      // Since observations are contiguous, assign entire category at once
      this_cat_obs_start = eta_obs_offset[this_cat_obs_start];
      this_cat_obs_end = eta_obs_offset[this_cat_obs_end] + Km1c;
      eta[this_cat_obs_start : this_cat_obs_end] = to_vector(eta_matrix_c');
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
  target += normal_lupdf(skill_thresholds | 0, 3.5);
  
  // priors for loadings
  target += student_t_lupdf(loadings_questions_m1 | 3, 0, 1);
}
generated quantities {
  vector[ordered_prob_length] ordered_prob_by_cat_qu_fit; // averaged probabilities per (cat_type, question), best fit
  vector[ordered_prob_length] ordered_prob_by_cat_qu_pr; // averaged probabilities per (cat_type, question), posterior predictive, potentially excluding X columns
  array[N_total] int<lower=0> ypred;
  array[N_total] real log_lik;
  array[C] vector[Q_max] ordinal_brier_score;
  matrix[R, Qcmp] participant_reduced_rho_by_cat_qu_fit; // per-participant paired reduced effect sizes (R x Qcmp)
  matrix[R, Qcmp] participant_reduced_rho_by_cat_qu_pr; // per-participant paired reduced effect sizes (R x Qcmp)
  
  {
    // Per-observation category probs, reused across the FIT then PR passes (reset to -1 between).
    array[C] matrix[N_max, K_max] this_cat_prob_obs;
    
    // ===== FIT pass: full eta (unit RE + all X) =====
    // 1. precompute this_cat_prob_obs + ypred + log_lik
    for (c in 1 : C) {
      int Nc = N[c];
      int Kc = K[c];
      int this_cat_obs_start = cat_obs_start[c];
      for (n in 1 : Nc) {
        int obs_idx = this_cat_obs_start + n - 1;
        vector[Kc] tmp_Kc = segment(eta, eta_obs_offset[obs_idx], Kc);
        this_cat_prob_obs[c][n, 1 : Kc] = softmax(tmp_Kc)';
        ypred[obs_idx] = categorical_logit_rng(tmp_Kc);
        log_lik[obs_idx] = categorical_logit_lpmf(y[obs_idx] | tmp_Kc);
      }
    }
    // compute ordered_prob_by_cat_qu_fit from this_cat_prob_obs
    for (c in 1 : C) {
      int Nc = N[c];
      int Qc = Q[c];
      int Kc = K[c];
      int Km1c = Kc - 1;
      int this_cat_obs_start = cat_obs_start[c];
      matrix[Qc, Kc] this_cat_ordered_prob_by_question = rep_matrix(0., Qc,
                                                                    Kc);
      array[Qc] int count_by_question = rep_array(0, Qc);
      ordinal_brier_score[c] = rep_vector(0., Q_max);
      for (n in 1 : Nc) {
        int obs_idx = this_cat_obs_start + n - 1;
        this_cat_ordered_prob_by_question[question_of_obs[obs_idx],  : ] += this_cat_prob_obs[c][n, 1 : Kc];
        count_by_question[question_of_obs[obs_idx]] += 1;
      }
      for (q in 1 : Qc) {
        int start = ordered_prob_offset[c] + (q - 1) * Kc + 1;
        ordered_prob_by_cat_qu_fit[start : (start + Kc - 1)] = to_vector
                                                               (
                                                                this_cat_ordered_prob_by_question[q,  : ])
                                                               / count_by_question[q];
      }
      // compute ordinal Brier from this_cat_prob_obs
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
                                           this_cat_prob_obs[c][1 : Nc, 1 : Kc],
                                           segment(question_of_obs,
                                                   this_cat_obs_start, Nc),
                                           cat_cumsum_matrix_brier[c][1 : Kc, 1 : Km1c],
                                           cat_cum_obs[c][1 : Nc, 1 : Km1c],
                                           q_start_c, q_end_c, sorted_idx_c);
      }
    }
    // 3. participant_reduced_rho_by_cat_qu_fit: compute effect sizes on per-participant paired joint questions
    participant_reduced_rho_by_cat_qu_fit = compute_reduced_rho(
                                              this_cat_prob_obs, R, Qcmp, U,
                                              question_cmp, cmp_cat_type,
                                              red_fn, K, N, cat_obs_start,
                                              unit_of_obs, question_of_obs);
    // 4. free the buffer before reusing it for the _pr pass
    for (c in 1 : C) {
      this_cat_prob_obs[c] = rep_matrix(-1., N_max, K_max);
    }
    
    // ===== PR pass: confounding-adjusted eta (only Xpr_id columns of X) =====
    if (P_pr == P) {
      // all predictors kept -> PR == FIT
      ordered_prob_by_cat_qu_pr = ordered_prob_by_cat_qu_fit;
      participant_reduced_rho_by_cat_qu_pr = participant_reduced_rho_by_cat_qu_fit;
    } else {
      // 5. precompute this_cat_prob_obs from eta_pr (only Xpr_id columns of X)
      for (c in 1 : C) {
        int Nc = N[c];
        int Qc = Q[c];
        int Kc = K[c];
        int Km1c = Kc - 1;
        int q_off = cat_question_offset[c];
        int this_cat_obs_start = cat_obs_start[c];
        int this_cat_obs_end = cat_obs_end[c];
        matrix[Qc, Km1c] skill_thresholds_c = to_matrix(
                                                        segment(
                                                                skill_thresholds,
                                                                skill_threshold_offset[c]
                                                                + 1,
                                                                Qc * Km1c),
                                                        Km1c, Qc)';
        matrix[Nc, Kc] eta_matrix_pr = pcm_get_etas(
                                                    append_row(1.0,
                                                               segment
                                                               (
                                                                loadings_questions_m1,
                                                                q_off
                                                                - (c - 1) + 1,
                                                                Qc - 1)),
                                                    skill_thresholds_c,
                                                    X[this_cat_obs_start : this_cat_obs_end, Xpr_id]
                                                    * latent_factor_beta[Xpr_id],
                                                    question_of_obs[this_cat_obs_start : this_cat_obs_end],
                                                    unit_of_obs[this_cat_obs_start : this_cat_obs_end],
                                                    cat_cumsum_matrix[c][1 : Km1c, 1 : Km1c]);
        for (n in 1 : Nc) {
          this_cat_prob_obs[c][n, 1 : Kc] = softmax(eta_matrix_pr[n,  : ]')';
        }
      }
      // 6. ordered_prob_by_cat_qu_pr from this_cat_prob_obs
      for (c in 1 : C) {
        int Nc = N[c];
        int Qc = Q[c];
        int Kc = K[c];
        int this_cat_obs_start = cat_obs_start[c];
        matrix[Qc, Kc] this_cat_ordered_prob_by_question_pr = rep_matrix(0.,
                                                                    Qc, Kc);
        array[Qc] int count_by_question = rep_array(0, Qc);
        for (n in 1 : Nc) {
          int obs_idx = this_cat_obs_start + n - 1;
          this_cat_ordered_prob_by_question_pr[question_of_obs[obs_idx],  : ] += this_cat_prob_obs[c][n, 1 : Kc];
          count_by_question[question_of_obs[obs_idx]] += 1;
        }
        for (q in 1 : Qc) {
          int start = ordered_prob_offset[c] + (q - 1) * Kc + 1;
          ordered_prob_by_cat_qu_pr[start : (start + Kc - 1)] = to_vector
                                                                (
                                                                 this_cat_ordered_prob_by_question_pr[q,  : ])
                                                                / count_by_question[q];
        }
      }
      // 7. participant_reduced_rho_by_cat_qu_pr from the same confounding-adjusted probs
      participant_reduced_rho_by_cat_qu_pr = compute_reduced_rho(
                                               this_cat_prob_obs, R, Qcmp, U,
                                               question_cmp, cmp_cat_type,
                                               red_fn, K, N, cat_obs_start,
                                               unit_of_obs, question_of_obs);
    }
  }
}
