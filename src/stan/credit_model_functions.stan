matrix cm_get_skill_thresholds(vector skill_thresholds_1,
                               matrix skill_thresholds_incs)
{
  int Q = size(skill_thresholds_1);
  int Km1 = cols(skill_thresholds_incs) + 1;
  matrix [Q, Km1] skill_thresholds;
  skill_thresholds = rep_matrix(skill_thresholds_1, Km1);
  for(j in 1:Q)
  {
    skill_thresholds[j,2:Km1] += cumulative_sum( skill_thresholds_incs[j,:] );
  }
  return skill_thresholds;
}

matrix cm_get_skill_thresholds2(vector skill_thresholds_1,
                               matrix skill_thresholds_incs)
{
  int Q = size(skill_thresholds_1);
  int Km1 = cols(skill_thresholds_incs) + 1;
  matrix [Q, Km1] skill_thresholds;
  skill_thresholds = rep_matrix(skill_thresholds_1, Km1);
  for(j in 1:Q)
  {
    skill_thresholds[j,2:Km1] += skill_thresholds_incs[j,:];
  }
  return skill_thresholds;
}

matrix cm_get_etas(vector loadings_questions,
                   matrix skill_thresholds,
                   vector latent_factor_obs,
                   array [] int question_of_obs,
                   array [] int unit_of_obs)
{
  int N = num_elements(unit_of_obs);
  int Q = rows(skill_thresholds);
  int K = cols(skill_thresholds) + 1;

  matrix[N, K] etas;
  etas[:,1] = rep_vector(0., N);
  etas[:,2:K] = rep_matrix(loadings_questions[question_of_obs], K-1) .*
                  (rep_matrix(latent_factor_obs, K-1) -
                   skill_thresholds[question_of_obs,:]);
  return etas;
}

matrix cm_get_etas(vector loadings_questions,
                   array[] vector skill_thresholds,
                   vector latent_factor_obs,
                   array [] int question_of_obs,
                   array [] int unit_of_obs)
{
  int N = num_elements(unit_of_obs);
  int Q = size(skill_thresholds);
  int K = num_elements(skill_thresholds[1]);
  matrix[Q, K] skill_thresholds_mat;
  matrix[N, K] etas;
  for(q in 1:Q) {
    skill_thresholds_mat[q,:] = skill_thresholds[q]';
  }  
  etas = rep_matrix(loadings_questions[question_of_obs], K) .*
                    (rep_matrix(latent_factor_obs, K) -
                    skill_thresholds_mat[question_of_obs,:]);
  return etas;
}

matrix cm_get_etas(matrix loadings_questions,
                   matrix skill_thresholds,
                   vector latent_factor_obs,
                   array [] int question_of_obs,
                   array [] int unit_of_obs)
{
  int N = num_elements(unit_of_obs);
  int Q = rows(skill_thresholds);
  int K = cols(skill_thresholds) + 1;

  matrix[N, K] etas;
  etas[:,1] = rep_vector(0., N);
  etas[:,2:K] = loadings_questions[question_of_obs, :] .*
                  (rep_matrix(latent_factor_obs, K-1) -
                   skill_thresholds[question_of_obs,:]);
  return etas;
}

// Partial credit model get_etas function
matrix pcm_get_etas(vector loadings_questions,
                   matrix skill_thresholds,
                   vector latent_factor_obs,
                   array [] int question_of_obs,
                   array [] int unit_of_obs,
                   matrix cumsum_matrix)
{
  int N = num_elements(unit_of_obs);
  int Q = rows(skill_thresholds);
  int Km1 = cols(skill_thresholds);
  int K = Km1 + 1;
  matrix[N, K] etas;
  matrix[N, Km1] eta_increments;
  eta_increments = rep_matrix(loadings_questions[question_of_obs], Km1) .*
                   (rep_matrix(latent_factor_obs, Km1) -
                    skill_thresholds[question_of_obs,:]);
  etas[:,1] = rep_vector(0., N);
  etas[:,2:K] = eta_increments * cumsum_matrix;
  return etas;
}

// Partial sum function for categorical logit likelihood (for threading)
real partial_categorical_logit_lpmf(array[] int slice_y, int start, int end,
                                    matrix eta) {
  real lp = 0;
  int slice_idx = 1;
  for (n in start : end) {
    lp += categorical_logit_lpmf(slice_y[slice_idx] | eta[n]');
    slice_idx += 1;
  }
  return lp;
}

// Compute ordinal Brier score (cumulative Brier score) by question
// Returns a vector of length Q with the average ordinal Brier score for each question
// cumsum_matrix should be a K x (K-1) matrix where entry [i,j] = 1 if i <= j, else 0
// cum_obs_matrix should be a N x (K-1) matrix precomputed in transformed data
// question_start/end and sorted_indices enable efficient aggregation by question
vector pcm_get_ordered_brier_score(matrix ordered_prob_by_obs,
                                    array [] int question_of_obs,
                                    matrix cumsum_matrix,
                                    matrix cum_obs_matrix,
                                    array [] int question_start,
                                    array [] int question_end,
                                    array [] int sorted_indices)
{
  int start_idx;
  int end_idx;
  int N = rows(ordered_prob_by_obs);
  int Km1 = cols(cumsum_matrix);
  int Q = max(question_of_obs);
  vector[Q] score_by_question = rep_vector(0.0, Q);
  
  // Compute cumulative probabilities via matrix multiplication: N x Km1
  matrix[N, Km1] cum_pred = ordered_prob_by_obs * cumsum_matrix;
  
  // Compute squared differences and sum over categories for each observation
  matrix[N, Km1] squared_diffs = square(cum_pred - cum_obs_matrix);
  vector[N] score_by_obs = to_vector(squared_diffs * rep_vector(1.0, Km1)) / Km1;
  
  // Aggregate by question using precomputed indices
  for (q in 1:Q) {
    start_idx = question_start[q];
    end_idx = question_end[q];
    score_by_question[q] = sum(score_by_obs[sorted_indices[start_idx:end_idx]]) / (end_idx - start_idx + 1);
  }
  
  return score_by_question;
}

// Compute ordinal Brier score (cumulative Brier score) by question for credit model
// Returns a vector of length Q with the average ordinal Brier score for each question
// cumsum_matrix should be a K x (K-1) matrix where entry [i,j] = 1 if i <= j, else 0
// cum_obs_matrix should be a N x (K-1) matrix precomputed in transformed data
// question_start/end and sorted_indices enable efficient aggregation by question
vector cm_get_ordered_brier_score(matrix ordered_prob_by_obs,
                                   array [] int question_of_obs,
                                   matrix cumsum_matrix,
                                   matrix cum_obs_matrix,
                                   array [] int question_start,
                                   array [] int question_end,
                                   array [] int sorted_indices)
{
  int start_idx;
  int end_idx;
  int N = rows(ordered_prob_by_obs);
  int Km1 = cols(cumsum_matrix);
  int Q = max(question_of_obs);
  vector[Q] score_by_question = rep_vector(0.0, Q);
  
  // Compute cumulative probabilities via matrix multiplication: N x Km1
  matrix[N, Km1] cum_pred = ordered_prob_by_obs * cumsum_matrix;
  
  // Compute squared differences and sum over categories for each observation
  matrix[N, Km1] squared_diffs = square(cum_pred - cum_obs_matrix);
  vector[N] score_by_obs = to_vector(squared_diffs * rep_vector(1.0, Km1)) / Km1;
  
  // Aggregate by question using precomputed indices
  for (q in 1:Q) {
    start_idx = question_start[q];
    end_idx = question_end[q];
    score_by_question[q] = sum(score_by_obs[sorted_indices[start_idx:end_idx]]) / (end_idx - start_idx + 1);
  }
  
  return score_by_question;
}
