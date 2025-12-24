matrix cm_get_cutpoints_noncentered(vector skill_thresholds_1,
                                    matrix skill_thresholds_incs,
                                    vector loadings_questions)
{
  int Q = size(skill_thresholds_1);
  int Km1 = cols(skill_thresholds_incs) + 1;
  matrix [Q, Km1] cutpoints;
  cutpoints = rep_matrix(skill_thresholds_1, Km1);
  for(j in 1:Q)
  {
    cutpoints[j,2:Km1] += cumulative_sum( skill_thresholds_incs[j,:] );
    cutpoints[j,] *= loadings_questions[j];
  }
  return cutpoints;
}

vector ol_get_etas(vector loadings_questions,
                   vector latent_factor_obs,
                   array [] int question_of_obs)
{
  return loadings_questions[question_of_obs] .* latent_factor_obs;
}

vector ol_get_etas(matrix loadings_questions,
                   matrix latent_factor_obs,
                   array [] int question_of_obs)
{
  return rows_dot_product(loadings_questions[question_of_obs, :], latent_factor_obs');
}

vector ol_get_etas(vector question_baselines,
                   vector loadings_questions,
                   vector latent_factor_obs,
                   array [] int question_of_obs)
{
  return question_baselines[question_of_obs] + 
         loadings_questions[question_of_obs] .* latent_factor_obs;
}

vector ol_get_etas(vector question_baselines,
                   matrix loadings_questions,
                   matrix latent_factor_obs,
                   array [] int question_of_obs)
{
  return question_baselines[question_of_obs] + 
         rows_dot_product(loadings_questions[question_of_obs, :], latent_factor_obs');
}

matrix ol_get_ordered_prob_by_obs(vector eta,
                                  matrix cut_points,
                                  array [] int question_of_obs)
{
  int N = num_elements(question_of_obs);
  int Km1 = cols(cut_points);
  int Q = rows(cut_points);  
  matrix[N,Km1 + 1] ordered_prob_by_obs;
  for(n in 1:N)
  {
    ordered_prob_by_obs[n,1]       = 1 - inv_logit( eta[n] - cut_points[question_of_obs[n],1] );
    ordered_prob_by_obs[n,2:(Km1)] = inv_logit( eta[n] - cut_points[question_of_obs[n],1:(Km1 - 1)] ) -
                                      inv_logit( eta[n] - cut_points[question_of_obs[n],2:Km1] );
    ordered_prob_by_obs[n,Km1 + 1] = inv_logit( eta[n] - cut_points[question_of_obs[n],Km1] );
  }
  return ordered_prob_by_obs;
}

real partial_ordered_logistic_lpmf(array[] int slice_y, int start, int end,
                                     vector eta, matrix cutpoints,
                                     array[] int question_of_obs) {
    real lp = 0;
    int slice_idx = 1;
    for (n in start : end) {
      lp += ordered_logistic_lpmf(slice_y[slice_idx] | eta[n],
              cutpoints[question_of_obs[n],  : ]');
      slice_idx += 1;
    }
    return lp;
  }