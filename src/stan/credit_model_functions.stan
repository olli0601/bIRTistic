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
  for(n in 1:N)
  {
    etas[n,2:K] = loadings_questions[question_of_obs[n]] *
                  (latent_factor_obs[n] -
                   skill_thresholds[question_of_obs[n],:]);
  }
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
  for(n in 1:N)
  {
    etas[n,2:K] = loadings_questions[question_of_obs[n], :] .*
                  (latent_factor_obs[n] -
                   skill_thresholds[question_of_obs[n],:]);
  }
  return etas;
}
