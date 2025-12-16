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
  s2z_sd_unit = inv(sqrt(1. - inv(U)));
}
parameters {
  sum_to_zero_vector[U] latent_factor_unit;
  vector[P] latent_factor_beta;
  
  vector[Qcat1] cat1_skill_thresholds_1;
  matrix[Qcat1, Kcat1 - 2] cat1_skill_thresholds_incs;
  vector<lower=0>[Qcat1] cat1_loadings_questions;
  
  vector[Qcat2] cat2_skill_thresholds_1;
  matrix[Qcat2, Kcat2 - 2] cat2_skill_thresholds_incs;
  vector<lower=0>[Qcat2] cat2_loadings_questions;
}
transformed parameters {
  matrix[Ncat1, Kcat1] cat1_eta;
  matrix[Ncat2, Kcat2] cat2_eta;
  
  {
    matrix[Qcat1, Kcat1 - 1] cat1_skill_thresholds;
    matrix[Qcat2, Kcat2 - 1] cat2_skill_thresholds;
    
    cat1_skill_thresholds = cm_get_skill_thresholds(cat1_skill_thresholds_1,
                              cat1_skill_thresholds_incs);
    cat2_skill_thresholds = cm_get_skill_thresholds(cat2_skill_thresholds_1,
                              cat2_skill_thresholds_incs);
    
    cat1_eta = cm_get_etas(cat1_loadings_questions, cat1_skill_thresholds,
                           latent_factor_unit[cat1_unit_of_obs]
                           + cat1_X * latent_factor_beta,
                           cat1_question_of_obs, cat1_unit_of_obs);
    cat2_eta = cm_get_etas(cat2_loadings_questions, cat2_skill_thresholds,
                           latent_factor_unit[cat2_unit_of_obs]
                           + cat2_X * latent_factor_beta,
                           cat2_question_of_obs, cat2_unit_of_obs);
  }
}
model {
  // likelihood under credit model
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
  target += normal_lupdf(cat1_skill_thresholds_1 | 0, 3.5);
  target += normal_lupdf(to_vector(cat1_skill_thresholds_incs) | 0, 2.5);
  target += normal_lupdf(cat2_skill_thresholds_1 | 0, 3.5);
  target += normal_lupdf(to_vector(cat2_skill_thresholds_incs) | 0, 2.5);
  
  // priors for loadings
  target += student_t_lupdf(cat1_loadings_questions | 3, 0, 1);
  target += student_t_lupdf(cat2_loadings_questions | 3, 0, 1);
}
generated quantities {
  matrix[Ncat1, Kcat1] cat1_ordered_prob_by_obs;
  matrix[Ncat2, Kcat2] cat2_ordered_prob_by_obs;
  array[Ncat1] int<lower=0> cat1_ypred;
  array[Ncat2] int<lower=0> cat2_ypred;
  array[Ncat1 + Ncat2] real log_lik;
  
  for (n in 1 : Ncat1) {
    cat1_ordered_prob_by_obs[n,  : ] = softmax(cat1_eta[n]')';
    cat1_ypred[n] = categorical_logit_rng(cat1_eta[n]');
  }
  for (n in 1 : Ncat2) {
    cat2_ordered_prob_by_obs[n,  : ] = softmax(cat2_eta[n]')';
    cat2_ypred[n] = categorical_logit_rng(cat2_eta[n]');
  }
  for (n in 1 : Ncat1) {
    log_lik[n] = categorical_logit_lpmf(cat1_y[n] | cat1_eta[n]');
  }
  for (n in 1 : Ncat2) {
    log_lik[Ncat1 + n] = categorical_logit_lpmf(cat2_y[n] | cat2_eta[n]');
  }
}
