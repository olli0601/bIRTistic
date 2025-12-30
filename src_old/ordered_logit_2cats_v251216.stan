functions {
  #include ordered_logit_functions.stan
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
  real s2z_sd_questions;
  s2z_sd_unit = inv(sqrt(1. - inv(U)));
  s2z_sd_questions = inv(sqrt(1. - inv(Qcat1 + Qcat2)));
}
parameters {
  sum_to_zero_vector[U] latent_factor_unit;
  sum_to_zero_vector[Qcat1 + Qcat2] questions_baselines;
  vector[P] latent_factor_beta;
  vector[Qcat1] cat1_skill_thresholds_1;
  matrix<lower=0>[Qcat1, Kcat1 - 2] cat1_skill_thresholds_incs;
  vector<lower=0>[Qcat1 - 1] cat1_loadings_questions_m1;
  vector[Qcat2] cat2_skill_thresholds_1;
  matrix<lower=0>[Qcat2, Kcat2 - 2] cat2_skill_thresholds_incs;
  vector<lower=0>[Qcat2 - 1] cat2_loadings_questions_m1;
}
transformed parameters {
  vector[Ncat1] cat1_eta;
  vector[Ncat2] cat2_eta;
  matrix[Qcat1, Kcat1 - 1] cat1_cutpoints;
  matrix[Qcat2, Kcat2 - 1] cat2_cutpoints;
  cat1_cutpoints = cm_get_cutpoints_noncentered(cat1_skill_thresholds_1,
                     cat1_skill_thresholds_incs,
                     append_row(1.0, cat1_loadings_questions_m1));
  cat2_cutpoints = cm_get_cutpoints_noncentered(cat2_skill_thresholds_1,
                     cat2_skill_thresholds_incs,
                     append_row(1.0, cat2_loadings_questions_m1));
  cat1_eta = ol_get_etas(questions_baselines[1 : Qcat1],
                         append_row(1.0, cat1_loadings_questions_m1),
                         latent_factor_unit[cat1_unit_of_obs]
                         + cat1_X * latent_factor_beta, cat1_question_of_obs);
  cat2_eta = ol_get_etas(questions_baselines[(Qcat1 + 1) : (Qcat1 + Qcat2)],
                         append_row(1.0, cat2_loadings_questions_m1),
                         latent_factor_unit[cat2_unit_of_obs]
                         + cat2_X * latent_factor_beta, cat2_question_of_obs);
}
model {
  // likelihood under credit model
  //print("\n cat1_cutpoints", cat1_cutpoints);
  for (n in 1 : Ncat1) {
    target += ordered_logistic_lupmf(cat1_y[n] | cat1_eta[n],
                cat1_cutpoints[cat1_question_of_obs[n],  : ]');
  }
  //print("\n cat2_cutpoints", cat1_cutpoints);
  for (n in 1 : Ncat2) {
    target += ordered_logistic_lupmf(cat2_y[n] | cat2_eta[n],
                cat2_cutpoints[cat2_question_of_obs[n],  : ]');
  }
  // priors for latent factors
  target += normal_lupdf(latent_factor_unit | 0, s2z_sd_unit);
  target += std_normal_lupdf(latent_factor_beta);
  // priors for skill thresholds
  target += normal_lupdf(cat1_skill_thresholds_1 | 0, 3.5);
  target += normal_lupdf(to_vector(cat1_skill_thresholds_incs) | 0, 2.5);
  target += normal_lupdf(cat2_skill_thresholds_1 | 0, 3.5);
  target += normal_lupdf(to_vector(cat2_skill_thresholds_incs) | 0, 2.5);
  // priors for question baselines
  target += normal_lupdf(questions_baselines | 0, s2z_sd_questions);
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
  cat1_ordered_prob_by_obs = ol_get_ordered_prob_by_obs(cat1_eta,
                               cat1_cutpoints, cat1_question_of_obs);
  cat2_ordered_prob_by_obs = ol_get_ordered_prob_by_obs(cat2_eta,
                               cat2_cutpoints, cat2_question_of_obs);
  for (n in 1 : Ncat1) {
    cat1_ypred[n] = ordered_logistic_rng(cat1_eta[n],
                      cat1_cutpoints[cat1_question_of_obs[n],  : ]');
  }
  for (n in 1 : Ncat2) {
    cat2_ypred[n] = ordered_logistic_rng(cat2_eta[n],
                      cat2_cutpoints[cat2_question_of_obs[n],  : ]');
  }
  for (n in 1 : Ncat1) {
    log_lik[n] = ordered_logistic_lpmf(cat1_y[n] | cat1_eta[n],
                   cat1_cutpoints[cat1_question_of_obs[n],  : ]');
  }
  for (n in 1 : Ncat2) {
    log_lik[Ncat1 + n] = ordered_logistic_lpmf(cat2_y[n] | cat2_eta[n],
                           cat2_cutpoints[cat2_question_of_obs[n],  : ]');
  }
}
