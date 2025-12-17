functions
{
  vector lf_get_eta( array [] vector beta_questions,
                     array [] vector loadings_questions,
                     array [] int group_of_obs,
                     array [] int question_of_obs,
                     array [] int unit_of_obs,
                     vector latent_factor_unit
                   )
  {
    int N = num_elements(unit_of_obs);
    vector [N] eta;
    eta = latent_factor_unit[unit_of_obs];
    for(n in 1:N)
    {
      eta[n]
        =
        beta_questions[group_of_obs[n]][question_of_obs[n]] +
        loadings_questions[group_of_obs[n]][question_of_obs[n]] * eta[n];
    }
    return eta;
  }
}
data
{
  int<lower=1> Ncat; // number of ordered categorial observations
  int<lower=2> Qcat; // number of ordered categorial questions
  int<lower=3> Kcat; // number of categories
  int<lower=1> Ncat2;
  int<lower=2> Qcat2;
  int<lower=3> Kcat2;
  int<lower=2> P; // number of groups
  int<lower=1> U; // number of units
  array [Ncat] int<lower=1, upper=Kcat> ycat; // observations
  array [Ncat] int<lower=1,upper=U> unit_of_cat_obs;
  array [Ncat] int<lower=1,upper=P> group_of_cat_obs;
  array [Ncat] int<lower=1,upper=Qcat> question_of_cat_obs;
  array [Ncat2] int<lower=1, upper=Kcat2> ycat2; // observations
  array [Ncat2] int<lower=1,upper=U> unit_of_cat2_obs;
  array [Ncat2] int<lower=1,upper=P> group_of_cat2_obs;
  array [Ncat2] int<lower=1,upper=Qcat2> question_of_cat2_obs;
}
transformed data
{
    real s2z_sd_unit;
    real s2z_sd_cat_questions;
    real s2z_sd_cat2_questions;
    s2z_sd_unit = inv(sqrt(1. - inv(U)));
    s2z_sd_cat_questions = inv(sqrt(1. - inv(Qcat)));
    s2z_sd_cat2_questions = inv(sqrt(1. - inv(Qcat2)));
}
parameters
{
  sum_to_zero_vector[U] latent_factor_unit;
  array [P] sum_to_zero_vector[Qcat] beta_cat_questions;
  array [P] vector<lower=0>[Qcat-1] loadings_cat_questions_m1;
  real cat_cut_point_1;
  matrix<lower=0>[Kcat - 2,P] cat_cut_point_gaps_groups;
  array [P] sum_to_zero_vector[Qcat2] beta_cat2_questions;
  array [P] vector<lower=0>[Qcat2-1] loadings_cat2_questions_m1;
  real cat2_cut_point_1;
  matrix<lower=0>[Kcat2 - 2,P] cat2_cut_point_gaps_groups;
}
transformed parameters
{
  array [P] ordered[Kcat-1] cat_cut_points;
  array [P] vector<lower=0>[Qcat] loadings_cat_questions;
  array [P] ordered[Kcat2-1] cat2_cut_points;
  array [P] vector<lower=0>[Qcat2] loadings_cat2_questions;
  vector [Ncat] cat_eta;
  vector [Ncat2] cat2_eta;

  for (p in 1:P)
  {
    cat_cut_points[p] =
      append_row(
        rep_vector(cat_cut_point_1, 1),
        rep_vector(cat_cut_point_1, Kcat-2) + cumulative_sum(cat_cut_point_gaps_groups[,p])
        );
    loadings_cat_questions[p] = append_row(1.0, loadings_cat_questions_m1[p]);
    cat2_cut_points[p] =
      append_row(
        rep_vector(cat2_cut_point_1, 1),
        rep_vector(cat2_cut_point_1, Kcat2-2) + cumulative_sum(cat2_cut_point_gaps_groups[,p])
        );
    loadings_cat2_questions[p] = append_row(1.0, loadings_cat2_questions_m1[p]);
  }
  cat_eta
        =
        lf_get_eta( beta_cat_questions,
                    loadings_cat_questions,
                    group_of_cat_obs,
                    question_of_cat_obs,
                    unit_of_cat_obs,
                    latent_factor_unit
                    );
  cat2_eta
        =
        lf_get_eta( beta_cat2_questions,
                    loadings_cat2_questions,
                    group_of_cat2_obs,
                    question_of_cat2_obs,
                    unit_of_cat2_obs,
                    latent_factor_unit
                    );
}
model
{
  // ordered cat likelihood, cannot be vectorized, for cat obs
  for (n in 1:Ncat)
  {
    ycat[n] ~ ordered_logistic(cat_eta[n], cat_cut_points[group_of_cat_obs[n]]);
  }
  for (n in 1:Ncat2)
  {
    ycat2[n] ~ ordered_logistic(cat2_eta[n], cat2_cut_points[group_of_cat2_obs[n]]);
  }

  // priors for latent factors
  latent_factor_unit ~ normal(0, s2z_sd_unit);

  // priors for cat loadings
  cat_cut_point_1 ~ normal(0, 3.5);
  to_vector(cat_cut_point_gaps_groups) ~ normal(0, 2.5);
  for (p in 1:P)
  {
    beta_cat_questions[p] ~ normal(0, s2z_sd_cat_questions);
    loadings_cat_questions_m1[p] ~ lognormal(0,1);
  }
  cat2_cut_point_1 ~ normal(0, 3.5);
  to_vector(cat2_cut_point_gaps_groups) ~ normal(0, 2.5);
  for (p in 1:P)
  {
    beta_cat2_questions[p] ~ normal(0, s2z_sd_cat2_questions);
    loadings_cat2_questions_m1[p] ~ lognormal(0,1);
  }
}
generated quantities
{
  array [Kcat,Ncat] real<lower=0, upper=1> cat_ordered_prob_by_obs;
  array [Kcat2,Ncat2] real<lower=0, upper=1> cat2_ordered_prob_by_obs;
  array [Ncat] int<lower=0> cat_ypred;
  array [Ncat2] int<lower=0> cat2_ypred;
  array [Ncat+Ncat2] real log_lik;

  for (n in 1:Ncat)
  {
    cat_ordered_prob_by_obs[1,n]          = 1 - inv_logit( cat_eta[n] - cat_cut_points[group_of_cat_obs[n]][1] );
    cat_ordered_prob_by_obs[2:(Kcat-1),n] = to_array_1d( inv_logit( cat_eta[n] - cat_cut_points[group_of_cat_obs[n]][1:(Kcat-2)] )
                                                         -
                                                         inv_logit( cat_eta[n] - cat_cut_points[group_of_cat_obs[n]][2:(Kcat-1)] )
                                                         );
    cat_ordered_prob_by_obs[Kcat,n]       = inv_logit( cat_eta[n] - cat_cut_points[group_of_cat_obs[n]][Kcat-1] );
  }
  for (n in 1:Ncat2)
  {
    cat2_ordered_prob_by_obs[1,n]           = 1 - inv_logit( cat2_eta[n] - cat2_cut_points[group_of_cat2_obs[n]][1] );
    cat2_ordered_prob_by_obs[2:(Kcat2-1),n] = to_array_1d( inv_logit( cat2_eta[n] - cat2_cut_points[group_of_cat2_obs[n]][1:(Kcat2-2)] )
                                                           -
                                                           inv_logit( cat2_eta[n] - cat2_cut_points[group_of_cat2_obs[n]][2:(Kcat2-1)] )
                                                           );
    cat2_ordered_prob_by_obs[Kcat2,n]       = inv_logit( cat2_eta[n] - cat2_cut_points[group_of_cat2_obs[n]][Kcat2-1] );
  }
  for( n in 1:Ncat )
  {
    cat_ypred[n] = ordered_logistic_rng( cat_eta[n], cat_cut_points[group_of_cat_obs[n]] );
  }
  for( n in 1:Ncat2 )
  {
    cat2_ypred[n] = ordered_logistic_rng( cat2_eta[n], cat2_cut_points[group_of_cat2_obs[n]] );
  }
  for(n in 1:Ncat)
  {
    log_lik[n] =  ordered_logistic_lpmf( ycat[n] | cat_eta[n], cat_cut_points[group_of_cat_obs[n]] );
  }
  for(n in 1:Ncat2)
  {
    log_lik[ Ncat+n ] = ordered_logistic_lpmf( ycat2[n] | cat2_eta[n], cat2_cut_points[group_of_cat2_obs[n]] );
  }
}
