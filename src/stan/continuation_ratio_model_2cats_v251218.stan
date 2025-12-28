functions {
  // Continuation ratio log probability mass function with log link
  // P(Y = k) computed from sequential conditional probabilities
  real continuation_ratio_lpmf(int y, real eta, vector cutpoints) {
    int K = num_elements(cutpoints) + 1;
    real lp = 0;
    
    // Continuation ratio with log link:
    // P(Y > k | Y >= k) = exp(-(eta + cutpoints[k]))
    // P(Y = k | Y >= k) = 1 - exp(-(eta + cutpoints[k]))
    
    // For categories 1 to y-1: must continue (Y > k)
    for (k in 1 : (y - 1)) {
      lp += -(eta + cutpoints[k]); // log P(Y > k | Y >= k)
    }
    
    // For category y: stop here
    if (y < K) {
      lp += log1m_exp(-(eta + cutpoints[y])); // log P(Y = y | Y >= y)
    } else {
      // For last category K, must have continued through all K-1
      // No additional term needed as it's already accumulated
    }
    
    return lp;
  }
  
  // Generate random sample from continuation ratio model
  int continuation_ratio_rng(real eta, vector cutpoints) {
    int K = num_elements(cutpoints) + 1;
    
    for (k in 1 : (K - 1)) {
      real p_continue = exp(-(eta + cutpoints[k]));
      if (bernoulli_rng(1 - p_continue)) {
        return k;
      }
    }
    return K;
  }
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
  matrix[Qcat1, Kcat1 - 1] cat1_cutpoints;
  matrix[Qcat2, Kcat2 - 1] cat2_cutpoints;
}
transformed parameters {
  vector[Ncat1] cat1_eta;
  vector[Ncat2] cat2_eta;
  
  cat1_eta = latent_factor_unit[cat1_unit_of_obs]
             + cat1_X * latent_factor_beta;
  cat2_eta = latent_factor_unit[cat2_unit_of_obs]
             + cat2_X * latent_factor_beta;
}
model {
  // likelihood under continuation ratio model
  for (n in 1 : Ncat1) {
    target += continuation_ratio_lpmf(cat1_y[n] | cat1_eta[n],
                cat1_cutpoints[cat1_question_of_obs[n],  : ]');
  }
  for (n in 1 : Ncat2) {
    target += continuation_ratio_lpmf(cat2_y[n] | cat2_eta[n],
                cat2_cutpoints[cat2_question_of_obs[n],  : ]');
  }
  
  // priors for latent factors
  target += normal_lupdf(latent_factor_unit | 0, s2z_sd_unit);
  target += std_normal_lupdf(latent_factor_beta);
  
  // priors for cutpoints
  target += normal_lupdf(to_vector(cat1_cutpoints) | 0, 3.5);
  target += normal_lupdf(to_vector(cat2_cutpoints) | 0, 3.5);
}
generated quantities {
  array[Ncat1] int<lower=0> cat1_ypred;
  array[Ncat2] int<lower=0> cat2_ypred;
  array[Ncat1 + Ncat2] real log_lik;
  
  for (n in 1 : Ncat1) {
    cat1_ypred[n] = continuation_ratio_rng(cat1_eta[n],
                      cat1_cutpoints[cat1_question_of_obs[n],  : ]');
  }
  for (n in 1 : Ncat2) {
    cat2_ypred[n] = continuation_ratio_rng(cat2_eta[n],
                      cat2_cutpoints[cat2_question_of_obs[n],  : ]');
  }
  for (n in 1 : Ncat1) {
    log_lik[n] = continuation_ratio_lpmf(cat1_y[n] | cat1_eta[n],
                   cat1_cutpoints[cat1_question_of_obs[n],  : ]');
  }
  for (n in 1 : Ncat2) {
    log_lik[Ncat1 + n] = continuation_ratio_lpmf(cat2_y[n] | cat2_eta[n],
                           cat2_cutpoints[cat2_question_of_obs[n],  : ]');
  }
}
