// Beta-Bernoulli model for the BinomialModel subclass in model_binomial.py.
// Single proportion p with a Beta(a, b) prior; observed y_i in {0, 1}.

data {
  int<lower=1> N;
  array[N] int<lower=0, upper=1> y;
  real<lower=0> a;
  real<lower=0> b;
}
transformed data {
  int<lower=0, upper=N> successes = sum(y);
}
parameters {
  real<lower=0, upper=1> p;
}
model {
  target += beta_lpdf(p | a, b);
  target += binomial_lpmf(successes | N, p);
}
generated quantities {
  array[N] int<lower=0, upper=1> y_rep;
  for (n in 1:N) {
    y_rep[n] = bernoulli_rng(p);
  }
}
