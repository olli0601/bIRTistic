// Beta-Bernoulli model for the BinomialModel subclass in model_binomial.py.
// Single proportion p with a Beta(a, b) prior; observed y_i in {0, 1}.

data {
  int<lower=1> N;
  array[N] int<lower=0, upper=1> y;
  real<lower=0> a;
  real<lower=0> b;
}
parameters {
  real<lower=0, upper=1> p;
}
model {
  p ~ beta(a, b);
  y ~ bernoulli(p);
}
