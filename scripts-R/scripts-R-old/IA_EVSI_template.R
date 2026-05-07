# Bayesian EVSI / ENBS interim stopping workflow (Stan + R, pseudo-real code)
# File: EVSI_ENBS_interim_workflow.R
# Purpose: Template workflow to implement EVSI / ENBS-based interim stopping for a
# trial with a binary efficacy endpoint and binary toxicity endpoint.
# Technologies: cmdstanr (Stan), posterior, data.table, foreach/parallel
# NOTES:
#  - This is a template. Tweak priors, utility function, and computational knobs.
#  - For production use: add rigorous input checks, logging, seed control, and unit tests.

# ---- Dependencies ----
library(cmdstanr)    # for Stan fitting (use cmdstanr; ensure CmdStan installed)
library(posterior)   # for draws handling
library(data.table)
library(foreach)
library(doParallel)

# ---- Stan model (joint binary efficacy & toxicity) ----
# Simple logistic model where treatment effect on efficacy and toxicity are
# modeled with log-odds differences. Control arm log-odds are parameters.
stan_code <- 
  "data {
  int<lower=1> N;               // number of observed patients
  int<lower=0,upper=1> y_eff[N]; // efficacy outcomes (0/1)
  int<lower=0,upper=1> y_tox[N]; // toxicity outcomes (0/1)
  int<lower=0,upper=1> trt[N];   // treatment indicator (0=control,1=active)
}
parameters {
  real alpha_eff;          // log-odds baseline efficacy (control)
  real beta_eff;           // treatment log-odds effect on efficacy
  real alpha_tox;          // log-odds baseline toxicity (control)
  real beta_tox;           // treatment log-odds effect on toxicity
}
model {
  // Priors -- weakly informative; replace/elicitate as needed
  alpha_eff ~ normal(0, 2);
  beta_eff  ~ normal(0, 1);
  alpha_tox ~ normal(-2, 2);
  beta_tox  ~ normal(0, 1);

  for (i in 1:N) {
    real p_eff = inv_logit(alpha_eff + beta_eff * trt[i]);
    real p_tox = inv_logit(alpha_tox + beta_tox * trt[i]);
    y_eff[i] ~ bernoulli(p_eff);
    y_tox[i] ~ bernoulli(p_tox);
  }
}
generated quantities {
  // nothing here for now
}
"

# Write Stan file
stan_file <- file.path(tempdir(), "joint_binomial.stan")
writeLines(stan_code, con = stan_file)
mod <- cmdstan_model(stan_file)

# ---- Utility function ----
# Define utility U(action, theta). Here we compare two actions at interim:
#  - STOP_NOW: adopt best current estimate? For simplicity we define stopping to
#    mean 'no further sampling' and decision to deploy treatment or not is
#    determined by posterior expected utility.
#  - CONTINUE: run additional n_future patients per arm, incur sampling cost,
#    then decide at final analysis (adopt or not) based on posterior.

# We express utility in "QALY-equivalent" units for illustration.
# Parameters to elicit / set:
w_benefit <- 1.0    # QALY per successful efficacy outcome (scale)
w_harm    <- 0.5    # QALY lost per toxicity (serious AE)
unit_cost <- 1000   # monetary cost per patient (for ENBS; same units as utility?)
# In practice convert monetary cost to QALY-equivalent (or keep separate and subtract)

utility_from_outcomes <- function(n_eff_successes_trt, n_tox_trt,
                                  n_eff_successes_ctrl, n_tox_ctrl) {
  # net population-level utility difference (treatment minus control)
  ben_trt <- w_benefit * n_eff_successes_trt
  harm_trt <- w_harm * n_tox_trt
  ben_ctrl <- w_benefit * n_eff_successes_ctrl
  harm_ctrl <- w_harm * n_tox_ctrl
  # return net incremental utility of adopting treatment over control
  return((ben_trt - harm_trt) - (ben_ctrl - harm_ctrl))
}

# ---- Posterior sampling helper ----
fit_posterior <- function(data_list, chains = 4, iter_warmup = 1000, iter_sampling = 1000) {
  fit <- mod$sample(data = data_list,
                    chains = chains,
                    parallel_chains = min(chains, 4),
                    iter_warmup = iter_warmup,
                    iter_sampling = iter_sampling,
                    refresh = 0)
  return(fit)
}

# Extract posterior draws into a data.frame/matrix
extract_posterior_draws <- function(fit) {
  draws <- as_draws_df(fit$draws())
  # select relevant params
  draws <- as.data.table(draws)[, .(alpha_eff, beta_eff, alpha_tox, beta_tox)]
  return(draws)
}

# ---- Predictive simulation for future data ----
# Given a parameter draw theta, simulate future outcomes for n_future per arm
simulate_future <- function(theta, n_future_per_arm) {
  # theta: list or vector with named elements
  p_eff_ctrl <- plogis(theta$alpha_eff)
  p_eff_trt  <- plogis(theta$alpha_eff + theta$beta_eff)
  p_tox_ctrl <- plogis(theta$alpha_tox)
  p_tox_trt  <- plogis(theta$alpha_tox + theta$beta_tox)
  
  # simulate counts
  eff_ctrl <- rbinom(1, n_future_per_arm, p_eff_ctrl)
  eff_trt  <- rbinom(1, n_future_per_arm, p_eff_trt)
  tox_ctrl <- rbinom(1, n_future_per_arm, p_tox_ctrl)
  tox_trt  <- rbinom(1, n_future_per_arm, p_tox_trt)
  
  return(list(eff_ctrl = eff_ctrl, eff_trt = eff_trt,
              tox_ctrl = tox_ctrl, tox_trt = tox_trt,
              p_eff_ctrl = p_eff_ctrl, p_eff_trt = p_eff_trt,
              p_tox_ctrl = p_tox_ctrl, p_tox_trt = p_tox_trt))
}

# ---- Evaluate utility of "STOP_NOW" (adopt-or-not based on current posterior) ----
# For this template, assume if posterior mean incremental utility > 0 we adopt treatment
# and utility equals expected incremental utility (scaled to future population); else 0.

utility_stop_now <- function(posterior_draws, future_population = 1000) {
  # approximate posterior expected incremental utility if decision made now
  # compute expected probabilities under each draw, scale to population
  draws <- posterior_draws
  draws[, p_eff_ctrl := plogis(alpha_eff)]
  draws[, p_eff_trt  := plogis(alpha_eff + beta_eff)]
  draws[, p_tox_ctrl := plogis(alpha_tox)]
  draws[, p_tox_trt  := plogis(alpha_tox + beta_tox)]
  
  # expected counts in a future population if adopted
  draws[, inc_util := utility_from_outcomes(
    n_eff_successes_trt = future_population * p_eff_trt,
    n_tox_trt = future_population * p_tox_trt,
    n_eff_successes_ctrl = future_population * p_eff_ctrl,
    n_tox_ctrl = future_population * p_tox_ctrl)]
  
  # posterior expected incremental utility if adopting now
  E_inc_util <- mean(draws$inc_util)
  # decision: adopt if E_inc_util > 0
  return(max(E_inc_util, 0))
}

# ---- Evaluate expected utility of CONTINUE (one-step lookahead) ----
# Monte Carlo nested: for each posterior draw, simulate K possible future datasets,
# compute posterior after augmenting observed + simulated, then compute decision utility.

utility_continue_one_step <- function(posterior_draws, observed_data, n_future_per_arm = 50,
                                      K = 200, sample_cost_per_patient = unit_cost) {
  # posterior_draws: data.table of posterior draws
  M <- nrow(posterior_draws)
  # For each posterior draw (theta), we simulate K future datasets, compute the post-update
  # decision (adopt or not) and utility. Average across K to approximate expected utility
  # conditional on theta, then average across posterior draws.
  
  # set up parallel backend
  ncores <- max(1, parallel::detectCores() - 1)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  res_theta <- foreach(i = 1:M, .combine = c, .packages = c('data.table','cmdstanr','posterior')) %dopar% {
    theta <- as.list(posterior_draws[i])
    u_given_theta <- numeric(K)
    for (k in 1:K) {
      fut <- simulate_future(theta, n_future_per_arm)
      # build augmented dataset (observed + future)
      # For computational speed, instead of refitting full Stan for each simulated dataset,
      # you could use Laplace or normal approximation; here we show the full refit template (costly)
      
      # Construct data list for Stan
      # observed_data should include vectors y_eff, y_tox, trt for N observed
      N_obs <- observed_data$N
      y_eff_obs <- observed_data$y_eff
      y_tox_obs <- observed_data$y_tox
      trt_obs   <- observed_data$trt
      
      # expand with simulated future: half control, half treated
      n_new <- 2 * n_future_per_arm
      y_eff_new <- c(rep(0, n_new))
      y_tox_new <- c(rep(0, n_new))
      trt_new   <- c(rep(0, n_future_per_arm), rep(1, n_future_per_arm))
      
      # assign simulated counts to vectors
      # control successes
      y_eff_new[1:n_future_per_arm] <- rbinom(n_future_per_arm, 1, fut$p_eff_ctrl)
      y_eff_new[(n_future_per_arm+1):n_new] <- rbinom(n_future_per_arm, 1, fut$p_eff_trt)
      y_tox_new[1:n_future_per_arm] <- rbinom(n_future_per_arm, 1, fut$p_tox_ctrl)
      y_tox_new[(n_future_per_arm+1):n_new] <- rbinom(n_future_per_arm, 1, fut$p_tox_trt)
      
      y_eff_aug <- c(y_eff_obs, y_eff_new)
      y_tox_aug <- c(y_tox_obs, y_tox_new)
      trt_aug   <- c(trt_obs, trt_new)
      
      data_aug <- list(N = N_obs + n_new,
                       y_eff = y_eff_aug,
                       y_tox = y_tox_aug,
                       trt = trt_aug)
      
      # Refit posterior on augmented data (in practice use fast approx or reuse draws to update)
      fit_aug <- mod$sample(data = data_aug, chains = 2, iter_warmup = 400, iter_sampling = 400, refresh = 0)
      draws_aug <- as_draws_df(fit_aug$draws())
      draws_aug <- as.data.table(draws_aug)[, .(alpha_eff, beta_eff, alpha_tox, beta_tox)]
      
      # Decision after data: compute posterior expected incremental utility for future_population
      E_inc_util_post <- mean(sapply(1:nrow(draws_aug), function(j) {
        da <- draws_aug[j]
        p_eff_ctrl <- plogis(da$alpha_eff)
        p_eff_trt  <- plogis(da$alpha_eff + da$beta_eff)
        p_tox_ctrl <- plogis(da$alpha_tox)
        p_tox_trt  <- plogis(da$alpha_tox + da$beta_tox)
        utility_from_outcomes(future_population * p_eff_trt, future_population * p_tox_trt,
                              future_population * p_eff_ctrl, future_population * p_tox_ctrl)
      }))
      
      # If posterior expected incremental utility > 0 -> adopt, utility = that; else 0
      u_final_decision <- max(E_inc_util_post, 0)
      # subtract sampling cost
      total_sampling_cost <- n_new * sample_cost_per_patient
      u_given_theta[k] <- u_final_decision - total_sampling_cost
    }
    mean(u_given_theta) # expected utility conditional on theta
  }
  
  stopCluster(cl)
  # average over posterior draws
  EU_continue <- mean(res_theta)
  return(EU_continue)
}

# ---- EVSI and ENBS calculation ----
# EVSI = EU(continue) - EU(stop now)
# ENBS = EVSI - cost_of_sampling (where cost_of_sampling is monetary or QALY-equivalent)

compute_evsi_enbs <- function(posterior_draws, observed_data, n_future_per_arm = 50, K = 200,
                              sample_cost_per_patient = unit_cost, future_population = 1000) {
  EU_stop <- utility_stop_now(posterior_draws, future_population)
  EU_cont <- utility_continue_one_step(posterior_draws, observed_data, n_future_per_arm, K, sample_cost_per_patient)
  EVSI <- EU_cont - EU_stop
  # cost_of_sampling already subtracted inside EU_cont in this implementation; if not, subtract here
  ENBS <- EVSI # if costs included, this equals ENBS already
  return(list(EU_stop = EU_stop, EU_continue = EU_cont, EVSI = EVSI, ENBS = ENBS))
}

# ---- Calibration by simulation ----
# To calibrate decision thresholds or sanity-check operating characteristics, run many
# simulated trials under null and alternative scenarios, applying the EVSI/ENBS rule at each interim.

calibrate_policy <- function(sim_scenarios, n_sim = 1000, ...) {
  # sim_scenarios: list of scenario definitions (true params, accrual patterns, etc.)
  # For each simulated trial, run the interim decision process and record outcomes.
  # This is computationally heavy; parallelize.
  message("Calibration routine: run sims, collect type I error, power, average sample size, bias")
}

# ---- Example usage (pseudo) ----
# observed_data <- list(N = 100, y_eff = y_eff_vector, y_tox = y_tox_vector, trt = trt_vector)
# fit <- fit_posterior(observed_data)
# posterior_draws <- extract_posterior_draws(fit)
# evsi_res <- compute_evsi_enbs(posterior_draws, observed_data, n_future_per_arm = 50, K = 200)
# print(evsi_res)

# ---- Notes on practical improvements & speedups ----
# 1. Re-fitting Stan for every simulated future dataset (nested Monte Carlo) is very slow.
#    Use approximations: Laplace approximation, normal approximation to the posterior, or
#    use importance sampling to update posterior draws instead of refitting. See literature
#    on fast EVSI estimators (regression-based, moment-matching).
# 2. Emulators: train a regression / GP mapping from summary statistics of augmented data to
#    posterior expected utility; then replace inner refit with emulator predictions.
# 3. Pre-compute look-up tables across a grid of sufficient statistics and use interpolation
#    during real-time interim decisions.
# 4. For regulatory submissions provide full simulation code, seeds, and extreme-scenario checks.

# End of template
