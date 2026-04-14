# Tests for bIRTistic

This directory contains tests for the Stan models and R code in the bIRTistic project.

## Unit Tests

### test_pcm_ncats_correctness.R

**Purpose:** Verify computational correctness of `partial_credit_model_ncats_v260413.stan`

**What it tests:**
- Both `partial_credit_model_2cats_v251224.stan` (reference) and `partial_credit_model_ncats_v260413.stan` (new flexible version) produce **identical** log-likelihood values when run with identical fixed parameter values
- Uses simple synthetic data with 2 category types, 2 questions per type, 4 total observations
- Runs both models with `fixed_param=TRUE` to evaluate computational logic without MCMC sampling variation

**How to run:**
```r
# Using testthat
library(testthat)
test_file("test/test_pcm_ncats_correctness.R")

# Or directly
Rscript test/test_pcm_ncats_correctness.R
```

**Expected output:**
```
══ Testing test_pcm_ncats_correctness.R ═══
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 7 ] Done!
```

**Framework:** Uses `testthat` package for structured assertions and better diagnostics

**Rationale:**
This test validates that the ncats model's more complex implementation (concatenated arrays, long vector storage, dynamic indexing) produces mathematically equivalent results to the simpler hardcoded 2cats model. The test uses fixed parameters to eliminate MCMC stochasticity and focus purely on computational correctness.

---

## Test Strategy

The unit test (`test_pcm_ncats_correctness.R`) provides fast, deterministic verification of computational correctness. Run it after any code changes to the Stan models to ensure numerical equivalence is maintained.

### Running All Tests

To run all unit tests at once:
```r
# Option 1: Using the helper script
Rscript test/run_all_tests.R

# Option 2: Using testthat directly
library(testthat)
test_dir("test")
```

### Individual Test Execution

See individual test sections above for how to run specific tests.

---

## Related Files

- `R/fit_partial_credit_model.R` - Original fitting function for 2cats model
- `R/fit_partial_credit_model_ncats.R` - New fitting function for ncats model with proper handling of long vector format for `ordered_prob_by_obs`
