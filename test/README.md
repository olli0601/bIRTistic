# Test Files

This folder contains test files for validating model implementations.

## Files

### test_pcm_2cats_vs_ncats.Rmd

Validates that `partial_credit_model_2cats_v251224.stan` and `partial_credit_model_ncats_v260413.stan` produce identical `log_lik` values when run on the same Colombia data with the same seed.

**Purpose**: Ensure the flexible ncats implementation is numerically equivalent to the hardcoded 2cats version.

**What it does**:
- Loads Colombia data
- Prepares data in both formats (old 2cats and new ncats)
- Runs both models for 100 iterations with the same seed
- Compares log-likelihood values
- Visualizes any differences

**Expected result**: Maximum absolute difference in log_lik values should be < 1e-8 (numerical tolerance).

## Related Files

- `R/fit_partial_credit_model.R` - Original fitting function for 2cats model
- `R/fit_partial_credit_model_ncats.R` - New fitting function for ncats model with proper handling of long vector format for `ordered_prob_by_obs`
