"""
Test that NumPyro partial_credit_model_2cats produces same log likelihood as Stan version.

This test compares log likelihood values between:
- Stan: test/stan/partial_credit_model_2cats_v251224.stan
- NumPyro: src/numpyro/partial_credit_model_2cats_v251224.pyro

Expected: max absolute difference < 1e-6 (accounting for numerical precision differences)

Date: 2026-05-08
"""

import pytest
import numpy as np
import jax
import jax.numpy as jnp
from cmdstanpy import CmdStanModel
import sys
from pathlib import Path

# Add src/numpyro to path
project_root = Path(__file__).parent.parent.parent
numpyro_path = str(project_root / 'src' / 'numpyro')
if numpyro_path not in sys.path:
    sys.path.insert(0, numpyro_path)

# Import the NumPyro model
import importlib.util
from importlib.machinery import SourceFileLoader

pyro_file = project_root / 'src' / 'numpyro' / 'partial_credit_model_2cats_v251224.pyro'
loader = SourceFileLoader("pcm_2cats", str(pyro_file))
spec = importlib.util.spec_from_loader("pcm_2cats", loader)
pcm_2cats = importlib.util.module_from_spec(spec)
loader.exec_module(pcm_2cats)


def test_numpyro_vs_stan_log_likelihood():
    """Test that NumPyro and Stan produce identical log likelihoods with fixed parameters."""
    
    # Set random seed for reproducibility
    np.random.seed(42)
    jax.config.update('jax_platform_name', 'cpu')
    
    # Create test data matching the R test
    stan_data = {
        'U': 2,
        'P': 1,
        'Ncat1': 2,
        'Qcat1': 2,
        'Kcat1': 3,
        'cat1_y': np.array([1, 2], dtype=np.int32),
        'cat1_question_of_obs': np.array([1, 2], dtype=np.int32),
        'cat1_unit_of_obs': np.array([1, 2], dtype=np.int32),
        'cat1_X': np.array([[0.5], [-0.3]]),
        'Ncat2': 2,
        'Qcat2': 2,
        'Kcat2': 4,
        'cat2_y': np.array([2, 3], dtype=np.int32),
        'cat2_question_of_obs': np.array([1, 2], dtype=np.int32),
        'cat2_unit_of_obs': np.array([1, 2], dtype=np.int32),
        'cat2_X': np.array([[0.1], [-0.2]])
    }
    
    # Fixed parameter values matching R test
    stan_params = {
        'latent_factor_unit': np.array([-0.5, 0.5]),
        'latent_factor_beta': np.array([1.0]),
        'cat1_skill_thresholds': np.array([
            [0.5, 1.0],
            [-0.5, 0.5]
        ]),
        'cat1_loadings_questions_m1': np.array([1.5]),
        'cat2_skill_thresholds': np.array([
            [0.3, 0.6, 0.9],
            [-0.3, -0.6, -0.9]
        ]),
        'cat2_loadings_questions_m1': np.array([0.8])
    }
    
    # =========================================================================
    # Run Stan model with fixed parameters
    # =========================================================================
    
    print("\nCompiling and running Stan model...")
    stan_model = CmdStanModel(
        stan_file=str(project_root / 'test' / 'stan' / 'partial_credit_model_2cats_v251224.stan'),
        stanc_options={'include-paths': str(project_root / 'src' / 'stan')}
    )
    
    # Run with fixed parameters (0 warmup, 1 sample, fixed_param=True)
    stan_fit = stan_model.sample(
        data=stan_data,
        inits=stan_params,
        chains=1,
        iter_warmup=0,
        iter_sampling=1,
        fixed_param=True,
        adapt_engaged=False,
        show_console=False
    )
    
    # Extract Stan log likelihoods
    stan_log_lik = stan_fit.stan_variable('log_lik').flatten()
    
    print(f"Stan log_lik shape: {stan_log_lik.shape}")
    print(f"Stan log_lik values:\n{stan_log_lik}")
    
    # =========================================================================
    # Run NumPyro model with same parameters
    # =========================================================================
    
    print("\nComputing NumPyro log likelihood...")
    
    # Convert data to JAX arrays
    jax_data = {k: jnp.array(v) if isinstance(v, np.ndarray) else v 
                for k, v in stan_data.items()}
    
    # Convert parameters to JAX arrays
    jax_params = {
        'latent_factor_unit': jnp.array(stan_params['latent_factor_unit']),
        'latent_factor_beta': jnp.array(stan_params['latent_factor_beta']),
        'cat1_skill_thresholds': jnp.array(stan_params['cat1_skill_thresholds']),
        'cat2_skill_thresholds': jnp.array(stan_params['cat2_skill_thresholds']),
        'cat1_loadings_questions_m1': jnp.array(stan_params['cat1_loadings_questions_m1']),
        'cat2_loadings_questions_m1': jnp.array(stan_params['cat2_loadings_questions_m1'])
    }
    
    # Compute log likelihood using NumPyro model
    numpyro_log_lik = pcm_2cats.get_log_likelihood(jax_data, jax_params)
    numpyro_log_lik = np.array(numpyro_log_lik)
    
    print(f"NumPyro log_lik shape: {numpyro_log_lik.shape}")
    print(f"NumPyro log_lik values:\n{numpyro_log_lik}")
    
    # =========================================================================
    # Compare results
    # =========================================================================
    
    print("\n" + "="*70)
    print("COMPARISON RESULTS")
    print("="*70)
    
    # Compute differences
    abs_diff = np.abs(stan_log_lik - numpyro_log_lik)
    max_abs_diff = np.max(abs_diff)
    mean_abs_diff = np.mean(abs_diff)
    
    print(f"\nNumber of observations: {len(stan_log_lik)}")
    print(f"Max absolute difference: {max_abs_diff:.2e}")
    print(f"Mean absolute difference: {mean_abs_diff:.2e}")
    
    print(f"\nPer-observation differences:")
    for i, (s, n, d) in enumerate(zip(stan_log_lik, numpyro_log_lik, abs_diff)):
        print(f"  Obs {i+1}: Stan={s:.8f}, NumPyro={n:.8f}, diff={d:.2e}")
    
    # Test assertions
    # Allow for small numerical differences due to different implementations
    tolerance = 1e-6
    
    assert len(stan_log_lik) == len(numpyro_log_lik), \
        f"Length mismatch: Stan={len(stan_log_lik)}, NumPyro={len(numpyro_log_lik)}"
    
    assert max_abs_diff < tolerance, \
        f"Max absolute difference {max_abs_diff:.2e} exceeds tolerance {tolerance:.2e}"
    
    print(f"\n✓ TEST PASSED: Max difference {max_abs_diff:.2e} < tolerance {tolerance:.2e}")
    print("="*70)


if __name__ == '__main__':
    # Run tests
    print("="*70)
    print("TESTING NUMPYRO PARTIAL CREDIT MODEL 2CATS")
    print("="*70)

    print("\n1. Testing log likelihood vs Stan...")
    test_numpyro_vs_stan_log_likelihood()
    
    print("\n" + "="*70)
    print("ALL TESTS PASSED!")
    print("="*70)
