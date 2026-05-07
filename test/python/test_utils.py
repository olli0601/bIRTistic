"""Test the testing utilities themselves."""
import pandas as pd
import numpy as np
import pytest
from utils.comparison import (
    compare_dataframes,
    compare_scalars,
    compare_posteriors,
    compare_arrays
)


class TestCompareDataFrames:
    """Tests for DataFrame comparison utility."""
    
    def test_identical_dataframes(self):
        """Test that identical DataFrames are detected as equal."""
        df = pd.DataFrame({
            'a': [1, 2, 3],
            'b': [4.0, 5.0, 6.0],
            'c': ['x', 'y', 'z']
        })
        matches, msg = compare_dataframes(df, df.copy())
        assert matches, msg
    
    def test_different_shape(self):
        """Test that different shapes are detected."""
        df1 = pd.DataFrame({'a': [1, 2, 3]})
        df2 = pd.DataFrame({'a': [1, 2]})
        matches, msg = compare_dataframes(df1, df2)
        assert not matches
        assert "Shape mismatch" in msg
    
    def test_missing_column(self):
        """Test that missing columns are detected."""
        df1 = pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]})
        df2 = pd.DataFrame({'a': [1, 2, 3]})
        matches, msg = compare_dataframes(df1, df2)
        assert not matches
        # Shape check catches this first (different number of columns)
        assert "mismatch" in msg.lower()
    
    def test_numerical_tolerance(self):
        """Test numerical comparison with tolerance."""
        df1 = pd.DataFrame({'a': [1.0, 2.0, 3.0]})
        df2 = pd.DataFrame({'a': [1.0001, 2.0001, 3.0001]})
        # Should fail with tight tolerance
        matches, msg = compare_dataframes(df1, df2, rtol=1e-5, atol=1e-8)
        assert not matches
        # Should pass with looser tolerance
        matches, msg = compare_dataframes(df1, df2, rtol=1e-3, atol=1e-3)
        assert matches
    
    def test_column_order_ignored(self):
        """Test that column order can be ignored."""
        df1 = pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]})
        df2 = pd.DataFrame({'b': [4, 5, 6], 'a': [1, 2, 3]})
        # Should pass when order doesn't matter
        matches, msg = compare_dataframes(df1, df2, check_column_order=False)
        assert matches


class TestCompareScalars:
    """Tests for scalar comparison utility."""
    
    def test_identical_values(self):
        """Test that identical values match."""
        matches, msg = compare_scalars(1.0, 1.0)
        assert matches
    
    def test_within_tolerance(self):
        """Test scalar comparison with tolerance."""
        matches, msg = compare_scalars(1.0, 1.001, rtol=0.01)
        assert matches
    
    def test_outside_tolerance(self):
        """Test that values outside tolerance are detected."""
        matches, msg = compare_scalars(1.0, 1.1, rtol=0.01)
        assert not matches
        assert "mismatch" in msg
    
    def test_named_value(self):
        """Test that value name appears in message."""
        matches, msg = compare_scalars(1.0, 2.0, name="test_param")
        assert not matches
        assert "test_param" in msg


class TestComparePosteriors:
    """Tests for posterior comparison utility."""
    
    def test_matching_posteriors(self):
        """Test posterior comparison with user-specified tolerances."""
        posterior_r = {
            'means': {'param1': 1.0, 'param2': 2.0},
            'rhat': {'param1': 1.01, 'param2': 1.02}
        }
        # Within 3% for means, within 0.005 for Rhat
        posterior_py = {
            'means': {'param1': 1.02, 'param2': 2.04},  # 2% and 2% diff
            'rhat': {'param1': 1.012, 'param2': 1.024}  # 0.002 and 0.004 diff
        }
        matches, msg = compare_posteriors(posterior_py, posterior_r)
        assert matches, msg
    
    def test_means_exceed_tolerance(self):
        """Test that means outside 3% tolerance are detected."""
        posterior_r = {
            'means': {'param1': 1.0},
            'rhat': {'param1': 1.01}
        }
        posterior_py = {
            'means': {'param1': 1.05},  # 5% difference
            'rhat': {'param1': 1.01}
        }
        matches, msg = compare_posteriors(
            posterior_py, posterior_r, mean_rtol=0.03
        )
        assert not matches
        assert "param1 mean" in msg
    
    def test_rhat_exceed_tolerance(self):
        """Test that Rhat outside 0.005 tolerance are detected."""
        posterior_r = {
            'means': {'param1': 1.0},
            'rhat': {'param1': 1.01}
        }
        posterior_py = {
            'means': {'param1': 1.0},
            'rhat': {'param1': 1.02}  # 0.01 difference
        }
        matches, msg = compare_posteriors(
            posterior_py, posterior_r, rhat_atol=0.005
        )
        assert not matches
        assert "param1 Rhat" in msg
    
    def test_missing_parameter(self):
        """Test that missing parameters are detected."""
        posterior_r = {
            'means': {'param1': 1.0, 'param2': 2.0},
            'rhat': {'param1': 1.01, 'param2': 1.02}
        }
        posterior_py = {
            'means': {'param1': 1.0},  # Missing param2
            'rhat': {'param1': 1.01}
        }
        matches, msg = compare_posteriors(posterior_py, posterior_r)
        assert not matches
        assert "Missing parameter" in msg


class TestCompareArrays:
    """Tests for array comparison utility."""
    
    def test_identical_arrays(self):
        """Test that identical arrays match."""
        arr = np.array([1.0, 2.0, 3.0])
        matches, msg = compare_arrays(arr, arr.copy())
        assert matches
    
    def test_different_shapes(self):
        """Test that different shapes are detected."""
        arr1 = np.array([1, 2, 3])
        arr2 = np.array([[1, 2, 3]])
        matches, msg = compare_arrays(arr1, arr2)
        assert not matches
        assert "shape mismatch" in msg
    
    def test_within_tolerance(self):
        """Test array comparison with tolerance."""
        arr1 = np.array([1.0, 2.0, 3.0])
        arr2 = np.array([1.001, 2.001, 3.001])
        matches, msg = compare_arrays(arr1, arr2, rtol=0.01)
        assert matches
    
    def test_nan_handling(self):
        """Test that NaN values are handled correctly."""
        arr1 = np.array([1.0, np.nan, 3.0])
        arr2 = np.array([1.0, np.nan, 3.0])
        matches, msg = compare_arrays(arr1, arr2)
        assert matches  # NaNs in same positions should match
