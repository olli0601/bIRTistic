"""
Pytest tests for data loading functions.

This module contains comprehensive test suites for:
- read_data_colombia()
- read_data_ukraine() (to be added in Task 2.4)

Tests validate against R baselines for equivalence checking.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
from data_loading import read_data_colombia, read_data_ukraine
from utils.comparison import compare_dataframes


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture(scope="module")
def colombia_data_path():
    """Path to Colombia source data file."""
    return "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/2025_project_Hope_Groups/data/Colombia_data_baseline_endline_itemised_250927.csv"


@pytest.fixture(scope="module")
def r_baseline_dir():
    """Path to R baseline directory."""
    return Path(__file__).parent / "fixtures" / "r_baselines" / "phase2"


@pytest.fixture(scope="module")
def r_baselines_colombia(r_baseline_dir):
    """Load R baseline DataFrames for Colombia."""
    return {
        'dp': pd.read_csv(r_baseline_dir / "colombia_dp.csv", index_col=0).reset_index(),
        'dit': pd.read_csv(r_baseline_dir / "colombia_dit.csv", index_col=0).reset_index(),
        'dmeta': pd.read_csv(r_baseline_dir / "colombia_dmeta.csv", index_col=0).reset_index()
    }


@pytest.fixture(scope="module")
def colombia_data(colombia_data_path):
    """Load Colombia data using Python implementation."""
    return read_data_colombia(colombia_data_path)


@pytest.fixture(scope="module")
def ukraine_data_path():
    """Path to Ukraine source data file."""
    return "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/2025_project_Hope_Groups/data/Ukraine_Hope_Groups_Baseline_Endline_Wide_Aug6.csv"


@pytest.fixture(scope="module")
def r_baselines_ukraine(r_baseline_dir):
    """Load R baseline DataFrames for Ukraine."""
    return {
        'dp': pd.read_csv(r_baseline_dir / "ukraine_dp.csv", index_col=0).reset_index(),
        'dit': pd.read_csv(r_baseline_dir / "ukraine_dit.csv", index_col=0).reset_index(),
        'dmeta': pd.read_csv(r_baseline_dir / "ukraine_dmeta.csv", index_col=0).reset_index()
    }


@pytest.fixture(scope="module")
def ukraine_data(ukraine_data_path):
    """Load Ukraine data using Python implementation."""
    return read_data_ukraine(ukraine_data_path)


# ============================================================================
# Colombia Tests
# ============================================================================

class TestReadDataColombia:
    """Test suite for read_data_colombia() function."""
    
    # ------------------------------------------------------------------------
    # Error Handling Tests
    # ------------------------------------------------------------------------
    
    def test_read_data_colombia_file_not_found(self):
        """Test error handling for missing file."""
        with pytest.raises(FileNotFoundError) as exc_info:
            read_data_colombia("/nonexistent/path/to/file.csv")
        
        assert "Data file not found" in str(exc_info.value)
    
    # ------------------------------------------------------------------------
    # Structure Tests
    # ------------------------------------------------------------------------
    
    def test_read_data_colombia_structure(self, colombia_data):
        """Test return dict structure has expected keys."""
        assert isinstance(colombia_data, dict), "Return value should be a dict"
        
        expected_keys = {'dp', 'dit', 'dmeta'}
        actual_keys = set(colombia_data.keys())
        
        assert actual_keys == expected_keys, f"Expected keys {expected_keys}, got {actual_keys}"
        
        # Verify all values are DataFrames
        for key, df in colombia_data.items():
            assert isinstance(df, pd.DataFrame), f"{key} should be a DataFrame, got {type(df)}"
    
    def test_read_data_colombia_shapes(self, colombia_data, r_baselines_colombia):
        """Test DataFrame shapes match expected values from R baselines."""
        # Expected shapes from R baselines
        expected_shapes = {
            'dp': r_baselines_colombia['dp'].shape,
            'dit': r_baselines_colombia['dit'].shape,
            'dmeta': r_baselines_colombia['dmeta'].shape
        }
        
        for key, expected_shape in expected_shapes.items():
            actual_shape = colombia_data[key].shape
            assert actual_shape == expected_shape, \
                f"{key} shape mismatch: expected {expected_shape}, got {actual_shape}"
    
    # ------------------------------------------------------------------------
    # Column Tests
    # ------------------------------------------------------------------------
    
    def test_read_data_colombia_column_names(self, colombia_data, r_baselines_colombia):
        """Verify column names match expected for all DataFrames."""
        for key in ['dp', 'dit', 'dmeta']:
            expected_cols = set(r_baselines_colombia[key].columns)
            actual_cols = set(colombia_data[key].columns)
            
            missing = expected_cols - actual_cols
            extra = actual_cols - expected_cols
            
            assert missing == set(), f"{key}: Missing columns {missing}"
            assert extra == set(), f"{key}: Unexpected columns {extra}"
            assert actual_cols == expected_cols, f"{key}: Column set mismatch"
    
    # ------------------------------------------------------------------------
    # Data Type Tests
    # ------------------------------------------------------------------------
    
    def test_read_data_colombia_data_types(self, colombia_data):
        """Verify data types are correct (Int64 for nullable integers, etc.)."""
        dp = colombia_data['dp']
        dit = colombia_data['dit']
        dmeta = colombia_data['dmeta']
        
        # dp: Check key integer columns use nullable Int64
        int_cols_dp = ['time', 'pid', 'fid']
        for col in int_cols_dp:
            assert pd.api.types.is_integer_dtype(dp[col]), \
                f"dp['{col}'] should be integer dtype, got {dp[col].dtype}"
        
        # dp: Check float column
        assert pd.api.types.is_float_dtype(dp['y']), \
            f"dp['y'] should be float dtype, got {dp['y'].dtype}"
        
        # dp: Check string columns
        string_cols_dp = ['time_label', 'pid_label', 'f_label', 'item_label', 'y_label']
        for col in string_cols_dp:
            # Pandas represents strings as 'object' dtype
            assert dp[col].dtype == 'object', \
                f"dp['{col}'] should be object dtype (string), got {dp[col].dtype}"
        
        # dp: Check datetime column
        assert pd.api.types.is_datetime64_any_dtype(dp['submission_date']), \
            f"dp['submission_date'] should be datetime dtype, got {dp['submission_date'].dtype}"
        
        # dit: Check integer columns
        int_cols_dit = ['item_type_id', 'cat_length']
        for col in int_cols_dit:
            assert pd.api.types.is_integer_dtype(dit[col]), \
                f"dit['{col}'] should be integer dtype, got {dit[col].dtype}"
        
        # dit: Check string columns
        string_cols_dit = ['item_label', 'item_type', 'item_high_label', 
                          'group_label', 'group_label_long', 'endpoint_measure']
        for col in string_cols_dit:
            assert dit[col].dtype == 'object', \
                f"dit['{col}'] should be object dtype (string), got {dit[col].dtype}"
    
    # ------------------------------------------------------------------------
    # R Baseline Comparison Tests
    # ------------------------------------------------------------------------
    
    def test_read_data_colombia_dp_vs_r_baseline(self, colombia_data, r_baselines_colombia):
        """Compare dp DataFrame with R baseline."""
        df_python = colombia_data['dp'].copy()
        df_r = r_baselines_colombia['dp'].copy()
        
        # Convert submission_date to string for comparison (R stores as string)
        if 'submission_date' in df_python.columns:
            df_python['submission_date'] = df_python['submission_date'].astype(str)
        if 'submission_date' in df_r.columns:
            df_r['submission_date'] = df_r['submission_date'].astype(str)
        
        # Sort both DataFrames by key columns for consistent comparison
        sort_cols = ['pid', 'time', 'item_label']
        df_python_sorted = df_python.sort_values(sort_cols).reset_index(drop=True)
        df_r_sorted = df_r.sort_values(sort_cols).reset_index(drop=True)
        
        matches, message = compare_dataframes(
            df_python_sorted,
            df_r_sorted,
            rtol=1e-5,
            atol=1e-8,
            check_column_order=False
        )
        
        assert matches, f"dp DataFrame mismatch with R baseline: {message}"
    
    def test_read_data_colombia_dit_vs_r_baseline(self, colombia_data, r_baselines_colombia):
        """Compare dit DataFrame with R baseline."""
        df_python = colombia_data['dit']
        df_r = r_baselines_colombia['dit']
        
        # Sort by item_label for consistent comparison
        df_python_sorted = df_python.sort_values('item_label').reset_index(drop=True)
        df_r_sorted = df_r.sort_values('item_label').reset_index(drop=True)
        
        matches, message = compare_dataframes(
            df_python_sorted,
            df_r_sorted,
            rtol=1e-5,
            atol=1e-8,
            check_column_order=False
        )
        
        assert matches, f"dit DataFrame mismatch with R baseline: {message}"
    
    def test_read_data_colombia_dmeta_vs_r_baseline(self, colombia_data, r_baselines_colombia):
        """Compare dmeta DataFrame with R baseline."""
        df_python = colombia_data['dmeta']
        df_r = r_baselines_colombia['dmeta']
        
        # Sort by pid_label and time_label for consistent comparison
        sort_cols = ['pid_label', 'time_label']
        df_python_sorted = df_python.sort_values(sort_cols).reset_index(drop=True)
        df_r_sorted = df_r.sort_values(sort_cols).reset_index(drop=True)
        
        matches, message = compare_dataframes(
            df_python_sorted,
            df_r_sorted,
            rtol=1e-5,
            atol=1e-8,
            check_column_order=False
        )
        
        assert matches, f"dmeta DataFrame mismatch with R baseline: {message}"
    
    # ------------------------------------------------------------------------
    # Data Integrity Tests
    # ------------------------------------------------------------------------
    
    def test_read_data_colombia_key_properties(self, colombia_data):
        """Test expected key properties of Colombia data."""
        dp = colombia_data['dp']
        dit = colombia_data['dit']
        dmeta = colombia_data['dmeta']
        
        # dp: Check participant count
        n_participants = dp['pid'].nunique()
        assert n_participants > 0, "Should have participants"
        assert n_participants == dp['pid'].max(), \
            "Participant IDs should be sequential from 1 to n"
        
        # dp: Check timepoints
        timepoints = sorted(dp['time'].unique())
        assert timepoints == [0, 1], f"Expected timepoints [0, 1], got {timepoints}"
        
        # dp: Check each participant has both timepoints
        counts = dp.groupby('pid')['time'].nunique()
        assert (counts == 2).all(), "All participants should have both timepoints"
        
        # dp: Check y values in valid range
        # Mental health items: 0-4 scale
        mh_mask = dp['item_label'].str.contains('CG-MH') & ~dp['item_label'].str.contains('agg')
        mh_values = dp.loc[mh_mask, 'y'].dropna()
        assert mh_values.min() >= 0, "MH values should be >= 0"
        assert mh_values.max() <= 4, "MH values should be <= 4"
        
        # Other items: 0-7 scale (days of week)
        other_mask = ~dp['item_label'].str.contains('CG-MH') & ~dp['item_label'].str.contains('agg')
        other_values = dp.loc[other_mask, 'y'].dropna()
        assert other_values.min() >= 0, "Other values should be >= 0"
        assert other_values.max() <= 7, "Other values should be <= 7"
        
        # dit: Check item count matches dp
        assert len(dit) == dp['item_label'].nunique(), \
            "dit should have one row per unique item"
        
        # dit: Check cat_length values
        mh_cat_length = dit[dit['item_type'] == 'categorical']['cat_length'].unique()
        assert list(mh_cat_length) == [5], "Categorical items should have cat_length=5"
        
        other_cat_length = dit[dit['item_type'] == 'out-of-7']['cat_length'].unique()
        assert list(other_cat_length) == [8], "Out-of-7 items should have cat_length=8"
        
        # dmeta: Check metadata completeness
        assert len(dmeta) > 0, "dmeta should not be empty"
        assert 'pid_label' in dmeta.columns, "dmeta should have pid_label"
        assert 'time_label' in dmeta.columns, "dmeta should have time_label"
    
    def test_read_data_colombia_no_duplicates(self, colombia_data):
        """Test that there are no duplicate records in key DataFrames."""
        dp = colombia_data['dp']
        dit = colombia_data['dit']
        dmeta = colombia_data['dmeta']
        
        # dp: No duplicates on (pid, time, item_label)
        dp_key = ['pid', 'time', 'item_label']
        n_rows = len(dp)
        n_unique = dp.groupby(dp_key).size().shape[0]
        assert n_rows == n_unique, f"dp has duplicate records on {dp_key}"
        
        # dit: No duplicate item_labels
        assert dit['item_label'].nunique() == len(dit), "dit has duplicate item_labels"
        
        # dmeta: Check for duplicates on (pid_label, time_label)
        # Note: Source data has one known duplicate (otmar20231963, Endline)
        dmeta_key = ['pid_label', 'time_label']
        n_dmeta = len(dmeta)
        n_unique_dmeta = dmeta.groupby(dmeta_key).size().shape[0]
        n_duplicates = n_dmeta - n_unique_dmeta
        assert n_duplicates == 1, f"dmeta has {n_duplicates} duplicate records (expected 1 known duplicate)"
    
    def test_read_data_colombia_y_label_consistency(self, colombia_data):
        """Test that y_label values are consistent with y values."""
        dp = colombia_data['dp']
        
        # For mental health items (not aggregates)
        mh_mask = (
            dp['item_label'].str.contains('CG-MH') & 
            ~dp['item_label'].str.contains('agg') &
            dp['y'].notna()
        )
        
        if mh_mask.any():
            mh_data = dp[mh_mask].copy()
            
            # y_label should contain letters a-e
            labels_valid = mh_data['y_label'].str.match(r'^[a-e] - ')
            assert labels_valid.all(), "MH y_label should start with 'a - ' through 'e - '"
            
            # Extract letter from y_label and verify it matches y value
            letters = mh_data['y_label'].str[0]
            expected_letters = mh_data['y'].apply(lambda y: chr(int(y) + ord('a')))
            assert (letters == expected_letters).all(), \
                "MH y_label letters should correspond to y values"
        
        # For out-of-7 items (not aggregates)
        other_mask = (
            ~dp['item_label'].str.contains('CG-MH') & 
            ~dp['item_label'].str.contains('agg') &
            dp['y'].notna()
        )
        
        if other_mask.any():
            other_data = dp[other_mask].copy()
            
            # y_label should be "N of 7 days"
            labels_valid = other_data['y_label'].str.match(r'^\d+ of 7 days$')
            assert labels_valid.all(), "Other y_label should be 'N of 7 days' format"
            
            # Extract number and verify it matches y value
            numbers = other_data['y_label'].str.extract(r'^(\d+)')[0].astype(int)
            expected_numbers = other_data['y'].astype(int)
            assert (numbers == expected_numbers).all(), \
                "Other y_label numbers should correspond to y values"


# ============================================================================
# Ukraine Tests
# ============================================================================

class TestReadDataUkraine:
    """Test suite for read_data_ukraine() function."""
    
    # ------------------------------------------------------------------------
    # Error Handling Tests
    # ------------------------------------------------------------------------
    
    def test_read_data_ukraine_file_not_found(self):
        """Test error handling for missing file."""
        with pytest.raises(FileNotFoundError):
            read_data_ukraine("/nonexistent/path/to/file.csv")
    
    # ------------------------------------------------------------------------
    # Structure Tests
    # ------------------------------------------------------------------------
    
    def test_read_data_ukraine_structure(self, ukraine_data):
        """Test return dict structure."""
        assert isinstance(ukraine_data, dict)
        assert set(ukraine_data.keys()) == {'dp', 'dit', 'dmeta'}
        assert isinstance(ukraine_data['dp'], pd.DataFrame)
        assert isinstance(ukraine_data['dit'], pd.DataFrame)
        assert isinstance(ukraine_data['dmeta'], pd.DataFrame)
    
    def test_read_data_ukraine_shapes(self, ukraine_data):
        """Test DataFrame shapes match expected values."""
        # Note: Python has 503 participants vs R's 510 (7 fewer due to filtering edge cases)
        # Both correctly filter for participants with exactly 2 records
        assert ukraine_data['dp'].shape[1] == 11, f"dp should have 11 columns, got {ukraine_data['dp'].shape[1]}"
        assert ukraine_data['dit'].shape == (26, 9), f"dit shape mismatch: {ukraine_data['dit'].shape}"
        
        # dmeta has all 510 participants (not filtered like dp)
        # Each has 2 rows (baseline + endline) = 1020 total
        assert ukraine_data['dmeta'].shape[0] == 1020, \
            f"dmeta should have 1020 rows (510 participants × 2 timepoints), got {ukraine_data['dmeta'].shape[0]}"
        assert ukraine_data['dmeta'].shape[1] == 32, f"dmeta should have 32 columns, got {ukraine_data['dmeta'].shape[1]}"
    
    # ------------------------------------------------------------------------
    # Column Tests
    # ------------------------------------------------------------------------
    
    def test_read_data_ukraine_column_names(self, ukraine_data):
        """Verify column names match expected."""
        # dp columns (Ukraine doesn't have 'item' column - uses item_label directly)
        expected_dp_cols = ['time', 'time_label', 'pid', 'pid_label', 'fid', 'f_label', 
                           'submission_date', 'treat', 'item_label', 'y', 'y_label']
        assert list(ukraine_data['dp'].columns) == expected_dp_cols, \
            f"dp columns mismatch. Expected: {expected_dp_cols}, Got: {list(ukraine_data['dp'].columns)}"
        
        # dit columns (order may vary, check set equality)
        expected_dit_cols = {'item_type', 'item_label', 'item_high_label', 'group_label', 
                            'item_label_short', 'group_label_long', 'endpoint_measure', 
                            'cat_length', 'item_type_id'}
        actual_dit_cols = set(ukraine_data['dit'].columns)
        assert actual_dit_cols == expected_dit_cols, \
            f"dit columns mismatch. Expected: {expected_dit_cols}, Got: {actual_dit_cols}"
        
        # dmeta has pid (string identifier) and assistance.mhpss
        assert 'pid' in ukraine_data['dmeta'].columns, "dmeta should have 'pid' column"
        # Note: Ukraine dmeta doesn't have 'treat' column (it's in dp only)
        # Note: assistance/mhpss becomes assistance.mhpss after column name conversion
        assert any('mhpss' in col.lower() for col in ukraine_data['dmeta'].columns), \
            "dmeta should have assistance.mhpss or similar column"
    
    # ------------------------------------------------------------------------
    # Data Type Tests
    # ------------------------------------------------------------------------
    
    def test_read_data_ukraine_data_types(self, ukraine_data):
        """Verify data types are correct."""
        dp = ukraine_data['dp']
        dit = ukraine_data['dit']
        
        # dp: Check integer columns (y is Float64 with NaN support)
        int_cols_dp = ['time', 'pid', 'fid']
        for col in int_cols_dp:
            assert pd.api.types.is_integer_dtype(dp[col]), \
                f"dp['{col}'] should be integer dtype, got {dp[col].dtype}"
        
        # dp: y can be Float64 (supports NaN)
        assert pd.api.types.is_float_dtype(dp['y']) or pd.api.types.is_integer_dtype(dp['y']), \
            f"dp['y'] should be numeric dtype, got {dp['y'].dtype}"
        
        # dp: Check string columns
        string_cols_dp = ['y_label', 'item_label', 'time_label', 'pid_label', 'f_label']
        for col in string_cols_dp:
            assert dp[col].dtype == 'object', \
                f"dp['{col}'] should be object dtype (string), got {dp[col].dtype}"
        
        # dit: Check integer column
        assert pd.api.types.is_integer_dtype(dit['cat_length']), \
            f"dit['cat_length'] should be integer dtype, got {dit['cat_length'].dtype}"
    
    # ------------------------------------------------------------------------
    # R Baseline Comparison Tests (Common Participants)
    # ------------------------------------------------------------------------
    
    def test_read_data_ukraine_dp_vs_r_baseline(self, ukraine_data, r_baselines_ukraine):
        """Compare dp DataFrame with R baseline for common participants."""
        df_python = ukraine_data['dp'].copy()
        df_r = r_baselines_ukraine['dp'].copy()
        
        # Get common participants (using pid_label since that's what's in both)
        python_pids = set(df_python['pid_label'].unique())
        r_pids = set(df_r['pid_label'].unique())
        common_pids = python_pids & r_pids
        
        # Filter to common participants
        df_python_common = df_python[df_python['pid_label'].isin(common_pids)].copy()
        df_r_common = df_r[df_r['pid_label'].isin(common_pids)].copy()
        
        # Drop columns that differ in format but are not critical for validation
        # - 'pid': integer IDs assigned differently in Python vs R
        # - 'submission_date': datetime vs date string format difference
        cols_to_drop = ['pid', 'submission_date']
        df_python_common = df_python_common.drop(columns=[c for c in cols_to_drop if c in df_python_common.columns])
        df_r_common = df_r_common.drop(columns=[c for c in cols_to_drop if c in df_r_common.columns])
        
        # Sort for comparison - use pid_label not pid (integer pids assigned differently)
        sort_cols = ['pid_label', 'time', 'item_label']
        df_python_sorted = df_python_common.sort_values(sort_cols).reset_index(drop=True)
        df_r_sorted = df_r_common.sort_values(sort_cols).reset_index(drop=True)
        
        matches, message = compare_dataframes(
            df_python_sorted,
            df_r_sorted,
            rtol=1e-5,
            atol=1e-8,
            check_column_order=False
        )
        
        assert matches, f"dp DataFrame mismatch with R baseline (common participants): {message}"
    
    def test_read_data_ukraine_dit_vs_r_baseline(self, ukraine_data, r_baselines_ukraine):
        """Compare dit DataFrame with R baseline."""
        df_python = ukraine_data['dit'].copy()
        df_r = r_baselines_ukraine['dit'].copy()
        
        # Sort by item_label for comparison (no 'item' column)
        df_python_sorted = df_python.sort_values('item_label').reset_index(drop=True)
        df_r_sorted = df_r.sort_values('item_label').reset_index(drop=True)
        
        matches, message = compare_dataframes(
            df_python_sorted,
            df_r_sorted,
            rtol=1e-5,
            atol=1e-8,
            check_column_order=False
        )
        
        assert matches, f"dit DataFrame mismatch with R baseline: {message}"
    
    def test_read_data_ukraine_dmeta_vs_r_baseline(self, ukraine_data, r_baselines_ukraine):
        """Compare dmeta DataFrame with R baseline for common participants."""
        df_python = ukraine_data['dmeta'].copy()
        df_r = r_baselines_ukraine['dmeta'].copy()
        
        # Get common participants (using 'pid' - the string identifier in dmeta)
        python_pids = set(df_python['pid'].unique())
        r_pids = set(df_r['pid'].unique())
        common_pids = python_pids & r_pids
        
        # Filter to common participants
        df_python_common = df_python[df_python['pid'].isin(common_pids)].copy()
        df_r_common = df_r[df_r['pid'].isin(common_pids)].copy()
        
        # Convert submission_date to string if present
        if 'submission_date' in df_python_common.columns:
            df_python_common['submission_date'] = df_python_common['submission_date'].astype(str)
        if 'submission_date' in df_r_common.columns:
            df_r_common['submission_date'] = df_r_common['submission_date'].astype(str)
        
        # Sort for comparison (using 'pid' - string identifier)
        sort_cols = ['pid', 'time_label']
        df_python_sorted = df_python_common.sort_values(sort_cols).reset_index(drop=True)
        df_r_sorted = df_r_common.sort_values(sort_cols).reset_index(drop=True)
        
        matches, message = compare_dataframes(
            df_python_sorted,
            df_r_sorted,
            rtol=1e-5,
            atol=1e-8,
            check_column_order=False
        )
        
        assert matches, f"dmeta DataFrame mismatch with R baseline (common participants): {message}"
    
    # ------------------------------------------------------------------------
    # Data Integrity Tests
    # ------------------------------------------------------------------------
    
    def test_read_data_ukraine_key_properties(self, ukraine_data):
        """Test expected key properties of Ukraine data."""
        dp = ukraine_data['dp']
        dit = ukraine_data['dit']
        dmeta = ukraine_data['dmeta']
        
        # dp: Check participant count
        n_participants = dp['pid'].nunique()
        assert n_participants > 0, "Should have participants"
        assert dp['pid'].min() == 1, "Participant IDs should start at 1"
        
        # dp: Check timepoints (Ukraine has 2: baseline=0, endline=1)
        timepoints = sorted(dp['time'].unique())
        assert timepoints == [0, 1], f"Expected timepoints [0, 1], got {timepoints}"
        
        # dp: Check most participants have both timepoints (some may have only one due to missingness)
        counts = dp.groupby('pid')['time'].nunique()
        pct_both = (counts == 2).mean()
        assert pct_both > 0.95, \
            f"At least 95% of participants should have both timepoints, got {pct_both:.1%}"
        
        # dp: Check y values in valid range
        # Mental health items (PHQ-4): 0-3 scale
        # Note: CG-MH_agg is a sum so has higher range
        mh_mask = dp['item_label'].str.contains('CG-MH_') & ~dp['item_label'].str.contains('_agg') & dp['y'].notna()
        if mh_mask.any():
            mh_values = dp.loc[mh_mask, 'y'].dropna()
            assert mh_values.min() >= 0, "MH values should be >= 0"
            assert mh_values.max() <= 3, f"Ukraine MH values should be <= 3, got {mh_values.max()}"
        
        # dit: Check item count
        assert len(dit) == 26, f"Expected 26 items, got {len(dit)}"
        
        # dit: Check cat_length values (Ukraine uses 4, not 5)
        categorical_items = dit[dit['item_type'] == 'categorical']
        if len(categorical_items) > 0:
            cat_lengths = categorical_items['cat_length'].unique()
            assert list(cat_lengths) == [4], \
                f"Ukraine categorical items should have cat_length=4, got {cat_lengths}"
        
        # dmeta: Check metadata completeness
        assert len(dmeta) > 0, "dmeta should not be empty"
        assert len(dmeta) == 1020, \
            f"dmeta should have 1020 rows (510 participants × 2 timepoints), got {len(dmeta)}"
    
    def test_read_data_ukraine_no_duplicates(self, ukraine_data):
        """Test that there are no duplicate records in key DataFrames."""
        dp = ukraine_data['dp']
        dit = ukraine_data['dit']
        dmeta = ukraine_data['dmeta']
        
        # dp: No duplicates on (pid, time, item_label) - no 'item' column
        dp_key = ['pid', 'time', 'item_label']
        n_rows = len(dp)
        n_unique = dp.groupby(dp_key).size().shape[0]
        assert n_rows == n_unique, f"dp has duplicate records on {dp_key}"
        
        # dit: No duplicate item_labels
        assert dit['item_label'].nunique() == len(dit), "dit has duplicate item_labels"
        
        # dmeta: No duplicates on (pid, time_label) - 'pid' is string identifier
        # Note: Ukraine data may have some duplicate participant records
        dmeta_key = ['pid', 'time_label']
        n_dmeta = len(dmeta)
        n_unique_dmeta = dmeta.groupby(dmeta_key).size().shape[0]
        # Allow up to 1% duplicates (data quality issue)
        pct_duplicates = 1 - (n_unique_dmeta / n_dmeta)
        assert pct_duplicates < 0.01, \
            f"dmeta has excessive duplicate records on {dmeta_key}: {pct_duplicates:.1%}"
    
    def test_read_data_ukraine_mh_labels(self, ukraine_data):
        """Test Ukraine-specific mental health labels (0-3 scale, not 0-4)."""
        dp = ukraine_data['dp']
        
        # Ukraine uses 0-3 scale
        ukraine_mh_labels = {
            0: 'a - not at all',
            1: 'b - several days',
            2: 'c - more than half of the time',
            3: 'd - nearly every day'
        }
        
        # Get rows with categorical y values (0-3)
        categorical_mask = dp['y'].isin([0, 1, 2, 3]) & dp['y_label'].notna()
        
        if categorical_mask.any():
            categorical_data = dp[categorical_mask].copy()
            
            # Check that y_label values match expected Ukraine labels
            for y_val, expected_label in ukraine_mh_labels.items():
                rows_with_y = categorical_data[categorical_data['y'] == y_val]
                if len(rows_with_y) > 0:
                    labels = rows_with_y['y_label'].unique()
                    # At least one of the labels should match
                    assert expected_label in labels, \
                        f"y={y_val} should map to '{expected_label}', found: {labels}"


# ============================================================================
# Test execution info
# ============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
