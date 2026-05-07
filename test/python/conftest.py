"""Shared pytest configuration and fixtures."""
import pytest
import pandas as pd
from pathlib import Path


@pytest.fixture
def test_data_dir():
    """Path to test fixtures directory."""
    return Path(__file__).parent / "fixtures"


@pytest.fixture
def r_baseline_dir(test_data_dir):
    """Path to R baseline outputs."""
    # Will store R outputs for comparison
    baseline_dir = test_data_dir / "r_baselines"
    baseline_dir.mkdir(exist_ok=True, parents=True)
    return baseline_dir


@pytest.fixture
def temp_output_dir(tmp_path):
    """Temporary directory for test outputs."""
    return tmp_path


@pytest.fixture
def sample_dataframe():
    """Sample DataFrame for testing."""
    return pd.DataFrame({
        'id': [1, 2, 3, 4, 5],
        'value': [10.5, 20.3, 15.7, 30.2, 25.1],
        'category': ['A', 'B', 'A', 'C', 'B']
    })
