# bIRTistic Test Suite

Comprehensive testing framework for validating Python implementations against R baselines.

## Overview

This test suite validates that Python ports of R/Stan code produce equivalent outputs within specified tolerances:

- **Data Loading:** DataFrames must match exactly (or within numerical tolerances)
- **Model Fitting:** Posterior means within 3% relative difference, Rhat within 0.005 absolute difference
- **Visualizations:** Numerical data + layout + colors must match

## Directory Structure

```
tests/
├── __init__.py
├── conftest.py              # pytest configuration and shared fixtures
├── test_utils.py            # Tests for comparison utilities
├── generate_r_baselines.py  # Script to generate R baseline outputs
├── fixtures/                # Test data and R baseline outputs
│   ├── README.md
│   └── r_baselines/         # Saved R outputs for comparison
├── utils/                   # Testing utilities
│   ├── __init__.py
│   └── comparison.py        # DataFrame/array/posterior comparison
└── test_*.py                # Test files (added as implementations progress)
```

## Running Tests

### Run all tests
```bash
cd ~/git/bIRTistic
pixi shell
pytest
```

### Run specific test file
```bash
pytest tests/test_utils.py
```

### Run with coverage
```bash
pytest --cov=birtistic --cov-report=html
```

### Run only fast tests (skip slow/integration)
```bash
pytest -m "not slow"
```

## Comparison Utilities

Located in `tests/utils/comparison.py`:

### `compare_dataframes()`

Compares Python and R DataFrames with configurable tolerances:

```python
from tests.utils.comparison import compare_dataframes

matches, message = compare_dataframes(
    df_python,
    df_r,
    rtol=1e-5,  # Relative tolerance for numerical columns
    atol=1e-8,  # Absolute tolerance
    check_column_order=False  # Column order doesn't need to match
)
assert matches, message
```

### `compare_scalars()`

Compares individual scalar values:

```python
from tests.utils.comparison import compare_scalars

matches, message = compare_scalars(
    value_python,
    value_r,
    rtol=1e-5,
    name="log_likelihood"
)
assert matches, message
```

### `compare_posteriors()`

Compares MCMC posterior outputs with user-specified tolerances:

```python
from tests.utils.comparison import compare_posteriors

posterior_py = {
    'means': {'alpha': 0.5, 'beta': 1.2},
    'rhat': {'alpha': 1.01, 'beta': 1.02}
}

posterior_r = {
    'means': {'alpha': 0.51, 'beta': 1.21},
    'rhat': {'alpha': 1.012, 'beta': 1.021}
}

matches, message = compare_posteriors(
    posterior_py,
    posterior_r,
    mean_rtol=0.03,   # 3% relative tolerance for means
    rhat_atol=0.005   # 0.005 absolute tolerance for Rhat
)
assert matches, message
```

### `compare_arrays()`

Compares NumPy arrays:

```python
from tests.utils.comparison import compare_arrays

matches, message = compare_arrays(
    arr_python,
    arr_r,
    rtol=1e-5,
    name="predictions"
)
assert matches, message
```

## Generating R Baselines

Before testing a ported function, generate R baseline outputs:

```bash
# List available baseline generators
python tests/generate_r_baselines.py --list

# Generate baselines for specific phase
python tests/generate_r_baselines.py --phase 2  # Data loading
python tests/generate_r_baselines.py --phase 3  # Model fitting
```

Baselines are saved to `tests/fixtures/r_baselines/`.

## Writing New Tests

### Step 1: Generate R Baseline

Add baseline generation code to `generate_r_baselines.py` for the function you're testing.

### Step 2: Create Test File

Create `test_<module>.py`:

```python
"""Tests for birtistic.data_loading module."""
import pandas as pd
import pytest
from birtistic.data_loading import read_data_colombia
from tests.utils.comparison import compare_dataframes


def test_read_data_colombia(r_baseline_dir):
    """Test Colombia data loading against R baseline."""
    # Load Python implementation
    data_py = read_data_colombia("path/to/data.csv")
    
    # Load R baseline
    dp_r = pd.read_csv(r_baseline_dir / "colombia_dp.csv")
    dit_r = pd.read_csv(r_baseline_dir / "colombia_dit.csv")
    
    # Compare
    matches, msg = compare_dataframes(data_py['dp'], dp_r)
    assert matches, msg
    
    matches, msg = compare_dataframes(data_py['dit'], dit_r)
    assert matches, msg
```

### Step 3: Run Test

```bash
pytest tests/test_data_loading.py -v
```

## Fixtures

### Available Fixtures (from `conftest.py`)

- **`test_data_dir`**: Path to `tests/fixtures/`
- **`r_baseline_dir`**: Path to `tests/fixtures/r_baselines/`
- **`temp_output_dir`**: Temporary directory for test outputs
- **`sample_dataframe`**: Sample DataFrame for utility testing

### Adding New Fixtures

Add to `conftest.py`:

```python
@pytest.fixture
def sample_posterior():
    """Sample posterior for testing."""
    return {
        'means': {'alpha': 0.5, 'beta': 1.2},
        'rhat': {'alpha': 1.01, 'beta': 1.02}
    }
```

## Validation Tolerances

Per user specification:

| Metric                          | Tolerance                     | Type                   |
| ------------------------------- | ----------------------------- | ---------------------- |
| DataFrame numerical columns     | Default: rtol=1e-5, atol=1e-8 | Relative + Absolute    |
| DataFrame non-numerical columns | Exact match                   | Exact                  |
| Posterior means                 | 3% (0.03)                     | Relative               |
| Posterior Rhat                  | 0.005                         | Absolute               |
| Same random seed                | Required for MCMC             | Exact seed match       |
| Plot numerical data             | Default: rtol=1e-5            | Relative               |
| Plot layout/colors              | Exact match                   | Qualitative inspection |

## Markers

Tests can be marked for selective execution:

```python
@pytest.mark.slow
def test_full_model_fit():
    """Long-running integration test."""
    pass

@pytest.mark.integration
def test_end_to_end_workflow():
    """Integration test requiring full pipeline."""
    pass

@pytest.mark.baseline
def test_requires_r_baseline():
    """Test requiring R baseline generation."""
    pass
```

Run specific markers:
```bash
pytest -m "not slow"  # Skip slow tests
pytest -m integration  # Run only integration tests
```

## Coverage Reports

After running tests with `--cov`:

```bash
# Open HTML coverage report
open htmlcov/index.html  # macOS
xdg-open htmlcov/index.html  # Linux
```

Target: >80% coverage for production code.

## Best Practices

1. **Use same random seeds** for MCMC/SVI comparisons
2. **Generate R baselines once** at the start of each phase
3. **Document tolerance choices** in test docstrings when deviating from defaults
4. **Test intermediate outputs**, not just final results
5. **Use descriptive assertion messages** with `assert matches, message`
6. **Mark slow tests** appropriately with `@pytest.mark.slow`

## Troubleshooting

### Test failures due to R baseline missing

Generate baselines first:
```bash
python tests/generate_r_baselines.py --phase <N>
```

### Numerical tolerance too tight

Adjust tolerances in comparison function calls:
```python
compare_dataframes(df_py, df_r, rtol=1e-3, atol=1e-6)
```

### Import errors for `birtistic` module

Ensure package structure exists:
```bash
mkdir -p python/birtistic
touch python/birtistic/__init__.py
```

## Next Steps

As implementation progresses:

1. **Phase 2 (Data Loading):** Add `test_data_loading.py`
2. **Phase 3 (Model Fitting):** Add `test_model_fitting.py`
3. **Phase 4 (Visualization):** Add `test_plotting.py`
4. **Phase 5 (NumPyro):** Add `test_numpyro_models.py`

Each test file should follow the pattern: generate R baseline → implement Python → test → validate.
