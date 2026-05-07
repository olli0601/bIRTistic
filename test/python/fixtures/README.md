# Test Fixtures

This directory contains test data and R baseline outputs for validation.

## Structure

```
fixtures/
├── README.md           # This file
└── r_baselines/        # R baseline outputs for comparison
    ├── phase2/         # Data loading baselines (to be created)
    ├── phase3/         # Model fitting baselines (to be created)
    ├── phase4/         # Visualization baselines (to be created)
    └── phase5/         # NumPyro baselines (to be created)
```

## Purpose

R baseline outputs serve as ground truth for validating Python implementations. Each baseline is generated once from the original R code and stored here for repeatable testing.

## Baseline Generation

Baselines are generated using:

```bash
python tests/generate_r_baselines.py --phase <N>
```

See `../generate_r_baselines.py` for implementation details.

## Baseline Types

### Phase 2: Data Loading Baselines

Generated from:
- `read_data_colombia.R` → `phase2/colombia_*.csv`
- `read_data_ukraine.R` → `phase2/ukraine_*.csv`

Contents:
- `dp.csv`: Participant data (long format)
- `dit.csv`: Item metadata
- `dmeta.csv`: Participant metadata

### Phase 3: Model Fitting Baselines

Generated from:
- `fit_partial_credit_model_ncats.R` → `phase3/model_fit_*.json`

Contents:
- Posterior means (JSON)
- Rhat values (JSON)
- Model metadata

### Phase 4: Visualization Baselines

Generated from:
- Notebook plots → `phase4/plot_data_*.csv` + `phase4/plot_meta_*.json`

Contents:
- Plot data layers (CSV)
- Plot aesthetics (colors, labels, etc. as JSON)

### Phase 5: NumPyro Baselines

Uses outputs from Phase 3 (Stan) as reference for NumPyro implementations.

## File Naming Convention

```
<phase>/<function_name>_<output_type>_<variant>.ext

Examples:
- phase2/read_data_colombia_dp.csv
- phase3/fit_partial_credit_model_posterior_means.json
- phase4/plot_item_characteristics_data.csv
```

## Version Control

⚠️ **Large baseline files (>10MB) should be excluded from git** and regenerated locally.

Add to `.gitignore` if baselines become too large:
```
python/tests/fixtures/r_baselines/*.csv
python/tests/fixtures/r_baselines/*.json
```

Current approach: Small baselines (<1MB) committed to git for reproducibility.

## Updating Baselines

Baselines should be regenerated when:
1. R code is modified (bug fixes, improvements)
2. Input data changes
3. Random seeds change
4. R package versions update significantly

**Important:** Document baseline regeneration in git commit message.

## Accessing in Tests

Use the `r_baseline_dir` fixture:

```python
def test_function(r_baseline_dir):
    baseline_path = r_baseline_dir / "phase2" / "colombia_dp.csv"
    df_r = pd.read_csv(baseline_path)
    # Compare to Python output
```

## Baseline Validation

Each baseline should be validated immediately after generation:

1. Check file exists and is non-empty
2. Verify expected columns/structure
3. Check for NaN/missing values (if unexpected)
4. Verify data types
5. Document any anomalies

## Maintenance

- Review baselines quarterly or when major R package updates occur
- Clean up obsolete baselines from previous phases
- Compress large baselines if necessary (e.g., parquet instead of CSV)
