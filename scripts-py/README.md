# Python Scripts for bIRTistic

This directory contains Python analysis scripts for the bIRTistic project.

## Main Analysis Script

### Colombia_analysis_ordered_logit_ADVI_vs_HMC.py

**The primary comprehensive analysis script** that reproduces the complete R markdown analysis for Colombia data.

This script performs:
1. Data loading and preprocessing
2. HMC model fitting (using zarr format for fast I/O)
3. ADVI model fitting (using zarr format for fast I/O)
4. Endpoint computation and comparison
5. Parameter comparison (ordered_prob, model params, ypred)
6. All comparison visualizations including ordered probability multi-panel plots
7. Automatically calls endpoint plotting utility

**Outputs:** All plots saved as PDF format for publication quality.

**Usage:**
```bash
cd /Users/or105/git/bIRTistic
pixi run python scripts-py/Colombia_analysis_ordered_logit_ADVI_vs_HMC.py
```

**Outputs:** All results saved to `/Users/or105/sandbox/bIRTistic/py-colombia-ordered_logit-260430-vanilla-advi-vs-hmc/`

**Performance:** Uses zarr format (8.9x faster than NetCDF) for optimal I/O performance.

---

## Utility Plotting Scripts

These scripts are called automatically by the main script but can also be run independently to regenerate specific plots:

### Colombia_endpoints_plots.py

Generates detailed endpoint comparison bar plots with faceting by category groups (matches R markdown format).

**Outputs:**
- `endpoints_comparison_barplot_diff.pdf` - Difference (Baseline - Endline) comparison
- `endpoints_comparison_barplot_ratio.pdf` - Ratio comparison

**Usage:**
```bash
pixi run python scripts-py/Colombia_endpoints_plots.py
```

---

## Interactive Development

### Option 1: Run Interactively in VS Code

**The easiest way** - Start a "Python (bIRTistic pixi)" terminal in VS Code:

1. Click the + dropdown in the terminal panel
2. Select "Python (bIRTistic pixi)"
3. The Python environment will start with paths automatically configured!

The `.pythonrc.py` startup script runs automatically and makes modules available:

```python
# You can immediately import without setup!
from data_loading import read_data_colombia
from model_fitting import fit_ordered_logit_model_ncats
from analysis import get_endpoints
```

### Option 2: Run from Command Line

```bash
cd /Users/or105/git/bIRTistic
pixi run python
```

Then in Python:
```python
# The startup script should run automatically
# If not, manually run:
exec(open('scripts-py/__init__.py').read())

# Now you can import modules
from data_loading import read_data_colombia
```

---

## Technical Notes

- All scripts use **zarr format only** for data I/O (8.9x faster than NetCDF, 4.2x overall speedup)
- All plots saved as **PDF format** for publication quality
- All plots match the R markdown output exactly
- Model fitting uses cmdstanpy with zarr InferenceData storage
- Use `pixi run python` to ensure correct environment
- Python 3.10.20 environment managed by pixi

---

## File Organization

```
scripts-py/
├── Colombia_analysis_ordered_logit_ADVI_vs_HMC.py  # Main comprehensive script (includes ordered_prob plots)
├── Colombia_endpoints_plots.py                      # Utility: endpoint plots (called by main)
├── README.md                                        # This file
└── __init__.py                                      # Path setup for interactive use
```
