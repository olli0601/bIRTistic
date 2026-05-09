# bIRTistic

Making Bayesian Item Response Theory models faster with amortised inference

## Overview

This repository contains tools for running Bayesian Item Response Theory (IRT) analyses using Stan, with support for parallel execution on HPC systems.

## Installation

### Setup Instructions

Please install our prerequisites:

- [Pixi](https://pixi.sh/) (recommended - fast, cross-platform package manager)
  - Install: `curl -fsSL https://pixi.sh/install.sh | bash`
  - Or with Homebrew: `brew install pixi`
- Alternatively: [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/)
- Git

We strongly recommend to use Pixi for installation:

### Environment Options

This repository supports two pixi environments:

1. **`default`** (Python + R): Full stack with Python packages (JAX, NumPyro, etc.) and R packages
2. **`r-only`**: Minimal R-only environment with CmdStan (no Python dependencies)

Choose the environment that fits your needs.

Option 1: Install default environment (Python + R)

```bash
# Clone the repository
git clone https://github.com/olli0601/bIRTistic.git
cd bIRTistic

# Install all dependencies (Python + R)
pixi install

# Setup cmdstanr etc from CRAN/R-universe
pixi run setup

# Verify R installation (includes 8 schools model test)
pixi run verify-r

# Verify Python installation (includes CmdStanPy and 8 schools test)
pixi run verify-python

# Activate the environment
pixi shell
```

Option 2: Install R environment only

```bash
# Clone the repository
git clone https://github.com/olli0601/bIRTistic.git
cd bIRTistic

# Install R-only dependencies (no Python)
pixi install -e r-only

# Setup cmdstanr etc from CRAN/R-universe
pixi run -e r-only setup

# Verify installation (includes 8 schools model test)
pixi run -e r-only verify-r

# Activate the R-only environment
pixi shell -e r-only
```

<details>
<summary>An alternative though slower installation routine with conda/mamba can be found here</summary>


#### On Imperial HPC:
```bash
module load miniforge/3 tools/prod zlib libxml2/2.11.5-GCCcore-13.2.0
eval "$(~/miniforge3/bin/conda shell.bash hook)"
eval "$(mamba shell hook --shell bash)"

# Install cmdstan only once and expose to other environment files
mamba create -n cmdstan-build -c conda-forge cmdstan=2.37.0 -y

# Minimum conda environment to avoid installation from source
mamba env create -f bIRTistic.yml -y
mamba activate birtistic
```

#### Install python packages:

```bash
# Install core packages
pip install numpy scipy pandas matplotlib
pip install jax jaxlib flax optax
pip install numpyro cmdstanpy
pip install plotnine
pip install pytest pytest-cov jupyterlab ipykernel

# Verify installation
python python/test_environment.py
```

#### Install R packages through R:

After activating the environment, install R packages through R (possibly from source):

```bash
# Installing data.table (requires zlib, expose via module load)
Rscript -e "packages <- c('data.table'); install.packages(packages[!packages %in% installed.packages()[,'Package']], repos = 'https://cloud.r-project.org')"

# Installing everything else (xml2 requires libxml2/2.11.5-GCCcore-13.2.0, expose via module load)
Rscript -e "packages <- c('argparse', 'bayesplot', 'bookdown', 'formatR', 'GGally', 'ggpattern', 'ggplot2', 'ggsci', 'gridExtra', 'ggpubr', 'here', 'hexbin', 'jsonlite', 'kableExtra', 'knitr', 'loo', 'MatchIt', 'matrixStats', 'patchwork', 'polycor', 'posterior', 'purrr','rmarkdown', 'scales'); install.packages(packages[!packages %in% installed.packages()[,'Package']], repos = 'https://cloud.r-project.org', ask = false)"
```

#### Verify installation:

The conda/mamba installation can be verified using the same tests that `pixi run verify-r` uses:

```bash
# Check CmdStan version
Rscript -e "library(cmdstanr); cmdstanr::cmdstan_version()"

# Check Bernoulli model compiles and runs
Rscript -e "library(cmdstanr); file <- file.path(cmdstan_path(), 'examples', 'bernoulli', 'bernoulli.stan'); mod <- cmdstan_model(file, force_recompile = TRUE); data_list <- list(N = 10, y = c(0,1,0,0,0,0,0,0,0,1)); fit <- mod\$sample(data = data_list, seed = 123, chains = 2, parallel_chains = 1, refresh = 500);"
```

#### Running Examples with Conda/Mamba on HPC:

```bash
# Load modules and activate environment
module load miniforge/3 tools/prod zlib libxml2/2.11.5-GCCcore-13.2.0
eval "$(~/miniforge3/bin/conda shell.bash hook)"
eval "$(mamba shell hook --shell bash)"
mamba activate birtistic
```

</details>

## Quick Start

### Data Loading

```python
from birtistic.data_loading import read_data_colombia, read_data_ukraine

# Load Colombia data
data = read_data_colombia("path/to/Colombia_data_baseline_endline_itemised.csv")
dp = data['dp']    # Participant data (long format)
dit = data['dit']  # Item metadata
dmeta = data['dmeta']  # Participant metadata
```

### Model Fitting

```python
from birtistic.model_fitting import fit_partial_credit_model_ncats

# Fit model with HMC
results = fit_partial_credit_model_ncats(
    data=dp,
    item_metadata=dit,
    method='hmc',  # or 'advi' for variational inference
    seed=12345
)
```

### Plotting

```python
from birtistic.plotting import plot_item_characteristics

# Create plots
fig = plot_item_characteristics(results)
fig.save("item_characteristics.png")
```


### Running R Analysis Examples

```bash
# Activate the environment
pixi shell

# Run our unit tests
Rscript test/run_all_tests.R

# then e.g. render a vignette with
Rscript -e "rmarkdown::render('vignettes/Colombia_analysis_ordered_logit_ADVI_vs_HMC.Rmd')"
```

**Example analyses** (see `vignettes/` directory):
- `Colombia_analysis_for_HGpaper_ordered_logit_v251224.Rmd` - Hope Groups paper analysis using ordered logit model 
- `Colombia_analysis_for_HGpaper_pcm_v251224.Rmd` - Hope Groups paper analysis using partial credit model 
- `Colombia_analysis_ordered_logit_ADVI_vs_HMC.Rmd` - Compare ADVI vs HMC inference methods
- `Colombia_compare_models.Rmd` - Compare IRT models on Colombia data (local)
- `hpc_Colombia_compare_models.Rmd` - Compare IRT models on Colombia data (HPC parallel)
- `Colombia_interim_analyses_v251007.Rmd` - Interim analyses on Colombia data (local)
- `Colombia_interim_analyses_parallel.R` - Interim analyses on Colombia data (HPC parallel)
- `Colombia_validate_credit_model.Rmd` - Train/test validation withholding participant responses (local)
- `hpc_Colombia_train_test_item_response_models.Rmd` - Train/test validation (HPC parallel)


## Testing

Run the test suite:

```bash
# Run all tests
pixi shell
pytest

# Run with coverage
pytest --cov=birtistic --cov-report=html

# Run specific test file
pytest tests/test_data_loading.py
```


## Project Structure

<details>
<summary>Click to expand</summary>

```
bIRTistic/
├── .vscode/                                 # VS Code workspace configuration
│   ├── birtistic_pixi_r.sh                  # R terminal launcher for pixi environment
│   ├── icacuk_hpc_r_birtistic.sh            # R terminal launcher for HPC mamba environment
│   ├── settings.json                        # Editor settings and R configuration
│   └── tasks.json                           # Build tasks (Rtex, LaTeX, Markdown)
├── R/                                       # R scripts and functions
│   ├── fit_credit_model_ncats.R             # Fit credit model (ncats)
│   ├── fit_ordered_logit_model_ncats.R      # Fit ordered logit model (ncats)
│   ├── fit_ordered_logit_model_ncats_advi.R # Fit ordered logit model with ADVI
│   ├── fit_partial_credit_model_ncats.R     # Fit partial credit model (ncats)
│   ├── get_endpoints.R                      # Compute endpoint statistics from model fits
│   ├── fit_item_response_models.Rscript     # Generic IRT model fitting script (for HPC)
│   ├── hpc_submit_generic_job.R             # HPC job submission helper
│   ├── read_data_colombia.R                 # Read and preprocess Colombia data
│   ├── read_data_ukraine.R                  # Read and preprocess Ukraine data
│   ├── train_test_split_data.R              # Train/test data splitting
│   ├── train_test_item_response_models.Rscript       # Train/test IRT models
│   └── train_test_item_response_models_template.json # Configuration template
├── src/
│   ├── stan/                                # Current Stan models (ncats versions)
│   │   ├── credit_model_ncats_v260413.stan
│   │   ├── credit_model_functions.stan
│   │   ├── ordered_logit_ncats_v260413.stan
│   │   ├── ordered_logit_functions.stan
│   │   ├── partial_credit_model_ncats_v260413.stan
│   │   └── [compiled model binaries]
│   └── hpc/                                 # HPC submission scripts
│       ├── submit_hpc_job.sh
│       └── README_HPC.md
├── test/                                    # Test suite and legacy 2cats models
│   ├── R/                                   # Legacy 2cats R functions (for validation)
│   │   ├── fit_credit_model_2cats.R
│   │   ├── fit_ordered_logit_model_2cats.R
│   │   └── fit_partial_credit_model_2cats.R
│   ├── stan/                                # Legacy 2cats Stan models (for validation)
│   │   ├── continuation_ratio_model_2cats_v251218.stan
│   │   ├── credit_model_2cats_v251223.stan
│   │   ├── credit_model_2cats_v251224.stan
│   │   ├── ordered_logit_2cats_v251222.stan
│   │   ├── partial_credit_model_2cats_v251224.stan
│   │   └── [compiled model binaries]
│   ├── python/                              # Python test suite
│   │   ├── test_data_loading.py             # Test data loading functions
│   │   ├── test_utils.py                    # Test utility functions
│   │   ├── conftest.py                      # Pytest configuration
│   │   ├── fixtures/                        # Test fixtures
│   │   ├── utils/                           # Test utilities
│   │   ├── generate_baselines_colombia.py   # Generate baseline test data
│   │   ├── generate_r_baselines.py          # Generate R baseline data
│   │   └── pytest.ini                       # Pytest settings
│   ├── test_credit_ncats_correctness.R      # Verify ncats matches 2cats
│   ├── test_ordered_logit_ncats_correctness.R
│   ├── test_pcm_ncats_correctness.R
│   ├── run_all_tests.R                      # Run all R test suites
│   └── README.md                            # Test documentation
├── scripts-R/                               # R analysis scripts (formerly vignettes/)
│   ├── Colombia_analysis_for_HGpaper_ordered_logit_v251224.Rmd  # Hope Groups paper (ordered logit)
│   ├── Colombia_analysis_for_HGpaper_pcm_v251224.Rmd            # Hope Groups paper (PCM)
│   ├── Colombia_analysis_ordered_logit_ADVI_vs_HMC.Rmd          # Compare ADVI vs HMC
│   ├── Colombia_compare_models.Rmd          # Compare IRT models (local)
│   ├── Colombia_interim_analyses_parallel.R # Interim analyses (HPC parallel)
│   ├── Colombia_interim_analyses_v251007.Rmd
│   ├── Colombia_latent-factor-model_v250929.Rmd
│   ├── Colombia_validate_credit_model.Rmd
│   ├── hpc_Colombia_compare_models.Rmd      # Compare IRT models (HPC)
│   └── hpc_Colombia_train_test_item_response_models.Rmd  # Train/test (HPC)
├── scripts-py/                              # Python analysis scripts
│   └── Colombia_analysis_ordered_logit_ADVI_vs_HMC.py  # ADVI vs HMC comparison
├── dev/                                     # Development documentation
│   ├── IMPLEMENTATION_COMPLETE.md           # Implementation completion notes
│   ├── IMPLEMENTATION_SUMMARY.md            # Implementation summary
│   ├── PYTHON_PORT_STATUS.md                # Python port status
│   └── amortised_decision_making.md         # Amortised decision making documentation
├── old/                                     # Archived Stan models
│   └── stan/                                # Very old 2cats models (pre-reorganization)
├── python/                                  # Python implementation (JAX/NumPyro/CmdStanPy)
│   ├── __init__.py
│   ├── data_loading.py                      # Data loading functions
│   ├── model_fitting.py                     # Model fitting functions
│   ├── analysis.py                          # Analysis utilities
│   ├── test_environment.py                  # Python environment verification script
│   └── requirements.txt                     # Python package requirements
├── pixi.toml                                # Pixi environment specification (default + r-only)
├── pixi.lock                                # Pixi lock file (exact versions)
├── bIRTistic.yml                            # Conda environment specification (alternative)
├── job_config.csv.example                   # Example HPC job configuration
├── ENVIRONMENTS.md                          # Guide to pixi environments (default vs r-only)
├── README-python-port.md                    # Python implementation documentation
└── README.md                                # This file
```

</details>


## Troubleshooting

<details>
<summary>Click to expand</summary>

### CmdStan Installation Issues

If you encounter issues with CmdStan:

```r
# Check toolchain
cmdstanr::check_cmdstan_toolchain(fix = TRUE)

# Reinstall CmdStan
cmdstanr::install_cmdstan(overwrite = TRUE)
```

### Memory Issues on HPC

Adjust memory allocation in `src/hpc/submit_hpc_job.sh`:
```bash
#PBS -l mem=32gb  # Increase as needed
```

</details>

### Compilation Errors

Ensure C++ compiler is available:
```bash
# On Linux
which g++

# On macOS
which clang++
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

See [LICENSE](LICENSE) file for details.

## Citation

If you use this code in your research, please cite:

```
[Citation to be added]
```

## Contact

For questions or issues, please open an issue on GitHub or contact the maintainers.
