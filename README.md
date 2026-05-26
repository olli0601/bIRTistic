# bIRTistic

Making Bayesian Item Response Theory models faster with amortised inference

## Overview

This repository contains tools for running Bayesian Item Response Theory (IRT)
analyses end-to-end in both R and Python. The project supports three IRT
models — partial credit, credit, and ordered logit — implemented twice for
parity: once in Stan (HMC + ADVI via `cmdstanr` / `cmdstanpy`) and once in
NumPyro (SVI with a choice of variational families). Helpers cover data
loading, model fitting, posterior summarisation, endpoint estimation, and
plotting, with optional parallel execution on HPC systems.

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

The Python helpers live in `python/` as plain modules (no `pip install`).
Either run one of the end-to-end analysis scripts in `scripts-py/`, or load
the helpers ad-hoc and call them yourself.

### End-to-end script (one command)

```bash
pixi shell
pixi run python scripts-py/Colombia_analysis_partial_credit_ADVI-HMC-pyro.py
```

This loads the Colombia CSV, fits the partial credit model with Stan HMC,
Stan ADVI, and five NumPyro SVI variants (`AutoDiagonalNormal`,
`AutoLaplaceApproximation`, `AutoMultivariateNormal`,
`AutoLowRankMultivariateNormal`, `AutoIAFNormal`), then writes endpoint
tables, comparison CSVs, posterior-predictive plots, and a timing chart to
`/Users/or105/sandbox/bIRTistic/py-colombia-...`. Pass `resume=True` is
already on, so re-runs only repeat post-fit comparisons.

### Ad-hoc API (interactive)

```python
import os
import sys
from pathlib import Path

# python/ is on sys.path so we can import the helpers directly.
sys.path.insert(0, str(Path("python").resolve()))

import numpy as np
from data_loading import read_data_colombia
from fit_partial_credit_model_ncats_stanhmc import (
    fit_partial_credit_model_ncats_stanhmc,
)

np.random.seed(42)

data_file = (
    "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/"
    "OR_Work/2025/2025_project_Hope_Groups/data/"
    "Colombia_data_baseline_endline_itemised_250927.csv"
)
raw = read_data_colombia(data_file)
dp_col, dit_col = raw["dp"].copy(), raw["dit"].copy()

# Preprocess: drop aggregate items, build item_time_id + oid + oidt.
dp1 = dp_col[~dp_col["item_label"].str.contains("agg")].copy()
dp1["y_stan"] = dp1["y"] + 1
dp1 = dp1.merge(dit_col[["item_label", "item_type"]], on="item_label", how="left")

item_time = (
    dp1[["item_type", "item_label", "time"]]
    .drop_duplicates()
    .sort_values(["item_type", "time", "item_label"])
    .reset_index(drop=True)
)
item_time["item_time_id"] = item_time.groupby("item_type").cumcount() + 1
dp1 = dp1.merge(item_time, on=["item_label", "time", "item_type"], how="left")
dp1 = dp1.merge(
    dit_col[["item_type", "item_type_id"]].drop_duplicates(), on="item_type", how="left"
)
dp1 = dp1.sort_values(["item_type_id", "pid", "time", "item_label"]).reset_index(drop=True)
dp1["oid"] = range(1, len(dp1) + 1)
dp1["oidt"] = dp1.groupby("item_type").cumcount() + 1

out_dir = "/tmp/birtistic_quickstart"
os.makedirs(out_dir, exist_ok=True)

result_hmc = fit_partial_credit_model_ncats_stanhmc(
    dit_col,
    dp1,
    output_file_prefix=os.path.join(out_dir, "pcm_1_hmc"),
    stan_file="src/stan/partial_credit_model_ncats_v260413.stan",
    chains=4,
    parallel_chains=4,
    threads_per_chain=1,
    iter_warmup=500,
    iter_sampling=1000,
    seed=123,
    x_formula="~ time - 1",
    resume=True,
    with_core_analyses=True,
    with_additional_analyses=True,
    show_messages=True,
)

print("Draws:", f"{out_dir}/pcm_1_hmc_draws.zarr")
print("Good chains:", result_hmc["good_chains"])
```

Swap `fit_partial_credit_model_ncats_stanhmc` for any of the other helpers
in `python/` (Stan ADVI, NumPyro SVI, credit model, ordered logit model) —
the call signature is the same. Endpoints and per-item probabilities can
then be computed with `get_endpoints.get_endpoints(...)`.

## Running R Analysis Examples

```bash
pixi shell

# Run the R unit test suite (R + Stan correctness checks).
Rscript test/run_all_tests.R

# Render any Rmd in scripts-R/. Example:
Rscript -e "rmarkdown::render('scripts-R/Colombia_analysis_for_HGpaper_pcm_v251224.Rmd')"
```

Available analyses (`scripts-R/`):

- `Colombia_analysis_for_HGpaper_ordered_logit_v251224.Rmd` — Hope Groups paper analysis using ordered logit model
- `Colombia_analysis_for_HGpaper_pcm_v251224.Rmd` — Hope Groups paper analysis using partial credit model
- `Colombia_analysis_ordered_logit_ADVI_vs_HMC.Rmd` — Compare ADVI vs HMC for the ordered logit model
- `Colombia_compare_models.Rmd` — Compare IRT models on the Colombia data (local)
- `hpc_Colombia_compare_models.Rmd` — Same comparison, HPC-parallel version
- `Colombia_interim_analyses_v251007.Rmd` — Interim analyses on the Colombia data
- `Colombia_latent-factor-model_v250929.Rmd` — Latent factor model exploration
- `Colombia_validate_credit_model.Rmd` — Train/test validation withholding participants
- `hpc_Colombia_train_test_item_response_models.Rmd` — Train/test validation, HPC-parallel

## Running Python Analysis Examples

```bash
pixi shell

# Run the Python test suite (numpyro/Stan parity + Laplace Hessian checks).
pixi run python -m pytest test/python -v

# Run any analysis script end-to-end. Example:
pixi run python scripts-py/Ukraine_analysis_partial_credit_ADVI-HMC-pyro.py
```

Available analyses (`scripts-py/`):

- `Colombia_analysis_partial_credit_ADVI-HMC-pyro.py` — Partial credit model on Colombia data; fits Stan HMC + Stan ADVI + 5 NumPyro SVI variants; saves endpoint, probability, ypred, timing and parameter comparisons.
- `Colombia_analysis_ordered_logit_ADVI-HMC-pyro.py` — Same harness for the ordered logit model on Colombia data.
- `Ukraine_analysis_partial_credit_ADVI-HMC-pyro.py` — Partial credit model on Ukraine data (uses `categorical_threshold=2`).
- `Ukraine_analysis_ordered_logit_ADVI-HMC-pyro.py` — Ordered logit model on Ukraine data.
- `Colombia_analysis_ordered_logit_ADVI_vs_HMC.py` — Legacy ADVI vs HMC comparison (single comparison, no NumPyro).

Each script writes a complete result directory under
`/Users/or105/sandbox/bIRTistic/...` containing per-method `*_draws.zarr`,
`*_timing.csv`, and `comparison_*.{csv,pdf}` artifacts.

## Testing

```bash
pixi shell

# Python (NumPyro / Stan parity + Laplace Hessian eigenvalue checks).
pixi run python -m pytest test/python -v

# R (ncats correctness against legacy 2cats baselines).
Rscript test/run_all_tests.R
```


## Project Structure

<details>
<summary>Click to expand</summary>

```
bIRTistic/
├── .vscode/                                 # VS Code workspace configuration
├── R/                                       # R scripts and functions (ncats only)
│   ├── fit_credit_model_ncats.R
│   ├── fit_ordered_logit_model_ncats.R
│   ├── fit_ordered_logit_model_ncats_advi.R
│   ├── fit_partial_credit_model_ncats.R
│   ├── get_endpoints.R                       # Compute endpoint stats from posteriors
│   ├── read_data_colombia.R
│   ├── read_data_ukraine.R
│   ├── train_test_split_data.R
│   ├── fit_item_response_models.Rscript      # Generic HPC fitting entrypoint
│   ├── hpc_submit_generic_job.R              # HPC job submission helper
│   ├── train_test_item_response_models.Rscript
│   └── train_test_item_response_models_template.json
├── python/                                   # Python implementation (CmdStanPy + NumPyro)
│   ├── data_loading.py                       # read_data_colombia / read_data_ukraine
│   ├── utils.py                              # Stan-data builders, summaries, plotting
│   ├── get_endpoints.py                      # Endpoint computation from zarr draws
│   ├── fit_partial_credit_model_ncats_stanhmc.py
│   ├── fit_partial_credit_model_ncats_stanadvi.py
│   ├── fit_partial_credit_model_ncats_pyrosvi.py
│   ├── fit_credit_model_ncats_stanhmc.py
│   ├── fit_credit_model_ncats_stanadvi.py
│   ├── fit_credit_model_ncats_pyrosvi.py
│   ├── fit_ordered_logit_model_ncats_stanhmc.py
│   ├── fit_ordered_logit_model_ncats_stanadvi.py
│   ├── fit_ordered_logit_model_ncats_pyrosvi.py
│   └── requirements.txt
├── src/
│   ├── stan/                                 # Stan ncats models + shared functions
│   │   ├── credit_model_functions.stan
│   │   ├── credit_model_ncats_v260413.stan
│   │   ├── ordered_logit_functions.stan
│   │   ├── ordered_logit_ncats_v260413.stan
│   │   └── partial_credit_model_ncats_v260413.stan
│   ├── numpyro/                              # NumPyro implementations of the ncats models
│   │   ├── credit_model_ncats_v260413.pyro
│   │   ├── ordered_logit_ncats_v260413.pyro
│   │   ├── partial_credit_model_ncats_v260413.pyro
│   │   └── partial_credit_model_2cats_v251224.pyro  # legacy 2cats baseline
│   └── hpc/                                  # HPC submission scripts
├── test/                                     # Test suite and legacy 2cats models
│   ├── R/                                    # Legacy 2cats R helpers (for validation)
│   ├── stan/                                 # Legacy 2cats Stan models (for validation)
│   ├── python/                               # Python test suite
│   │   ├── test_data_loading.py
│   │   ├── test_utils.py
│   │   ├── test_numpyro_pcm_ncats_stanhmc.py            # Stan parity (partial credit)
│   │   ├── test_numpyro_credit_ncats_stanhmc.py         # Stan parity (credit)
│   │   ├── test_numpyro_ordered_logit_ncats_stanhmc.py  # Stan parity (ordered logit)
│   │   ├── test_numpyro_laplace_hessian.py              # AutoLaplace Hessian rank
│   │   ├── test_numpyro_pcm_2cats.py
│   │   ├── test_environment.py
│   │   ├── conftest.py, pytest.ini, fixtures/, utils/
│   │   ├── generate_baselines_colombia.py
│   │   └── generate_r_baselines.py
│   ├── test_credit_ncats_correctness.R
│   ├── test_ordered_logit_ncats_correctness.R
│   ├── test_pcm_ncats_correctness.R
│   ├── run_all_tests.R
│   └── README.md
├── scripts-R/                                # R analysis scripts (Rmd vignettes)
│   ├── Colombia_analysis_for_HGpaper_ordered_logit_v251224.Rmd
│   ├── Colombia_analysis_for_HGpaper_pcm_v251224.Rmd
│   ├── Colombia_analysis_ordered_logit_ADVI_vs_HMC.Rmd
│   ├── Colombia_compare_models.Rmd
│   ├── Colombia_interim_analyses_v251007.Rmd
│   ├── Colombia_latent-factor-model_v250929.Rmd
│   ├── Colombia_validate_credit_model.Rmd
│   ├── hpc_Colombia_compare_models.Rmd
│   ├── hpc_Colombia_train_test_item_response_models.Rmd
│   └── scripts-R-old/                        # Archived early prototypes
├── scripts-py/                               # Python analysis scripts
│   ├── Colombia_analysis_partial_credit_ADVI-HMC-pyro.py
│   ├── Colombia_analysis_ordered_logit_ADVI-HMC-pyro.py
│   ├── Ukraine_analysis_partial_credit_ADVI-HMC-pyro.py
│   ├── Ukraine_analysis_ordered_logit_ADVI-HMC-pyro.py
│   ├── Colombia_analysis_ordered_logit_ADVI_vs_HMC.py  # legacy ADVI-vs-HMC
│   └── README.md
├── dev/                                      # Development documentation
├── old/                                      # Archived Stan models (pre-ncats)
├── pixi.toml                                 # Pixi environment specification
├── pixi.lock                                 # Pixi lock file
├── job_config.csv.example                    # Example HPC job configuration
└── README.md                                 # This file
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

### Compilation Errors

Ensure a C++ compiler is available (pixi installs `compilers` automatically):
```bash
# On Linux
which g++

# On macOS
which clang++
```

</details>

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
