# bIRTistic

Making Bayesian Item Response Theory models faster with amortised inference

## Overview

This repository contains tools for running Bayesian Item Response Theory (IRT) analyses using Stan, with support for parallel execution on HPC systems.

## Installation

Please install our prerequisites:

- [Pixi](https://pixi.sh/) (recommended - fast, cross-platform package manager)
  - Install: `curl -fsSL https://pixi.sh/install.sh | bash`
  - Or with Homebrew: `brew install pixi`
- Alternatively: [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/)
- Git

We strongly recommend to use Pixi for installation:

```bash
# Clone the repository
git clone https://github.com/olli0601/bIRTistic.git
cd bIRTistic

# Install all dependencies
pixi install

# Setup cmdstanr etc from CRAN/R-universe
pixi run setup

# Verify installation (includes model compilation and sampling test)
pixi run verify
```

Then activate the environment:

```bash
# Activate the environment
pixi shell
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

#### Install R packages through R:

After activating the environment, install R packages through R (possibly from source):

```bash
# Installing data.table (requires zlib, expose via module load)
Rscript -e "packages <- c('data.table'); install.packages(packages[!packages %in% installed.packages()[,'Package']], repos = 'https://cloud.r-project.org')"

# Installing everything else (xml2 requires libxml2/2.11.5-GCCcore-13.2.0, expose via module load)
Rscript -e "packages <- c('argparse', 'bayesplot', 'bookdown', 'formatR', 'GGally', 'ggpattern', 'ggplot2', 'ggsci', 'gridExtra', 'ggpubr', 'here', 'hexbin', 'jsonlite', 'kableExtra', 'knitr', 'loo', 'MatchIt', 'matrixStats', 'patchwork', 'polycor', 'posterior', 'purrr','rmarkdown', 'scales'); install.packages(packages[!packages %in% installed.packages()[,'Package']], repos = 'https://cloud.r-project.org', ask = false)"
```

#### Verify installation:

The conda/mamba installation can be verified using the same tests that `pixi run verify` uses:

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

### Running Analysis Examples

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
- `Colombia_compare_models_2cat.Rmd` - Compare IRT models on Colombia data, on laptop
- `hpc_Colombia_compare_models_2cat.Rmd` - Compare IRT models on Colombia data, run jobs in parallel on HPC
- `Colombia_interim_analyses_v251007.Rmd` - Interim analyses on Colombia data, on laptop
- `Colombia_interim_analyses_parallel.R` - Interim analyses on Colombia data, run jobs in parallel on HPC
- `Colombia_validate_credit_model.Rmd` - Additional train/test IRT models that withholds all responses from some participants, on laptop
- `hpc_Colombia_train_test_item_response_models.Rmd` - Additional train/test IRT models, run jobs in parallel on HPC



## Project Structure

<details>
<summary>Click to expand</summary>
```
bIRTistic/
├── R/                                      # R scripts
│   ├── fit_credit_model_ncats.R             # Fit credit model 
│   ├── fit_ordered_logit_model_ncats.R      # Fit ordered logit model 
│   ├── fit_ordered_logit_model_ncats_advi.R # Fit ordered logit model with ADVI
│   ├── fit_partial_credit_model_ncats.R     # Fit partial credit model
│   ├── get_endpoints.R                      # Compute endpoint statistics from model fits
│   ├── fit_item_response_models.Rscript     # Generic IRT model fitting script 
│   ├── hpc_submit_generic_job.R             # HPC job submission helper
│   ├── read_data_colombia.R                 # Read and preprocess Colombia data
│   ├── read_data_ukraine.R                  # Read and preprocess Ukraine data
│   ├── train_test_split_data.R              # Train/test data splitting
│   ├── train_test_item_response_models.Rscript  # Train/test IRT models (uses R/fit_credit_model_ncats.R)
│   └── train_test_item_response_models_template.json  # Configuration template
├── src/
│   ├── stan/                                # Stan model files (ncats versions)
│   │   ├── credit_model_ncats_v260413.stan
│   │   ├── credit_model_functions.stan
│   │   ├── ordered_logit_ncats_v260413.stan
│   │   ├── ordered_logit_functions.stan
│   │   ├── partial_credit_model_ncats_v260413.stan
│   │   └── [compiled model files]
│   └── hpc/                                 # HPC submission scripts
│       ├── submit_hpc_job.sh
│       └── README_HPC.md
├── test/                                    # Test files and legacy 2cats models
│   ├── R/                                   # Legacy 2cats R functions (for testing)
│   │   ├── fit_credit_model_2cats.R
│   │   ├── fit_ordered_logit_model_2cats.R
│   │   └── fit_partial_credit_model_2cats.R
│   ├── stan/                                # Legacy 2cats Stan models (for testing)
│   │   ├── continuation_ratio_model_2cats_v251218.stan
│   │   ├── credit_model_2cats_v251223.stan
│   │   ├── credit_model_2cats_v251224.stan
│   │   ├── ordered_logit_2cats_v251222.stan
│   │   └── partial_credit_model_2cats_v251224.stan
│   ├── test_credit_ncats_correctness.R      # Verify ncats matches 2cats
│   ├── test_ordered_logit_ncats_correctness.R
│   ├── test_pcm_ncats_correctness.R
│   ├── run_all_tests.R                      # Run all test suites
│   └── README.md                            # Test documentation
├── vignettes/                               # R Markdown analysis examples
│   ├── Colombia_analysis_for_HGpaper_ordered_logit_v251224.Rmd  # Hope Groups paper (ordered logit)
│   ├── Colombia_analysis_for_HGpaper_pcm_v251224.Rmd            # Hope Groups paper (PCM)
│   ├── Colombia_analysis_ordered_logit_ADVI_vs_HMC.Rmd          # Compare ADVI vs HMC
│   ├── Colombia_compare_models_2cat.Rmd     # Compare IRT models 
│   ├── Colombia_interim_analyses_parallel.R
│   ├── Colombia_interim_analyses_v251007.Rmd
│   ├── Colombia_latent-factor-model_v250929.Rmd
│   ├── Colombia_validate_credit_model.Rmd
│   ├── hpc_Colombia_compare_models_2cat.Rmd # HPC model comparison 
│   └── hpc_Colombia_train_test_item_response_models.Rmd  # HPC train/test 
└── bIRTistic.yml                            # Conda environment specification
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
