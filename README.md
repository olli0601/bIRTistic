# bIRTistic

Making Bayesian Item Response Theory models faster with amortised inference

## Overview

This repository contains tools for running Bayesian Item Response Theory (IRT) analyses using Stan, with support for parallel execution on HPC systems.

## Installation

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/) (recommended for faster installation)
- Git

### Step 1: Clone the repository

```bash
git clone https://github.com/olli0601/bIRTistic.git
cd bIRTistic
```

### Step 2: Create conda environment

Using Imperial HPC:
```bash
module load miniforge/3 tools/prod zlib libxml2/2.11.5-GCCcore-13.2.0
eval "$(~/miniforge3/bin/conda shell.bash hook)"
eval "$(mamba shell hook --shell bash)"
# install cmdstan only once and expose to other environment files
mamba create -n cmdstan-build -c conda-forge cmdstan=2.37.0 -y
# minimum conda environment to avoid installation from source
mamba env create -f bIRTistic.yml -y
mamba activate birtistic
```

### Step 3: Install R packages through R

After activating the environment, install R packages through R (possibly from source):

```bash
# installing data.table
# requires zlib, expose via module load
Rscript -e "packages <- c('data.table'); install.packages(packages[!packages %in% installed.packages()[,'Package']], repos = 'https://cloud.r-project.org')"
# installing cmdstanr
#Rscript -e "install.packages('cmdstanr', repos = c('https://stan-dev.r-universe.dev', 'https://cloud.r-project.org'), ask = false)"
# installing everything else 
# xml2 requires libxml2/2.11.5-GCCcore-13.2.0, expose via module load
Rscript -e "packages <- c('argparse', 'bayesplot', 'bookdown', 'formatR', 'GGally', 'ggpattern', 'ggplot2', 'ggsci', 'gridExtra', 'ggpubr', 'here', 'hexbin', 'jsonlite', 'kableExtra', 'knitr', 'loo', 'MatchIt', 'matrixStats', 'patchwork', 'polycor', 'posterior', 'purrr','rmarkdown', 'scales'); install.packages(packages[!packages %in% installed.packages()[,'Package']], repos = 'https://cloud.r-project.org', ask = false)"
```

### Step 4: Verify installation

Test that everything is working:

```bash
# check using the latest version of cmdstan
Rscript -e "library(cmdstanr); cmdstanr::cmdstan_version()"

# check Hello World Bernoulli model compiles and runs
Rscript -e "library(cmdstanr); file <- file.path(cmdstan_path(), 'examples', 'bernoulli', 'bernoulli.stan'); mod <- cmdstan_model(file, force_recompile = TRUE); data_list <- list(N = 10, y = c(0,1,0,0,0,0,0,0,0,1)); fit <- mod\$sample(data = data_list, seed = 123, chains = 2, parallel_chains = 1, refresh = 500);"
```

## Project Structure

```
bIRTistic/
├── R/                                      # R scripts
│   ├── fit_credit_model.R                   # Fit credit model
│   ├── fit_ordered_logit_model.R            # Fit ordered logit model
│   ├── fit_partial_credit_model.R           # Fit partial credit model
│   ├── fit_item_response_models.Rscript     # Generic IRT model fitting script
│   ├── hpc_submit_generic_job.R             # HPC job submission helper
│   ├── read_data_colombia.R                 # Read and preprocess Colombia data
│   ├── read_data_ukraine.R                  # Read and preprocess Ukraine data
│   ├── train_test_split_data.R              # Train/test data splitting
│   ├── train_test_item_response_models.Rscript  # Train/test IRT models
│   └── train_test_item_response_models_template.json  # Configuration template
├── src/
│   ├── stan/                                # Stan model files
│   │   ├── continuation_ratio_model_2cats_v251218.stan
│   │   ├── credit_model_2cats_v251223.stan
│   │   ├── credit_model_2cats_v251224.stan
│   │   ├── credit_model_functions.stan
│   │   ├── ordered_logit_2cats_v251222.stan
│   │   ├── ordered_logit_functions.stan
│   │   ├── partial_credit_model_2cats_v251224.stan
│   │   └── [compiled model files]
│   └── hpc/                                 # HPC submission scripts
│       ├── submit_hpc_job.sh
│       └── README_HPC.md
├── vignettes/                               # R Markdown analysis examples
│   ├── Colombia_analysis_for_HGpaper_pcm_v251224.Rmd
│   ├── Colombia_compare_models_2cat.Rmd
│   ├── Colombia_interim_analyses_parallel.R
│   ├── Colombia_interim_analyses_v251007.Rmd
│   ├── Colombia_latent-factor-model_v250929.Rmd
│   ├── Colombia_validate_credit_model.Rmd
│   ├── hpc_Colombia_compare_models_2cat.Rmd
│   └── hpc_Colombia_train_test_item_response_models.Rmd
└── bIRTistic.yml                            # Conda environment specification
```

## Quick Start

### Running Analysis Examples

Load precompiled libraries and activate shell hooks if needed:

```bash
module load miniforge/3 tools/prod zlib libxml2/2.11.5-GCCcore-13.2.0
eval "$(~/miniforge3/bin/conda shell.bash hook)"
eval "$(mamba shell hook --shell bash)"
```

Activate environment:

```bash
mamba activate birtistic
```

Examples:

See the `vignettes/` directory for example analyses:
- `Colombia_analysis_for_HGpaper_pcm_v251224.Rmd` - Hope Groups paper analysis using partial credit model
- `Colombia_compare_models_2cat.Rmd` - Compare IRT models on Colombia data, on laptop
- `hpc_Colombia_compare_models_2cat.Rmd` - Compare IRT models on Colombia data, run jobs in parallel on HPC
- - `Colombia_interim_analyses_v251007.Rmd` - Interim analyses on Colombia data, on laptop
- `Colombia_interim_analyses_parallel.R` - Interim analyses on Colombia data, run jobs in parallel on HPC
- `Colombia_validate_credit_model.Rmd` - Additional train/test IRT models that withholds all responses from some participants, on laptop
- `hpc_Colombia_train_test_item_response_models.Rmd` - Additional train/test IRT models that withholds all responses from some participants, run jobs in parallel on HPC

## Troubleshooting

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
