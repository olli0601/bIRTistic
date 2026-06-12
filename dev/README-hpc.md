# bIRTistic

Making Bayesian Item Response Theory models faster with amortised inference

## Overview

This README includes earlier notes to install and run on HPC systems

## Installation

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

### Memory Issues on HPC

Adjust memory allocation in `src/hpc/submit_hpc_job.sh`:
```bash
#PBS -l mem=32gb  # Increase as needed
```

