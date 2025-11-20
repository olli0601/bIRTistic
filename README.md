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

Using conda:
```bash
conda env create -f bIRTistic.yml
conda activate birtistic
```

Or using mamba (faster):
```bash
mamba env create -f bIRTistic.yml
mamba activate birtistic
```

### Step 3: Install R packages not available in conda

After activating the environment, open R and install the remaining packages:

```r
# Install cmdstanr from Stan's R-universe
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))

# Verify CmdStan installation
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan()

# Install any other missing packages
packages <- c("data.table", "ggplot2", "ggsci", "bayesplot", 
              "argparse", "here", "knitr", "posterior")
install.packages(packages[!packages %in% installed.packages()[,"Package"]])
```

### Step 4: Verify installation

Test that everything is working:

```bash
Rscript -e "library(cmdstanr); cmdstanr::cmdstan_version()"
```

## Project Structure

```
bIRTistic/
├── R/                           # R scripts
│   └── credit_model_run_analysis.Rscript  # Main analysis script
├── src/
│   ├── stan/                    # Stan model files
│   │   ├── credit_model_2cats_v251120.stan
│   │   └── credit_model_functions.stan
│   └── hpc/                     # HPC submission scripts
│       ├── submit_hpc_job.sh
│       └── README_HPC.md
├── vignettes/                   # R Markdown analysis examples
└── environment.yml              # Conda environment specification
```

## Quick Start

### Running a Single Analysis

```bash
Rscript R/credit_model_run_analysis.Rscript \
  --stan_file src/stan/credit_model_2cats_v251120.stan \
  --stan_include_dir src/stan \
  --dit_file data/dit.rds \
  --dcati_file data/dcati.rds \
  --job_id 1 \
  --output_dir output/results \
  --output_prefix analysis \
  --chains 4 \
  --parallel_chains 4 \
  --iter_warmup 500 \
  --iter_sampling 1500 \
  --seed 123
```

### Running on HPC

See [src/hpc/README_HPC.md](src/hpc/README_HPC.md) for detailed instructions on:
- Preparing data for HPC execution
- Submitting batch jobs with PBS/SLURM
- Running multiple analyses in parallel

## Usage

### Command Line Arguments

**Required:**
- `--stan_file`: Path to Stan model file
- `--stan_include_dir`: Directory containing Stan include files
- `--dit_file`: Path to item definitions (dit.rds)
- `--dcati_file`: Path to processed data (dcati.rds)
- `--job_id`: Integer ID for this analysis job
- `--output_dir`: Directory for output files
- `--output_prefix`: Prefix for output filenames

**Optional:**
- `--chains`: Number of MCMC chains (default: 2)
- `--parallel_chains`: Number of parallel chains (default: 2)
- `--iter_warmup`: Warmup iterations (default: 500)
- `--iter_sampling`: Sampling iterations (default: 1500)
- `--seed`: Random seed for MCMC (default: 123)
- `--with_additional_analyses`: Include diagnostic plots (default: false)

## Examples

See the `vignettes/` directory for example analyses:
- `Colombia_interim_analyses_v251007.Rmd` - Complete interim analysis workflow

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
