# bIRTistic

Making Bayesian Item Response Theory models faster with amortised inference.

This branch (`Colombia-analysis-pcm`) ships an **R + Stan only** stack focused on
the partial credit model (PCM) analyses for the Hope Groups Colombia data. The
Python / NumPyro / JAX implementation lives on the `main` branch.

## Overview

This repository contains tools for running Bayesian Item Response Theory (IRT)
analyses using Stan, with support for parallel execution on HPC systems.

## Installation

### Prerequisites

- [Pixi](https://pixi.sh/) (recommended — fast, cross-platform package manager)
  - Install: `curl -fsSL https://pixi.sh/install.sh | bash`
  - Or with Homebrew: `brew install pixi`
- Git
- C++ toolchain (auto-installed via the pixi `compilers` package)

### Setup

```bash
# Clone the repository
git clone https://github.com/olli0601/bIRTistic.git
cd bIRTistic
git checkout Colombia-analysis-pcm

# Install R + CmdStan dependencies via pixi
pixi install

# Install cmdstanr from R-universe, ggpattern from CRAN, and verify
pixi run setup

# (Equivalent to: pixi run install-cmdstanr && pixi run install-extra && pixi run verify-r)

# Activate the environment
pixi shell
```

`pixi run setup` runs three subtasks: `install-cmdstanr`, `install-extra`,
`verify-r`. The verify step compiles and samples the CmdStan `bernoulli`
example, confirming the toolchain is wired up.

## Quick Start

### Running R analyses

```bash
# Activate the environment
pixi shell

# Run our R unit tests
Rscript test/run_all_tests.R

# Render the Hope Groups PCM analysis
Rscript -e "rmarkdown::render('scripts-R/Colombia_analysis_for_HGpaper_pcm_v251224.Rmd')"
```

**Example analyses** (see `scripts-R/`):

- `Colombia_analysis_for_HGpaper_pcm_v251224.Rmd` — Hope Groups paper analysis using the partial credit model (focus of this branch)
- `Colombia_analysis_for_HGpaper_ordered_logit_v251224.Rmd` — Hope Groups paper analysis using the ordered logit model
- `Colombia_analysis_ordered_logit_ADVI_vs_HMC.Rmd` — Compare ADVI vs HMC inference methods
- `Colombia_compare_models.Rmd` — Compare IRT models on Colombia data (local)
- `hpc_Colombia_compare_models.Rmd` — Compare IRT models on Colombia data (HPC parallel)
- `Colombia_interim_analyses_v251007.Rmd` — Interim analyses (local)
- `Colombia_validate_credit_model.Rmd` — Train/test validation withholding participant responses
- `hpc_Colombia_train_test_item_response_models.Rmd` — Train/test validation (HPC parallel)

## Testing

```bash
# Run all R correctness tests (validates ncats models against 2cats baselines)
pixi shell
Rscript test/run_all_tests.R
```

## Project Structure

<details>
<summary>Click to expand</summary>

```
bIRTistic/
├── .vscode/                                 # VS Code workspace configuration
├── R/                                       # R scripts and functions
│   ├── fit_credit_model_ncats.R
│   ├── fit_ordered_logit_model_ncats.R
│   ├── fit_ordered_logit_model_ncats_advi.R
│   ├── fit_partial_credit_model_ncats.R
│   ├── get_endpoints.R
│   ├── read_data_colombia.R
│   ├── read_data_ukraine.R
│   ├── train_test_split_data.R
│   └── ...                                  # HPC helpers + templates
├── src/
│   ├── stan/                                # Current Stan models (ncats versions)
│   │   ├── credit_model_ncats_v260413.stan
│   │   ├── credit_model_functions.stan
│   │   ├── ordered_logit_ncats_v260413.stan
│   │   ├── ordered_logit_functions.stan
│   │   └── partial_credit_model_ncats_v260413.stan
│   └── hpc/                                 # HPC submission scripts
├── test/                                    # R test suite and legacy 2cats models
│   ├── R/                                   # Legacy 2cats R functions (for validation)
│   ├── stan/                                # Legacy 2cats Stan models (for validation)
│   ├── test_credit_ncats_correctness.R
│   ├── test_ordered_logit_ncats_correctness.R
│   ├── test_pcm_ncats_correctness.R
│   └── run_all_tests.R
├── scripts-R/                               # R analysis scripts (Rmd vignettes)
├── dev/                                     # Development notes
├── old/                                     # Archived Stan models
├── pixi.toml                                # Pixi environment specification (R + CmdStan)
├── pixi.lock
└── README.md
```

</details>

## Troubleshooting

<details>
<summary>Click to expand</summary>

### CmdStan Installation Issues

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

Ensure a C++ compiler is available:

```bash
# On Linux
which g++

# On macOS
which clang++
```

</details>

## Contributing

Contributions are welcome — please open a Pull Request.

## License

See [LICENSE](LICENSE) for details.

## Citation

```
[Citation to be added]
```

## Contact

For questions or issues, open an issue on GitHub or contact the maintainers.
