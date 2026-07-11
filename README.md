# bIRTistic

Making Bayesian Item Response Theory models faster with amortised inference

## Overview

This repository contains tools for running Bayesian Item Response Theory (IRT)
analyses end-to-end in both R and Python. The project supports three IRT
models — partial credit, credit, and ordered logit — and Binomial and Multivariatenormal models 
that are analytically tractable. Each model has methods implementing NUTS and SVI inference.
For interim analyses, nested HMC is available; as
well as importance sampling, importance sampling with moment matching, and SMC sampling of the 
posterior over future data; and regression-based estimation of endpoints and the PPS.

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

This repository's two analysis stacks live on separate branches:

1. **`main`** (default branch) — full Python stack (JAX, NumPyro,
   cmdstanpy) for interim-analysis workflows. The
   `default` pixi environment installs everything needed; the
   `mps-experimental` environment adds the JAX MPS backend for Apple
   Silicon.
2. **`Colombia-analysis-pcm`** — R + Stan code for the Villaveces et al
   preprint *'Parenting with Hope' program among bereaved families in
   Colombia: A pre-post and quasi-experimental evaluation*. This branch
   has its own R-only `pixi.toml`; check it out separately if you want
   to reproduce the preprint analyses.

Option 1: Install default environment (Python full stack)

```bash
# Clone the repository
git clone https://github.com/olli0601/bIRTistic.git
cd bIRTistic

# Install Python + cmdstan (default branch is main).
pixi install

# Verify Python installation (includes CmdStanPy + 8 schools test).
pixi run verify-python

# Activate the environment.
pixi shell
```

Option 2: Reproduce the Villaveces et al Colombia PCM preprint (R + Stan)

The preprint code lives on the `Colombia-analysis-pcm` branch and ships
its own R-only `pixi.toml`. Check out the branch before installing so
the `scripts-R/` Rmd vignettes and `R/` helpers are on disk.

```bash
# Clone the repository.
git clone https://github.com/olli0601/bIRTistic.git
cd bIRTistic

# Switch to the preprint branch (R-only).
git checkout Colombia-analysis-pcm

# Install R + CmdStan dependencies.
pixi install

# Setup cmdstanr etc from CRAN/R-universe.
pixi run setup

# Verify installation (includes Bernoulli model test).
pixi run verify-r

# Activate the environment.
pixi shell
```

## Quick Start

The Python helpers live in `python/` as plain modules (no `pip install`).
Either run one of the end-to-end analysis scripts in `scripts-py/`, or load
the helpers ad-hoc and call them yourself.


### Ad-hoc API (interactive)

Given a current data set ``dp_x`` (one row per Bernoulli trial), fit the
Binomial x-posterior with NumPyro NUTS, then estimate the predictive
probability of success (PPS) by nested-MC HMC over ``S`` hypothetical
future data sets.

```python
import os
import sys
from pathlib import Path

# python/ is on sys.path so we can import the helpers directly.
sys.path.insert(0, str(Path("python").resolve()))

import numpy as np
import pandas as pd
from model_binomial import BinomialModel
from fit_interim import fit_interim_posterior_xz_with_nested_monte_carlo

out_dir = "/tmp/birtistic_quickstart"
os.makedirs(out_dir, exist_ok=True)
seed = 123
rng = np.random.default_rng(seed)

# Synthetic current cohort xi: 80 Bernoulli trials with true p = 0.4.
n_obs = 80
m_future = 200                            # remaining trials at end of trial
dp_x = pd.DataFrame({
    'pid': np.arange(n_obs),
    'y':   rng.binomial(1, 0.4, size=n_obs),
    'oid': np.arange(1, n_obs + 1),
})
dit = pd.DataFrame({
    'item_label':       ['intervention'],
    'item_type':        ['binomial'],
    'item_high_label':  ['higher_is_better'],
})

# 1) Fit x-posterior via NumPyro NUTS.
model = BinomialModel(
    dit=dit, dcati=dp_x,
    prior_a=1.0, prior_b=1.0, p_0=0.5, seed=seed,
)
prefix = os.path.join(out_dir, 'binomial_x')
fit_x = model.fit_pyro_hmc(
    output_file_prefix=prefix,
    chains=2, iter_warmup=500, iter_sampling=1000,
    save_to_file=True, resume=False, verbose=False,
)

# 2) Expand to future data: draw S posterior-predictive samples zi_s of
#    the m_future missing trials.
S = 200
zi = model.get_interim_z_from_ypredi(
    f"{prefix}_draws.zarr", m_future,
    pps_z_total=S, seed=seed,
)

# 3) Nested-MC PPS: for each s in 1..S refit the posterior on
#    (xi + zi_s) and record P(H_1 | x, z_s).
p_h1_xz = fit_interim_posterior_xz_with_nested_monte_carlo(
    model, zi,
    interim_method_args={
        'pps_z_total':        S,
        'pps_H1_min_effect_size_thresh':         0.25,        # H_1: 1 - p/p_0 > 0.25
        'pps_ProbH1_target_lwr_quantile':  0.89,
        'fit_method':         'fit_pyro_hmc',
        'seed':               seed,
        'save_to_file':       False,
        'verbose':            False,
        'cpu_n':              1,
        'output_file_prefix': os.path.join(out_dir, 'binomial_pps'),
    },
    fit_method_args={
        'x_formula':     '~ 1',
        'chains':        2,
        'iter_warmup':   500,
        'iter_sampling': 1000,
    },
)

# 4) Aggregate: PPS = fraction of samples crossing the decision threshold.
pps = float((p_h1_xz['p_h1_xz'] > 0.89).mean())
print(f"PPS = {pps:.3f}")
```

Swap `fit_pyro_hmc` for `fit_pyro_svi` (NumPyro SVI) or `fit_stan_hmc` /
`fit_stan_svi` (cmdstan) for a different inference backend. The same
pattern works with the IRT models (`PartialCreditModel`, `CreditModel`,
`OrderedLogit`).

## Running Python Analysis Examples

```bash
pixi shell

# Run the Python test suite (numpyro/Stan parity + Laplace Hessian checks).
pixi run python -m pytest test/python -v

# Run any analysis script end-to-end. Example:
pixi run python scripts-py/Binomial_interim_analyses_with_nested_MC_hmc.py
```

### Binomial benchmark

Beta-Bernoulli simulation: ``N = 500`` Bernoulli outcomes with true
``p = 0.4`` over one year of monthly interim cutoffs. The closed-form
Beta-Binomial tail sum is the analytic PPS reference; each script
estimates the same PPS by a different algorithm.

```bash
pixi shell

# Nested-MC HMC: per interim, refit the posterior on (xi + z_s) for S
# future-data samples and average the decision indicator. ~30 min.
pixi run python scripts-py/Binomial_interim_analyses_with_nested_MC_hmc.py

# Same nested-MC pattern but with NumPyro SVI (AutoDiagonalNormal)
# per refit (faster, slightly less accurate). ~10 min.
pixi run python scripts-py/Binomial_interim_analyses_with_nested_MC_svi_autodiagnormal.py

# Self-normalised IS: reweight the x-posterior draws by
# softmax_k log p(z_s | theta_k); no refits. <1 min.
pixi run python scripts-py/Binomial_interim_analyses_with_IS_from_x.py

# Strong-Oakley binomial-GLM regression of 1{theta in H_1} on a
# per-(item, draw) summary w(z). ~1 min.
pixi run python scripts-py/Binomial_interim_analysis_regression_H1x_on_wz.py

# Strong-Oakley Gaussian-GLM regression of the continuous endpoint
# 1 - p/p_0 on w(z); Phi-tail conversion to label probability. ~1 min.
pixi run python scripts-py/Binomial_interim_analysis_regression_endptx_on_wz.py

# Cross-method comparison: PPS boxplot, PPS bars with bootstrap CI,
# timing, ESS / particle for IS. <1 min.
pixi run python scripts-py/Binomial_interim_analyses_compare_methods.py
```

### Multivariate normal benchmark 

The MVN benchmark is split into a sequence of standalone scripts so each
phase (simulation, interim-data construction, inference, plotting) can
be run independently.

```bash
pixi shell

# 1) Simulate one cohort per J in {20, 60, 100} under the block-
#    equicorrelation K. Persists dp, closed-form per-component
#    P(H_1j | x) and closed-form PPS, plus diagnostic + heatmap plots.
#    <1 min.
pixi run python scripts-py/MVN_interim_analyses_make_sim_data.py

# 2) Fit the outer x-posterior per (J, interim) (resume from cached
#    zarrs), expand to S future-data samples zi, persist the long-form
#    zi + the first S posterior mu draws in mvn_J{J}_interim_data.pkl.
#    First run: ~15 min (HMC across J=20/60/100, 7 interims). Resume:
#    a few seconds.
pixi run python scripts-py/MVN_interim_analyses_make_interim_data.py

# 3) Nested-MC HMC PPS using the cached zi; emits per-J boxplot /
#    bars / scatter / _all / analytic-vs-HMC plots. ~10 hours full
#    deploy (J=20/60/100, S=200, 7 interims); reads cached pkl + only
#    regenerates plots if mvn_pps_nested_mc.pkl exists.
pixi run python scripts-py/MVN_interim_analyses_with_nested_MC_hmc.py

# 4) Self-normalised IS using a fresh NumPyro SVI
#    (AutoMultivariateNormal) variational posterior; vectorised
#    sufficient-statistic implementation of log p(z_s | theta_k).
#    ~1 min total (SVI fit + IS reweight per interim).
pixi run python scripts-py/MVN_interim_analyses_with_IS_from_x.py

# 5) Strong-Oakley Gaussian-GLM regression endpt-x on w(z); uses the
#    cached zi + mu_draws. <1 min.
pixi run python scripts-py/MVN_interim_analysis_regression_endptx_on_wz.py

# 6) Cross-J / cross-method comparison: per-J PPS bars + boxplot +
#    timing + ESS + cross-J MSE-vs-J plot. <1 min.
pixi run python scripts-py/MVN_interim_analyses_compare_methods.py
```

Each script writes its artifacts to a dedicated directory under
`/Users/or105/sandbox/bIRTistic/...` (sim / nested-MC / IS / RGE /
compare). The Binomial and MVN scripts share the same
`fit_interim_posterior_xz_with_nested_monte_carlo`,
`fit_interim_regress_endptx_on_wz`, etc. drivers from `fit_interim.py`.

### IRT analyses on Colombia Hope Groups data

```bash
pixi shell

# Partial credit model: Stan HMC + Stan ADVI + 5 NumPyro SVI variants
# (AutoDiagonalNormal, AutoLaplaceApproximation, AutoMultivariateNormal,
# AutoLowRankMultivariateNormal, AutoIAFNormal). ~2 hours full sweep.
pixi run python scripts-py/Colombia_analysis_partial_credit_ADVI-HMC-pyro.py

# Same harness for the ordered logit model. ~2 hours.
pixi run python scripts-py/Colombia_analysis_ordered_logit_ADVI-HMC-pyro.py

# ADVI vs HMC comparison (single comparison, no NumPyro). ~15 min.
pixi run python scripts-py/Colombia_analysis_ordered_logit_ADVI_vs_HMC.py

# Interim analyses on the Colombia partial-credit case study. ~30 min.
pixi run python scripts-py/Colombia_interim_analyses.py
```

### IRT analyses on Ukraine Hope Groups data

```bash
pixi shell

# Partial credit fit harness (uses categorical_threshold=2). ~2 hours.
pixi run python scripts-py/Ukraine_analysis_partial_credit_ADVI-HMC-pyro.py

# Ordered logit fit harness. ~2 hours.
pixi run python scripts-py/Ukraine_analysis_ordered_logit_ADVI-HMC-pyro.py

# Interim PPS pipeline (full case study):

# Nested-MC HMC: fit on (xi + z_s). ~6 hours.
pixi run python scripts-py/Ukraine_interim_analyses_with_nested_MC.py

# Self-normalised IS reweighting. <5 min.
pixi run python scripts-py/Ukraine_interim_analysis_with_IS_from_x.py

# IS with moment matching (Paananen 2021). <10 min.
pixi run python scripts-py/Ukraine_interim_analysis_with_IS_moment_matching_from_x.py

# SMC sampler with resample-move (annealed). ~3 hours.
pixi run python scripts-py/Ukraine_interim_analysis_with_SMC_resample_from_x.py

# Strong-Oakley regression of 1{theta in H_1} on w(z). ~10 min.
pixi run python scripts-py/Ukraine_interim_analysis_regression_H1x_on_wz.py

# Strong-Oakley regression of endpoint on w(z). ~10 min.
pixi run python scripts-py/Ukraine_interim_analysis_regression_endptx_on_wz.py

# Cross-method comparison: PPS / boxplot / timing / ESS plots. <1 min.
pixi run python scripts-py/Ukraine_interim_analysis_compare_methods.py
```

## Testing

```bash
pixi shell

# Python (NumPyro / Stan parity + Laplace Hessian eigenvalue checks).
pixi run python -m pytest test/python -v
```

The R + Stan test suite for the PCM preprint lives on the
`Colombia-analysis-pcm` branch; check out that branch and run
``Rscript test/run_all_tests.R`` there.


## Project Structure

<details>
<summary>Click to expand</summary>

```
bIRTistic/
├── .vscode/                                 # VS Code workspace configuration
├── python/                                   # Model classes + fitters + interim drivers
│   ├── model.py                              # Abstract Model + IRTModel base
│   ├── model_binomial.py                     # Beta-Bernoulli (§3.1) closed-form + HMC/SVI
│   ├── model_mvn.py                          # Multivariate normal (§3.3) closed-form + HMC/SVI
│   ├── model_irt.py                          # IRT mixin (per-time endpoints, w(z))
│   ├── model_pcm.py                          # PartialCreditModel
│   ├── model_credit.py                       # CreditModel
│   ├── model_ordered_logit.py                # OrderedLogit
│   ├── fit_interim.py                        # Interim drivers: nested-MC, IS, MM, SMC,
│   │                                         # regression-H1x, regression-endptx
│   ├── data_loading.py                       # read_data_colombia / read_data_ukraine
│   ├── utils.py                              # Stan-data builders, autoguide factory,
│   │                                         # palettes (futurama, material)
│   └── requirements.txt
├── src/
│   ├── stan/                                 # Stan model files
│   │   ├── binomial_v260603.stan
│   │   ├── credit_model_functions.stan
│   │   ├── credit_model_ncats_v260413.stan
│   │   ├── ordered_logit_functions.stan
│   │   ├── ordered_logit_ncats_v260413.stan
│   │   └── partial_credit_model_ncats_v260413.stan
│   └── numpyro/                              # NumPyro model programs
│       ├── binomial_v260603.pyro
│       ├── mvn_v260609.pyro
│       ├── credit_model_ncats_v260413.pyro
│       ├── ordered_logit_ncats_v260413.pyro
│       └── partial_credit_model_ncats_v260413.pyro
├── test/python/                              # Python test suite
│   ├── test_model_abc.py                     # Model ABC contract
│   ├── test_model_binomial.py                # BinomialModel + closed-form PPS
│   ├── test_model_pcm.py                     # PartialCreditModel
│   ├── test_model_credit_ordered_logit.py    # Credit + OrderedLogit
│   ├── test_interim_alignment.py
│   ├── test_interim_determinism.py
│   ├── test_interim_helpers.py
│   ├── test_get_endpoints.py
│   └── test_no_legacy_imports.py
├── scripts-py/                               # Python analysis scripts
│   ├── Binomial_interim_analyses_with_nested_MC_hmc.py
│   ├── Binomial_interim_analyses_with_nested_MC_svi_autodiagnormal.py
│   ├── Binomial_interim_analyses_with_IS_from_x.py
│   ├── Binomial_interim_analysis_regression_H1x_on_wz.py
│   ├── Binomial_interim_analysis_regression_endptx_on_wz.py
│   ├── Binomial_interim_analyses_compare_methods.py
│   ├── MVN_interim_analyses_make_sim_data.py
│   ├── MVN_interim_analyses_make_interim_data.py
│   ├── MVN_interim_analyses_with_nested_MC_hmc.py
│   ├── MVN_interim_analyses_with_IS_from_x.py
│   ├── MVN_interim_analysis_regression_endptx_on_wz.py
│   ├── MVN_interim_analyses_compare_methods.py
│   ├── Ukraine_interim_analyses_with_nested_MC.py
│   ├── Ukraine_interim_analysis_with_IS_from_x.py
│   ├── Ukraine_interim_analysis_with_IS_moment_matching_from_x.py
│   ├── Ukraine_interim_analysis_with_SMC_resample_from_x.py
│   ├── Ukraine_interim_analysis_regression_H1x_on_wz.py
│   ├── Ukraine_interim_analysis_regression_endptx_on_wz.py
│   ├── Ukraine_interim_analysis_compare_methods.py
│   ├── Colombia_analysis_partial_credit_ADVI-HMC-pyro.py
│   ├── Colombia_analysis_ordered_logit_ADVI-HMC-pyro.py
│   ├── Colombia_analysis_ordered_logit_ADVI_vs_HMC.py
│   ├── Colombia_interim_analyses.py
│   ├── Ukraine_analysis_partial_credit_ADVI-HMC-pyro.py
│   ├── Ukraine_analysis_ordered_logit_ADVI-HMC-pyro.py
│   └── README.md
├── dev/                                      # Development documentation, including
│   │                                         # `amortised_decision_making.md` (§3.x specs)
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
