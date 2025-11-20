# HPC Analysis Scripts

This directory contains scripts for running analyses on arbitrary datasets on HPC/cloud systems.

## Files

1. **`credit_model_run_analysis.Rscript`** (in `R/` directory) - Self-contained R script that runs a single analysis on any dataset
2. **`prepare_data_for_hpc.R`** - Helper script to export data and create job configuration (to be created in `R/` directory)
3. **`submit_hpc_job.sh`** - PBS batch submission script (in `src/hpc/` directory)

## Workflow

### 1. Prepare Data

First, run the main R Markdown document up to the data preparation section, then export the data:

```r
# In R, after running the Rmd up to the data preparation:
source("R/prepare_data_for_hpc.R")
```

This creates:
- A `data_for_hpc/` directory with:
  - `dit.rds` - Item definitions
  - `dcati_job_*.rds` - Pre-processed data for each analysis job
- A `job_config.csv` file with columns:
  - `JOB_ID` - Integer identifier for each job
  - `DCATI_FILE` - Filename of the corresponding dcati RDS file

### 2. Transfer Files to HPC

Transfer these files to your HPC system:
```bash
# Example using rsync
rsync -avz data_for_hpc/ user@hpc.system:/path/to/project/data_for_hpc/
rsync -avz job_config.csv user@hpc.system:/path/to/project/
rsync -avz src/ user@hpc.system:/path/to/project/src/
rsync -avz R/ user@hpc.system:/path/to/project/R/
```

### 3. Submit Jobs

On the HPC system:

```bash
# Make scripts executable
chmod +x src/hpc/submit_hpc_job.sh

# Edit submit_hpc_job.sh to set correct paths and array size
nano src/hpc/submit_hpc_job.sh

# Update the array size to match number of jobs in job_config.csv
# For example, if you have 11 jobs: #PBS -t 1-11

# Submit array job for all analyses (PBS)
qsub src/hpc/submit_hpc_job.sh
```

Or run a single analysis manually:

```bash
Rscript R/credit_model_run_analysis.Rscript \
  --stan_file src/stan/credit_model_2cats_v251120.stan \
  --stan_include_dir src/stan \
  --dit_file data_for_hpc/dit.rds \
  --dcati_file data_for_hpc/dcati_job_1.rds \
  --job_id 1 \
  --output_dir output/results \
  --output_prefix analysis \
  --chains 4 \
  --parallel_chains 4 \
  --iter_warmup 500 \
  --iter_sampling 1500 \
  --seed 123
```

## Command Line Arguments

### Required Arguments:
- `--stan_file`: Path to Stan model file
- `--stan_include_dir`: Directory containing Stan include files
- `--dit_file`: Path to dit.rds file
- `--dcati_file`: Path to dcati RDS file for this job
- `--job_id`: Integer ID for this analysis job
- `--output_dir`: Directory for output files
- `--output_prefix`: Prefix for output filenames

### Optional Arguments (with defaults):
- `--chains`: Number of MCMC chains (default: 2)
- `--parallel_chains`: Number of parallel chains (default: 2)
- `--iter_warmup`: Warmup iterations (default: 500)
- `--iter_sampling`: Sampling iterations (default: 1500)
- `--seed`: Random seed for MCMC (default: 123)
- `--with_additional_analyses`: Include additional diagnostic plots (default: false)

## Output Files

For each analysis job, the script generates:
- `*_stan.rds` - Fitted Stan model object
- `*_convergence_mixing.csv` - Convergence diagnostics
- `*_worsttrace.png` - Trace plot for worst-mixing parameters (if --with_additional_analyses)
- `*_intervals.png` - Parameter interval plot (if --with_additional_analyses)
- `*_areas.png` - Parameter density plot (if --with_additional_analyses)
- `*_ppcheck.png` - Posterior predictive check (if --with_additional_analyses)
- `*_probs_barplot_v2.png/.pdf` - Probability plots

## Adapting for Different HPC Systems

The default `submit_hpc_job.sh` script uses PBS/Torque. For other schedulers:

### SLURM
Replace PBS directives with SLURM equivalents:
```bash
#SBATCH --job-name=interim_analysis
#SBATCH --output=logs/interim_%a.out
#SBATCH --error=logs/interim_%a.err
#SBATCH --array=1-11
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
```
And change `${PBS_ARRAYID}` to `${SLURM_ARRAY_TASK_ID}`, then submit with `sbatch`.

### SGE
```bash
#$ -N interim_analysis
#$ -t 1-11
#$ -l h_rt=24:00:00
#$ -l h_vmem=16G
#$ -pe smp 4
```

### Cloud (AWS Batch, Google Cloud, etc.)
The R script is self-contained and can be containerized with Docker:

```dockerfile
FROM rocker/r-ver:4.3.0

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev

RUN R -e "install.packages(c('data.table', 'ggplot2', 'ggsci', 'bayesplot', 'cmdstanr'))"

COPY src/run_interim_analysis.R /app/
WORKDIR /app

ENTRYPOINT ["Rscript", "run_interim_analysis.R"]
```
