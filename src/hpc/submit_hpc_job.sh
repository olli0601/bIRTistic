#!/bin/bash
#PBS -N analysis_job
#PBS -o logs/job_${PBS_ARRAYID}.out
#PBS -e logs/job_${PBS_ARRAYID}.err
#PBS -t 1-11
#PBS -l walltime=24:00:00
#PBS -l mem=16gb
#PBS -l ncpus=4

# Example PBS submission script for HPC
# Adjust parameters based on your HPC system requirements
# Submit with: qsub submit_hpc_job.sh
#
# Requires a CSV file (job_config.csv) with columns: JOB_ID,DCATI_FILE
# Example:
#   JOB_ID,DCATI_FILE
#   1,dcati_dataset1.rds
#   2,dcati_dataset2.rds
#   3,dcati_dataset3.rds

# Load R module (adjust version as needed)
module load R/4.3.0

# Set up paths
PROJECT_DIR=/path/to/bIRTistic
DATA_DIR=${PROJECT_DIR}/data_for_hpc
STAN_DIR=${PROJECT_DIR}/src/stan
OUTPUT_DIR=${PROJECT_DIR}/output/hpc_results
SCRIPT=${PROJECT_DIR}/R/credit_model_run_analysis.Rscript
CONFIG_CSV=${PROJECT_DIR}/job_config.csv

# Create output and log directories
mkdir -p ${OUTPUT_DIR}
mkdir -p logs

# Get job configuration from CSV file
# Skip header line, then get the row corresponding to PBS_ARRAYID
JOB_ID=$(awk -F',' -v row=${PBS_ARRAYID} 'NR==row+1 {print $1}' ${CONFIG_CSV})
DCATI_FILE=$(awk -F',' -v row=${PBS_ARRAYID} 'NR==row+1 {print $2}' ${CONFIG_CSV})

echo "========================================="
echo "Running analysis job ${JOB_ID}"
echo "Data file: ${DCATI_FILE}"
echo "========================================="

# Run the analysis
Rscript ${SCRIPT} \
  --stan_file ${STAN_DIR}/credit_model_2cats_v251120.stan \
  --stan_include_dir ${STAN_DIR} \
  --dit_file ${DATA_DIR}/dit.rds \
  --dcati_file ${DATA_DIR}/${DCATI_FILE} \
  --job_id ${JOB_ID} \
  --output_dir ${OUTPUT_DIR} \
  --output_prefix analysis \
  --chains 4 \
  --parallel_chains 4 \
  --iter_warmup 500 \
  --iter_sampling 1500 \
  --seed 123

echo "========================================="
echo "Analysis job ${JOB_ID} complete"
echo "========================================="
