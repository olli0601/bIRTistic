#!/bin/bash
# R wrapper script for birtistic environment
module load miniforge/3 tools/prod zlib libxml2/2.11.5-GCCcore-13.2.0
eval "$(~/miniforge3/bin/conda shell.bash hook)"
eval "$(mamba shell hook --shell bash)"
mamba activate birtistic
exec ~/miniforge3/envs/birtistic/bin/R "$@"
