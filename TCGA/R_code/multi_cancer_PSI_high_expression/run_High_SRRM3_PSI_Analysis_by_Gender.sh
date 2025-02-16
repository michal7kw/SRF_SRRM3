#!/bin/bash

###############################################################################
# High SRRM3 PSI Analysis SLURM Job Script
#
# This script performs PSI analysis on samples with high SRRM3 expression across
# multiple cancer types using SLURM array jobs.
###############################################################################

#SBATCH --job-name=High_SRRM3_PSI
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/High_SRRM3_PSI_Analysis_%a.err"
#SBATCH --output="logs/High_SRRM3_PSI_Analysis_%a.out"
#SBATCH --array=0-5

# Load required modules and activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3/TCGA/multi_cancer_PSI_high_expression

# Create necessary directories
mkdir -p logs
mkdir -p results
mkdir -p cache

# Set R environment variables
export R_MAX_NUM_DLLS=150
export R_GC_MEM_GROW=3
export R_ENABLE_JIT=3
export OMP_NUM_THREADS=$SLURM_NTASKS
export OPENBLAS_NUM_THREADS=$SLURM_NTASKS
export MKL_NUM_THREADS=$SLURM_NTASKS

# Define cancer types
declare -a CANCER_TYPES=("ACC" "UVM" "SKCM" "LGG" "GBM")

# Get current cancer type
CANCER_TYPE=${CANCER_TYPES[$SLURM_ARRAY_TASK_ID]}

echo "Starting high SRRM3 PSI analysis for ${CANCER_TYPE} at $(date)"

# Export cancer type for R script
export CANCER_TYPE

# Run R script directly
R --vanilla --max-ppsize=500000 -f "./High_SRRM3_PSI_Analysis_by_Gender.R"
R_EXIT_CODE=$?

# Check if R script succeeded
if [ $R_EXIT_CODE -ne 0 ]; then
    echo "Error: R script failed for ${CANCER_TYPE}"
    exit 1
fi

echo "Completed analysis for ${CANCER_TYPE} at $(date)"