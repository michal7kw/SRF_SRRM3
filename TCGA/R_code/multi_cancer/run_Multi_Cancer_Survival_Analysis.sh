#!/bin/bash

###############################################################################
# Multi-Cancer Survival Analysis SLURM Job Script
# 
# This script performs survival analysis across multiple cancer types in parallel
# using SLURM array jobs. For each cancer type, it runs R-based survival analysis
# comparing high vs low SRRM3 expression/PSI values.
#
# Input Files:
# - Multi_Cancer_Survival_Analysis.R: Main R analysis script
# - Clinical data files for each cancer type in TCGA format
# - Expression/PSI data files for each cancer type
#
# Output Files:
# - Individual results: results/<CANCER_TYPE>_results.rds 
# - Combined analysis: results/multi_cancer_analysis/*
# - Log files: logs/Multi_Cancer_Survival_Analysis_*.{err,out}
# - Resource monitoring: logs/resources_*.log
###############################################################################

#SBATCH --job-name=Multi_Cancer_Survival
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/Multi_Cancer_Survival_Analysis_%a.err"
#SBATCH --output="logs/Multi_Cancer_Survival_Analysis_%a.out"
#SBATCH --array=0-5  # One job per cancer type (5 cancer types)

# Load required modules and activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Set working directory
WORK_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3/TCGA/R_code/multi_cancer"
cd $WORK_DIR

# Create necessary directories for outputs and temporary files
mkdir -p logs
mkdir -p results
mkdir -p cache

# Set R environment variables for optimal performance
export R_MAX_NUM_DLLS=150          # Increase DLL limit
export R_GC_MEM_GROW=3             # Memory growth factor
export R_ENABLE_JIT=3              # Enable JIT compilation
export OMP_NUM_THREADS=$SLURM_NTASKS       # Set thread count
export OPENBLAS_NUM_THREADS=$SLURM_NTASKS  # OpenBLAS threads
export MKL_NUM_THREADS=$SLURM_NTASKS       # MKL threads

# Function to monitor system resource usage
monitor_resources() {
    while true; do
        echo "$(date): Memory usage: $(free -h)"
        echo "$(date): CPU usage: $(top -bn1 | head -n 3)"
        sleep 300  # Log every 5 minutes
    done
}

# Start resource monitoring in background
monitor_resources > "logs/resources_${SLURM_ARRAY_TASK_ID}.log" &
MONITOR_PID=$!

# Install required packages if needed
Rscript setup_packages.R

# Run the analysis with working directory
echo "Starting analysis for array task ${SLURM_ARRAY_TASK_ID} at $(date)"
R --vanilla --max-ppsize=500000 -e "
work_dir <- '$WORK_DIR'
setwd(work_dir)
source('Multi_Cancer_Survival_Analysis.R')
"
echo "Completed analysis for array task ${SLURM_ARRAY_TASK_ID} at $(date)"

# Kill the resource monitoring process
kill $MONITOR_PID

# If this is the last array job, combine all results
if [ $SLURM_ARRAY_TASK_ID -eq 4 ]; then
    echo "Combining results from all analyses..."
    Rscript combine_results.R
fi