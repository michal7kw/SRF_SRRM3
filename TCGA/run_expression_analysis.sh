#!/bin/bash

###############################################################################
# Expression Survival Analysis SLURM Job Script
#
# This script performs survival analysis based on gene expression data across
# multiple cancer types in parallel using SLURM array jobs. For each cancer type,
# it runs R-based survival analysis comparing high vs low expression values.
#
# Input Files:
# - Expression_Survival_Analysis.R: Main R analysis script
# - Clinical data files for each cancer type in TCGA format
# - Expression data files for each cancer type
#
# Output Files:
# - Individual results: results/<CANCER_TYPE>_expression_results.rds
# - Log files: logs/Expression_Survival_Analysis_*.{err,out}
# - Resource monitoring: logs/resources_*.log
###############################################################################

# SLURM job configuration
#SBATCH --job-name=Expression_Survival
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB                # Memory allocation
#SBATCH --time=24:00:00           # Maximum runtime
#SBATCH --nodes=1                 # Use one node
#SBATCH --ntasks=32               # Number of CPU cores
#SBATCH --mail-type=ALL           # Email notifications
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/Expression_Survival_Analysis_%a.err"  # Error log path
#SBATCH --output="logs/Expression_Survival_Analysis_%a.out" # Output log path
#SBATCH --array=0-5               # Array job for 6 cancer types

# Load required modules and activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake  # Activate the snakemake conda environment

# Set working directory to project location
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3/TCGA/

# Create necessary directories for outputs and temporary files
mkdir -p logs    # Directory for log files
mkdir -p results # Directory for analysis results
mkdir -p cache   # Directory for cached data

# Configure R environment for optimal performance
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
MONITOR_PID=$!  # Store process ID for later termination

# Define array of cancer types to analyze
declare -a CANCER_TYPES=("ACC" "UVM" "SKCM" "LLG" "GBM")

# Get the cancer type for this specific array job
CANCER_TYPE=${CANCER_TYPES[$SLURM_ARRAY_TASK_ID]}

# Log analysis start information
echo "Starting analysis for ${CANCER_TYPE} at $(date)"
echo "Running on node: $(hostname)"
echo "Number of cores: $SLURM_NTASKS"
echo "Memory allocated: 64GB"

# Execute the main R analysis script
Rscript Expression_Survival_Analysis.R

# Calculate and log execution time
END_TIME=$(date +%s)
START_TIME=$(date +%s -r "logs/Expression_Survival_Analysis_${SLURM_ARRAY_TASK_ID}.out")
DURATION=$((END_TIME - START_TIME))
echo "Analysis completed in ${DURATION} seconds"

# Terminate the resource monitoring process
kill $MONITOR_PID

# Optional: Compress log files to save space
# gzip "logs/resources_${SLURM_ARRAY_TASK_ID}.log"
# gzip "logs/Expression_Survival_Analysis_${SLURM_ARRAY_TASK_ID}.out"
# gzip "logs/Expression_Survival_Analysis_${SLURM_ARRAY_TASK_ID}.err"

# Final job completion message
echo "Job completed at $(date)"
