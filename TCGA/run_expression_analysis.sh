#!/bin/bash
#SBATCH --job-name=Expression_Survival
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/Expression_Survival_Analysis_%a.err"
#SBATCH --output="logs/Expression_Survival_Analysis_%a.out"
#SBATCH --array=0-5  # One job per cancer type

# Load required modules and activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3/TCGA/

# Create necessary directories
mkdir -p logs
mkdir -p results
mkdir -p cache

# Set R environment variables for better performance
export R_MAX_NUM_DLLS=150
export R_GC_MEM_GROW=3
export R_ENABLE_JIT=3
export OMP_NUM_THREADS=$SLURM_NTASKS
export OPENBLAS_NUM_THREADS=$SLURM_NTASKS
export MKL_NUM_THREADS=$SLURM_NTASKS

# Function to monitor resource usage
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

# Define cancer types array
declare -a CANCER_TYPES=("ACC" "UVM" "SKCM" "LLG" "GBM")

# Get the cancer type for this array job
CANCER_TYPE=${CANCER_TYPES[$SLURM_ARRAY_TASK_ID]}

# Log start time and cancer type
echo "Starting analysis for ${CANCER_TYPE} at $(date)"
echo "Running on node: $(hostname)"
echo "Number of cores: $SLURM_NTASKS"
echo "Memory allocated: 64GB"

# Run the R script
Rscript Expression_Survival_Analysis.R

# Calculate and log execution time
END_TIME=$(date +%s)
START_TIME=$(date +%s -r "logs/Expression_Survival_Analysis_${SLURM_ARRAY_TASK_ID}.out")
DURATION=$((END_TIME - START_TIME))
echo "Analysis completed in ${DURATION} seconds"

# Kill the resource monitoring process
kill $MONITOR_PID

# # Compress log files
# gzip "logs/resources_${SLURM_ARRAY_TASK_ID}.log"
# gzip "logs/Expression_Survival_Analysis_${SLURM_ARRAY_TASK_ID}.out"
# gzip "logs/Expression_Survival_Analysis_${SLURM_ARRAY_TASK_ID}.err"

echo "Job completed at $(date)"
