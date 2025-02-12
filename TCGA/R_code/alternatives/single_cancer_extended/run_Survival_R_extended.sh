#!/bin/bash
#SBATCH --job-name=Survival_R_Extended
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/Survival_R_extended.err"
#SBATCH --output="logs/Survival_R_extended.out"

# Load required modules and activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3/TCGA/

# Create necessary directories
mkdir -p logs
mkdir -p results_ext

# Function to monitor resource usage
monitor_resources() {
    while true; do
        echo "$(date): Memory usage: $(free -h)"
        echo "$(date): CPU usage: $(top -bn1 | head -n 3)"
        sleep 300  # Log every 5 minutes
    done
}

# Start resource monitoring in background
monitor_resources > "logs/resources_${SLURM_JOB_ID}.log" &
MONITOR_PID=$!

# Run the R script
echo "Starting Survival R Extended Analysis at $(date)"
Rscript Survival_R_extended.R

# Kill the resource monitor
kill $MONITOR_PID

echo "Analysis completed at $(date)"

# Print summary of results
echo "Results summary:"
ls -l results_ext/