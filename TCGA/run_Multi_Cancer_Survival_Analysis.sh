#!/bin/bash
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
declare -a CANCER_TYPES=("BRCA" "ACC" "UVM" "SKCM" "LGG" "GBM")

# Get the cancer type for this array job
CANCER_TYPE=${CANCER_TYPES[$SLURM_ARRAY_TASK_ID]}

echo "Starting analysis for ${CANCER_TYPE} at $(date)"

# Create a temporary R script for this specific array job
TEMP_SCRIPT="temp_analysis_${SLURM_ARRAY_TASK_ID}.R"

# Create the R script for this specific cancer type
cat << EOF > temp_analysis_${SLURM_ARRAY_TASK_ID}.R
# Load required libraries and set options
options(run.main=FALSE)
options(verbose = TRUE)
options(error = function() traceback(2))
options(future.globals.maxSize = 8000 * 1024^2)
options(mc.cores = parallel::detectCores() - 1)

# Source the analysis script
source("Multi_Cancer_Survival_Analysis.R")

# Create directories
if (!dir.exists("results")) dir.create("results", recursive = TRUE)
if (!dir.exists("cache")) dir.create("cache", recursive = TRUE)

# Print analysis info
cat(sprintf("\nStarting analysis for %s\n", "$CANCER_TYPE"))

# Run analysis with error handling
tryCatch({
  results <- perform_multi_cancer_analysis(
    cancer_types = c("$CANCER_TYPE"),
    analysis_type = "PSI",
    gene = "SRRM3",
    grouping_method = "quartile"
  )
  
  if (!is.null(results)) {
    saveRDS(results, file.path("results", paste0("${CANCER_TYPE}_results.rds")))
    cat(sprintf("\nSuccessfully saved results for %s\n", "$CANCER_TYPE"))
  } else {
    cat(sprintf("\nNo results generated for %s\n", "$CANCER_TYPE"))
  }
}, error = function(e) {
  cat(sprintf("\nError in analysis for %s: %s\n", "$CANCER_TYPE", e$message))
  quit(status = 1)
})

cat(sprintf("\nCompleted analysis for %s\n", "$CANCER_TYPE"))
EOF

# Add this before running the R script
cat << EOF > setup.R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("sparseMatrixStats")
EOF

R --vanilla < setup.R

# Run the R script with increased memory limit
R --vanilla --max-ppsize=500000 < temp_analysis_${SLURM_ARRAY_TASK_ID}.R

echo "Completed analysis for ${CANCER_TYPE} at $(date)"

# Clean up temporary script
rm temp_analysis_${SLURM_ARRAY_TASK_ID}.R

# Kill the resource monitor
kill $MONITOR_PID

# If this is the last array job, combine all results
if [ $SLURM_ARRAY_TASK_ID -eq 4 ]; then
    cat << EOF > combine_results.R
    # Load required libraries
    library(tidyverse)
    library(survival)
    library(survminer)
    
    # Source the analysis script
    source("Multi_Cancer_Survival_Analysis.R")
    
    # Create output directory
    dir.create("results/multi_cancer_analysis", recursive = TRUE, showWarnings = FALSE)
    
    # List all individual result files
    result_files <- list.files("results", pattern = "*_results.rds", full.names = TRUE)
    
    if (length(result_files) > 0) {
        # Load and combine all results
        all_results <- list()
        for (file in result_files) {
            cancer_type <- gsub("_results.rds", "", basename(file))
            all_results[[cancer_type]] <- readRDS(file)
        }
        
        # Generate and save comparative analyses
        comparative_results <- generate_comparative_analyses(all_results)
        save_combined_results(all_results, comparative_results, "results/multi_cancer_analysis")
        
        message("Successfully generated comparative analyses")
    } else {
        message("No result files found to combine")
    }
EOF

    # Run the combination script
    Rscript combine_results.R
    rm combine_results.R
fi