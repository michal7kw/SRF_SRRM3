#!/bin/bash

###############################################################################
# High SRRM3 PSI Analysis SLURM Job Script
#
# This script performs PSI (Percent Spliced In) analysis on samples with high
# SRRM3 expression across multiple cancer types using SLURM array jobs. For each
# cancer type, it identifies samples in the top 75th percentile of SRRM3 expression
# and analyzes their splicing patterns.
#
# Input Files:
# - High_SRRM3_PSI_Analysis.R: Main R analysis script
# - Clinical data files for each cancer type in TCGA format
# - Expression data files for each cancer type
# - PSI data files for each cancer type
#
# Output Files:
# - Individual results: results/high_srrm3_analysis/<CANCER_TYPE>_summary.csv
# - Combined results: results/high_srrm3_analysis/combined_summary.csv
# - Log files: logs/High_SRRM3_PSI_Analysis_*.{err,out}
###############################################################################

# SLURM job configuration
#SBATCH --job-name=High_SRRM3_PSI
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB                # Memory allocation
#SBATCH --time=24:00:00           # Maximum runtime
#SBATCH --nodes=1                 # Use one node
#SBATCH --ntasks=32               # Number of CPU cores
#SBATCH --mail-type=ALL           # Email notifications
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/High_SRRM3_PSI_Analysis_%a.err"  # Error log path
#SBATCH --output="logs/High_SRRM3_PSI_Analysis_%a.out" # Output log path
#SBATCH --array=0-5               # Array job for 6 cancer types

# Load required modules and activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake  # Activate the snakemake conda environment

# Set working directory to project location
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3/TCGA/

# Create necessary directories for outputs and temporary files
mkdir -p logs    # Directory for log files
mkdir -p results/high_srrm3_analysis # Directory for analysis results
mkdir -p cache   # Directory for cached data

# Configure R environment for optimal performance
export R_MAX_NUM_DLLS=150          # Increase DLL limit
export R_GC_MEM_GROW=3             # Memory growth factor
export R_ENABLE_JIT=3              # Enable JIT compilation
export OMP_NUM_THREADS=$SLURM_NTASKS       # Set thread count
export OPENBLAS_NUM_THREADS=$SLURM_NTASKS  # OpenBLAS threads
export MKL_NUM_THREADS=$SLURM_NTASKS       # MKL threads

# Define array of cancer types to analyze
declare -a CANCER_TYPES=("BRCA" "ACC" "UVM" "SKCM" "LGG" "GBM")

# Get the cancer type for this specific array job
CANCER_TYPE=${CANCER_TYPES[$SLURM_ARRAY_TASK_ID]}

echo "Starting high SRRM3 PSI analysis for ${CANCER_TYPE} at $(date)"

# Create temporary R script for this specific array job
cat << EOF > temp_high_srrm3_analysis_${SLURM_ARRAY_TASK_ID}.R
# Set R options for better error handling and performance
options(run.main=FALSE)
options(verbose = TRUE)
options(error = function() traceback(2))
options(future.globals.maxSize = 8000 * 1024^2)  # Increase memory limit
options(mc.cores = parallel::detectCores() - 1)  # Use available cores

# Source the main analysis script containing the analysis functions
source("High_SRRM3_PSI_Analysis.R")

# Run analysis with error handling
tryCatch({
  # Perform PSI analysis on samples with high SRRM3 expression
  results <- perform_high_srrm3_psi_analysis(
    cancer_type = "$CANCER_TYPE",  # Current cancer type
    expression_threshold = 0.75    # Analyze top 25% of SRRM3 expression
  )
  
  # Check if results were generated
  if (!is.null(results)) {
    message("Analysis completed successfully")
  } else {
    message("No results generated")
  }
}, error = function(e) {
  # Handle any errors that occur during analysis
  message(sprintf("Error in analysis: %s", e$message))
  quit(status = 1)  # Exit with error status
})
EOF

# Execute the R analysis script with increased memory limit
R --vanilla --max-ppsize=500000 < temp_high_srrm3_analysis_${SLURM_ARRAY_TASK_ID}.R

echo "Completed analysis for ${CANCER_TYPE} at $(date)"

# Clean up temporary R script
rm temp_high_srrm3_analysis_${SLURM_ARRAY_TASK_ID}.R

# If this is the last array job (index 5), combine all results
if [ $SLURM_ARRAY_TASK_ID -eq 5 ]; then
    echo "Combining results from all cancer types..."
    
    # Create R script to combine individual results
    cat << EOF > combine_results.R
    library(tidyverse)
    
    # Define results directory
    results_dir <- "results/high_srrm3_analysis"
    
    # Combine all summary files
    summaries <- list.files(results_dir, 
                           pattern = "*_summary.csv", 
                           full.names = TRUE) %>%
      lapply(read.csv) %>%  # Read each summary file
      bind_rows()           # Combine into single data frame
    
    # Save combined summary
    write.csv(summaries, 
              file.path(results_dir, "combined_summary.csv"), 
              row.names = FALSE)
    
    message("Successfully combined results")
EOF

    # Run the combination script
    R --vanilla < combine_results.R
    
    # Clean up combination script
    rm combine_results.R
fi 