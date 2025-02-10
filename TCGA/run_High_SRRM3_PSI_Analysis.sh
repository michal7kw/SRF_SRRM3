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
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/High_SRRM3_PSI_Analysis_%a.err"
#SBATCH --output="logs/High_SRRM3_PSI_Analysis_%a.out"
#SBATCH --array=0-5

# Load required modules and activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3/TCGA/

# Create necessary directories
mkdir -p logs
mkdir -p results/high_srrm3_analysis
mkdir -p cache

# Set R environment variables
export R_MAX_NUM_DLLS=150
export R_GC_MEM_GROW=3
export R_ENABLE_JIT=3
export OMP_NUM_THREADS=$SLURM_NTASKS
export OPENBLAS_NUM_THREADS=$SLURM_NTASKS
export MKL_NUM_THREADS=$SLURM_NTASKS

# Function to monitor system resource usage
monitor_resources() {
    while true; do
        echo "$(date): Memory usage: $(free -h)"
        echo "$(date): CPU usage: $(top -bn1 | head -n 3)"
        sleep 300
    done
}

# Start resource monitoring
monitor_resources > "logs/resources_high_srrm3_${SLURM_ARRAY_TASK_ID}.log" &
MONITOR_PID=$!

# Define cancer types
declare -a CANCER_TYPES=("BRCA" "ACC" "UVM" "SKCM" "LGG" "GBM")

# Get current cancer type
CANCER_TYPE=${CANCER_TYPES[$SLURM_ARRAY_TASK_ID]}

echo "Starting high SRRM3 PSI analysis for ${CANCER_TYPE} at $(date)"

# Create temporary R script
cat << 'EOF' > "temp_high_srrm3_analysis_${SLURM_ARRAY_TASK_ID}.R"
# Set options
options(run.main=FALSE)
options(verbose = TRUE)
options(error = function() traceback(2))
options(future.globals.maxSize = 8000 * 1024^2)
options(mc.cores = parallel::detectCores() - 1)

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(survminer)
  library(recount3)
  library(biomaRt)
  library(parallel)
  library(BiocParallel)
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(tidyverse)
  library(DESeq2)
  library(httr)
  library(retry)
  library(futile.logger)
  library(GenomicFeatures)
  library(rtracklayer)
  library(matrixStats)
  library(sparseMatrixStats)
  library(viridis)
  library(data.table)
})

# Source analysis scripts
source("Multi_Cancer_Survival_Analysis.R")  # Source this first as it contains base functions
source("High_SRRM3_PSI_Analysis.R")        # Source this second as it depends on the first

# Run analysis
tryCatch({
    cancer_type <- Sys.getenv("CANCER_TYPE")
    message(sprintf("\nAnalyzing %s", cancer_type))
    
    results <- perform_high_srrm3_psi_analysis(
        cancer_type = cancer_type,
        expression_threshold = 0.75
    )
    
    if (!is.null(results)) {
        message("Analysis completed successfully")
    } else {
        message("No results generated")
        quit(status = 1)
    }
}, error = function(e) {
    message(sprintf("Error in analysis: %s", e$message))
    quit(status = 1)
})
EOF

# Export cancer type for R script
export CANCER_TYPE

# Run R script
R --vanilla --max-ppsize=500000 < "temp_high_srrm3_analysis_${SLURM_ARRAY_TASK_ID}.R"
R_EXIT_CODE=$?

# Check if R script succeeded
if [ $R_EXIT_CODE -ne 0 ]; then
    echo "Error: R script failed for ${CANCER_TYPE}"
    exit 1
fi

echo "Completed analysis for ${CANCER_TYPE} at $(date)"

# Clean up
rm "temp_high_srrm3_analysis_${SLURM_ARRAY_TASK_ID}.R"
kill $MONITOR_PID

# Combine results if this is the last job
if [ $SLURM_ARRAY_TASK_ID -eq 5 ]; then
    echo "Combining results from all cancer types..."
    
    cat << 'EOF' > combine_results.R
    library(tidyverse)
    
    results_dir <- "results/high_srrm3_analysis"
    
    # Combine summaries
    summaries <- list.files(results_dir, 
                           pattern = "*_summary.csv",
                           full.names = TRUE) %>%
        lapply(read.csv) %>%
        bind_rows()
    
    # Add additional statistics
    summaries <- summaries %>%
        mutate(
            psi_analysis_rate = analyzed_samples / high_expr_samples,
            survival_ratio = median_survival_high / median_survival_low
        )
    
    # Save combined results
    write.csv(summaries,
              file.path(results_dir, "combined_summary.csv"),
              row.names = FALSE)
    
    # Create summary visualization
    p <- ggplot(summaries, aes(x = cancer_type)) +
        geom_bar(aes(y = high_expr_samples), stat = "identity") +
        geom_bar(aes(y = analyzed_samples), stat = "identity", alpha = 0.5) +
        theme_minimal() +
        labs(title = "Sample counts by cancer type",
             y = "Number of samples",
             x = "Cancer type") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(results_dir, "sample_counts.pdf"), p, width = 10, height = 6)
EOF

    R --vanilla < combine_results.R
    rm combine_results.R
fi 