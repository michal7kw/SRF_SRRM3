# Load required libraries
library(tidyverse)
library(survival)
library(survminer)

# Get working directory
work_dir <- getwd()
RESULTS_DIR <- file.path(work_dir, "results")

# Source the analysis script
source("Multi_Cancer_Survival_Analysis.R")

# List all individual result files using absolute path
result_files <- list.files(RESULTS_DIR, pattern = "*_results.rds", full.names = TRUE)

if (length(result_files) > 0) {
    # Load and combine all results
    all_results <- list()
    for (file in result_files) {
        cancer_type <- gsub("_results.rds", "", basename(file))
        all_results[[cancer_type]] <- readRDS(file)
    }
    
    # Generate and save comparative analyses
    comparative_results <- generate_comparative_analyses(all_results)
    save_combined_results(all_results, comparative_results, RESULTS_DIR)
    
    message("Successfully generated comparative analyses")
} else {
    message("No result files found to combine")
} 