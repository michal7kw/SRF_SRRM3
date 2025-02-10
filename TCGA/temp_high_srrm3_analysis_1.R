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
