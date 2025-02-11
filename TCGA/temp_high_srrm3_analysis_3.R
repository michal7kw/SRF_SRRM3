# Set error handling to get more information
options(error = function() {
  traceback(3)
  if (!interactive()) quit(status = 1)
})

# Set options
options(run.main=FALSE)
options(verbose = TRUE)
options(future.globals.maxSize = 8000 * 1024^2)
options(mc.cores = parallel::detectCores() - 1)

# Load required libraries with error checking
required_packages <- c(
  "dplyr", "survival", "survminer", "recount3", "biomaRt",
  "parallel", "BiocParallel", "TCGAbiolinks", "SummarizedExperiment",
  "tidyverse", "DESeq2", "httr", "retry", "futile.logger",
  "GenomicFeatures", "rtracklayer", "matrixStats", "sparseMatrixStats",
  "viridis", "data.table", "R.utils"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    stop(sprintf("Required package '%s' is not installed", pkg))
  }
}

# Source analysis scripts with error checking
for (script in c("./R_code/multi_cancer/Multi_Cancer_Survival_Analysis.R", "./R_code/multi_cancer_PSI_high_expression/High_SRRM3_PSI_Analysis.R")) {
  message(sprintf("Sourcing %s...", script))
  if (!file.exists(script)) {
    stop(sprintf("Script file '%s' not found", script))
  }
  source(script)
  message(sprintf("%s sourced successfully", script))
}

# Verify function exists
if (!exists("perform_high_srrm3_psi_analysis")) {
  stop("perform_high_srrm3_psi_analysis function not defined after sourcing scripts")
}

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
    message("Traceback:")
    print(sys.calls())
    quit(status = 1)
})
