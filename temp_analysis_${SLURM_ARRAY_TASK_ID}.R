# Load required libraries and set options
options(run.main=FALSE)
options(verbose = TRUE)
options(error = function() traceback(2))
options(future.globals.maxSize = 8000 * 1024^2)
options(mc.cores = parallel::detectCores() - 1)

# Source the analysis script with error handling
tryCatch({
  source("Multi_Cancer_Survival_Analysis.R")
}, error = function(e) {
  cat(sprintf("\nError loading analysis script: %s\n", e$message))
  quit(status = 1)
})

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