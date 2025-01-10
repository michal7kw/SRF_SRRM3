# Source the optimized code
source("TCGA/Survival_R_Alternative_extended_optimized.R")

# Run a single analysis
single_analysis <- function(cancer_type = "BRCA", 
                          analysis_type = "SRRM3_expression",
                          survival_type = "OS") {
  
  result <- perform_survival_analysis_enhanced(
    cancer_type = cancer_type,
    analysis_type = analysis_type,
    survival_type = survival_type,
    grouping_method = "quartile"
  )
  
  # Save the results
  if (!is.null(result)) {
    save_analysis_results(result, cancer_type)
    return(result)
  }
}

# Run all analyses for a cancer type
run_all_analyses <- function(cancer_type = "BRCA") {
  results <- run_survival_analysis(
    cancer_type = cancer_type,
    debug = TRUE,
    timeout = 600  # 10 minutes timeout
  )
  
  if (!is.null(results)) {
    save_analysis_results(results, cancer_type)
  }
  
  return(results)
}

# Example usage:
# 1. For a single analysis:
result_single <- single_analysis(
  cancer_type = "BRCA",
  analysis_type = "SRRM3_expression",
  survival_type = "OS"
)

# 2. For all analyses:
results_all <- run_all_analyses("BRCA")

# 3. To view a specific plot:
if (!is.null(results_all)) {
  # View SRRM3 expression OS plot
  print(results_all$SRRM3_expression_OS$plot)
  
  # View SRRM3 PSI OS plot
  print(results_all$SRRM3_PSI_OS$plot)
  
  # View SRRM4 expression OS plot
  print(results_all$SRRM4_expression_OS$plot)
} 