# Source the main analysis script
source("Multi_Cancer_Survival_Analysis.R")

# Create results directory
if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

# Run analysis for single cancer type
results <- perform_multi_cancer_analysis(
  cancer_types = c("LGG"),
  analysis_type = "PSI",
  gene = "SRRM3",
  grouping_method = "quartile"
)

# Save results for this cancer type
saveRDS(results, file.path("results", paste0("LGG_results.rds")))
