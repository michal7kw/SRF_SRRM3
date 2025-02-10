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

# Source the original analysis script for shared functions
source("Multi_Cancer_Survival_Analysis.R")

# New function to perform PSI analysis on high SRRM3 expression samples
perform_high_srrm3_psi_analysis <- function(cancer_type, expression_threshold = 0.75) {
  message(sprintf("\nAnalyzing %s for high SRRM3 expression samples", cancer_type))
  
  # Get clinical data
  clinical_data <- get_clinical_data(cancer_type)
  
  # Get SRRM3 expression data
  expression_data <- get_expression_data(cancer_type, "SRRM3")
  
  # Identify high SRRM3 expression samples
  expression_cutoff <- quantile(expression_data$expression, expression_threshold, na.rm = TRUE)
  high_expression_samples <- expression_data %>%
    filter(expression >= expression_cutoff) %>%
    pull(case_id)
  
  message(sprintf("Found %d samples with high SRRM3 expression (>%.2f percentile)",
                 length(high_expression_samples), expression_threshold * 100))
  
  # Filter clinical data for high expression samples
  filtered_clinical <- clinical_data %>% filter(case_id %in% high_expression_samples)
  
  # Retrieve PSI data only for high expression samples
  psi_data <- get_psi_data(cancer_type, "SRRM3", sample_ids = high_expression_samples)
  
  # Merge filtered clinical data and PSI data
  merged_data <- inner_join(filtered_clinical, psi_data, by = "case_id")
  
  if (nrow(merged_data) < 10) {
    stop("Insufficient samples for analysis after filtering")
  }
  
  # Split into high/low PSI groups using median
  median_psi <- median(merged_data$psi, na.rm = TRUE)
  merged_data$group <- if_else(merged_data$psi > median_psi, "High PSI", "Low PSI")
  
  # Fit survival model
  fit <- survfit(Surv(time = overall_survival, 
                     event = deceased) ~ group, 
                data = merged_data)
  
  # Create survival plot
  plot <- ggsurvplot(
    fit,
    data = merged_data,
    pval = TRUE,
    risk.table = TRUE,
    tables.height = 0.3,
    title = sprintf("%s survival by PSI in high SRRM3 samples", cancer_type)
  )
  
  # Create results directory
  results_dir <- file.path("results", "high_srrm3_analysis")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save plots
  plot_base <- file.path(results_dir, sprintf("%s_high_srrm3_psi_survival", cancer_type))
  ggsave(paste0(plot_base, ".pdf"), plot = plot$plot, width = 10, height = 8)
  ggsave(paste0(plot_base, ".png"), plot = plot$plot, width = 10, height = 8, dpi = 300)
  
  # Prepare results
  results <- list(
    fit = fit,
    data = merged_data,
    summary = list(
      cancer_type = cancer_type,
      total_samples = nrow(merged_data),
      expression_cutoff = expression_cutoff,
      median_psi = median_psi,
      group_sizes = table(merged_data$group)
    )
  )
  
  # Save results
  saveRDS(results, file.path(results_dir, sprintf("%s_results.rds", cancer_type)))
  
  # Save summary
  summary_df <- data.frame(
    cancer_type = cancer_type,
    total_samples = nrow(merged_data),
    expression_cutoff = expression_cutoff,
    median_psi = median_psi,
    high_psi_samples = sum(merged_data$group == "High PSI"),
    low_psi_samples = sum(merged_data$group == "Low PSI")
  )
  write.csv(summary_df, 
            file.path(results_dir, sprintf("%s_summary.csv", cancer_type)), 
            row.names = FALSE)
  
  return(results)
} 