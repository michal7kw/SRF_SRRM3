#####################################################################
# High SRRM3 PSI Analysis Script
#####################################################################
# This script performs PSI analysis specifically for samples with high
# SRRM3 expression. It reuses functions from Multi_Cancer_Survival_Analysis.R
# and adds specific functionality for high SRRM3 expression analysis.

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

# Source the main analysis script to reuse common functions
source("Multi_Cancer_Survival_Analysis.R")

#####################################################################
# High SRRM3 Sample Selection
#####################################################################
# Function to identify samples with high SRRM3 expression
get_high_srrm3_samples <- function(cancer_type, threshold = 0.75) {
  message(sprintf("Getting high SRRM3 expression samples for %s", cancer_type))
  
  # Get SRRM3 expression data
  expression_data <- get_expression_data(cancer_type, "SRRM3")
  
  if (nrow(expression_data) == 0) {
    stop("No expression data available")
  }
  
  # Calculate expression threshold
  cutoff <- quantile(expression_data$expression, probs = threshold, na.rm = TRUE)
  
  # Select high expression samples
  high_expr_samples <- expression_data %>%
    filter(expression >= cutoff) %>%
    pull(case_id)
  
  message(sprintf("Found %d samples with high SRRM3 expression (>%.2f percentile)",
                 length(high_expr_samples), threshold * 100))
  
  return(high_expr_samples)
}

#####################################################################
# High SRRM3 PSI Analysis
#####################################################################
# Function to perform PSI analysis on high SRRM3 expression samples
perform_high_srrm3_psi_analysis <- function(cancer_type, 
                                          expression_threshold = 0.75,
                                          grouping_method = "quartile") {
  tryCatch({
    # Get high SRRM3 expression samples
    high_expr_samples <- get_high_srrm3_samples(cancer_type, expression_threshold)
    
    if (length(high_expr_samples) < 10) {
      stop("Insufficient samples with high SRRM3 expression")
    }
    
    # Get clinical data
    clinical_data <- get_clinical_data(cancer_type)
    
    # Get PSI data specifically for high expression samples
    psi_data <- get_psi_data(cancer_type, "SRRM3", high_expr_samples)
    
    # Merge clinical and PSI data
    merged_data <- inner_join(
      clinical_data %>% mutate(merge_id = substr(case_id, 1, 12)),
      psi_data %>% mutate(merge_id = substr(case_id, 1, 12)),
      by = "merge_id"
    )
    
    if (nrow(merged_data) == 0) {
      stop("No samples after merging clinical and PSI data")
    }
    
    # Group samples based on PSI values
    if (grouping_method == "quartile") {
      quartiles <- quantile(merged_data$psi, probs = c(0.25, 0.75), na.rm = TRUE)
      merged_data$group <- case_when(
        merged_data$psi <= quartiles[1] ~ "Low",
        merged_data$psi >= quartiles[2] ~ "High",
        TRUE ~ "Medium"
      )
    }
    
    # Remove any NA groups
    merged_data <- merged_data %>% filter(!is.na(group))
    
    # Print group sizes
    message("Sample sizes per group:")
    print(table(merged_data$group))
    
    # Perform survival analysis
    fit <- survfit(
      Surv(time = overall_survival, event = deceased) ~ group,
      data = merged_data
    )
    
    # Create survival plot
    plot <- ggsurvplot(
      fit,
      data = merged_data,
      pval = TRUE,
      risk.table = TRUE,
      tables.height = 0.3,
      title = sprintf("%s survival by PSI (High SRRM3 samples)", cancer_type)
    )
    
    # Create results directory
    results_dir <- file.path("results", "high_srrm3_analysis")
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Save plots
    plot_base <- file.path(results_dir, sprintf("%s_high_srrm3_psi_survival", cancer_type))
    ggsave(paste0(plot_base, ".pdf"), plot = plot$plot, width = 10, height = 8)
    ggsave(paste0(plot_base, ".png"), plot = plot$plot, width = 10, height = 8, dpi = 300)
    
    # Create summary
    summary_df <- data.frame(
      cancer_type = cancer_type,
      total_samples = nrow(merged_data),
      high_expr_samples = length(high_expr_samples),
      analyzed_samples = nrow(merged_data),
      high_psi_group = sum(merged_data$group == "High"),
      low_psi_group = sum(merged_data$group == "Low"),
      medium_psi_group = sum(merged_data$group == "Medium"),
      median_survival_high = summary(fit)$table["median", "group=High"],
      median_survival_low = summary(fit)$table["median", "group=Low"]
    )
    
    # Save summary
    write.csv(
      summary_df,
      file.path(results_dir, sprintf("%s_summary.csv", cancer_type)),
      row.names = FALSE
    )
    
    # Save full results object
    saveRDS(
      list(fit = fit, data = merged_data, summary = summary_df),
      file.path(results_dir, sprintf("%s_results.rds", cancer_type))
    )
    
    return(list(fit = fit, data = merged_data, summary = summary_df))
    
  }, error = function(e) {
    message(sprintf("Error in high SRRM3 PSI analysis for %s: %s", 
                   cancer_type, e$message))
    return(NULL)
  })
} 