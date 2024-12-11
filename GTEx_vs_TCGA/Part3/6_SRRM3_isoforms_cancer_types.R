# Analysis of SRRM3 isoform expression across cancer types
# Load required libraries
library(recount3)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(SummarizedExperiment)
library(progress)
library(futile.logger)
library(R.utils)
library(viridis)

# Set up logging
flog.appender(appender.file("srrm3_cancer_analysis.log"))
flog.threshold(DEBUG)

# Source the functions from the previous script
source("5_PSI_values_opt_checkpoints_examination.R")

# Function to get cancer type from TCGA project name
get_cancer_type <- function(project) {
  # Extract cancer type from project name (e.g., "TCGA-BRCA" -> "BRCA")
  gsub("TCGA-", "", project)
}

# Function to analyze SRRM3 isoforms across cancer types
analyze_cancer_types <- function() {
  # Get available TCGA projects
  projects <- available_projects()
  tcga_projects <- projects[grep("^TCGA-", projects$project),]
  
  # Store results
  results_list <- list()
  
  # Progress bar
  pb <- progress_bar$new(
    format = "Processing :project [:bar] :percent eta: :eta",
    total = nrow(tcga_projects)
  )
  
  for(i in 1:nrow(tcga_projects)) {
    project_info <- tcga_projects[i,]
    project_name <- project_info$project
    cancer_type <- get_cancer_type(project_name)
    
    pb$tick(tokens = list(project = project_name))
    
    # Create RSE object
    rse <- create_rse_safe(project_info)
    if(is.null(rse)) {
      flog.warn("Skipping %s - failed to create RSE", project_name)
      next
    }
    
    # Calculate PSI values
    psi_values <- calculate_exon15_psi(rse)
    if(is.null(psi_values)) {
      flog.warn("Skipping %s - failed to calculate PSI values", project_name)
      next
    }
    
    # Store results with cancer type
    results_list[[project_name]] <- data.frame(
      cancer_type = cancer_type,
      sample_id = names(psi_values),
      psi = psi_values,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine all results
  all_results <- do.call(rbind, results_list)
  return(all_results)
}

# Function to plot cancer type distributions
plot_cancer_distributions <- function(results) {
  # Calculate summary statistics
  summary_stats <- results %>%
    group_by(cancer_type) %>%
    summarise(
      median_psi = median(psi, na.rm = TRUE),
      mean_psi = mean(psi, na.rm = TRUE),
      sd_psi = sd(psi, na.rm = TRUE),
      n_samples = n(),
      .groups = "drop"
    ) %>%
    filter(n_samples >= 10)  # Filter for cancer types with at least 10 samples
  
  # Create violin plot
  p1 <- ggplot(results %>% filter(cancer_type %in% summary_stats$cancer_type), 
               aes(x = reorder(cancer_type, psi, FUN = median), y = psi)) +
    geom_violin(aes(fill = cancer_type), alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    scale_fill_viridis(discrete = TRUE) +
    labs(
      title = "SRRM3 Exon 15 PSI Distribution Across Cancer Types",
      x = "Cancer Type",
      y = "PSI Value"
    )
  
  # Create summary statistics plot
  p2 <- ggplot(summary_stats, 
               aes(x = reorder(cancer_type, median_psi), y = median_psi)) +
    geom_point(aes(size = n_samples), alpha = 0.7) +
    geom_errorbar(aes(ymin = median_psi - sd_psi, 
                     ymax = median_psi + sd_psi), 
                 width = 0.2) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = "Median PSI Values by Cancer Type",
      x = "Cancer Type",
      y = "Median PSI Value",
      size = "Number of Samples"
    )
  
  return(list(distribution_plot = p1, summary_plot = p2))
}

# Main execution
main <- function() {
  flog.info("Starting SRRM3 cancer type analysis")
  
  # Analyze cancer types
  results <- analyze_cancer_types()
  
  # Save results
  write.csv(results, "srrm3_cancer_type_results.csv", row.names = FALSE)
  flog.info("Results saved to srrm3_cancer_type_results.csv")
  
  # Create plots
  plots <- plot_cancer_distributions(results)
  
  # Save plots
  ggsave("srrm3_distribution_plot.pdf", plots$distribution_plot, 
         width = 12, height = 8)
  ggsave("srrm3_summary_plot.pdf", plots$summary_plot, 
         width = 12, height = 8)
  
  flog.info("Analysis complete. Plots saved.")
  
  return(results)
}

# Run the analysis if this script is being run directly
if (!interactive()) {
  main()
}
