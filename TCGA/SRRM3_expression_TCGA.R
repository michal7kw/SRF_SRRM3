### Output #################################
# "SRRM3_expression_TCGA_short_names.pdf"
# "SRRM3_expression_TCGA_full_names.pdf"
# "SRRM3_expression_TCGA_short_names.csv"
# "SRRM3_expression_TCGA_full_names.csv"
###########################################

library(recount3)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(SummarizedExperiment)
library(futile.logger)
library(viridis)
library(scales)

# Set up logging
flog.appender(appender.file("./logs/SRRM3_expression_TCGA.log"))
flog.threshold(DEBUG)

# Function to get expression data
get_expression_data <- function(project_info) {
  project_name <- project_info$project
  cache_file <- file.path(CACHE_DIR, paste0("srrm3_expr_", project_name, ".rds"))
  
  # Check cache first
  cached_data <- get_cached_data(cache_file)
  if (!is.null(cached_data)) {
    return(cached_data)
  }
  
  flog.info("\nProcessing project: %s", project_name)
  
  tryCatch({
    # Create RSE object
    rse_gene <- create_rse(subset(available_projects(), 
                                 project == project_info$project & 
                                   file_source == "tcga"), 
                          type = "gene")
    
    # Get SRRM3 expression
    srrm3_idx <- which(rowData(rse_gene)$gene_name == "SRRM3")
    
    if(length(srrm3_idx) == 0) {
      flog.warn("SRRM3 not found in %s", project_name)
      return(NULL)
    }
    
    # Extract expression data and clinical info
    data <- data.frame(
      sample_id = colnames(rse_gene),
      srrm3_expr = assay(rse_gene)[srrm3_idx, ],
      stringsAsFactors = FALSE
    )
    
    # Add cancer type information
    data$cancer_type_short <- sub("^TCGA-", "", project_info$project)
    data$cancer_type_full <- NA
    
    # Get histological diagnosis from clinical data
    clinical <- as.data.frame(colData(rse_gene))
    if ("tcga.cgc_case_histological_diagnosis" %in% colnames(clinical)) {
      hist_diag <- clinical$tcga.cgc_case_histological_diagnosis
      if (!all(is.na(hist_diag))) {
        data$cancer_type_full <- hist_diag
      }
    }
    
    # If no histological diagnosis available, use short name
    data$cancer_type_full[is.na(data$cancer_type_full)] <- 
      data$cancer_type_short[is.na(data$cancer_type_full)]
    
    # Cache the results
    save_cached_data(data, cache_file)
    
    return(data)
  }, error = function(e) {
    flog.error("Error processing project %s: %s", project_name, e$message)
    return(NULL)
  })
}

# Function to analyze cancer types
analyze_cancer_types <- function() {
  # Get TCGA projects
  projects <- get_tcga_projects()
  flog.info("Found %d TCGA projects", nrow(projects))
  
  # Process each project
  results <- lapply(seq_len(nrow(projects)), function(i) {
    flog.info("Processing project %d/%d", i, nrow(projects))
    get_expression_data(projects[i,])
  })
  
  # Combine all results
  results <- do.call(rbind, results[!sapply(results, is.null)])
  
  if(is.null(results) || nrow(results) == 0) {
    stop("No valid results to analyze")
  }
  
  # Calculate statistics for both naming schemes
  analysis_data <- list(
    short_names = results %>%
      group_by(cancer_type_short) %>%
      summarise(
        mean_expr = mean(srrm3_expr, na.rm = TRUE),
        sd_expr = sd(srrm3_expr, na.rm = TRUE),
        sample_size = n(),
        .groups = "drop"
      ) %>%
      rename(cancer_type = cancer_type_short),
    
    full_names = results %>%
      group_by(cancer_type_full) %>%
      summarise(
        mean_expr = mean(srrm3_expr, na.rm = TRUE),
        sd_expr = sd(srrm3_expr, na.rm = TRUE),
        sample_size = n(),
        .groups = "drop"
      ) %>%
      rename(cancer_type = cancer_type_full)
  )
  
  return(analysis_data)
}

# Function to plot cancer distributions
plot_cancer_distributions <- function(analysis_data, use_full_names = FALSE) {
  # Select appropriate data
  data_to_plot <- if(use_full_names) analysis_data$full_names else analysis_data$short_names
  name_type <- if(use_full_names) "full_names" else "short_names"
  
  if(is.null(data_to_plot) || nrow(data_to_plot) == 0) {
    stop(sprintf("No data available for %s", name_type))
  }
  
  # Calculate standard error
  plot_data <- data_to_plot %>%
    mutate(se = sd_expr / sqrt(sample_size))
  
  # Create plot
  p1 <- ggplot(plot_data, 
               aes(x = reorder(cancer_type, mean_expr), 
                   y = mean_expr)) +
    geom_bar(stat = "identity", fill = "#2166AC", alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_expr - se, 
                      ymax = mean_expr + se),
                  width = 0.25) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = "SRRM3 Expression Across Cancer Types",
      x = "Cancer Type",
      y = "Mean Expression (normalized counts)"
    )
  
  # Save plot
  ggsave(sprintf("./output/SRRM3_expression_TCGA_%s.pdf", name_type), 
         p1, width = 15, height = 10, dpi = 300)
  ggsave(sprintf("./output/SRRM3_expression_TCGA_%s.png", name_type), 
         p1, width = 15, height = 10, dpi = 300)
  
  return(p1)
}

# Main execution
results <- analyze_cancer_types()
if(!is.null(results)) {
  # Create plots with both naming schemes
  plot_short <- plot_cancer_distributions(results, use_full_names = FALSE)
  plot_full <- plot_cancer_distributions(results, use_full_names = TRUE)
  
  # Save results
  write.csv(results$short_names, 
            "./output/SRRM3_expression_TCGA_short_names.csv", 
            row.names = FALSE)
  write.csv(results$full_names, 
            "./output/SRRM3_expression_TCGA_full_names.csv", 
            row.names = FALSE)
}