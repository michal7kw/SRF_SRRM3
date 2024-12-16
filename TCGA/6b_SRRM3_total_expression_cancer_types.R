# Analysis of total SRRM3 expression across cancer types
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
flog.appender(appender.file("srrm3_total_expression_analysis.log"))
flog.threshold(DEBUG)

# Define cache directory
CACHE_DIR <- "cache"
if (!dir.exists(CACHE_DIR)) {
  dir.create(CACHE_DIR)
}

# Function to get cached data
get_cached_data <- function(cache_file) {
  if (file.exists(cache_file)) {
    flog.info("Loading cached data from %s", cache_file)
    return(readRDS(cache_file))
  }
  return(NULL)
}

# Function to save cached data
save_cached_data <- function(data, cache_file) {
  flog.info("Saving data to cache: %s", cache_file)
  saveRDS(data, cache_file)
}

# Define SRRM3 information
SRRM3_INFO <- list(
  gene = list(
    name = "SRRM3",
    chr = "chr7",
    start = 76201896,
    end = 76287287
  )
)

# Function to get cancer type from project name
get_cancer_type <- function(project_name) {
  # Extract the cancer type code (e.g., "BRCA" from "TCGA-BRCA")
  cancer_type <- sub("^TCGA-", "", project_name)
  return(cancer_type)
}

# Function to create RSE object for gene expression
create_rse_gene <- function(project_info) {
  tryCatch({
    # Create RSE object
    rse <- create_rse(
      project_info,
      type = "gene",
      verbose = TRUE
    )
    
    # Filter to SRRM3
    rse_filtered <- rse[rowRanges(rse)$gene_name == "SRRM3", ]
    
    # Print available assay names for debugging
    flog.info("Available assays: %s", paste(assayNames(rse_filtered), collapse = ", "))
    
    # Get the first available assay if 'counts' is not available
    if (length(assayNames(rse_filtered)) > 0) {
      assay_name <- assayNames(rse_filtered)[1]
      flog.info("Using assay: %s", assay_name)
      
      # Get the data
      assay_data <- assay(rse_filtered, assay_name)
      
      rm(rse)
      gc()
      
      return(list(
        rse = rse_filtered,
        data = assay_data
      ))
    } else {
      flog.error("No assays available in RSE object")
      return(NULL)
    }
  }, error = function(e) {
    flog.error("Error creating RSE: %s", e$message)
    return(NULL)
  })
}

# Main function to analyze cancer types
analyze_cancer_types <- function() {
  # Check for cached results
  cache_file <- file.path(CACHE_DIR, "srrm3_expression_results.rds")
  cached_results <- get_cached_data(cache_file)
  
  if (!is.null(cached_results)) {
    flog.info("Using cached results")
    return(cached_results)
  }
  
  # Get available TCGA projects for human
  projects <- available_projects(organism = "human")
  
  # Filter for TCGA projects only
  tcga_projects <- projects[projects$file_source == "tcga", ]
  
  # Print available projects for debugging
  flog.info("Found %d TCGA projects", nrow(tcga_projects))
  
  # Store results
  results_list <- list()
  
  # Progress tracking
  total_projects <- nrow(tcga_projects)
  flog.info("Starting analysis of %d TCGA projects", total_projects)
  
  if(total_projects == 0) {
    flog.error("No TCGA projects found")
    return(NULL)
  }
  
  for(i in 1:nrow(tcga_projects)) {
    project_info <- tcga_projects[i,]
    project_name <- project_info$project
    cancer_type <- get_cancer_type(project_name)
    
    # Check for cached project data
    project_cache_file <- file.path(CACHE_DIR, paste0("srrm3_", project_name, ".rds"))
    project_data <- get_cached_data(project_cache_file)
    
    if (!is.null(project_data)) {
      results_list[[project_name]] <- project_data
      flog.info("Loaded cached data for %s: %d samples", project_name, nrow(project_data))
      next
    }
    
    flog.info("==========================================")
    flog.info("Processing project %d/%d: %s", i, total_projects, project_name)
    
    tryCatch({
      rse_result <- create_rse_gene(project_info)
      if(is.null(rse_result) || ncol(rse_result$rse) == 0 || nrow(rse_result$rse) == 0) {
        flog.warn("Invalid RSE object for %s", project_name)
        next
      }
      
      # Get expression values and normalize to TPM
      expression_values <- as.numeric(rse_result$data)
      gene_length <- width(rowRanges(rse_result$rse))
      tpm <- exp(log(expression_values) - log(gene_length) + log(1e6))
      
      # Store results
      project_results <- data.frame(
        cancer_type = cancer_type,
        sample_id = colnames(rse_result$rse),
        expression = tpm,
        stringsAsFactors = FALSE
      )
      
      # Cache project results
      save_cached_data(project_results, project_cache_file)
      
      results_list[[project_name]] <- project_results
      flog.info("Processed %s: %d samples", project_name, length(tpm))
      
    }, error = function(e) {
      flog.error("Error processing %s: %s", project_name, e$message)
    })
  }
  
  if(length(results_list) == 0) {
    flog.error("No results were generated - all projects failed processing")
    return(NULL)
  }
  
  all_results <- do.call(rbind, results_list)
  
  # Cache complete results
  save_cached_data(all_results, cache_file)
  
  flog.info("Analysis complete: %d samples across %d cancer types", 
            nrow(all_results), length(unique(all_results$cancer_type)))
  
  return(all_results)
}

# Function to plot cancer type distributions
plot_cancer_distributions <- function(results) {
  # Calculate summary statistics for cancer types with at least 10 samples
  summary_stats <- results %>%
    group_by(cancer_type) %>%
    summarise(
      median_expression = median(expression, na.rm = TRUE),
      mean_expression = mean(expression, na.rm = TRUE),
      sd_expression = sd(expression, na.rm = TRUE),
      n_samples = n(),
      .groups = "drop"
    ) %>%
    filter(n_samples >= 10)
  
  # Add sample size to cancer type labels
  results$cancer_type_label <- paste0(
    results$cancer_type, "\n(n=", 
    summary_stats$n_samples[match(results$cancer_type, summary_stats$cancer_type)], 
    ")"
  )
  
  # Create violin plot with improved aesthetics
  p1 <- ggplot(
    results %>% filter(cancer_type %in% summary_stats$cancer_type),
    aes(x = reorder(cancer_type_label, expression, FUN = median),
        y = expression)
  ) +
    geom_violin(aes(fill = cancer_type), alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.2, fill = "white", alpha = 0.7, outlier.size = 0.5) +
    scale_fill_viridis(discrete = TRUE) +
    scale_y_continuous(
      trans = "log10",
      labels = trans_format("log10", math_format(10^.x))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    ) +
    labs(
      title = "SRRM3 Expression Across Cancer Types",
      subtitle = paste("Analysis of", length(unique(results$cancer_type)), "cancer types"),
      x = "Cancer Type",
      y = "Expression (TPM)"
    )
  
  # Save plot in high resolution
  ggsave("SRRM3_expression_cancer_types.pdf", p1, width = 15, height = 10, dpi = 300)
  ggsave("SRRM3_expression_cancer_types.png", p1, width = 15, height = 10, dpi = 300)
  
  return(p1)
}

# Run analysis
results <- analyze_cancer_types()
if(!is.null(results)) {
  plot <- plot_cancer_distributions(results)
  print(plot)
  
  # Save results
  write.csv(results, "SRRM3_expression_results.csv", row.names = FALSE)
}
