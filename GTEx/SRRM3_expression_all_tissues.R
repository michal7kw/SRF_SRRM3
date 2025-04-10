### Output #################################
# "SRRM3_expression_all_tissues.pdf"
# "SRRM3_expression_all_tissues.png"
# "SRRM3_expression_all_tissues.csv"
###########################################

# Analysis of total SRRM3 expression across GTEx tissues
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
flog.appender(appender.file("./logs/SRRM3_expression_all_tissues.log"))
flog.threshold(DEBUG)

# Define cache directory
CACHE_DIR <- "../cache"
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

# Function to get tissue type from GTEx metadata
get_tissue_type <- function(rse) {
  # Extract metadata
  metadata <- colData(rse) %>% 
    as.data.frame()
  
  # Try different possible metadata fields
  tissue_type <- NULL
  
  # Check different possible metadata fields
  if("gtex.smtsd" %in% colnames(metadata)) {
    tissue_type <- metadata$gtex.smtsd
  } else if("gtex.smts" %in% colnames(metadata)) {
    tissue_type <- metadata$gtex.smts
  } else if("tissue" %in% colnames(metadata)) {
    tissue_type <- metadata$tissue
  }
  
  if(is.null(tissue_type)) {
    flog.warn("Could not find tissue type information in metadata")
    flog.info("Available metadata columns: %s", 
              paste(colnames(metadata), collapse=", "))
    return(rep("Unknown", nrow(metadata)))
  }
  
  # Log unique tissue types found
  unique_tissues <- unique(tissue_type)
  flog.info("Found %d unique tissue types:", length(unique_tissues))
  for(tissue in unique_tissues) {
    flog.info("  - %s", tissue)
  }
  
  return(tissue_type)
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

# Main function to analyze GTEx tissues
analyze_gtex_tissues <- function() {
  # Check for cached results
  cache_file <- file.path(CACHE_DIR, "SRRM3_expression_all_tissues.rds")
  cached_results <- get_cached_data(cache_file)
  
  if (!is.null(cached_results)) {
    flog.info("Using cached results")
    return(cached_results)
  }
  
  # Get all available projects
  projects <- available_projects()
  # Filter for GTEx projects - note they use different project names
  gtex_projects <- subset(projects, file_source == "gtex")
  
  # Print available GTEx projects for debugging
  flog.info("Found %d GTEx projects", nrow(gtex_projects))
  for(i in 1:nrow(gtex_projects)) {
    flog.info("Project %d: %s", i, gtex_projects$project[i])
  }
  
  # Store results for all GTEx projects
  results_list <- list()
  
  # Process each GTEx project
  for(i in 1:nrow(gtex_projects)) {
    project_info <- gtex_projects[i,]
    project_name <- project_info$project
    
    flog.info("Processing GTEx project: %s", project_name)
    
    rse_result <- create_rse_gene(project_info)
    if(is.null(rse_result) || ncol(rse_result$rse) == 0 || nrow(rse_result$rse) == 0) {
      flog.warn("Invalid RSE object for %s", project_name)
      next
    }
    
    # Get expression values and normalize to TPM
    expression_values <- as.numeric(rse_result$data)
    gene_length <- width(rowRanges(rse_result$rse))
    tpm <- exp(log(expression_values) - log(gene_length) + log(1e6))
    
    # Get tissue types
    tissue_types <- get_tissue_type(rse_result$rse)
    
    # Store results
    results_list[[project_name]] <- data.frame(
      tissue_type = tissue_types,
      sample_id = colnames(rse_result$rse),
      expression = tpm,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine all results
  if(length(results_list) == 0) {
    flog.error("No valid results found from any GTEx project")
    return(NULL)
  }
  
  all_results <- do.call(rbind, results_list)
  
  # Cache results
  save_cached_data(all_results, cache_file)
  
  flog.info("Processed GTEx: %d samples across %d tissue types", 
            nrow(all_results), length(unique(all_results$tissue_type)))
  
  return(all_results)
}

# Function to plot tissue distributions
plot_tissue_distributions <- function(results) {
  # Calculate summary statistics for tissues with at least 10 samples
  summary_stats <- results %>%
    group_by(tissue_type) %>%
    summarise(
      median_expression = median(expression, na.rm = TRUE),
      mean_expression = mean(expression, na.rm = TRUE),
      sd_expression = sd(expression, na.rm = TRUE),
      n_samples = n(),
      .groups = "drop"
    ) %>%
    filter(n_samples >= 10)
  
  # Add sample size to tissue labels
  results$tissue_label <- paste0(
    results$tissue_type, "\n(n=", 
    summary_stats$n_samples[match(results$tissue_type, summary_stats$tissue_type)], 
    ")"
  )
  
  # Create violin plot
  p1 <- ggplot(
    results %>% filter(tissue_type %in% summary_stats$tissue_type),
    aes(x = reorder(tissue_label, expression, FUN = median),
        y = expression,
        fill = tissue_type)
  ) +
    geom_violin(alpha = 0.7, scale = "width") +
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
      title = "SRRM3 Expression Across GTEx Tissues",
      subtitle = paste("Analysis of", length(unique(results$tissue_type)), "tissue types"),
      x = "Tissue Type",
      y = "Expression (TPM)"
    )
  
  # Create means plot
  p2 <- ggplot(summary_stats,
               aes(x = reorder(tissue_type, mean_expression),
                   y = mean_expression,
                   color = tissue_type)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_expression - sd_expression,
                      ymax = mean_expression + sd_expression),
                  width = 0.2) +
    scale_color_viridis(discrete = TRUE) +
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
      title = "Mean SRRM3 Expression by Tissue Type",
      x = "Tissue Type",
      y = "Mean Expression (TPM)"
    )
  
  # Save outputs with consistent dimensions
  tryCatch({
    ggsave("./output/SRRM3_expression_all_tissues_violin.pdf", p1, width = 15, height = 8)
    ggsave("./output/SRRM3_expression_all_tissues_violin.png", p1, width = 15, height = 8)
    ggsave("./output/SRRM3_expression_all_tissues_means.pdf", p2, width = 15, height = 8)
    ggsave("./output/SRRM3_expression_all_tissues_means.png", p2, width = 15, height = 8)
    write.csv(results, "./output/SRRM3_expression_all_tissues.csv", row.names = FALSE)
  }, error = function(e) {
    flog.error("Error saving outputs: %s", e$message)
  })
  
  return(list(violin = p1, means = p2))
}

# Run analysis
results <- analyze_gtex_tissues()
if(!is.null(results)) {
  plot <- plot_tissue_distributions(results)
  print(plot)
  
  # Save results
  write.csv(results, "./output/SRRM3_expression_all_tissues.csv", row.names = FALSE)
} 