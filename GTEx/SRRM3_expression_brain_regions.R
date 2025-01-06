### Output #################################
# "SRRM3_expression_brain_regions.pdf"
# "SRRM3_expression_brain_regions.pdf"
# "SRRM3_expression_brain_regions.csv"
###########################################

# Load required libraries
library(recount3)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(SummarizedExperiment)
library(futile.logger)
library(digest)
library(viridis)
library(scales)

# Set up logging
flog.appender(appender.file("./logs/SRRM3_expression_brain_regions.log"))
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

# Modified RSE creation function with caching
create_rse_brain <- function(project_info) {
  # Create cache filename based on project info
  cache_file <- file.path(CACHE_DIR, paste0("rse_brain_expr_", project_info$project, ".rds"))
  
  # Check cache first
  cached_data <- get_cached_data(cache_file)
  if (!is.null(cached_data)) {
    return(cached_data)
  }
  
  tryCatch({
    flog.info("Creating RSE for project: %s", project_info$project)
    
    rse <- create_rse(
      project_info,
      type = "gene",
      verbose = TRUE
    )
    
    # Filter to SRRM3
    rse_filtered <- rse[rowRanges(rse)$gene_name == "SRRM3", ]
    
    if(nrow(rse_filtered) == 0) {
      flog.warn("No SRRM3 gene found in the data")
      return(NULL)
    }
    
    # Cache the filtered RSE object
    save_cached_data(rse_filtered, cache_file)
    
    return(rse_filtered)
  }, error = function(e) {
    flog.error("Error creating RSE: %s", e$message)
    return(NULL)
  })
}

# Function to process brain metadata (similar to PSI_values_GTEx.R)
process_brain_metadata <- function(rse) {
  # Create cache filename based on metadata characteristics
  cache_key <- digest::digest(as.data.frame(colData(rse)))
  cache_file <- file.path(CACHE_DIR, paste0("metadata_brain_expr_", cache_key, ".rds"))
  
  # Check cache first
  cached_data <- get_cached_data(cache_file)
  if (!is.null(cached_data)) {
    return(cached_data)
  }
  
  # Extract metadata
  metadata <- colData(rse) %>% 
    as.data.frame()
  
  # Clean up tissue type information from GTEx metadata
  metadata$brain_region <- metadata$gtex.smtsd
  
  # Remove "Brain -" prefix if present
  metadata$brain_region <- gsub("^Brain - ", "", metadata$brain_region)
  
  # Log unique brain regions found
  unique_regions <- unique(metadata$brain_region)
  flog.info("Found %d unique brain regions:", length(unique_regions))
  for(region in unique_regions) {
    flog.info("  - %s", region)
  }
  
  # Cache the processed metadata
  save_cached_data(metadata, cache_file)
  
  return(metadata)
}

# Main analysis function
main_analysis <- function() {
  # Check for cached final results
  final_cache_file <- file.path(CACHE_DIR, "SRRM3_expression_brain_regions.rds")
  cached_results <- get_cached_data(final_cache_file)
  if (!is.null(cached_results)) {
    return(cached_results)
  }
  
  # Get projects
  projects <- available_projects()
  
  # Process GTEx brain data
  gtex_project <- subset(projects, project == "BRAIN" & file_source == "gtex")
  gtex_rse <- create_rse_brain(gtex_project)
  
  if(is.null(gtex_rse)) {
    flog.error("Failed to create RSE for GTEx brain data")
    stop("RSE creation failed")
  }
  
  # Get expression values and normalize to TPM
  expression_values <- as.numeric(assay(gtex_rse))
  gene_length <- width(rowRanges(gtex_rse))
  tpm <- exp(log(expression_values) - log(gene_length) + log(1e6))
  
  # Process metadata
  metadata <- process_brain_metadata(gtex_rse)
  
  # Combine results
  final_results <- list(
    rse = gtex_rse,
    tpm = tpm,
    metadata = metadata
  )
  
  # Cache final results
  save_cached_data(final_results, final_cache_file)
  
  return(final_results)
}

# Run analysis
results <- main_analysis()

# Create plot data
plot_data <- data.frame(
  sample_id = colnames(results$rse),
  expression = results$tpm,
  brain_region = results$metadata$brain_region
) %>%
  filter(!is.na(expression)) %>%  # Remove NA values
  filter(brain_region != "") # Remove empty brain regions

# Calculate summary statistics
summary_stats <- plot_data %>%
  group_by(brain_region) %>%
  summarize(
    mean_expression = mean(expression),
    median_expression = median(expression),
    sd_expression = sd(expression),
    n_samples = n(),
    .groups = 'drop'
  )

# Add small pseudocount to expression values for log transformation
plot_data$expression_log <- plot_data$expression + 1  # Add pseudocount of 1

# Create violin plot
p1 <- ggplot(plot_data, aes(x = reorder(brain_region, expression_log, FUN = median), 
                           y = expression_log, 
                           fill = brain_region)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.7, outlier.size = 0.5) +
  scale_y_continuous(
    trans = "log10",
    labels = trans_format("log10", math_format(10^.x)),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    limits = c(NA, NA)  # Allow automatic limit detection
  ) +
  scale_fill_viridis(discrete = TRUE) +
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
    title = "SRRM3 Expression Across Brain Regions",
    subtitle = "Log10(TPM + 1) transformation",
    y = "Expression (TPM + 1)",
    x = "Brain Region"
  )

# Update means plot to use pseudocount too
summary_stats$mean_expression_log <- log10(summary_stats$mean_expression + 1)
summary_stats$sd_expression_log <- summary_stats$sd_expression / (log(10) * (summary_stats$mean_expression + 1))

p2 <- ggplot(summary_stats, 
             aes(x = reorder(brain_region, mean_expression_log), 
                 y = mean_expression_log, 
                 color = brain_region)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_expression_log - sd_expression_log, 
                    ymax = mean_expression_log + sd_expression_log), 
                width = 0.2) +
  scale_color_viridis(discrete = TRUE) +
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
    title = "Mean SRRM3 Expression by Brain Region",
    subtitle = "Log10(TPM + 1) transformation",
    y = "Mean Expression (log10(TPM + 1))",
    x = "Brain Region"
  )

# Save outputs
tryCatch({
  ggsave("./output/SRRM3_expression_brain_regions_violin.pdf", p1, width = 12, height = 7)
  ggsave("./output/SRRM3_expression_brain_regions_means.pdf", p2, width = 12, height = 7)
  write.csv(summary_stats, "./output/SRRM3_expression_brain_regions.csv", row.names = FALSE)
}, error = function(e) {
  flog.error("Error saving outputs: %s", e$message)
}) 