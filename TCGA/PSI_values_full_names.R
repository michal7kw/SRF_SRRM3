### Output #################################
# "PSI_values_full_names.pdf" - Violin plot of PSI distributions
# "PSI_values_full_names.pdf" - Point plot of mean PSI values
# "PSI_values_full_names.csv" - Summary statistics of PSI values
###########################################

# Load required libraries for data processing and visualization
library(recount3)        # For accessing TCGA data
library(GenomicRanges)   # For genomic coordinate operations
library(dplyr)           # For data manipulation
library(tidyr)           # For data tidying
library(ggplot2)         # For creating visualizations
library(ggrepel)         # For improved text labels in plots
library(SummarizedExperiment) # For working with RNA-seq data
library(progress)        # For progress tracking
library(futile.logger)   # For logging operations
library(R.utils)         # For utility functions
library(viridis)         # For color scales in plots

# Source common theme for consistent plot styling
source("common_theme.R")

# Set up logging to track script execution
flog.appender(appender.file("./logs/PSI_values_full_names.log"))
flog.threshold(DEBUG)

# Define cache directory to store intermediate results
CACHE_DIR <- "../cache"
if (!dir.exists(CACHE_DIR)) {
  dir.create(CACHE_DIR)
}

# Function to get cached data to avoid redundant computations
get_cached_data <- function(cache_file) {
  if (file.exists(cache_file)) {
    flog.info("Loading cached data from %s", cache_file)
    return(readRDS(cache_file))
  }
  return(NULL)
}

# Function to save processed data to cache for future use
save_cached_data <- function(data, cache_file) {
  flog.info("Saving data to cache: %s", cache_file)
  saveRDS(data, cache_file)
}

# Define SRRM3 gene and exon 15 information for analysis
SRRM3_INFO <- list(
  gene = list(
    name = "SRRM3",
    chr = "chr7",
    start = 76201896,
    end = 76287287
  ),
  exon15 = list(
    start = 76283524,
    end = 76283602,
    length = 79
  ),
  transcripts = list(
    with_exon15 = "NM_001291831.2",    # Transcript including exon 15
    without_exon15 = "NM_001110199.3"   # Transcript excluding exon 15
  )
)

# Function to identify relevant splicing junctions around exon 15
find_exon15_junctions <- function(jxn_coords) {
  exon_start <- SRRM3_INFO$exon15$start
  exon_end <- SRRM3_INFO$exon15$end
  
  # Debug coordinates for troubleshooting
  flog.debug("Looking for junctions around exon 15: %d-%d", exon_start, exon_end)
  
  # Find upstream junctions (ending at exon start)
  upstream_jxns <- which(abs(end(jxn_coords) - exon_start) <= 5)  # Allow 5bp flexibility
  
  # Find downstream junctions (starting at exon end)
  downstream_jxns <- which(abs(start(jxn_coords) - exon_end) <= 5)  # Allow 5bp flexibility
  
  # Find exclusion junctions (those that skip exon 15)
  exclusion_jxns <- which(
    start(jxn_coords) < (exon_start - 5) & 
    end(jxn_coords) > (exon_end + 5)
  )
  
  # Log junction coordinates for debugging
  if(length(upstream_jxns) > 0) {
    flog.debug("Upstream junction coordinates:")
    for(i in upstream_jxns) {
      flog.debug("Junction %d: %d-%d", i, start(jxn_coords[i]), end(jxn_coords[i]))
    }
  }
  
  if(length(downstream_jxns) > 0) {
    flog.debug("Downstream junction coordinates:")
    for(i in downstream_jxns) {
      flog.debug("Junction %d: %d-%d", i, start(jxn_coords[i]), end(jxn_coords[i]))
    }
  }
  
  if(length(exclusion_jxns) > 0) {
    flog.debug("Exclusion junction coordinates:")
    for(i in exclusion_jxns) {
      flog.debug("Junction %d: %d-%d", i, start(jxn_coords[i]), end(jxn_coords[i]))
    }
  }
  
  # Combine inclusion junctions
  inclusion_jxns <- unique(c(upstream_jxns, downstream_jxns))
  
  flog.debug("Upstream junctions found: %d", length(upstream_jxns))
  flog.debug("Downstream junctions found: %d", length(downstream_jxns))
  flog.debug("Total inclusion junctions: %d", length(inclusion_jxns))
  flog.debug("Exclusion junctions found: %d", length(exclusion_jxns))
  
  return(list(
    inclusion = inclusion_jxns,
    exclusion = exclusion_jxns,
    details = list(
      upstream = upstream_jxns,
      downstream = downstream_jxns
    )
  ))
}

# Function to create RangedSummarizedExperiment (RSE) objects with caching
create_rse_safe <- function(project_info) {  
  # Create cache filename based on project info
  cache_file <- file.path(CACHE_DIR, paste0("rse_", project_info$project, ".rds"))
  
  # Check cache first to avoid redundant processing
  cached_data <- get_cached_data(cache_file)
  if (!is.null(cached_data)) {
    return(cached_data)
  }
  
  tryCatch({
    flog.info("Creating RSE for project: %s", project_info$project)
    
    # Create genomic region around exon 15
    region <- GRanges(
      seqnames = SRRM3_INFO$gene$chr,
      ranges = IRanges(
        start = SRRM3_INFO$exon15$start - 5000,
        end = SRRM3_INFO$exon15$end + 5000
      )
    )
    
    # Use UNIQUE junction format for GTEx and TCGA data
    jxn_format <- if(project_info$file_source %in% c("gtex", "tcga")) "UNIQUE" else "ALL"
    
    # Create RSE object for junction data
    rse <- create_rse(
      project_info,
      type = "jxn",
      jxn_format = jxn_format,
      verbose = TRUE
    )
    
    # Filter to relevant genomic region
    relevant_rows <- which(
      as.character(seqnames(rowRanges(rse))) == as.character(seqnames(region)) &
        start(rowRanges(rse)) >= (start(region)) &
        end(rowRanges(rse)) <= (end(region))
    )
    
    if(length(relevant_rows) == 0) {
      flog.warn("No features found in the specified region")
      return(NULL)
    }
    
    # Create filtered RSE object
    rse_filtered <- rse[relevant_rows, ]
    rm(rse)
    gc()
    
    # Cache the filtered RSE object
    save_cached_data(rse_filtered, cache_file)
    
    return(rse_filtered)
  }, error = function(e) {
    flog.error("Error creating RSE: %s", e$message)
    return(NULL)
  })
}

# Function to calculate Percent Spliced In (PSI) values for exon 15
calculate_exon15_psi <- function(rse) {
  if(is.null(rse)) return(NULL)
  
  # Create cache key based on RSE characteristics
  cache_key <- digest::digest(list(
    coords = as.data.frame(rowRanges(rse)),
    counts = assay(rse)
  ))
  cache_file <- file.path(CACHE_DIR, paste0("psi_", cache_key, ".rds"))
  
  # Check cache first
  cached_data <- get_cached_data(cache_file)
  if (!is.null(cached_data)) {
    return(cached_data)
  }
  
  # Get junction counts from RSE object
  junction_counts <- assay(rse)
  jxn_coords <- rowRanges(rse)
  
  # Find relevant junctions for PSI calculation
  junctions <- find_exon15_junctions(jxn_coords)
  
  if(length(junctions$inclusion) == 0 || length(junctions$exclusion) == 0) {
    flog.warn("Missing junctions - Inclusion: %d, Exclusion: %d", 
              length(junctions$inclusion), length(junctions$exclusion))
    return(NULL)
  }
  
  # Calculate PSI for each sample
  psi_values <- sapply(seq_len(ncol(junction_counts)), function(i) {
    inclusion_reads <- sum(junction_counts[junctions$inclusion, i])
    exclusion_reads <- sum(junction_counts[junctions$exclusion, i])
    total_reads <- inclusion_reads + exclusion_reads
    
    if(total_reads >= 10) {  # Apply minimum coverage threshold
      psi <- (inclusion_reads / total_reads) * 100
      return(psi)
    } else {
      return(NA)
    }
  })
  
  # Create results object
  results <- list(
    psi = psi_values,
    junctions = list(
      inclusion = rownames(junction_counts)[junctions$inclusion],
      exclusion = rownames(junction_counts)[junctions$exclusion]
    )
  )
  
  # Cache results
  save_cached_data(results, cache_file)
  
  # Log summary statistics
  flog.info("PSI calculation complete:")
  flog.info("Total samples: %d", length(psi_values))
  flog.info("Non-NA samples: %d", sum(!is.na(psi_values)))
  flog.info("Mean PSI: %.2f", mean(psi_values, na.rm = TRUE))
  
  return(results)
}

# Modified process_tcga_metadata function with caching
process_tcga_metadata <- function(rse) {
  # Create cache filename based on metadata characteristics
  cache_key <- digest::digest(colData(rse))
  cache_file <- file.path(CACHE_DIR, paste0("metadata_", cache_key, ".rds"))
  
  # Check cache first
  cached_data <- get_cached_data(cache_file)
  if (!is.null(cached_data)) {
    return(cached_data)
  }
  
  # Extract metadata
  metadata <- colData(rse) %>% 
    as.data.frame()
  
  # Extract cancer type from TCGA metadata
  metadata$cancer_type <- metadata$tcga.cgc_case_histological_diagnosis
  
  # Clean up cancer type names
  metadata$cancer_type <- gsub("_", " ", metadata$cancer_type)
  
  # Log unique cancer types
  unique_types <- unique(metadata$cancer_type)
  flog.info("Found %d unique cancer types:", length(unique_types))
  for(type in unique_types) {
    flog.info("  - %s", type)
  }
  
  # Cache metadata
  save_cached_data(metadata, cache_file)
  
  return(metadata)
}

# Main analysis with caching
main_analysis <- function() {
  # Check for cached final results
  final_cache_file <- file.path(CACHE_DIR, "PSI_values_full_names.rds")
  cached_results <- get_cached_data(final_cache_file)
  if (!is.null(cached_results)) {
    return(cached_results)
  }
  
  # Get available TCGA projects
  projects <- available_projects()
  tcga_projects <- subset(projects, file_source == "tcga")
  
  # Initialize lists to store results
  all_psi_data <- list()
  all_metadata <- list()
  
  # Process each TCGA project
  for(i in 1:nrow(tcga_projects)) {
    project_name <- tcga_projects$project[i]
    flog.info("Processing TCGA project %d/%d: %s", i, nrow(tcga_projects), project_name)
    
    # Create RSE object
    rse <- create_rse_safe(tcga_projects[i,])
    
    if(is.null(rse)) {
      flog.error("Failed to create RSE for project %s", project_name)
      next
    }
    
    # Calculate PSI values
    psi_result <- calculate_exon15_psi(rse)
    
    if(is.null(psi_result)) {
      flog.error("Failed to calculate PSI values for project %s", project_name)
      next
    }
    
    # Process metadata
    metadata <- process_tcga_metadata(rse)
    
    # Store results
    all_psi_data[[project_name]] <- psi_result$psi
    all_metadata[[project_name]] <- metadata
  }
  
  # Combine all data
  final_results <- list(
    psi_data = all_psi_data,
    metadata = all_metadata
  )
  
  # Cache final results
  save_cached_data(final_results, final_cache_file)
  
  return(final_results)
}
# Run analysis and create plots
results <- main_analysis()

# Create plot data
plot_data <- data.frame(
  sample_id = unlist(lapply(names(results$psi_data), function(proj) {
    paste0(proj, "_", seq_along(results$psi_data[[proj]]))
  })),
  PSI = unlist(results$psi_data),
  cancer_type = unlist(lapply(names(results$metadata), function(proj) {
    results$metadata[[proj]]$cancer_type
  }))
) %>%
  filter(!is.na(PSI)) %>%
  filter(cancer_type != "")

# Calculate summary statistics
summary_stats <- plot_data %>%
  group_by(cancer_type) %>%
  summarize(
    mean_PSI = mean(PSI),
    median_PSI = median(PSI),
    sd_PSI = sd(PSI),
    n_samples = n(),
    .groups = 'drop'
  )

# Create and save plots
p1 <- ggplot(plot_data, aes(x = reorder(cancer_type, PSI, FUN = median), 
                           y = PSI, 
                           fill = cancer_type)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
  scale_fill_viridis_d() +
  get_tcga_theme() +
  labs(
    title = "SRRM3 Exon 15 PSI Values Across Cancer Types",
    subtitle = "Distribution of Percent Spliced In values",
    x = "Cancer Type",
    y = "Percent Spliced In (PSI)",
    fill = "Cancer Type"
  )

p2 <- ggplot(summary_stats, 
             aes(x = reorder(cancer_type, mean_PSI), 
                 y = mean_PSI)) +
  geom_point(aes(size = n_samples), alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_PSI - sd_PSI, 
                    ymax = mean_PSI + sd_PSI), 
                width = 0.2) +
  scale_size_continuous(name = "Number of Samples") +
  get_tcga_theme() +
  labs(
    title = "Mean SRRM3 Exon 15 PSI Values by Cancer Type",
    subtitle = "Error bars show standard deviation",
    x = "Cancer Type",
    y = "Mean PSI"
  )

# Save plots and results
tryCatch({
  ggsave("./output/PSI_values_full_namess_violin.pdf", p1, width = 15, height = 8)
  ggsave("./output/PSI_values_full_names_means.pdf", p2, width = 15, height = 8)
  write.csv(summary_stats, "./output/PSI_values_full_names_summary.csv", row.names = FALSE)
}, error = function(e) {
  flog.error("Error saving outputs: %s", e$message)
})
