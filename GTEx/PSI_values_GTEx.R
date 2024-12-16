### Output #################################
# "PSI_values_GTEx_brain_regions_violin.pdf"
# "PSI_values_GTEx_brain_regions_means.pdf"
# "PSI_values_GTEx_brain_regions_summary.csv"
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

# Set up logging
flog.appender(appender.file("./logs/psi_values_gtex.log"))
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
  ),
  exon15 = list(
    start = 76283524,
    end = 76283602,
    length = 79
  )
)

# Function to find relevant junctions for exon 15
find_exon15_junctions <- function(jxn_coords) {
  exon_start <- SRRM3_INFO$exon15$start
  exon_end <- SRRM3_INFO$exon15$end
  
  # Debug coordinates
  flog.debug("Looking for junctions around exon 15: %d-%d", exon_start, exon_end)
  
  # Find upstream junction (ending at exon start)
  upstream_jxns <- which(abs(end(jxn_coords) - exon_start) <= 5)  # Allow 5bp flexibility
  
  # Find downstream junction (starting at exon end)
  downstream_jxns <- which(abs(start(jxn_coords) - exon_end) <= 5)  # Allow 5bp flexibility
  
  # Find exclusion junctions (those that skip exon 15)
  exclusion_jxns <- which(
    start(jxn_coords) < (exon_start - 5) & 
    end(jxn_coords) > (exon_end + 5)
  )
  
  # Debug junction coordinates
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

# Modified RSE creation function with caching
create_rse_safe <- function(project_info) {  
  # Create cache filename based on project info
  cache_file <- file.path(CACHE_DIR, paste0("rse_gtex_", project_info$project, ".rds"))
  
  # Check cache first
  cached_data <- get_cached_data(cache_file)
  if (!is.null(cached_data)) {
    return(cached_data)
  }
  
  tryCatch({
    flog.info("Creating RSE for project: %s", project_info$project)
    
    # Create region around exon 15
    region <- GRanges(
      seqnames = SRRM3_INFO$gene$chr,
      ranges = IRanges(
        start = SRRM3_INFO$exon15$start - 5000,  # Include surrounding region
        end = SRRM3_INFO$exon15$end + 5000
      )
    )
    
    # Use UNIQUE junction format for GTEx and TCGA
    jxn_format <- if(project_info$file_source %in% c("gtex", "tcga")) "UNIQUE" else "ALL"
    
    rse <- create_rse(
      project_info,
      type = "jxn",
      jxn_format = jxn_format,
      verbose = TRUE
    )
    
    # Filter to relevant region
    relevant_rows <- which(
      as.character(seqnames(rowRanges(rse))) == as.character(seqnames(region)) &
        start(rowRanges(rse)) >= (start(region)) &
        end(rowRanges(rse)) <= (end(region))
    )
    
    if(length(relevant_rows) == 0) {
      flog.warn("No features found in the specified region")
      return(NULL)
    }
    
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

# Modified PSI calculation function with caching
calculate_exon15_psi <- function(rse) {
  if(is.null(rse)) return(NULL)
  
  # Create cache filename based on RSE object characteristics
  cache_key <- digest::digest(list(
    coords = as.data.frame(rowRanges(rse)),
    counts = assay(rse)
  ))
  cache_file <- file.path(CACHE_DIR, paste0("psi_gtex_", cache_key, ".rds"))
  
  # Check cache first
  cached_data <- get_cached_data(cache_file)
  if (!is.null(cached_data)) {
    return(cached_data)
  }
  
  # Get junction counts
  junction_counts <- assay(rse)
  jxn_coords <- rowRanges(rse)
  
  # Find relevant junctions
  junctions <- find_exon15_junctions(jxn_coords)
  
  if(length(junctions$inclusion) == 0 || length(junctions$exclusion) == 0) {
    flog.warn("Missing junctions - Inclusion: %d, Exclusion: %d", 
              length(junctions$inclusion), length(junctions$exclusion))
    return(NULL)
  }
  
  # Calculate PSI for each sample
  psi_values <- sapply(seq_len(ncol(junction_counts)), function(i) {
    # Get read counts
    inclusion_reads <- sum(junction_counts[junctions$inclusion, i])
    exclusion_reads <- sum(junction_counts[junctions$exclusion, i])
    total_reads <- inclusion_reads + exclusion_reads
    
    if(total_reads >= 10) {  # Minimum coverage threshold
      psi <- (inclusion_reads / total_reads) * 100
      return(psi)
    } else {
      return(NA)
    }
  })
  
  result <- list(
    psi = psi_values,
    junctions = list(
      inclusion = rownames(junction_counts)[junctions$inclusion],
      exclusion = rownames(junction_counts)[junctions$exclusion]
    )
  )
  
  # Cache the results
  save_cached_data(result, cache_file)
  
  return(result)
}

# Modified metadata processing function with caching
process_brain_metadata <- function(rse) {
  # Create cache filename based on metadata characteristics
  cache_key <- digest::digest(as.data.frame(colData(rse)))
  cache_file <- file.path(CACHE_DIR, paste0("metadata_gtex_", cache_key, ".rds"))
  
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

# Main analysis with caching
main_analysis <- function() {
  # Check for cached final results
  final_cache_file <- file.path(CACHE_DIR, "final_gtex_analysis.rds")
  cached_results <- get_cached_data(final_cache_file)
  if (!is.null(cached_results)) {
    return(cached_results)
  }
  
  # Get projects
  projects <- available_projects()
  
  # Process GTEx brain data
  gtex_project <- subset(projects, project == "BRAIN" & file_source == "gtex")
  gtex_rse <- create_rse_safe(gtex_project)
  
  if(is.null(gtex_rse)) {
    flog.error("Failed to create RSE for GTEx brain data")
    stop("RSE creation failed")
  }
  
  # Calculate PSI values
  gtex_psi <- calculate_exon15_psi(gtex_rse)
  
  if(is.null(gtex_psi)) {
    flog.error("Failed to calculate PSI values for GTEx brain data")
    stop("PSI calculation failed")
  }
  
  # Process metadata
  metadata <- process_brain_metadata(gtex_rse)
  
  # Combine results
  final_results <- list(
    rse = gtex_rse,
    psi = gtex_psi,
    metadata = metadata
  )
  
  # Cache final results
  save_cached_data(final_results, final_cache_file)
  
  return(final_results)
}

# Run analysis and create plots
results <- main_analysis()

# Create plot data
plot_data <- data.frame(
  sample_id = colnames(results$rse),
  PSI = results$psi$psi,
  brain_region = results$metadata$brain_region
) %>%
  filter(!is.na(PSI)) %>%  # Remove NA values
  filter(brain_region != "") # Remove empty brain regions

# Calculate summary statistics
summary_stats <- plot_data %>%
  group_by(brain_region) %>%
  summarize(
    mean_PSI = mean(PSI),
    median_PSI = median(PSI),
    sd_PSI = sd(PSI),
    n_samples = n(),
    .groups = 'drop'
  )

# Create violin plot
p1 <- ggplot(plot_data, aes(x = reorder(brain_region, PSI, FUN = median), 
                           y = PSI, 
                           fill = brain_region)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "SRRM3 Exon 15 PSI Values Across Brain Regions",
    y = "Percent Spliced In (PSI)",
    x = "Brain Region"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(size = 12)
  )

# Create scatter plot with means and error bars
p2 <- ggplot(summary_stats, 
             aes(x = reorder(brain_region, mean_PSI), 
                 y = mean_PSI, 
                 color = brain_region)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_PSI - sd_PSI, 
                    ymax = mean_PSI + sd_PSI), 
                width = 0.2) +
  theme_minimal() +
  labs(
    title = "Mean SRRM3 Exon 15 PSI Values by Brain Region",
    y = "Mean PSI",
    x = "Brain Region"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(size = 12)
  )

# Save outputs
tryCatch({
  ggsave("./output/PSI_values_GTEx_brain_regions_violin.pdf", p1, width = 12, height = 7)
  ggsave("./output/PSI_values_GTEx_brain_regions_means.pdf", p2, width = 12, height = 7)
  write.csv(summary_stats, "./output/PSI_values_GTEx_brain_regions_summary.csv", row.names = FALSE)
}, error = function(e) {
  flog.error("Error saving outputs: %s", e$message)
})
