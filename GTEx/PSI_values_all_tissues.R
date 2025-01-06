### Output #################################
# "PSI_values_all_tissues.pdf"
# "PSI_values_all_tissues.pdf"
# "PSI_values_all_tissues.csv"
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
flog.appender(appender.file("./logs/PSI_values_all_tissues.log"))
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
  
  # Combine inclusion junctions
  inclusion_jxns <- unique(c(upstream_jxns, downstream_jxns))
  
  return(list(
    inclusion = inclusion_jxns,
    exclusion = exclusion_jxns,
    details = list(
      upstream = upstream_jxns,
      downstream = downstream_jxns
    )
  ))
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

# Modified main analysis function to process all GTEx tissues
main_analysis_all_tissues <- function() {
  # Check for cached final results
  final_cache_file <- file.path(CACHE_DIR, "PSI_values_all_tissues.rds")
  cached_results <- get_cached_data(final_cache_file)
  if (!is.null(cached_results)) {
    return(cached_results)
  }
  
  # Get all GTEx projects
  projects <- available_projects()
  gtex_projects <- subset(projects, file_source == "gtex")
  
  all_results <- list()
  
  # Process each GTEx project
  for(i in 1:nrow(gtex_projects)) {
    project_info <- gtex_projects[i,]
    project_name <- project_info$project
    
    flog.info("Processing GTEx project: %s", project_name)
    
    # Create RSE and calculate PSI
    rse <- create_rse_safe(project_info)
    if(is.null(rse)) {
      flog.warn("Skipping project %s - RSE creation failed", project_name)
      next
    }
    
    psi_values <- calculate_exon15_psi(rse)
    if(is.null(psi_values)) {
      flog.warn("Skipping project %s - PSI calculation failed", project_name)
      next
    }
    
    # Get tissue information
    metadata <- colData(rse) %>% 
      as.data.frame()
    
    # Use the same tissue type field as in SRRM3_total_expression_GTEx.R
    tissue_type <- if("gtex.smtsd" %in% colnames(metadata)) {
      metadata$gtex.smtsd
    } else if("gtex.smts" %in% colnames(metadata)) {
      metadata$gtex.smts
    } else {
      rep("Unknown", nrow(metadata))
    }
    
    # Store results
    all_results[[project_name]] <- data.frame(
      project = project_name,
      sample_id = colnames(rse),
      PSI = psi_values$psi,
      tissue_type = tissue_type,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine all results
  final_results <- do.call(rbind, all_results)
  
  # Cache final results
  save_cached_data(final_results, final_cache_file)
  
  return(final_results)
}

# Run analysis
results <- main_analysis_all_tissues()

# Process results for plotting
plot_data <- results %>%
  filter(!is.na(PSI)) %>%
  filter(tissue_type != "")

# Calculate summary statistics
summary_stats <- plot_data %>%
  group_by(tissue_type) %>%
  summarize(
    mean_PSI = mean(PSI),
    median_PSI = median(PSI),
    sd_PSI = sd(PSI),
    n_samples = n(),
    .groups = 'drop'
  ) %>%
  filter(n_samples >= 10)  # Filter for tissues with at least 10 samples

# Add sample size to tissue labels
plot_data$tissue_label <- paste0(
  plot_data$tissue_type, "\n(n=", 
  summary_stats$n_samples[match(plot_data$tissue_type, summary_stats$tissue_type)], 
  ")"
)

# Create violin plot
p1 <- ggplot(
  plot_data %>% filter(tissue_type %in% summary_stats$tissue_type),
  aes(x = reorder(tissue_label, PSI, FUN = median),
      y = PSI,
      fill = tissue_type)
) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.7, outlier.size = 0.5) +
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
    title = "SRRM3 Exon 15 PSI Values Across GTEx Tissues",
    subtitle = paste("Analysis of", length(unique(plot_data$tissue_type)), "tissue types"),
    x = "Tissue Type",
    y = "Percent Spliced In (PSI)"
  )

# Create means plot
p2 <- ggplot(summary_stats,
             aes(x = reorder(tissue_type, mean_PSI),
                 y = mean_PSI,
                 color = tissue_type)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_PSI - sd_PSI,
                    ymax = mean_PSI + sd_PSI),
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
    title = "Mean SRRM3 Exon 15 PSI Values by Tissue Type",
    x = "Tissue Type",
    y = "Mean PSI"
  )

# Save outputs with consistent dimensions
tryCatch({
  ggsave("./output/PSI_values_all_tissues_violin.pdf", p1, width = 15, height = 8)
  ggsave("./output/PSI_values_all_tissues_means.pdf", p2, width = 15, height = 8)
  write.csv(summary_stats, "./output/PSI_values_all_tissues_summary.csv", row.names = FALSE)
}, error = function(e) {
  flog.error("Error saving outputs: %s", e$message)
}) 