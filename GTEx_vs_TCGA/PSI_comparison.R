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
library(digest)

# Set up logging
flog.appender(appender.file("./logs/PSI_comparison.log"))
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
  ),
  transcripts = list(
    with_exon15 = "NM_001291831.2",    # 16 exons
    without_exon15 = "NM_001110199.3"   # 15 exons
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
  cache_file <- file.path(CACHE_DIR, paste0("rse_comparison_", 
                                          project_info$project, "_",
                                          project_info$file_source, ".rds"))
  
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
  cache_file <- file.path(CACHE_DIR, paste0("psi_comparison_", cache_key, ".rds"))
  
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

# Main analysis with caching
main_analysis <- function() {
  # Check for cached final results
  final_cache_file <- file.path(CACHE_DIR, "final_comparison_analysis.rds")
  cached_results <- get_cached_data(final_cache_file)
  if (!is.null(cached_results)) {
    return(cached_results)
  }
  
  # Get projects
  projects <- available_projects()
  
  # Process GTEx brain data
  gtex_project <- subset(projects, project == "BRAIN" & file_source == "gtex")
  gtex_rse <- create_rse_safe(gtex_project)
  
  # Process TCGA GBM data
  tcga_project <- subset(projects, project == "GBM" & file_source == "tcga")
  tcga_rse <- create_rse_safe(tcga_project)
  
  # Combine results
  final_results <- list(
    gtex = list(
      rse = gtex_rse,
      psi = calculate_exon15_psi(gtex_rse)
    ),
    tcga = list(
      GBM = list(
        rse = tcga_rse,
        psi = calculate_exon15_psi(tcga_rse)
      )
    )
  )
  
  # Cache final results
  save_cached_data(final_results, final_cache_file)
  
  return(final_results)
}

# Run the analysis
results <- main_analysis()

# Prepare data for plotting
plot_data <- data.frame(
  PSI = numeric(),
  dataset = character(),
  stringsAsFactors = FALSE
)

# Add GTEx data
if (!is.null(results$gtex$psi)) {
  plot_data <- rbind(plot_data, data.frame(
    PSI = results$gtex$psi$psi,
    dataset = "GTEx Brain",
    stringsAsFactors = FALSE
  ))
}

# Add TCGA data
for (project in names(results$tcga)) {
  if (!is.null(results$tcga[[project]]$psi)) {
    plot_data <- rbind(plot_data, data.frame(
      PSI = results$tcga[[project]]$psi$psi,
      dataset = project,
      stringsAsFactors = FALSE
    ))
  }
}

# After loading the data, clean it up
plot_data <- plot_data %>%
  # Remove NA values
  filter(!is.na(PSI)) %>%
  # Ensure PSI values are between 0 and 100
  mutate(PSI = pmax(0, pmin(100, PSI))) %>%
  # Combine all TCGA samples into one group
  mutate(dataset = ifelse(dataset == "GTEx Brain", "GTEx", "TCGA"))

# Create plot
p <- ggplot(plot_data, aes(x = dataset, y = PSI, fill = dataset)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Hide boxplot outliers since we're showing all points
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.5, color = "grey30") +
  scale_fill_manual(values = c("GTEx" = "#FFB6B6", "TCGA" = "#7FDBDB")) +
  theme_minimal() +
  labs(
    title = "SRRM3 Exon 15 PSI Values",
    y = "Percent Spliced In (PSI)",
    x = "Dataset"
  ) +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    axis.text = element_text(color = "grey30"),
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 12)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25)
  )

# Save results
output_dir <- "./output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

ggsave(file.path(output_dir, "PSI_comparison.pdf"), p, width = 8, height = 6)

# Perform statistical test
test_result <- wilcox.test(
  PSI ~ dataset, 
  data = plot_data,
  alternative = "two.sided"
)

# Save summary statistics with combined TCGA samples
summary_stats <- plot_data %>%
  group_by(dataset) %>%
  summarise(
    mean_PSI = mean(PSI, na.rm = TRUE),
    median_PSI = median(PSI, na.rm = TRUE),
    sd_PSI = sd(PSI, na.rm = TRUE),
    n_samples = n()
  )

write.csv(
  summary_stats,
  file.path(output_dir, "PSI_comparison_summary.csv"),
  row.names = FALSE
)

# Log results
flog.info("Analysis complete. Results saved to output directory.")
flog.info("Statistical test p-value: %g", test_result$p.value)
