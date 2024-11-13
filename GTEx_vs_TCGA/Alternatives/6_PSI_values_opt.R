# Load required libraries
library(recount3)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(SummarizedExperiment)
library(progress)  # For progress bars
library(futile.logger)  # For better logging

# Set up logging
flog.appender(appender.file("srrm3_analysis.log"))
flog.threshold(DEBUG)

# Set memory options
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB
if(.Platform$OS.type == "windows") {
  memory.limit(size = 32000)  # 32GB for Windows
}

# Function to safely create RSE with progress tracking
create_rse_safe <- function(project_info, type, region) {
  tryCatch({
    flog.info("Creating RSE for project: %s", project_info$project)
    
    # Create RSE without filtering
    rse <- create_rse(project_info, type = type)
    
    # Filter to relevant chromosome and region after creation
    chr <- as.character(seqnames(region))
    relevant_rows <- which(
      as.character(seqnames(rowRanges(rse))) == chr &
        start(rowRanges(rse)) >= (start(region) - 10000) &
        end(rowRanges(rse)) <= (end(region) + 10000)
    )
    
    # Subset the RSE object
    rse_filtered <- rse[relevant_rows, ]
    
    flog.info("Successfully created RSE with dimensions: %d x %d", 
              nrow(rse_filtered), ncol(rse_filtered))
    return(rse_filtered)
    
  }, error = function(e) {
    flog.error("Error creating RSE: %s", e$message)
    stop(e)
  })
}

# Improved PSI calculation function
calculate_psi <- function(project_info, gene_coords, chunk_size = 500) {
  flog.info("Starting PSI calculation for project: %s", project_info$project)
  
  # Create region GRanges
  region <- GRanges(
    seqnames = "chr7",
    ranges = IRanges(
      start = gene_coords$start,
      end = gene_coords$end
    )
  )
  
  # Get junction-level data
  rse_jxn <- create_rse_safe(project_info, "jxn", region)
  
  # Get overlapping junctions
  flog.info("Finding overlapping junctions")
  overlaps <- findOverlaps(rowRanges(rse_jxn), region)
  
  if(length(overlaps) == 0) {
    flog.warn("No junctions found in the specified region")
    return(NULL)
  }
  
  # Extract junction data
  jxn_data <- rse_jxn[queryHits(overlaps),]
  junction_counts <- assay(jxn_data)
  jxn_coords <- rowRanges(jxn_data)
  
  # Create junction identifiers
  jxn_ids <- paste0(
    seqnames(jxn_coords), ":",
    start(jxn_coords), "-",
    end(jxn_coords)
  )
  
  # Initialize results matrix
  psi_matrix <- matrix(NA, nrow = length(jxn_ids), ncol = ncol(junction_counts))
  rownames(psi_matrix) <- jxn_ids
  colnames(psi_matrix) <- colnames(junction_counts)
  
  # Process samples in chunks
  sample_chunks <- split(
    seq_len(ncol(junction_counts)), 
    ceiling(seq_len(ncol(junction_counts))/chunk_size)
  )
  
  # Create progress bar
  pb <- progress_bar$new(
    format = "Processing samples [:bar] :percent eta: :eta",
    total = length(sample_chunks)
  )
  
  # Process each chunk
  for(chunk_idx in seq_along(sample_chunks)) {
    chunk_samples <- sample_chunks[[chunk_idx]]
    
    flog.debug("Processing chunk %d/%d (samples %d-%d)", 
               chunk_idx, length(sample_chunks),
               min(chunk_samples), max(chunk_samples))
    
    # For each junction
    for(i in seq_along(jxn_ids)) {
      # Find alternative junctions
      alt_5prime <- which(start(jxn_coords) == start(jxn_coords[i]) & 
                            seq_along(jxn_coords) != i)
      alt_3prime <- which(end(jxn_coords) == end(jxn_coords[i]) & 
                            seq_along(jxn_coords) != i)
      
      # Calculate PSI for chunk samples
      for(j in chunk_samples) {
        inclusion_reads <- junction_counts[i, j]
        
        alt_5prime_reads <- if(length(alt_5prime) > 0) 
          sum(junction_counts[alt_5prime, j]) else 0
        alt_3prime_reads <- if(length(alt_3prime) > 0) 
          sum(junction_counts[alt_3prime, j]) else 0
        
        total_reads <- inclusion_reads + alt_5prime_reads + alt_3prime_reads
        if(total_reads > 0) {
          psi_matrix[i, j] <- inclusion_reads / total_reads * 100
        }
      }
    }
    
    # Update progress bar
    pb$tick()
    
    # Clean up memory
    gc()
  }
  
  flog.info("PSI calculation completed")
  
  return(list(
    psi_values = psi_matrix,
    junction_info = data.frame(
      junction_id = jxn_ids,
      start = start(jxn_coords),
      end = end(jxn_coords),
      stringsAsFactors = FALSE
    )
  ))
}

# Improved comparison function
create_psi_comparison <- function(tcga_psi, gtex_psi) {
  flog.info("Creating PSI comparison")
  
  # Find common junctions
  common_junctions <- intersect(
    rownames(tcga_psi$psi_values),
    rownames(gtex_psi$psi_values)
  )
  
  flog.info("Found %d common junctions", length(common_junctions))
  
  # Create progress bar for statistical tests
  pb <- progress_bar$new(
    format = "Calculating statistics [:bar] :percent eta: :eta",
    total = length(common_junctions)
  )
  
  # Calculate statistics
  comparison_stats <- lapply(common_junctions, function(j) {
    stats <- tryCatch({
      test_result <- wilcox.test(
        tcga_psi$psi_values[j,],
        gtex_psi$psi_values[j,],
        na.rm = TRUE
      )
      
      list(
        junction_id = j,
        tcga_psi = mean(tcga_psi$psi_values[j,], na.rm = TRUE),
        gtex_psi = mean(gtex_psi$psi_values[j,], na.rm = TRUE),
        p_value = test_result$p.value
      )
    }, error = function(e) {
      flog.warn("Error calculating stats for junction %s: %s", j, e$message)
      NULL
    })
    
    pb$tick()
    return(stats)
  })
  
  # Remove NULL results and convert to data frame
  comparison_df <- do.call(rbind, Filter(Negate(is.null), comparison_stats)) %>%
    as.data.frame() %>%
    mutate(
      delta_psi = tcga_psi - gtex_psi,
      significance = -log10(p_value),
      is_significant = p_value < 0.05 & abs(delta_psi) > 10
    )
  
  flog.info("Comparison completed. Found %d significant changes", 
            sum(comparison_df$is_significant))
  
  return(comparison_df)
}

# Main execution
main <- function() {
  flog.info("Starting SRRM3 analysis")
  
  # SRRM3 coordinates
  srrm3_coords <- list(
    start = 76201900,
    end = 76287288
  )
  
  # Get projects
  projects <- available_projects()
  
  # Calculate PSI values
  tcga_project <- subset(projects, project == "GBM" & file_source == "tcga")
  gtex_project <- subset(projects, project == "BRAIN" & file_source == "gtex")
  
  flog.info("Processing TCGA data")
  tcga_psi <- calculate_psi(tcga_project, srrm3_coords)
  
  flog.info("Processing GTEx data")
  gtex_psi <- calculate_psi(gtex_project, srrm3_coords)
  
  # Create comparison
  comparison_data <- create_psi_comparison(tcga_psi, gtex_psi)
  
  # Create visualization
  flog.info("Creating visualization")
  p <- ggplot(comparison_data, aes(x = gtex_psi, y = tcga_psi)) +
    geom_point(aes(color = significance), alpha = 0.7, size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    geom_label_repel(
      data = subset(comparison_data, is_significant),
      aes(label = junction_id),
      size = 3,
      box.padding = 0.5,
      point.padding = 0.5,
      max.overlaps = 20
    ) +
    scale_color_viridis_c(name = "-log10(p-value)") +
    scale_x_continuous(limits = c(0, 100)) +
    scale_y_continuous(limits = c(0, 100)) +
    theme_minimal() +
    labs(
      title = "SRRM3 PSI Comparison",
      subtitle = "TCGA Glioma vs GTEx Brain",
      x = "Mean PSI (GTEx)",
      y = "Mean PSI (TCGA)",
      caption = paste0(
        "Significant changes (p < 0.05, |ΔPSI| > 10): ",
        sum(comparison_data$is_significant)
      )
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  
  # Save results
  flog.info("Saving results")
  pdf("srrm3_psi_comparison_opt.pdf", width = 12, height = 10)
  print(p)
  dev.off()
  
  write.csv(comparison_data, "srrm3_psi_summary_opt.csv", row.names = FALSE)
  
  flog.info("Analysis completed successfully")
}

# Run the analysis with error handling
tryCatch({
  main()
}, error = function(e) {
  flog.error("Analysis failed: %s", e$message)
  stop(e)
})