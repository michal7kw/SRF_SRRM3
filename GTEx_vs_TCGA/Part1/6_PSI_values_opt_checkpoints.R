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
library(viridis)  # For color scales

# Set up logging
flog.appender(appender.file("srrm3_analysis.log"))
flog.threshold(DEBUG)

# Set memory options and use AWS mirror for potentially faster access
options(future.globals.maxSize = 8000 * 1024^2)
options(recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release")

# Function to save checkpoint
save_checkpoint <- function(data, checkpoint_name) {
  checkpoint_file <- paste0("checkpoint_", checkpoint_name, ".rds")
  flog.info("Saving checkpoint: %s", checkpoint_file)
  saveRDS(data, checkpoint_file)
}

# Function to load checkpoint
load_checkpoint <- function(checkpoint_name) {
  checkpoint_file <- paste0("checkpoint_", checkpoint_name, ".rds")
  if(file.exists(checkpoint_file)) {
    flog.info("Loading checkpoint: %s", checkpoint_file)
    return(readRDS(checkpoint_file))
  }
  return(NULL)
}

# Optimized RSE creation function
create_rse_safe <- function(project_info, type, region) {
  checkpoint_name <- paste0(project_info$project, "_rse")
  rse_filtered <- load_checkpoint(checkpoint_name)
  
  if(!is.null(rse_filtered)) {
    flog.info("Loaded RSE from checkpoint")
    return(rse_filtered)
  }
  
  tryCatch({
    flog.info("Creating RSE for project: %s", project_info$project)
    
    # Use UNIQUE junction format for GTEx and TCGA
    jxn_format <- if(project_info$file_source %in% c("gtex", "tcga")) "UNIQUE" else "ALL"
    
    # Create RSE with optimized parameters
    rse <- create_rse(
      project_info,
      type = type,
      jxn_format = jxn_format,
      verbose = TRUE
    )
    
    # Filter to SRRM3 region
    chr <- as.character(seqnames(region))
    relevant_rows <- which(
      as.character(seqnames(rowRanges(rse))) == chr &
        start(rowRanges(rse)) >= (start(region) - 10000) &
        end(rowRanges(rse)) <= (end(region) + 10000)
    )
    
    if(length(relevant_rows) == 0) {
      flog.warn("No features found in the specified region")
      return(NULL)
    }
    
    # Subset the RSE object
    rse_filtered <- rse[relevant_rows, ]
    
    # Clean up to free memory
    rm(rse)
    gc()
    
    flog.info("Successfully created RSE with dimensions: %d x %d", 
              nrow(rse_filtered), ncol(rse_filtered))
    
    # Save checkpoint
    save_checkpoint(rse_filtered, checkpoint_name)
    
    return(rse_filtered)
    
  }, error = function(e) {
    flog.error("Error creating RSE: %s", e$message)
    stop(e)
  })
}

# Helper function to calculate PSI for a chunk of samples
calculate_chunk_psi <- function(rse, chunk_samples) {
  junction_counts <- assay(rse)[, chunk_samples, drop = FALSE]
  jxn_coords <- rowRanges(rse)
  
  psi_chunk <- matrix(NA, nrow = nrow(junction_counts), 
                      ncol = length(chunk_samples))
  
  for(i in seq_len(nrow(junction_counts))) {
    # Find alternative junctions
    alt_5prime <- which(start(jxn_coords) == start(jxn_coords[i]) & 
                          seq_along(jxn_coords) != i)
    alt_3prime <- which(end(jxn_coords) == end(jxn_coords[i]) & 
                          seq_along(jxn_coords) != i)
    
    # Calculate PSI
    inclusion_reads <- junction_counts[i,]
    alt_5prime_reads <- if(length(alt_5prime) > 0) 
      colSums(junction_counts[alt_5prime,, drop = FALSE]) else 0
    alt_3prime_reads <- if(length(alt_3prime) > 0) 
      colSums(junction_counts[alt_3prime,, drop = FALSE]) else 0
    
    total_reads <- inclusion_reads + alt_5prime_reads + alt_3prime_reads
    psi_chunk[i,] <- ifelse(total_reads > 0,
                            inclusion_reads / total_reads * 100,
                            NA)
  }
  
  return(psi_chunk)
}

# Optimized PSI calculation function
calculate_psi <- function(project_info, gene_coords, chunk_size = 500) {
  project_name <- project_info$project
  
  # Try to load final results first
  final_results <- load_checkpoint(paste0(project_name, "_final_psi"))
  if(!is.null(final_results)) {
    flog.info("Loaded final PSI results from checkpoint")
    return(final_results)
  }
  
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
  if(is.null(rse_jxn)) {
    return(NULL)
  }
  
  # Try to load partial results
  partial_results <- load_checkpoint(paste0(project_name, "_partial_psi"))
  if(!is.null(partial_results)) {
    flog.info("Resuming from partial results")
    psi_matrix <- partial_results$psi_matrix
    last_processed_sample <- partial_results$last_processed_sample
  } else {
    # Initialize new analysis
    junction_counts <- assay(rse_jxn)
    jxn_coords <- rowRanges(rse_jxn)
    
    # Create junction identifiers
    jxn_ids <- paste0(
      seqnames(jxn_coords), ":",
      start(jxn_coords), "-",
      end(jxn_coords)
    )
    
    psi_matrix <- matrix(NA, nrow = length(jxn_ids), ncol = ncol(junction_counts))
    rownames(psi_matrix) <- jxn_ids
    colnames(psi_matrix) <- colnames(junction_counts)
    last_processed_sample <- 0
  }
  
  # Process remaining samples
  total_samples <- ncol(psi_matrix)
  remaining_samples <- (last_processed_sample + 1):total_samples
  
  if(length(remaining_samples) > 0) {
    # Create progress bar
    pb <- progress_bar$new(
      format = "Processing samples [:bar] :percent eta: :eta",
      total = length(remaining_samples)
    )
    
    # Process in chunks
    chunks <- split(remaining_samples, 
                    ceiling(seq_along(remaining_samples)/chunk_size))
    
    for(chunk in chunks) {
      # Calculate PSI for chunk
      chunk_matrix <- calculate_chunk_psi(rse_jxn, chunk)
      psi_matrix[, chunk] <- chunk_matrix
      
      # Save checkpoint
      save_checkpoint(
        list(
          psi_matrix = psi_matrix,
          last_processed_sample = max(chunk)
        ),
        paste0(project_name, "_partial_psi")
      )
      
      # Update progress
      pb$tick(length(chunk))
    }
  }
  
  # Create final results
  results <- list(
    psi_values = psi_matrix,
    junction_info = data.frame(
      junction_id = rownames(psi_matrix),
      start = start(rowRanges(rse_jxn)),
      end = end(rowRanges(rse_jxn)),
      stringsAsFactors = FALSE
    )
  )
  
  # Save final results
  save_checkpoint(results, paste0(project_name, "_final_psi"))
  
  return(results)
}

# PSI comparison function
create_psi_comparison <- function(tcga_psi, gtex_psi) {
  flog.info("Creating PSI comparison")
  
  # Find common junctions
  common_junctions <- intersect(
    rownames(tcga_psi$psi_values),
    rownames(gtex_psi$psi_values)
  )
  
  if(length(common_junctions) == 0) {
    flog.warn("No common junctions found between datasets")
    return(NULL)
  }
  
  flog.info("Found %d common junctions", length(common_junctions))
  
  # Create progress bar
  pb <- progress_bar$new(
    format = "Calculating statistics [:bar] :percent eta: :eta",
    total = length(common_junctions)
  )
  
  # Initialize results list
  comparison_stats <- list()
  
  # Calculate statistics for each junction
  for(junction in common_junctions) {
    tryCatch({
      # Extract PSI values
      tcga_vals <- as.numeric(tcga_psi$psi_values[junction,])
      gtex_vals <- as.numeric(gtex_psi$psi_values[junction,])
      
      # Remove NA values
      tcga_vals <- tcga_vals[!is.na(tcga_vals)]
      gtex_vals <- gtex_vals[!is.na(gtex_vals)]
      
      # Only proceed if we have enough data points
      if(length(tcga_vals) >= 3 && length(gtex_vals) >= 3) {
        # Perform Wilcoxon test
        test_result <- wilcox.test(tcga_vals, gtex_vals)
        
        # Calculate means
        tcga_mean <- mean(tcga_vals, na.rm = TRUE)
        gtex_mean <- mean(gtex_vals, na.rm = TRUE)
        
        comparison_stats[[junction]] <- data.frame(
          junction_id = junction,
          tcga_psi = tcga_mean,
          gtex_psi = gtex_mean,
          delta_psi = tcga_mean - gtex_mean,
          p_value = test_result$p.value,
          tcga_samples = length(tcga_vals),
          gtex_samples = length(gtex_vals),
          stringsAsFactors = FALSE
        )
      }
    }, error = function(e) {
      flog.warn("Error processing junction %s: %s", junction, e$message)
    })
    
    pb$tick()
  }
  
  # Combine all results into a single data frame
  comparison_df <- do.call(rbind, comparison_stats)
  
  # Add additional statistics
  comparison_df <- comparison_df %>%
    mutate(
      significance = -log10(p_value),
      is_significant = p_value < 0.05 & abs(delta_psi) > 10,
      direction = case_when(
        delta_psi > 0 ~ "Higher in TCGA",
        delta_psi < 0 ~ "Higher in GTEx",
        TRUE ~ "No change"
      )
    )
  
  # Add junction coordinates
  comparison_df <- comparison_df %>%
    left_join(
      tcga_psi$junction_info,
      by = c("junction_id")
    )
  
  flog.info("Comparison completed. Found %d significant changes", 
            sum(comparison_df$is_significant))
  
  return(comparison_df)
}

# Visualization function
create_visualization <- function(comparison_data) {
  flog.info("Creating visualization")
  
  p <- ggplot(comparison_data, 
              aes(x = gtex_psi, y = tcga_psi, color = significance)) +
    geom_point(alpha = 0.7, size = 3) +
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
      subtitle = "TCGA GBM vs GTEx Brain",
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
  
  # Save plot
  pdf("srrm3_psi_comparison.pdf", width = 12, height = 10)
  print(p)
  dev.off()
  
  return(p)
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
  
  if(!is.null(tcga_psi) && !is.null(gtex_psi)) {
    # Create comparison
    comparison_data <- create_psi_comparison(tcga_psi, gtex_psi)
    
    if(!is.null(comparison_data)) {
      # Save comparison results
      write.csv(comparison_data, "srrm3_psi_summary.csv", row.names = FALSE)
      
      # Create and save visualization
      # Create and save visualization
      create_visualization(comparison_data)
      
      flog.info("Analysis completed successfully")
    } else {
      flog.warn("No valid comparison could be generated")
    }
  } else {
    flog.warn("Missing data for comparison")
  }
}

# Run the analysis with error handling
tryCatch({
  main()
}, error = function(e) {
  flog.error("Analysis failed: %s", e$message)
  stop(e)
})