# Load required packages
library(dplyr)
library(SummarizedExperiment)
library(ggplot2)

inspect_rse <- function(rse, max_rows = 5) {
  # Main object information
  cat("=== RangedSummarizedExperiment Object Summary ===\n\n")
  
  # Basic dimensions
  cat(sprintf("Dimensions: %d rows x %d columns\n", nrow(rse), ncol(rse)))
  cat(sprintf("Class: %s\n\n", class(rse)[1]))
  
  # Metadata
  cat("=== Metadata ===\n")
  meta_names <- names(metadata(rse))
  if (length(meta_names) > 0) {
    cat(sprintf("Available metadata (%d): %s\n\n", 
                length(meta_names), 
                paste(meta_names, collapse = ", ")))
  } else {
    cat("No metadata available\n\n")
  }
  
  # Assays
  cat("=== Assays ===\n")
  assay_names <- assayNames(rse)
  cat(sprintf("Available assays (%d): %s\n\n", 
              length(assay_names), 
              paste(assay_names, collapse = ", ")))
  
  # Row Data
  cat("=== Row Data ===\n")
  rd_names <- names(rowData(rse))
  if (length(rd_names) > 0) {
    cat(sprintf("Row annotations (%d): %s\n", 
                length(rd_names), 
                paste(rd_names, collapse = ", ")))
    
    # Sample of first few rows
    cat("\nFirst few rows of rowData:\n")
    print(head(as.data.frame(rowData(rse)), max_rows))
    cat("\n")
  }
  
  # Column Data
  cat("=== Column Data ===\n")
  cd_names <- names(colData(rse))
  if (length(cd_names) > 0) {
    cat(sprintf("Column annotations (%d): %s\n", 
                length(cd_names), 
                paste(cd_names, collapse = ", ")))
    
    # Sample of first few columns
    cat("\nFirst few rows of colData:\n")
    print(head(as.data.frame(colData(rse)), max_rows))
  }
  
  # Assay data preview
  if (length(assay_names) > 0) {
    cat("\n=== Assay Data Preview ===\n")
    for (assay_name in assay_names) {
      cat(sprintf("\nFirst %d x %d values of %s:\n", max_rows, max_rows, assay_name))
      print(assay(rse, assay_name)[1:max_rows, 1:max_rows])
    }
  }
}

# Function to plot exon coverage distribution
plot_exon_coverage <- function(rse, target_exon_id) {  
  # Get library sizes (total reads per sample)
  lib_sizes <- colSums(assay(rse))
  
  # Extract coverage values and normalize to RPM
  coverage_data <- data.frame(
    coverage = as.numeric(assay(rse)[target_exon_id, ]) / lib_sizes * 1e6
  )
  
  # Create histogram to show the frequency distribution of coverage values
  ggplot(coverage_data, aes(x = coverage)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    theme_minimal() +
    labs(
      title = "Target Exon Coverage Distribution",
      x = "Coverage (RPM)",
      y = "Number of Samples"
    )
}

get_target_exon_id <- function(rse, target_info) {
  # Get the row ranges as a data frame
  exon_coords <- rowRanges(rse)
  
  # Find the exact match for our target exon
  target_idx <- which(
    seqnames(exon_coords) == target_info$gene$chr &
      start(exon_coords) == target_info$exon24$start &
      end(exon_coords) == target_info$exon24$end
  )
  
  if(length(target_idx) == 0) {
    message("Target exon not found!")
    return(NULL)
  }
  
  if(length(target_idx) > 1) {
    message("Multiple matches found, using the first one")
  }
  
  # Return the row name (exon ID) for the target exon
  target_exon_id <- rownames(rse)[target_idx[1]]
  
  # Print the found coordinates for verification
  message(sprintf("Found target exon: %s at %s:%d-%d", 
                  target_exon_id,
                  target_info$gene$chr,
                  target_info$exon24$start,
                  target_info$exon24$end))
  
  return(target_exon_id)
}