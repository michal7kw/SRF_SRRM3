

library(recount3)
library(SummarizedExperiment)
library(GenomicRanges)

# Function to create coverage plot for SRRM3 gene
create_srrm3_coverage_plot <- function() {
  # Define SRRM3 coordinates
  chr <- "chr7"
  start_pos <- 76282502
  end_pos <- 76287486
  
  # Create GRanges object
  gr <- GRanges(chr, IRanges(start=start_pos, end=end_pos))
  
  # Get coverage data using coverage_matrix function
  print('Downloading TCGA data...')
  SRRM3_tcga <- coverage_matrix('TCGA', chr, gr)
  
  print('Downloading GTEx data...')
  SRRM3_gtex <- coverage_matrix('SRP012682', chr, gr)
  
  # Normalize the coverage data
  normalize_matrix <- function(mat) {
    # Apply normalization to each sample (column)
    normalized <- apply(mat, 2, function(x) {
      max_val <- max(x, na.rm = TRUE)
      if(max_val > 0) x/max_val else x
    })
    return(normalized)
  }
  
  tcga_normalized <- normalize_matrix(SRRM3_tcga)
  gtex_normalized <- normalize_matrix(SRRM3_gtex)
  
  # Create plot
  positions <- start_pos:end_pos
  
  # Set up the plotting area
  plot(NULL, 
       xlim = c(start_pos, end_pos),
       ylim = c(0, 1),
       xlab = "Genomic Position on Chr7",
       ylab = "Normalized Coverage",
       main = "SRRM3 Coverage Plot")
  
  # Add TCGA samples in blue
  for(i in 1:ncol(tcga_normalized)) {
    lines(positions, 
          tcga_normalized[,i], 
          col = adjustcolor("blue", alpha = 0.3))
  }
  
  # Add GTEx samples in red
  for(i in 1:ncol(gtex_normalized)) {
    lines(positions, 
          gtex_normalized[,i], 
          col = adjustcolor("red", alpha = 0.3))
  }
  
  # Add legend
  legend("topright", 
         legend = c("TCGA samples", "GTEx samples"),
         col = c("blue", "red"),
         lty = 1,
         bty = "n")
  
  # Print sample counts
  print(paste("Number of TCGA samples:", ncol(tcga_normalized)))
  print(paste("Number of GTEx samples:", ncol(gtex_normalized)))
}

# Usage:
create_srrm3_coverage_plot()