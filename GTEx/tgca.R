# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("recount", quietly = TRUE))
  BiocManager::install("recount")
if (!requireNamespace("GenomicRanges", quietly = TRUE))
  BiocManager::install("GenomicRanges")

library(recount)
library(GenomicRanges)

# Function to safely download TCGA coverage data
get_tcga_coverage <- function(chromosome, start_pos, end_pos) {
  tryCatch({
    # First download the resource that we want to use
    url <- download_study('TCGA')
    
    # Load the RangedSummarizedExperiment object
    load(url)
    
    # Create GRanges object for the region of interest
    region <- GRanges(
      seqnames = chromosome,
      ranges = IRanges(start = start_pos, end = end_pos)
    )
    
    # Get coverage matrix for the specified region
    coverage_data <- coverage_matrix(
      project = 'TCGA',
      chr = chromosome,
      regions = region
    )
    
    return(coverage_data)
    
  }, error = function(e) {
    message("Error occurred while downloading data: ", e$message)
    message("Please check:\n",
            "1. Internet connection\n",
            "2. Valid chromosome name (e.g., 'chr7')\n",
            "3. Valid genomic coordinates\n",
            "4. Sufficient disk space\n",
            "5. Latest version of recount package")
    return(NULL)
  })
}

# Example usage for SRRM3 region
srrm3_coverage <- get_tcga_coverage(
  chromosome = "chr7",
  start_pos = 76282502,
  end_pos = 76287486
)

# If successful, you can work with the coverage data
if (!is.null(srrm3_coverage)) {
  # Basic summary statistics
  print(summary(srrm3_coverage))
  
  # Number of samples
  print(paste("Number of samples:", ncol(srrm3_coverage)))
  
  # Save the data if needed
  save(srrm3_coverage, file = "srrm3_tcga_coverage.RData")
}