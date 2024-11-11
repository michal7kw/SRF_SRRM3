# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("recount", quietly = TRUE))
    BiocManager::install("recount")

library(recount)
library(GenomicRanges)

setwd("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3/GTEx")

# Create a data directory if it doesn't exist
if (!dir.exists("DATA")) {
  dir.create("DATA")
}

# Define the genomic region of interest for SRRM3
srrm3_region <- GRanges(
  seqnames = "chr7",
  ranges = IRanges(start = 76282502, end = 76287486)
)

# Create a progress monitoring function
progress_monitor <- function(x) {
  counter <<- counter + 1
  if (counter %% 10 == 0) {
    cat(sprintf('\rProgress: %d/%d files (%.1f%%)', 
                counter, total_files, 
                counter/total_files*100))
  }
}

# Function to safely save progress
save_progress <- function(gtex_data = NULL, prefix = "CHECKPOINT") {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  filename <- sprintf("DATA/%s_SRRM3_coverage_%s.rda", prefix, timestamp)
  save(gtex_data, file = filename)
  cat(sprintf("\nProgress saved to: %s\n", filename))
}

# Set up error handling and interruption
tryCatch({
  # Download GTEx data
  print('Downloading GTEx data...')
  files_info <- download_study('SRP012682')
  total_files <- nrow(files_info)
  counter <- 0
  
  SRRM3_gtex <- coverage_matrix(
    project = 'SRP012682',
    chr = 'chr7',
    region = srrm3_region,
    verbose = TRUE,
    parallel = TRUE,
    cores = parallel::detectCores() - 1,
    downloadCallback = progress_monitor
  )
  
  # Save final results
  save(SRRM3_gtex, file = 'DATA/SRRM3_gtex_coverage.rda')
  
}, error = function(e) {
  cat("\nError occurred:", conditionMessage(e), "\n")
  save_progress(gtex_data = if(exists("SRRM3_gtex")) SRRM3_gtex else NULL,
               prefix = "ERROR_GTEX")
}, interrupt = function(i) {
  cat("\nInterrupted by user. Saving progress...\n")
  save_progress(gtex_data = if(exists("SRRM3_gtex")) SRRM3_gtex else NULL,
               prefix = "INTERRUPT_GTEX")
})

# Print dimensions of downloaded data to verify
print("Dimensions of GTEx data:")
print(dim(SRRM3_gtex)) 