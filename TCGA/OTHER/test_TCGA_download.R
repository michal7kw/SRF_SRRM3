# Test script for TCGA data download
library(TCGAbiolinks)
library(SummarizedExperiment)
library(futile.logger)
library(dplyr)

# Set up logging
if (!dir.exists("logs")) dir.create("logs")
flog.appender(appender.file("logs/test_download.log"))
flog.threshold(DEBUG)

# Create cache directory if it doesn't exist
CACHE_DIR <- "cache"
if (!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR)

test_TCGA_download <- function(cancer_type = "BRCA", gene = "SRRM3") {
  flog.info("Testing TCGA download for %s, gene %s", cancer_type, gene)
  
  # Cache file path
  cache_file <- file.path(CACHE_DIR, paste0(cancer_type, "_", gene, "_test.rds"))
  
  # Check cache first
  if (file.exists(cache_file)) {
    flog.info("Using cached test data")
    return(readRDS(cache_file))
  }
  
  tryCatch({
    # Create query
    query <- GDCquery(
      project = paste0("TCGA-", cancer_type),
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      sample.type = c("Primary Tumor")
    )
    
    # Get first sample only for testing
    samples <- getResults(query)
    if (nrow(samples) == 0) {
      flog.error("No samples found for %s", cancer_type)
      return(FALSE)
    }
    
    test_query <- GDCquery(
      project = paste0("TCGA-", cancer_type),
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      sample.type = c("Primary Tumor"),
      barcode = samples$cases[1]
    )
    
    flog.info("Downloading test data...")
    GDCdownload(test_query)
    
    flog.info("Preparing test data...")
    test_data <- GDCprepare(test_query)
    
    # Check if gene exists
    gene_idx <- which(rowData(test_data)$gene_name == gene)
    if (length(gene_idx) == 0) {
      flog.error("Gene %s not found in dataset", gene)
      return(FALSE)
    }
    
    # Get expression data for the gene
    expr_data <- assay(test_data)[gene_idx, ]
    
    # Save to cache
    saveRDS(list(success = TRUE, data = expr_data), cache_file)
    
    flog.info("Test successful!")
    return(TRUE)
    
  }, error = function(e) {
    flog.error("Test failed: %s", e$message)
    return(FALSE)
  })
}

# Run test
result <- test_TCGA_download("BRCA", "SRRM3")
if (result) {
  cat("Test completed successfully. Check logs/test_download.log for details.\n")
} else {
  cat("Test failed. Check logs/test_download.log for details.\n")
} 