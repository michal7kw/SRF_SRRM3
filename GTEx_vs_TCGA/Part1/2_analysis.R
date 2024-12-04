library(recount3)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(SummarizedExperiment)

# First get available projects
message("Getting available projects...")
projects <- available_projects()

# Function to analyze exon structure
analyze_srrm3_structure <- function(project_info) {
  message("\nAnalyzing SRRM3 structure for project: ", project_info$project)
  
  # Get exon-level data
  rse_exon <- create_rse(project_info, type = "exon")
  
  # Define SRRM3 region
  region <- GRanges(
    seqnames = "chr7",
    ranges = IRanges(start = 76201900, end = 76287288)
  )
  
  # Get overlapping exons
  overlaps <- findOverlaps(rowRanges(rse_exon), region)
  
  if(length(overlaps) == 0) {
    message("No exons found in SRRM3 region")
    return(NULL)
  }
  
  # Extract exon data
  exon_data <- rse_exon[queryHits(overlaps),]
  exon_ranges <- rowRanges(exon_data)
  
  # Create a data frame with exon information
  exon_df <- data.frame(
    exon_id = names(exon_ranges),
    seqnames = as.character(seqnames(exon_ranges)),
    start = start(exon_ranges),
    end = end(exon_ranges),
    width = width(exon_ranges),
    stringsAsFactors = FALSE
  )
  
  # Add metadata columns
  meta_columns <- c("transcript_id", "transcript_name", "exon_number", 
                    "transcript_type", "gene_name")
  
  for(col in meta_columns) {
    if(col %in% names(mcols(exon_ranges))) {
      exon_df[[col]] <- mcols(exon_ranges)[[col]]
    }
  }
  
  # Sort by start position and transcript
  exon_df <- exon_df[order(exon_df$start, exon_df$transcript_id),]
  
  # Print available metadata columns
  message("\nAvailable metadata columns:")
  print(names(mcols(exon_ranges)))
  
  # Print first few rows of exon data
  message("\nFirst few exons:")
  print(head(exon_df))
  
  # Get expression data for the exons
  exon_expr <- assay(exon_data)
  
  # Calculate mean expression per exon
  mean_expr <- rowMeans(exon_expr)
  exon_df$mean_expression <- mean_expr
  
  message("\nExpression summary per exon:")
  print(summary(exon_df$mean_expression))
  
  return(list(
    exon_details = exon_df,
    expression_data = exon_expr,
    raw_ranges = exon_ranges
  ))
}

# Print available projects
message("\nAvailable TCGA projects:")
print(subset(projects, file_source == "tcga")$project)

message("\nAvailable GTEx projects:")
print(subset(projects, file_source == "gtex")$project)

# Get TCGA brain data (GBM)
message("\n=== Analyzing TCGA GBM SRRM3 structure ===")
tcga_project_info <- subset(projects, 
                            project == "GBM" & 
                              file_source == "tcga")
tcga_structure <- analyze_srrm3_structure(tcga_project_info)

# Get GTEx brain data
message("\n=== Analyzing GTEx brain SRRM3 structure ===")
gtex_project_info <- subset(projects, 
                            project == "BRAIN" & 
                              file_source == "gtex")
gtex_structure <- analyze_srrm3_structure(gtex_project_info)

# Print detailed information about the transcripts
if(!is.null(tcga_structure)) {
  message("\nTCGA SRRM3 transcript details:")
  ranges <- tcga_structure$raw_ranges
  if("transcript_id" %in% names(mcols(ranges))) {
    unique_transcripts <- unique(mcols(ranges)$transcript_id)
    message("Found ", length(unique_transcripts), " unique transcripts:")
    print(unique_transcripts)
  } else {
    message("Available metadata:")
    print(names(mcols(ranges)))
  }
}

if(!is.null(gtex_structure)) {
  message("\nGTEx SRRM3 transcript details:")
  ranges <- gtex_structure$raw_ranges
  if("transcript_id" %in% names(mcols(ranges))) {
    unique_transcripts <- unique(mcols(ranges)$transcript_id)
    message("Found ", length(unique_transcripts), " unique transcripts:")
    print(unique_transcripts)
  } else {
    message("Available metadata:")
    print(names(mcols(ranges)))
  }
}

# Save the detailed structure information
save(tcga_structure, gtex_structure, 
     file = "2_srrm3_structure_analysis.RData")