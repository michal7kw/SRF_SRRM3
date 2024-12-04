# Load required libraries
library(recount3)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(SummarizedExperiment)
source("helpers.R")

# Define SRRM3 information with genomic coordinates for exon 24
SRRM3_INFO <- list(
  gene = list(
    name = "SRRM3",
    chr = "chr7",
    start = 76201896,  # GRCh38 coordinates
    end = 76287287
  ),
  exon24 = list(
    start = 76283524,  # Genomic coordinates from GTF
    end = 76283602,
    length = 79
  )
)

# Get projects
projects <- available_projects()
projects

################################################################################
######################## Get TCGA GBM data ##################################### 
################################################################################
tcga_project <- subset(projects, project == "GBM" & file_source == "tcga")
tcga_project
tcga_project$project

tcga_exon <- create_rse(tcga_project, type = "exon")
tcga_exon
tcga_exon@metadata

# Get the row ranges as a data frame
tcga_exon_coords <- rowRanges(tcga_exon)
tcga_exon_coords[1099686:1099691]

########## Get the target exon ID
tcga_target_exon_id <- get_target_exon_id(tcga_exon, SRRM3_INFO)
tcga_target_exon_id

########## Now you can use this ID to plot the coverage
if(!is.null(tcga_target_exon_id)) {
  plot_exon_coverage(tcga_exon, tcga_target_exon_id)
}

########## get_target_exon_id ##################################################
# Find the exact match for our target exon
tcga_target_idx <- which(
  seqnames(tcga_exon_coords) == SRRM3_INFO$gene$chr &
    start(tcga_exon_coords) == SRRM3_INFO$exon24$start &
    end(tcga_exon_coords) == SRRM3_INFO$exon24$end
)
tcga_target_idx
rownames(tcga_exon)[tcga_target_idx[1]]
################################################################################

########## plot_exon_coverage ##################################################
# To check raw counts
tcga_raw_counts <- as.numeric(assay(tcga_exon)[tcga_target_idx, ])
head(tcga_raw_counts)

# To check RPM values
lib_sizes <- colSums(assay(tcga_exon))
rpm_values <- tcga_raw_counts / lib_sizes * 1e6
head(rpm_values)

coverage_data <- data.frame(
  coverage = as.numeric(assay(tcga_exon)[tcga_target_idx, ])
)
################################################################################






################################################################################
######################## Get GTEx brain data ###################################
################################################################################
gtex_project <- subset(projects, project == "BRAIN" & file_source == "gtex")
gtex_project
gtex_project$project

gtex_exon <- create_rse(gtex_project, type = "exon")
gtex_exon
gtex_exon@metadata

# Get the row ranges as a data frame
gtex_exon_coords <- rowRanges(gtex_exon)
gtex_exon_coords[1099686:1099691]

########## Get the target exon ID
target_exon_id <- get_target_exon_id(gtex_exon, SRRM3_INFO)
target_exon_id

########## use this ID to plot the coverage
if(!is.null(target_exon_id)) {
  plot_exon_coverage(gtex_exon, target_exon_id)
}

##########  Extract coverage data
gtex_exon_counts <- assay(gtex_exon)
gtex_exon_counts[1099686:1099691, 1:3]

########## get_target_exon_id ##################################################
# Find the exact match for our target exon
gtex_target_idx <- which(
  seqnames(gtex_exon_coords) == SRRM3_INFO$gene$chr &
    start(gtex_exon_coords) == SRRM3_INFO$exon24$start &
    end(gtex_exon_coords) == SRRM3_INFO$exon24$end
)
gtex_target_idx
rownames(gtex_exon)[gtex_target_idx[1]]
################################################################################

########## Calculate normalized counts (RPM)
gtex_totals <- colSums(gtex_exon_counts)

plot(density(gtex_totals),
     main="Density Plot of Column Sums",
     xlab="Total Counts")

gtex_normalized_counts <- t(t(gtex_exon_counts) / gtex_totals) * 1e6
gtex_normalized_counts[1099686:1099691, 1:3]

gtex_exon_counts <- assay(gtex_exon)
gtex_exon_counts[1:3, 1:3]



# Save GTEx exon data
saveRDS(gtex_exon, "gtex_exon.rds")
gtex_exon <- readRDS("gtex_exon.rds")

saveRDS(tcga_exon, "tcga_exon.rds")
tcga_exon <- readRDS("tcga_exon.rds")
