#####################################################################
# Multi-Cancer Survival Analysis Script
#####################################################################
# This script performs survival analysis across multiple cancer types
# using both PSI (Percent Spliced In) values and gene expression data.
# It analyzes the relationship between SRRM3/SRRM4 molecular features
# and patient survival outcomes using TCGA data.

#####################################################################
# Load Required Libraries
#####################################################################
suppressPackageStartupMessages({
  # Data manipulation and analysis
  library(dplyr)        # Data manipulation and transformation
  library(tidyverse)    # Collection of data science packages
  library(data.table)   # Fast data manipulation
  
  # Survival analysis
  library(survival)     # Core survival analysis functions
  library(survminer)    # Survival analysis visualization
  
  # Genomic data handling
  library(recount3)     # Access to RNA-seq data
  library(biomaRt)      # Access to genomic annotations
  library(TCGAbiolinks) # TCGA data access
  library(SummarizedExperiment) # Container for genomic data
  library(DESeq2)      # RNA-seq analysis
  library(GenomicFeatures) # Genomic feature handling
  library(rtracklayer)  # Import/export genomic tracks
  
  # Matrix operations
  library(matrixStats)  # Matrix calculations
  library(sparseMatrixStats) # Sparse matrix operations
  
  # Parallel processing
  library(parallel)     # Base R parallel processing
  library(BiocParallel) # Bioconductor parallel processing
  
  # Utilities
  library(httr)         # HTTP requests
  library(retry)        # Retry failed operations
  library(futile.logger) # Logging functionality
  library(viridis)      # Color palettes
})

#####################################################################
# Setup and Configuration
#####################################################################
# Get working directory from environment or use current
work_dir <- getwd()

# Define paths
RESULTS_DIR <- file.path(work_dir, "results")
CACHE_DIR <- file.path(work_dir, "cache")
LOGS_DIR <- file.path(work_dir, "logs")
CLINICAL_DIR <- file.path(work_dir, "../../TCGA_clinical")

# Create necessary directories
for (dir in c(RESULTS_DIR, CACHE_DIR, LOGS_DIR)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    message(sprintf("Created directory: %s", dir))
  }
}

#####################################################################
# Set Up Parallel Processing
#####################################################################
# Detect number of cores from SLURM environment or default to 1
num_cores <- as.numeric(Sys.getenv("SLURM_NTASKS", unset = "1"))
if (num_cores > 1) {
  message(sprintf("Setting up parallel processing with %d cores", num_cores))
  BiocParallel::register(MulticoreParam(workers = num_cores))
} else {
  message("Running in single-core mode")
}

#####################################################################
# Ensembl Connection Management
#####################################################################
# Cache connection to Ensembl to avoid repeated connections
ensembl <- NULL
get_ensembl_connection <- function() {
  if (is.null(ensembl)) {
    ensembl <<- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  }
  return(ensembl)
}

#####################################################################
# Clinical Data Retrieval and Processing
#####################################################################
# Function to retrieve and clean clinical data from TCGA
get_clinical_data <- function(cancer_type) {
  message(sprintf("Getting clinical data for %s", cancer_type))
  
  # Construct path to clinical data directory using absolute path
  clinical_dir <- file.path(CLINICAL_DIR,
                           sprintf("clinical.project-tcga-%s.2025-02-12", 
                                 tolower(cancer_type)))
  message(sprintf("Looking for clinical data in directory: %s", clinical_dir))
  
  # Check if directory exists
  if (!dir.exists(clinical_dir)) {
    stop(sprintf("Clinical data directory not found for %s: %s", cancer_type, clinical_dir))
  }
  
  # List files in directory for debugging
  message("Files in clinical directory:")
  message(paste(list.files(clinical_dir), collapse="\n"))
  
  # Read clinical.tsv and follow_up.tsv files
  clinical <- read.delim(file.path(clinical_dir, "clinical.tsv"), 
                        sep="\t", stringsAsFactors=FALSE, na.strings="'--")
  follow_up <- read.delim(file.path(clinical_dir, "follow_up.tsv"), 
                         sep="\t", stringsAsFactors=FALSE, na.strings="'--")
  
  if (nrow(clinical) == 0) {
    stop(sprintf("No clinical data found for %s", cancer_type))
  }
  
  # Get the latest follow-up data for each patient
  latest_follow_up <- follow_up %>%
    group_by(case_submitter_id) %>%
    summarize(
      days_to_follow_up = max(days_to_follow_up, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(!is.infinite(days_to_follow_up))  # Remove cases where all values were NA
  
  # Clean and prepare clinical data
  clinical_clean <- clinical %>%
    left_join(latest_follow_up, by = "case_submitter_id") %>%
    dplyr::transmute(
      case_id = case_submitter_id,
      vital_status = vital_status,
      # Calculate survival time using either death or last follow-up
      overall_survival = case_when(
        !is.na(days_to_death) ~ as.numeric(days_to_death),
        !is.na(days_to_follow_up) ~ as.numeric(days_to_follow_up),
        TRUE ~ NA_real_
      ),
      # Standardize vital status to boolean
      deceased = case_when(
        tolower(vital_status) == "dead" ~ TRUE,
        tolower(vital_status) == "alive" ~ FALSE,
        TRUE ~ NA
      )
    ) %>%
    dplyr::filter(!is.na(overall_survival), !is.na(deceased))
  
  # Print summary statistics
  message(sprintf("Processed clinical data summary:"))
  message(sprintf("- Total patients: %d", nrow(clinical_clean)))
  message(sprintf("- Patients with death events: %d", sum(clinical_clean$deceased)))
  message(sprintf("- Median survival time: %.1f days", 
                 median(clinical_clean$overall_survival)))
  
  return(clinical_clean)
}

#####################################################################
# Expression Data Retrieval and Processing
#####################################################################
# Function to get gene expression data from TCGA via recount3
get_expression_data <- function(cancer_type, gene_name) {
  cache_file <- file.path(CACHE_DIR, 
                         paste0("expression_data_", cancer_type, "_", gene_name, ".rds"))
  
  # Check cache first
  if (file.exists(cache_file)) {
    message("Loading cached expression data...")
    return(readRDS(cache_file))
  }
  
  message("Getting expression data from recount3...")
  projects <- available_projects()
  
  # Find specific project in recount3

  project_info <- subset(projects, 
                        project == toupper(cancer_type) &
                        file_source == "tcga" &
                        project_type == "data_sources")
  
  if (nrow(project_info) != 1) {
    stop(sprintf("Could not find unique project for %s", cancer_type))
  }
  
  message("Creating RSE object...")
  rse_gene <- create_rse(
    project_info[1, ],
    type = "gene",
    annotation = "gencode_v26"
  )
  
  # Extract expression for specific gene
  gene_idx <- which(rowData(rse_gene)$gene_name == gene_name)[1]
  if (is.na(gene_idx)) {
    stop(sprintf("Gene %s not found in dataset", gene_name))
  }
  
  # Create data frame with expression values
  expression_data <- data.frame(
    case_id = substr(colData(rse_gene)$tcga.tcga_barcode, 1, 12),
    expression = assays(rse_gene)$counts[gene_idx, ]
  )
  
  # Handle duplicate samples by averaging
  if (any(duplicated(expression_data$case_id))) {
    expression_data <- expression_data %>%
      group_by(case_id) %>%
      summarize(expression = mean(expression, na.rm = TRUE)) %>%
      ungroup()
  }
  
  saveRDS(expression_data, cache_file)
  return(expression_data)
}

#####################################################################
# PSI (Percent Spliced In) Data Retrieval and Processing
#####################################################################
# Function to calculate PSI values for specific exons
get_psi_data <- function(cancer_type, gene_name, sample_ids = NULL) {
  cache_file <- file.path(CACHE_DIR, 
                         paste0("psi_data_", cancer_type, "_", gene_name, ".rds"))
  
  if (file.exists(cache_file)) {
    message("Loading cached PSI data...")
    return(readRDS(cache_file))
  }
  
    # Define genomic coordinates for SRRM3 exon 15
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
    )
  )
  
  message("Getting junction data from recount3...")
  projects <- available_projects()
  
  # Find specific project
  project_info <- subset(projects, 
                        project == toupper(cancer_type) &
                        file_source == "tcga" &
                        project_type == "data_sources")
  
  if (nrow(project_info) != 1) {
    message("\nAvailable TCGA projects:")
    print(subset(projects, file_source == "tcga" & project_type == "data_sources"))
    stop("Could not find unique project for ", cancer_type)
  }

    # Get junction data
  message("Creating RSE object...")
  rse_jxn <- create_rse(
    project_info[1, ],
    type = "jxn",
    jxn_format = "UNIQUE",
    verbose = TRUE
  )
  
  junction_counts <- assay(rse_jxn)
  jxn_coords <- rowRanges(rse_jxn)
  
  # Find junctions in SRRM3 region
  message("Finding relevant junctions...")
  srrm3_region <- which(
    start(jxn_coords) >= (SRRM3_INFO$exon15$start - 10000) &
    end(jxn_coords) <= (SRRM3_INFO$exon15$end + 10000) &
    seqnames(jxn_coords) == "chr7"
  )
  
  if (length(srrm3_region) == 0) {
    stop("No junctions found in SRRM3 region")
  }
  
  message(sprintf("Found %d junctions in SRRM3 region", length(srrm3_region)))
  
  # Extract relevant junctions
  jxn_coords <- jxn_coords[srrm3_region]
  junction_counts <- junction_counts[srrm3_region, ]
  
  # Identify inclusion and exclusion junctions
  inclusion_jxns <- which(
    (abs(end(jxn_coords) - SRRM3_INFO$exon15$start) <= 5) |
    (abs(start(jxn_coords) - SRRM3_INFO$exon15$end) <= 5)
  )
  
  exclusion_jxns <- which(
    start(jxn_coords) < (SRRM3_INFO$exon15$start - 5) &
    end(jxn_coords) > (SRRM3_INFO$exon15$end + 5)
  )
  
  message(sprintf("Found %d inclusion and %d exclusion junctions", 
                 length(inclusion_jxns), length(exclusion_jxns)))
  
  if (length(inclusion_jxns) == 0 || length(exclusion_jxns) == 0) {
    stop("Could not find both inclusion and exclusion junctions")
  }
  
  # Calculate PSI values for each sample
  psi_data <- data.frame(
    case_id = colData(rse_jxn)$tcga.tcga_barcode,
    psi = sapply(seq_len(ncol(junction_counts)), function(i) {
      inclusion_reads <- sum(junction_counts[inclusion_jxns, i])
      exclusion_reads <- sum(junction_counts[exclusion_jxns, i])
      total_reads <- inclusion_reads + exclusion_reads
      
      if(total_reads >= 10) {  # Minimum read coverage threshold
        return((inclusion_reads / total_reads) * 100)
      } else {
        return(NA)
      }
    })
  )
  
  # Clean and format data
  psi_data$case_id <- substr(psi_data$case_id, 1, 12)
  
  # Remove NA values and handle duplicates
  psi_data <- psi_data %>%
    filter(!is.na(psi)) %>%
    group_by(case_id) %>%
    summarize(psi = mean(psi, na.rm = TRUE)) %>%
    ungroup()
  
  message(sprintf("Found %d samples with valid PSI values", nrow(psi_data)))
  
  # If a subset of sample IDs is provided, filter the PSI data accordingly
  if (!is.null(sample_ids)) {
    message(sprintf("Filtering PSI data to %d high-expression samples", length(sample_ids)))
    psi_data <- psi_data %>% filter(case_id %in% sample_ids)
  }
  
  saveRDS(psi_data, cache_file)
  return(psi_data)
}

#####################################################################
# Gene Information Retrieval
#####################################################################
# Function to get gene coordinates and information from Ensembl
get_gene_info <- function(gene_name) {
  # Connect to Ensembl
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Get gene coordinates and information
  gene_info <- getBM(
    attributes = c(
      "hgnc_symbol",
      "chromosome_name",
      "start_position",
      "end_position",
      "strand"
    ),
    filters = "hgnc_symbol",
    values = gene_name,
    mart = ensembl
  )
  
  if (nrow(gene_info) == 0) {
    stop(paste("Could not find gene info for", gene_name))
  }
  
  # Handle multiple entries
  if (nrow(gene_info) > 1) {
    message("Multiple entries found for ", gene_name, ", using first entry")
    gene_info <- gene_info[1, ]
  }
  
  # Standardize chromosome format
  if (!grepl("^chr", gene_info$chromosome_name)) {
    gene_info$chromosome_name <- paste0("chr", gene_info$chromosome_name)
  }
  
  return(gene_info)
}

#####################################################################
# Multi-Cancer Analysis
#####################################################################
# Function to perform analysis across multiple cancer types
perform_multi_cancer_analysis <- function(cancer_types, analysis_type = "PSI", gene = "SRRM3", grouping_method = "quartile") {
  results_list <- list()
  
  for(cancer_type in cancer_types) {
    tryCatch({
      message(sprintf("\nProcessing %s", cancer_type))
      results <- perform_survival_analysis(cancer_type, analysis_type, gene, grouping_method)
      results_list[[cancer_type]] <- results
    }, error = function(e) {
      message(sprintf("Error analyzing %s: %s", cancer_type, e$message))
    })
  }
  
  return(results_list)
}

#####################################################################
# Survival Analysis
#####################################################################
# Function to perform survival analysis for a single cancer type
perform_survival_analysis <- function(cancer_type, data_type = "PSI", gene = "SRRM3", grouping_method = "quartile") {
  # Get clinical data
  clinical_data <- get_clinical_data(cancer_type)
  message(sprintf("Clinical data: %d samples", nrow(clinical_data)))
  
  # Get molecular data (PSI or expression)
  tryCatch({
    molecular_data <- if (data_type == "PSI") {
      get_psi_data(cancer_type, gene)
    } else {
      get_expression_data(cancer_type, gene)
    }
    message(sprintf("Molecular data: %d samples", nrow(molecular_data)))
    
    if (nrow(molecular_data) == 0) {
      stop("No samples in molecular data")
    }
    
    # Ensure IDs are in the same format before merging
    clinical_data$merge_id <- substr(clinical_data$case_id, 1, 12)
    molecular_data$merge_id <- substr(molecular_data$case_id, 1, 12)
    
    # Print sample IDs for debugging
    message("First few clinical IDs: ", paste(head(clinical_data$merge_id), collapse=", "))
    message("First few molecular IDs: ", paste(head(molecular_data$merge_id), collapse=", "))
    
    # Before merging clinical data, add this code:
    clinical_data <- clinical_data %>%
      filter(merge_id %in% molecular_data$merge_id)

    molecular_data <- molecular_data %>%
      filter(merge_id %in% clinical_data$merge_id)

    # Now both datasets should have matching IDs before merging
    
    # Merge the data
    merged_data <- inner_join(clinical_data, molecular_data, by = "merge_id")
    message(sprintf("Merged data samples: %d", nrow(merged_data)))
    
    if (nrow(merged_data) == 0) {
      stop("No samples after merging clinical and molecular data")
    }
    
    # Group samples based on molecular values
    value_col <- if(data_type == "PSI") "psi" else "expression"
    
    # Perform grouping based on specified method
    if (grouping_method == "quartile") {
      quartiles <- quantile(merged_data[[value_col]], probs = c(0.25, 0.75), na.rm = TRUE)
      if (quartiles[1] != quartiles[2]) {  # Check if quartiles are different
        merged_data$group <- case_when(
          merged_data[[value_col]] <= quartiles[1] ~ "Low",
          merged_data[[value_col]] >= quartiles[2] ~ "High",
          TRUE ~ "Medium"
        )
      } else {
        # If quartiles are the same, use median split
        median_val <- median(merged_data[[value_col]], na.rm = TRUE)
        merged_data$group <- if_else(merged_data[[value_col]] > median_val, "High", "Low")
      }
    }
    
    # Remove any NA groups
    merged_data <- merged_data %>% filter(!is.na(group))
    
    if (nrow(merged_data) == 0) {
      stop("No samples after grouping")
    }
    
    # Print group sizes
    message("Sample sizes per group:")
    print(table(merged_data$group))
    
    # Fit survival model
    fit <- survfit(Surv(time = overall_survival, 
                       event = deceased) ~ group, 
                  data = merged_data)
    
    # Create plot
    plot <- ggsurvplot(
      fit,
      data = merged_data,
      pval = TRUE,
      risk.table = TRUE,
      tables.height = 0.3,
      title = sprintf("%s survival by %s %s", cancer_type, gene, data_type)
    )
    
    # Save plots with absolute paths
    plot_base <- file.path(RESULTS_DIR, 
                          sprintf("%s_%s_%s_survival", cancer_type, gene, data_type))
    ggsave(paste0(plot_base, ".pdf"), plot = plot$plot, width = 10, height = 8)
    ggsave(paste0(plot_base, ".png"), plot = plot$plot, width = 10, height = 8, dpi = 300)
    
    message(sprintf("Saved plots to %s.{pdf,png}", plot_base))
    
    # Save results
    results <- list(
      fit = fit,
      data = merged_data,
      summary = list(
        cancer_type = cancer_type,
        gene = gene,
        data_type = data_type,
        total_samples = nrow(merged_data),
        group_sizes = table(merged_data$group)
      )
    )
    
    # Save results to RDS
    results_file <- file.path(RESULTS_DIR, sprintf("%s_results.rds", cancer_type))
    saveRDS(results, results_file)
    message(sprintf("Saved results to: %s", results_file))
    
    # Save summary to CSV
    summary_df <- data.frame(
      cancer_type = cancer_type,
      gene = gene,
      data_type = data_type,
      total_samples = nrow(merged_data),
      high_group = sum(merged_data$group == "High"),
      low_group = sum(merged_data$group == "Low"),
      medium_group = sum(merged_data$group == "Medium")
    )
    write.csv(summary_df, file.path(RESULTS_DIR, sprintf("%s_summary.csv", cancer_type)), row.names = FALSE)
    
    return(results)
    
  }, error = function(e) {
    message(sprintf("Error in %s %s analysis: %s", gene, data_type, e$message))
    return(NULL)
  })
}

#####################################################################
# Main Execution
#####################################################################

# Main analysis execution
main <- function() {
  # Initialize results list
  results_list <- list()
  
  # Get cancer type from SLURM array task ID
  array_task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "-1"))
  message(sprintf("Starting main function with array task ID: %d", array_task_id))
  
  if (array_task_id >= 0) {
    cancer_types <- c("ACC", "UVM", "SKCM", "LGG", "GBM")
    if (array_task_id < length(cancer_types)) {
      cancer_type <- cancer_types[array_task_id + 1]
      message(sprintf("\nProcessing cancer type: %s", cancer_type))
      
      # Check if clinical data directory exists
      clinical_dir <- file.path(work_dir, "../../TCGA_clinical")
      message(sprintf("Looking for clinical data in: %s", clinical_dir))
      if (!dir.exists(clinical_dir)) {
        stop(sprintf("Clinical data directory not found: %s", clinical_dir))
      }
      
      # Run each analysis type separately with its own error handling
      
      # 1. PSI Analysis for SRRM3
      tryCatch({
        message("\nRunning SRRM3 PSI analysis...")
        results_list[["SRRM3_PSI"]] <- perform_survival_analysis(cancer_type, "PSI", "SRRM3")
      }, error = function(e) {
        message(sprintf("Error in SRRM3 PSI analysis: %s", e$message))
      })
      
      # 2. Expression Analysis for SRRM3
      tryCatch({
        message("\nRunning SRRM3 expression analysis...")
        results_list[["SRRM3_expression"]] <- perform_survival_analysis(cancer_type, "expression", "SRRM3")
      }, error = function(e) {
        message(sprintf("Error in SRRM3 expression analysis: %s", e$message))
      })
      
      # 3. Expression Analysis for SRRM4
      tryCatch({
        message("\nRunning SRRM4 expression analysis...")
        results_list[["SRRM4_expression"]] <- perform_survival_analysis(cancer_type, "expression", "SRRM4")
      }, error = function(e) {
        message(sprintf("Error in SRRM4 expression analysis: %s", e$message))
      })
      
      # Save all results with absolute path
      results_file <- file.path(RESULTS_DIR, sprintf("%s_results.rds", cancer_type))
      saveRDS(results_list, results_file)
      message(sprintf("Saved results to: %s", results_file))
      
      message(sprintf("\nCompleted all analyses for %s at %s", 
                     cancer_type, 
                     format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    }
  }
}

# Execute main function if script is run directly
if (sys.nframe() == 0) {
  main()
}