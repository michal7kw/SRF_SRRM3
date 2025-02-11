#####################################################################
# High SRRM3 PSI Analysis Script
#####################################################################
# This script performs PSI analysis specifically for samples with high
# SRRM3 expression. It reuses functions from Multi_Cancer_Survival_Analysis.R
# and adds specific functionality for high SRRM3 expression analysis.

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(survminer)
  library(recount3)
  library(biomaRt)
  library(parallel)
  library(BiocParallel)
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(tidyverse)
  library(DESeq2)
  library(httr)
  library(retry)
  library(futile.logger)
  library(GenomicFeatures)
  library(rtracklayer)
  library(matrixStats)
  library(sparseMatrixStats)
  library(viridis)
  library(data.table)
  library(R.utils)
})

# At the start of the file, add options
options(run.main=FALSE)
options(verbose = TRUE)
options(error = function() traceback(2))
options(future.globals.maxSize = 8000 * 1024^2)
options(mc.cores = parallel::detectCores() - 1)

#####################################################################
# Clinical Data Retrieval and Processing
#####################################################################
# Function to retrieve and clean clinical data from TCGA
get_clinical_data <- function(cancer_type) {
  message(sprintf("Getting clinical data for %s", cancer_type))
  
  tryCatch({
    # Get clinical data
    clinical <- as.data.frame(GDCquery_clinic(project = paste0("TCGA-", cancer_type), type = "clinical"))
    
    if (nrow(clinical) == 0) {
      stop(sprintf("No clinical data found for %s", cancer_type))
    }
    
    # Clean and prepare clinical data
    clinical_clean <- clinical %>%
      dplyr::transmute(
        case_id = submitter_id,
        vital_status = vital_status,
        overall_survival = case_when(
          !is.na(days_to_death) ~ as.numeric(days_to_death),
          !is.na(days_to_last_follow_up) ~ as.numeric(days_to_last_follow_up),
          TRUE ~ NA_real_
        ),
        deceased = case_when(
          tolower(vital_status) == "dead" ~ TRUE,
          tolower(vital_status) == "alive" ~ FALSE,
          TRUE ~ NA
        )
      ) %>%
      dplyr::filter(!is.na(overall_survival), !is.na(deceased))
    
    return(clinical_clean)
    
  }, error = function(e) {
    stop(sprintf("Error retrieving clinical data for %s: %s", cancer_type, e$message))
  })
}

#####################################################################
# Expression Data Retrieval
#####################################################################
get_expression_data <- function(cancer_type, gene_name) {
  cache_file <- file.path("cache", paste0("expression_data_", cancer_type, "_", gene_name, ".rds"))
  
  if (file.exists(cache_file)) {
    message("Loading cached expression data...")
    return(readRDS(cache_file))
  }
  
  message("Getting expression data from recount3...")
  projects <- available_projects()
  
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
  
  # Debug information
  message("RSE object dimensions: ", paste(dim(rse_gene), collapse=" x "))
  message("Looking for gene ", gene_name, " in dataset...")
  gene_names <- rowData(rse_gene)$gene_name
  message("First few gene names: ", paste(head(gene_names), collapse=", "))
  
  # Extract expression for specific gene
  gene_idx <- which(rowData(rse_gene)$gene_name == gene_name)[1]
  if (is.na(gene_idx)) {
    # Try case-insensitive search
    gene_idx <- which(toupper(rowData(rse_gene)$gene_name) == toupper(gene_name))[1]
  }
  
  if (is.na(gene_idx)) {
    message("Available gene name patterns:")
    message(paste(grep(toupper(gene_name), toupper(gene_names), value=TRUE), collapse=", "))
    stop(sprintf("Gene %s not found in dataset", gene_name))
  }
  
  message("Found gene at index: ", gene_idx)
  message("Creating expression data frame...")
  
  expression_data <- data.frame(
    case_id = substr(colData(rse_gene)$tcga.tcga_barcode, 1, 12),
    expression = assays(rse_gene)$counts[gene_idx, ]
  )
  
  message("Initial expression data dimensions: ", paste(dim(expression_data), collapse=" x "))
  
  # Handle duplicates
  if (any(duplicated(expression_data$case_id))) {
    message("Found duplicated case IDs, aggregating...")
    expression_data <- expression_data %>%
      group_by(case_id) %>%
      summarize(expression = mean(expression, na.rm = TRUE)) %>%
      ungroup()
    message("After aggregation dimensions: ", paste(dim(expression_data), collapse=" x "))
  }
  
  # Ensure no NA or infinite values
  expression_data <- expression_data %>%
    filter(!is.na(expression), is.finite(expression))
  message("Final expression data dimensions: ", paste(dim(expression_data), collapse=" x "))
  
  saveRDS(expression_data, cache_file)
  return(expression_data)
}

#####################################################################
# High SRRM3 Sample Selection
#####################################################################
# Function to get high SRRM3 expression samples
get_high_srrm3_samples <- function(cancer_type, threshold = 0.75) {
  message(sprintf("Getting high SRRM3 expression samples for %s", cancer_type))
  
  expression_data <- get_expression_data(cancer_type, "SRRM3")
  
  if (nrow(expression_data) == 0) {
    stop("No expression data available")
  }
  
  # Ensure expression values are numeric
  expression_data$expression <- as.numeric(as.character(expression_data$expression))
  
  # Remove any NA or infinite values
  expression_data <- expression_data %>%
    filter(!is.na(expression), is.finite(expression))
  
  if (nrow(expression_data) == 0) {
    stop("No valid expression values after cleaning")
  }
  
  # Calculate cutoff and get high expression samples
  cutoff <- quantile(expression_data$expression, probs = threshold, na.rm = TRUE)
  message(sprintf("Expression cutoff: %.2f", cutoff))
  
  high_expr_samples <- expression_data %>%
    filter(!is.na(expression), is.finite(expression)) %>%
    filter(expression >= cutoff) %>%
    pull(case_id)
  
  message(sprintf("Found %d samples with high SRRM3 expression (>%.2f percentile)",
                 length(high_expr_samples), threshold * 100))
  
  return(high_expr_samples)
}

#####################################################################
# PSI Data Retrieval and Processing
#####################################################################
get_psi_data <- function(cancer_type, gene_name, sample_ids = NULL) {
  if (is.null(cancer_type) || is.null(gene_name)) {
    stop("Cancer type and gene name must be provided")
  }
  
  cache_file <- file.path("cache", paste0("psi_data_", cancer_type, "_", gene_name, ".rds"))
  
  # Check cache
  if (file.exists(cache_file)) {
    message("Loading cached PSI data...")
    return(readRDS(cache_file))
  }
  
  # Define SRRM3 info
  SRRM3_INFO <- list(
    gene = list(name = "SRRM3", chr = "chr7", start = 76201896, end = 76287287),
    exon15 = list(start = 76283524, end = 76283602, length = 79)
  )
  
  message("Getting junction data from recount3...")
  projects <- available_projects()
  project_info <- subset(projects, 
                        project == toupper(cancer_type) & 
                        file_source == "tcga" & 
                        project_type == "data_sources")
  
  if (nrow(project_info) != 1) {
    stop(sprintf("Could not find unique project for %s", cancer_type))
  }
  
  # Get junction data
  rse_jxn <- create_rse(project_info[1,], type = "jxn", jxn_format = "UNIQUE")
  junction_counts <- assay(rse_jxn)
  jxn_coords <- rowRanges(rse_jxn)
  
  # Find relevant junctions
  srrm3_region <- which(
    start(jxn_coords) >= (SRRM3_INFO$exon15$start - 10000) &
    end(jxn_coords) <= (SRRM3_INFO$exon15$end + 10000) &
    seqnames(jxn_coords) == "chr7"
  )
  
  if (length(srrm3_region) == 0) {
    stop("No junctions found in SRRM3 region")
  }
  
  # Process junction data and calculate PSI
  jxn_coords <- jxn_coords[srrm3_region]
  junction_counts <- junction_counts[srrm3_region,]
  
  # Find inclusion/exclusion junctions
  inclusion_jxns <- which(
    (abs(end(jxn_coords) - SRRM3_INFO$exon15$start) <= 5) |
    (abs(start(jxn_coords) - SRRM3_INFO$exon15$end) <= 5)
  )
  
  exclusion_jxns <- which(
    start(jxn_coords) < (SRRM3_INFO$exon15$start - 5) &
    end(jxn_coords) > (SRRM3_INFO$exon15$end + 5)
  )
  
  # Calculate PSI values
  psi_values <- sapply(seq_len(ncol(junction_counts)), function(i) {
    inclusion_reads <- sum(junction_counts[inclusion_jxns, i])
    exclusion_reads <- sum(junction_counts[exclusion_jxns, i])
    total_reads <- inclusion_reads + exclusion_reads
    
    if (total_reads >= 10) {
      return((inclusion_reads/total_reads) * 100)
    } else {
      return(NA_real_)
    }
  })
  
  # Create PSI data frame
  psi_data <- data.frame(
    case_id = substr(colData(rse_jxn)$tcga.tcga_barcode, 1, 12),
    psi = psi_values
  ) %>%
    filter(!is.na(psi)) %>%
    group_by(case_id) %>%
    summarize(psi = mean(psi, na.rm = TRUE)) %>%
    ungroup()
  
  # Filter by sample IDs if provided
  if (!is.null(sample_ids)) {
    psi_data <- psi_data %>% filter(case_id %in% sample_ids)
  }
  
  # Cache results
  saveRDS(psi_data, cache_file)
  
  return(psi_data)
}

#####################################################################
# High SRRM3 PSI Analysis
#####################################################################
perform_high_srrm3_psi_analysis <- function(cancer_type, expression_threshold = 0.75,
                                          grouping_method = "quartile") {
  tryCatch({
    high_expr_samples <- get_high_srrm3_samples(cancer_type, expression_threshold)
    
    message("Number of high expression samples: ", length(high_expr_samples))
    message("First few sample IDs: ", paste(head(high_expr_samples), collapse=", "))
    
    if (length(high_expr_samples) < 10) {
      stop("Insufficient samples with high SRRM3 expression")
    }
    
    clinical_data <- get_clinical_data(cancer_type)
    message("Number of clinical samples: ", nrow(clinical_data))
    
    psi_data <- get_psi_data(cancer_type, "SRRM3", high_expr_samples)
    message("Number of PSI samples: ", nrow(psi_data))
    
    # Merge data
    merged_data <- inner_join(
      clinical_data %>% mutate(merge_id = substr(case_id, 1, 12)),
      psi_data %>% mutate(merge_id = substr(case_id, 1, 12)),
      by = "merge_id"
    )
    
    message("Number of samples after merging: ", nrow(merged_data))
    
    if (nrow(merged_data) == 0) {
      stop("No samples after merging clinical and PSI data")
    }
    
    # Group samples
    quartiles <- quantile(merged_data$psi, probs = c(0.25, 0.75), na.rm = TRUE)
    merged_data$group <- case_when(
      merged_data$psi <= quartiles[1] ~ "Low",
      merged_data$psi >= quartiles[2] ~ "High",
      TRUE ~ "Medium"
    )
    
    merged_data <- merged_data %>% 
      filter(!is.na(group))
    
    message("Sample sizes per group:")
    print(table(merged_data$group))
    
    # Survival analysis
    fit <- survfit(Surv(time = overall_survival, event = deceased) ~ group,
                  data = merged_data)
    
    # Create plot
    plot <- ggsurvplot(
      fit, data = merged_data,
      pval = TRUE,
      risk.table = TRUE,
      tables.height = 0.3,
      title = sprintf("%s survival by PSI (High SRRM3 samples)", cancer_type)
    )
    
    # Save results
    results_dir <- file.path("results", "high_srrm3_analysis")
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    
    plot_base <- file.path(results_dir, sprintf("%s_high_srrm3_psi_survival", cancer_type))
    ggsave(paste0(plot_base, ".pdf"), plot = plot$plot, width = 10, height = 8)
    ggsave(paste0(plot_base, ".png"), plot = plot$plot, width = 10, height = 8, dpi = 300)
    
    # Create summary
    summary_df <- data.frame(
      cancer_type = cancer_type,
      total_samples = nrow(merged_data),
      high_expr_samples = length(high_expr_samples),
      analyzed_samples = nrow(merged_data),
      high_psi_group = sum(merged_data$group == "High"),
      low_psi_group = sum(merged_data$group == "Low"),
      medium_psi_group = sum(merged_data$group == "Medium")
    )
    
    write.csv(summary_df, 
              file.path(results_dir, sprintf("%s_summary.csv", cancer_type)), 
              row.names = FALSE)
    
    saveRDS(list(fit = fit, data = merged_data, summary = summary_df),
            file.path(results_dir, sprintf("%s_results.rds", cancer_type)))
    
    return(list(fit = fit, data = merged_data, summary = summary_df))
    
  }, error = function(e) {
    message(sprintf("Error in high SRRM3 PSI analysis for %s: %s", cancer_type, e$message))
    return(NULL)
  })
}

# Get cancer type from environment variable
cancer_type <- Sys.getenv("CANCER_TYPE")
if (nchar(cancer_type) > 0) {
  message(sprintf("\nAnalyzing %s", cancer_type))
  
  # Remove any corrupted cache files
  cache_files <- list.files("cache", pattern = paste0(cancer_type, ".*\\.rds$"), full.names = TRUE)
  for (file in cache_files) {
    tryCatch({
      data <- readRDS(file)
    }, error = function(e) {
      message("Removing corrupted cache file: ", file)
      unlink(file)
    })
  }
  
  results <- perform_high_srrm3_psi_analysis(cancer_type, expression_threshold = 0.75)
  
  if (is.null(results)) {
    message("No results generated")
    quit(status = 1)
  }
} 