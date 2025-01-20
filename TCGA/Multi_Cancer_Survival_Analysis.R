# Load required libraries
library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(recount3)
library(GenomicRanges)
library(viridis)
library(httr)
library(retry)
library(futile.logger)
library(GenomicFeatures)
library(rtracklayer)
library(matrixStats)  # for rowVars
library(biomaRt)
library(sparseMatrixStats)  # for rowVars with sparse matrices

# Create cache directory
if (!dir.exists("cache")) {
  dir.create("cache")
}

# Function to get clinical data
get_clinical_data <- function(cancer_type) {
  clinical <- GDCquery_clinic(paste0("TCGA-", cancer_type))
  
  # Process survival information and ensure consistent column names
  clinical <- clinical %>%
    dplyr::mutate(
      submitter_id = bcr_patient_barcode,  # Use bcr_patient_barcode instead of submitter_id
      deceased = ifelse(vital_status == "Alive", FALSE, TRUE),
      overall_survival = ifelse(vital_status == "Alive",
                              days_to_last_follow_up,
                              days_to_death)
    ) %>%
    dplyr::filter(!is.na(overall_survival)) %>%  # Remove samples with missing survival
    dplyr::select(submitter_id, deceased, overall_survival, dplyr::everything())  # Ensure key columns are present
  
  message(sprintf("Clinical data has %d patients with valid survival data", nrow(clinical)))
  return(clinical)
}

# Function to get expression data
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
    stop("Could not find unique project for ", cancer_type)
  }
  
  message("Creating RSE object...")
  rse_gene <- create_rse(
    project_info[1, ],
    type = "gene"
  )
  
  # Get gene expression
  gene_idx <- which(rowData(rse_gene)$gene_name == gene_name)
  if (length(gene_idx) == 0) {
    stop(paste("Gene", gene_name, "not found in dataset"))
  }
  
  expression_data <- data.frame(
    case_id = substr(colData(rse_gene)$tcga.tcga_barcode, 1, 12),  # Trim to first 12 characters
    expression = assays(rse_gene)$counts[gene_idx[1], ]
  )
  
  # Remove any duplicate IDs by taking the mean expression value
  if (any(duplicated(expression_data$case_id))) {
    message("Found duplicate sample IDs, taking mean expression values")
    expression_data <- expression_data %>%
      dplyr::group_by(case_id) %>%
      dplyr::summarize(expression = mean(expression, na.rm = TRUE)) %>%
      dplyr::ungroup()
  }
  
  saveRDS(expression_data, cache_file)
  return(expression_data)
}

# Update get_gene_info function
get_gene_info <- function(gene_name) {
  # Connect to Ensembl
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Get gene coordinates
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
  
  # If multiple entries exist, use the first one
  if (nrow(gene_info) > 1) {
    message("Multiple entries found for ", gene_name, ", using first entry")
    gene_info <- gene_info[1, ]
  }
  
  # Add chromosome prefix if missing
  if (!grepl("^chr", gene_info$chromosome_name)) {
    gene_info$chromosome_name <- paste0("chr", gene_info$chromosome_name)
  }
  
  return(gene_info)
}

# Update get_psi_data function
get_psi_data <- function(cancer_type, gene_name) {
  cache_file <- file.path("cache", paste0("psi_data_", cancer_type, "_", gene_name, ".rds"))
  
  if (file.exists(cache_file)) {
    message("Loading cached PSI data...")
    return(readRDS(cache_file))
  }
  
  message("Getting PSI data from recount3...")
  projects <- available_projects()
  
  project_info <- subset(projects, 
                        project == toupper(cancer_type) &
                        file_source == "tcga" &
                        project_type == "data_sources")
  
  message("Creating RSE object...")
  rse_jxn <- create_rse(
    project_info[1, ],
    type = "jxn"
  )
  
  # Get gene information including coordinates
  gene_info <- get_gene_info(gene_name)
  message(sprintf("Found gene %s at %s:%d-%d", 
                 gene_name, 
                 gene_info$chromosome_name, 
                 gene_info$start_position, 
                 gene_info$end_position))
  
  # Get junctions overlapping the gene region
  gene_junctions <- rowRanges(rse_jxn)
  
  # Fix chromosome names to match
  seqlevels(gene_junctions) <- sub("^chr", "", seqlevels(gene_junctions))
  gene_info$chromosome_name <- sub("^chr", "", gene_info$chromosome_name)
  
  # Create GRanges object for the gene
  gene_range <- GRanges(
    seqnames = gene_info$chromosome_name,
    ranges = IRanges(
      start = gene_info$start_position,
      end = gene_info$end_position
    )
  )
  
  # Find overlapping junctions
  overlaps <- findOverlaps(gene_junctions, gene_range)
  
  if (length(overlaps) == 0) {
    # Debug information
    message("Available chromosomes in junction data: ", 
            paste(seqlevels(gene_junctions), collapse=", "))
    message("Looking for chromosome: ", gene_info$chromosome_name)
    stop(paste("No junctions found in the region of", gene_name))
  }
  
  message(sprintf("Found %d junctions overlapping %s", length(overlaps), gene_name))
  
  # Get junction counts and calculate PSI
  junction_counts <- assays(rse_jxn)$counts[queryHits(overlaps), ]
  
  # Convert sparse matrix to dense if needed
  if (is(junction_counts, "dgTMatrix") || is(junction_counts, "dgCMatrix")) {
    junction_counts <- as.matrix(junction_counts)
  }
  
  junction_vars <- rowVars(junction_counts)
  
  # Get sample IDs and ensure they are in correct TCGA format
  sample_ids <- colData(rse_jxn)$tcga.tcga_barcode
  if (is.null(sample_ids)) {
    stop("Could not find TCGA barcodes in the data")
  }
  
  # Ensure consistent 12-character TCGA ID format
  sample_ids <- substr(sample_ids, 1, 12)
  
  # Calculate PSI for the most variable junction
  max_var_idx <- which.max(junction_vars)
  psi_values <- junction_counts[max_var_idx, ] / rowSums(junction_counts)
  
  # Create data frame with PSI values
  psi_data <- data.frame(
    case_id = sample_ids,
    psi = as.numeric(psi_values),
    stringsAsFactors = FALSE
  ) %>%
    # Remove any NA values
    filter(!is.na(psi))
    
  # Remove any duplicate IDs by taking the mean PSI value
  if (any(duplicated(psi_data$case_id))) {
    message("Found duplicate sample IDs, taking mean PSI values")
    psi_data <- psi_data %>%
      group_by(case_id) %>%
      summarize(psi = mean(psi, na.rm = TRUE)) %>%
      ungroup()
  }
  
  message(sprintf("Found %d samples with valid PSI values", nrow(psi_data)))
  
  saveRDS(psi_data, cache_file)
  return(psi_data)
}

# Add this helper function
check_sample_ids <- function(clinical_ids, molecular_ids) {
  # Get first few examples of each
  clin_examples <- head(clinical_ids)
  mol_examples <- head(molecular_ids)
  
  message("Sample ID format check:")
  message("Clinical ID examples: ", paste(clin_examples, collapse=", "))
  message("Molecular ID examples: ", paste(mol_examples, collapse=", "))
  message("Molecular IDs trimmed to 12 chars: ", 
          paste(head(substr(molecular_ids, 1, 12)), collapse=", "))
  
  # Check for common TCGA barcode lengths
  clin_lengths <- unique(nchar(clinical_ids))
  mol_lengths <- unique(nchar(molecular_ids))
  
  message("Clinical ID lengths: ", paste(clin_lengths, collapse=", "))
  message("Molecular ID lengths: ", paste(mol_lengths, collapse=", "))
  
  # Find matching samples using trimmed molecular IDs
  matches <- sum(clinical_ids %in% substr(molecular_ids, 1, 12))
  message(sprintf("Found %d matching samples between clinical and molecular data", matches))
}

# Function to perform survival analysis
perform_survival_analysis <- function(cancer_type, analysis_type = "PSI", gene = "SRRM3", grouping_method = "quartile") {
  message(sprintf("\nAnalyzing %s for %s using %s grouping", gene, cancer_type, grouping_method))
  
  # Get clinical data
  clinical_data <- get_clinical_data(cancer_type)
  message(sprintf("Clinical data: %d samples", nrow(clinical_data)))
  
  # Get molecular data based on analysis type
  if (analysis_type == "PSI") {
    message("Loading cached PSI data...")
    molecular_data <- get_psi_data(cancer_type, gene)
    message(sprintf("Molecular data: %d samples", nrow(molecular_data)))
  } else {
    molecular_data <- get_expression_data(cancer_type, gene)
    message(sprintf("Molecular data: %d samples", nrow(molecular_data)))
  }
  
  # Debug information for sample IDs
  message("Sample ID format check:")
  message("Clinical ID examples: ", paste(head(clinical_data$submitter_id), collapse=", "))
  message("Molecular ID examples: ", paste(head(molecular_data$case_id), collapse=", "))
  message("Molecular IDs trimmed to 12 chars: ", paste(head(substr(molecular_data$case_id, 1, 12)), collapse=", "))
  message("Clinical ID lengths: ", unique(nchar(clinical_data$submitter_id)))
  message("Molecular ID lengths: ", unique(nchar(molecular_data$case_id)))
  
  # Ensure IDs are in the same format before merging
  clinical_data$merge_id <- clinical_data$submitter_id
  molecular_data$merge_id <- substr(molecular_data$case_id, 1, 12)  # Trim molecular IDs to 12 chars
  
  # Find matching samples
  matching_samples <- intersect(clinical_data$merge_id, molecular_data$merge_id)
  message(sprintf("Found %d matching samples between clinical and molecular data", length(matching_samples)))
  
  # Additional debug info about sample matching
  message("\nSample matching details:")
  message("Number of unique clinical IDs: ", length(unique(clinical_data$merge_id)))
  message("Number of unique molecular IDs: ", length(unique(molecular_data$merge_id)))
  message("Number of duplicate clinical IDs: ", sum(duplicated(clinical_data$merge_id)))
  message("Number of duplicate molecular IDs: ", sum(duplicated(molecular_data$merge_id)))
  
  # Merge the data
  merged_data <- inner_join(
    clinical_data,
    molecular_data,
    by = "merge_id"
  )
  
  message(sprintf("After merging: %d samples", nrow(merged_data)))
  
  # Check for any unexpected duplicates in merged data
  if (nrow(merged_data) > length(matching_samples)) {
    message("\nWarning: More rows in merged data than matching samples")
    message("This may be due to duplicate IDs in either clinical or molecular data")
    message("Samples that appear multiple times:")
    duplicate_samples <- merged_data$merge_id[duplicated(merged_data$merge_id)]
    message(paste(duplicate_samples, collapse=", "))
  }
  
  # Add debug info
  message("Debug info:")
  message("First few clinical IDs: ", paste(head(clinical_data$submitter_id), collapse=", "))
  message("First few molecular IDs: ", paste(head(molecular_data$case_id), collapse=", "))
  message("First few clinical bcr_patient_barcodes: ", paste(head(clinical_data$bcr_patient_barcode), collapse=", "))
  message("Clinical data columns: ", paste(colnames(clinical_data), collapse=", "))
  
  if (nrow(merged_data) == 0) {
    stop("No valid samples after merging and filtering")
  }
  
  # Group samples based on molecular values
  value_col <- ifelse(analysis_type == "PSI", "psi", "expression")
  
  if (grouping_method == "median") {
    merged_data$group <- ifelse(
      merged_data[[value_col]] > median(merged_data[[value_col]], na.rm = TRUE), 
      "High", "Low"
    )
  } else if (grouping_method == "quartile") {
    quartiles <- quantile(merged_data[[value_col]], probs = c(0.25, 0.75), na.rm = TRUE)
    merged_data$group <- case_when(
      merged_data[[value_col]] <= quartiles[1] ~ "Low",
      merged_data[[value_col]] >= quartiles[2] ~ "High",
      TRUE ~ "Medium"
    )
  }
  
  # Clean up data
  merged_data <- merged_data %>%
    dplyr::mutate(
      overall_survival = as.numeric(overall_survival),
      deceased = as.logical(deceased)
    ) %>%
    dplyr::filter(!is.na(overall_survival), !is.na(deceased))
  
  message(sprintf("After cleaning: %d samples", nrow(merged_data)))
  
  if (nrow(merged_data) == 0) {
    stop("No valid samples after merging and filtering")
  }
  
  # Create survival object and fit
  message("Fitting survival model...")
  fit <- survfit(Surv(overall_survival, deceased) ~ group, data = merged_data)
  
  # Add summary statistics calculation
  message("Calculating group statistics...")
  group_stats <- merged_data %>%
    group_by(group) %>%
    summarise(
      mean_val = mean(!!sym(value_col), na.rm = TRUE),
      median_val = median(!!sym(value_col), na.rm = TRUE),
      sd_val = sd(!!sym(value_col), na.rm = TRUE),
      min_val = min(!!sym(value_col), na.rm = TRUE),
      max_val = max(!!sym(value_col), na.rm = TRUE),
      n = n(),
      .groups = 'drop'
    )
  
  # Create subtitle with ranges
  subtitle <- paste(
    sapply(1:nrow(group_stats), function(i) {
      group <- group_stats$group[i]
      if(analysis_type == "PSI") {
        sprintf("%s: PSI %.1f-%.1f%% (median: %.1f%%, n=%d)",
                group,
                group_stats$min_val[i] * 100,
                group_stats$max_val[i] * 100,
                group_stats$median_val[i] * 100,
                group_stats$n[i])
      } else {
        sprintf("%s: Expr %.1f-%.1f (median: %.1f, n=%d)",
                group,
                group_stats$min_val[i],
                group_stats$max_val[i],
                group_stats$median_val[i],
                group_stats$n[i])
      }
    }),
    collapse = "\n"
  )
  
  # Create plot
  message("Creating survival plot...")
  surv_plot <- ggsurvplot(
    fit,
    data = merged_data,
    pval = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    xlab = "Time (days)",
    ylab = "Overall Survival Probability",
    title = sprintf("%s %s Survival Analysis - %s", gene, analysis_type, cancer_type),
    subtitle = subtitle,
    legend.title = paste(analysis_type, "Level"),
    palette = c("High" = "#E7B800", "Low" = "#2E9FDF", "Medium" = "#868686"),
    font.subtitle = c(10, "plain"),
    break.time.by = 365,
    ggtheme = theme_bw() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
      )
  )
  
  # Create results directory if it doesn't exist
  if (!dir.exists("results")) {
    dir.create("results")
  }
  
  # Save plots
  plot_base_name <- file.path("results", 
                             sprintf("%s_%s_%s_survival", 
                                   cancer_type, gene, analysis_type))
  
  # Save as PDF
  pdf(paste0(plot_base_name, ".pdf"), width = 10, height = 8)
  print(surv_plot)
  dev.off()
  
  # Save as PNG
  png(paste0(plot_base_name, ".png"), 
      width = 10*100, height = 8*100, res = 100)
  print(surv_plot)
  dev.off()
  
  # Save results as RDS
  saveRDS(list(
    plot = surv_plot,
    fit = fit,
    data = merged_data,
    group_stats = group_stats
  ), paste0(plot_base_name, "_results.rds"))
  
  message(sprintf("Results saved to %s.{pdf,png,rds}", plot_base_name))
  
  return(list(
    plot = surv_plot,
    fit = fit,
    data = merged_data,
    group_stats = group_stats
  ))
}

# Create results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

# Get list of cancer types
cancer_types <- c("BRCA", "LUAD", "COAD", "GBM", "KIRC", "PRAD")

# Run analysis for each cancer type
results_list <- list()

for(cancer_type in cancer_types) {
  tryCatch({
    message(sprintf("\nProcessing %s", cancer_type))
    
    # Run PSI analysis for SRRM3
    message("\nRunning SRRM3 PSI analysis...")
    results_psi <- perform_survival_analysis(cancer_type, "PSI", "SRRM3", "quartile")
    results_list[[paste0(cancer_type, "_SRRM3_PSI")]] <- results_psi
    
    # Run expression analysis for SRRM3
    message("\nRunning SRRM3 expression analysis...")
    results_srrm3_expr <- perform_survival_analysis(cancer_type, "expression", "SRRM3", "quartile")
    results_list[[paste0(cancer_type, "_SRRM3_expression")]] <- results_srrm3_expr
    
    # Run expression analysis for SRRM4
    message("\nRunning SRRM4 expression analysis...")
    results_srrm4_expr <- perform_survival_analysis(cancer_type, "expression", "SRRM4", "quartile")
    results_list[[paste0(cancer_type, "_SRRM4_expression")]] <- results_srrm4_expr
    
  }, error = function(e) {
    message(sprintf("Error analyzing %s: %s", cancer_type, e$message))
  })
}

# Create a summary of all analyses
message("\nCreating analysis summary...")
summary_data <- data.frame()

for(result_name in names(results_list)) {
  result <- results_list[[result_name]]
  if (!is.null(result)) {
    # Extract cancer type and analysis type from result name
    parts <- strsplit(result_name, "_")[[1]]
    cancer_type <- parts[1]
    gene <- parts[2]
    analysis_type <- parts[3]
    
    # Get summary statistics
    summary_stats <- data.frame(
      cancer_type = cancer_type,
      gene = gene,
      analysis_type = analysis_type,
      n_samples = nrow(result$data),
      n_high = sum(result$data$group == "High"),
      n_medium = sum(result$data$group == "Medium"),
      n_low = sum(result$data$group == "Low"),
      logrank_pvalue = 1 - pchisq(result$fit$chisq, df = length(unique(result$data$group)) - 1)
    )
    
    summary_data <- rbind(summary_data, summary_stats)
  }
}

# Save summary to CSV
write.csv(summary_data, 
          file.path("results", "survival_analysis_summary.csv"), 
          row.names = FALSE)

message("\nAnalysis complete! Results have been saved to the results directory")
message("\nSummary of analyses:")
print(summary_data)