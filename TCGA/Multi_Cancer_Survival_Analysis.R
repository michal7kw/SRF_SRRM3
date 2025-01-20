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
})

# Set up parallel processing
num_cores <- as.numeric(Sys.getenv("SLURM_NTASKS", unset = "1"))
if (num_cores > 1) {
  message(sprintf("Setting up parallel processing with %d cores", num_cores))
  BiocParallel::register(MulticoreParam(workers = num_cores))
} else {
  message("Running in single-core mode")
}

# Create necessary directories
dir.create("results", showWarnings = FALSE)
dir.create("cache", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

# Cache connection to Ensembl
ensembl <- NULL
get_ensembl_connection <- function() {
  if (is.null(ensembl)) {
    ensembl <<- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  }
  return(ensembl)
}

# Function to get clinical data
get_clinical_data <- function(cancer_type) {
  message(sprintf("Getting clinical data for %s", cancer_type))
  clinical <- GDCquery_clinic(project = paste0("TCGA-", cancer_type), type = "clinical")
  
  if (nrow(clinical) == 0) {
    stop(sprintf("No clinical data found for %s", cancer_type))
  }
  
  # Print column names for debugging
  message("Available clinical columns:")
  message(paste(colnames(clinical), collapse=", "))
  
  # Clean and prepare clinical data
  clinical_clean <- clinical %>%
    dplyr::transmute(
      case_id = submitter_id,
      vital_status = vital_status,
      # Handle survival time with different column names
      overall_survival = case_when(
        !is.na(days_to_death) ~ as.numeric(days_to_death),
        !is.na(days_to_last_follow_up) ~ as.numeric(days_to_last_follow_up),
        TRUE ~ NA_real_
      ),
      # Handle different vital status formats
      deceased = case_when(
        tolower(vital_status) == "dead" ~ TRUE,
        tolower(vital_status) == "alive" ~ FALSE,
        TRUE ~ NA
      )
    ) %>%
    dplyr::filter(!is.na(overall_survival), !is.na(deceased))
  
  # Print summary of processed data
  message(sprintf("Processed clinical data summary:"))
  message(sprintf("- Total patients: %d", nrow(clinical_clean)))
  message(sprintf("- Patients with death events: %d", sum(clinical_clean$deceased)))
  message(sprintf("- Median survival time: %.1f days", median(clinical_clean$overall_survival)))
  
  return(clinical_clean)
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
    stop(sprintf("Could not find unique project for %s", cancer_type))
  }
  
  message("Creating RSE object...")
  rse_gene <- create_rse(
    project_info[1, ],
    type = "gene",
    annotation = "gencode_v26"
  )
  
  # Get gene expression
  gene_idx <- which(rowData(rse_gene)$gene_name == gene_name)[1]
  if (is.na(gene_idx)) {
    stop(sprintf("Gene %s not found in dataset", gene_name))
  }
  
  # Create data frame
  expression_data <- data.frame(
    case_id = substr(colData(rse_gene)$tcga.tcga_barcode, 1, 12),
    expression = assays(rse_gene)$counts[gene_idx, ]
  )
  
  # Handle duplicates using dplyr
  if (any(duplicated(expression_data$case_id))) {
    expression_data <- expression_data %>%
      group_by(case_id) %>%
      summarize(expression = mean(expression, na.rm = TRUE)) %>%
      ungroup()
  }
  
  saveRDS(expression_data, cache_file)
  return(expression_data)
}

# Function to get PSI data
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

# Function to perform multi-cancer analysis
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

# Function to perform survival analysis
perform_survival_analysis <- function(cancer_type, data_type = "PSI", gene = "SRRM3", grouping_method = "quartile") {
  # Get clinical data
  clinical_data <- get_clinical_data(cancer_type)
  message(sprintf("Clinical data: %d samples", nrow(clinical_data)))
  
  # Get molecular data based on data type
  molecular_data <- if (data_type == "PSI") {
    get_psi_data(cancer_type, gene)
  } else {
    get_expression_data(cancer_type, gene)
  }
  message(sprintf("Molecular data: %d samples", nrow(molecular_data)))
  
  # Ensure IDs are in the same format before merging
  clinical_data$merge_id <- clinical_data$case_id
  molecular_data$merge_id <- substr(molecular_data$case_id, 1, 12)
  
  # Merge the data
  merged_data <- inner_join(clinical_data, molecular_data, by = "merge_id")
  message(sprintf("Merged data samples: %d", nrow(merged_data)))
  
  if (nrow(merged_data) == 0) {
    stop("No samples after merging clinical and molecular data")
  }
  
  # Group samples based on molecular values
  value_col <- if(data_type == "PSI") "psi" else "expression"
  
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
  
  # Save plots
  plot_base <- file.path("results", sprintf("%s_%s_%s_survival", cancer_type, gene, data_type))
  ggsave(paste0(plot_base, ".pdf"), plot = plot$plot, width = 10, height = 8)
  ggsave(paste0(plot_base, ".png"), plot = plot$plot, width = 10, height = 8, dpi = 300)
  
  message(sprintf("Saved plots to %s.{pdf,png}", plot_base))
  
  return(list(fit = fit, data = merged_data))
}

# Main analysis
cancer_type <- NULL
results_list <- list()  # Initialize results list at the top level

array_task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "-1"))
if (array_task_id >= 0) {
  cancer_types <- c("BRCA", "LUAD", "COAD", "GBM", "KIRC", "PRAD")
  if (array_task_id < length(cancer_types)) {
    cancer_type <- cancer_types[array_task_id + 1]
    message(sprintf("\nProcessing cancer type: %s", cancer_type))
  }
}

if (!is.null(cancer_type)) {
  # Run each analysis type separately with its own error handling
  
  # 1. PSI Analysis for SRRM3
  tryCatch({
    message("\nRunning SRRM3 PSI analysis...")
    results_psi <- perform_survival_analysis(cancer_type, "PSI", "SRRM3", "quartile")
    results_list[[paste0(cancer_type, "_SRRM3_PSI")]] <- results_psi
    message("SRRM3 PSI analysis completed successfully")
  }, error = function(e) {
    message(sprintf("Error in SRRM3 PSI analysis: %s", conditionMessage(e)))
  })
  
  # 2. Expression Analysis for SRRM3
  tryCatch({
    message("\nRunning SRRM3 expression analysis...")
    results_expr_srrm3 <- perform_survival_analysis(cancer_type, "expression", "SRRM3", "quartile")
    results_list[[paste0(cancer_type, "_SRRM3_expression")]] <- results_expr_srrm3
    message("SRRM3 expression analysis completed successfully")
  }, error = function(e) {
    message(sprintf("Error in SRRM3 expression analysis: %s", conditionMessage(e)))
  })
  
  # 3. Expression Analysis for SRRM4
  tryCatch({
    message("\nRunning SRRM4 expression analysis...")
    results_expr_srrm4 <- perform_survival_analysis(cancer_type, "expression", "SRRM4", "quartile")
    results_list[[paste0(cancer_type, "_SRRM4_expression")]] <- results_expr_srrm4
    message("SRRM4 expression analysis completed successfully")
  }, error = function(e) {
    message(sprintf("Error in SRRM4 expression analysis: %s", conditionMessage(e)))
  })
  
  # Save results if any analyses were successful
  if (length(results_list) > 0) {
    results_file <- file.path("results", paste0(cancer_type, "_results.rds"))
    saveRDS(results_list, results_file)
    message(sprintf("\nSaved results to %s", results_file))
    
    # Create summary for this cancer type
    summary_data <- data.frame()
    
    for (result_name in names(results_list)) {
      result <- results_list[[result_name]]
      if (!is.null(result) && !is.null(result$data) && !is.null(result$fit)) {
        # Extract analysis components from result name
        parts <- strsplit(result_name, "_")[[1]]
        gene <- parts[2]
        analysis_type <- parts[3]
        
        # Calculate log-rank p-value safely
        p_value <- tryCatch({
          1 - pchisq(result$fit$chisq, df = length(unique(result$data$group)) - 1)
        }, error = function(e) {
          NA_real_
        })
        
        # Create summary row
        summary_stats <- data.frame(
          cancer_type = cancer_type,
          gene = gene,
          analysis_type = analysis_type,
          n_samples = nrow(result$data),
          n_high = sum(result$data$group == "High"),
          n_medium = sum(result$data$group == "Medium"),
          n_low = sum(result$data$group == "Low"),
          logrank_pvalue = p_value
        )
        
        summary_data <- rbind(summary_data, summary_stats)
      }
    }
    
    # Save summary if we have data
    if (nrow(summary_data) > 0) {
      summary_file <- file.path("results", paste0(cancer_type, "_summary.csv"))
      write.csv(summary_data, summary_file, row.names = FALSE)
      message(sprintf("Saved summary to %s", summary_file))
      message("\nAnalysis summary:")
      print(summary_data)
    } else {
      message("No successful analyses to summarize")
    }
  } else {
    message("\nNo successful analyses completed")
  }
  
  message(sprintf("\nCompleted all analyses for %s at %s", cancer_type, format(Sys.time())))
} else {
  message("No valid cancer type specified")
}