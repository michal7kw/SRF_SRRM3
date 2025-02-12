#!/usr/bin/env Rscript

# Set memory limit
memory.limit(size = 64000)  # Set to 64GB

# Configure parallel processing
options(future.globals.maxSize = 50000 * 1024^2)  # 50GB
options(future.rng.onMisuse = "ignore")

# Load required libraries with optimized settings
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
  library(viridis)
  library(data.table)
  library(pryr)  # For memory tracking
})

# Set up parallel processing
num_cores <- as.numeric(Sys.getenv("SLURM_NTASKS", unset = "1"))
if (num_cores > 1) {
  message(sprintf("Setting up parallel processing with %d cores", num_cores))
  BiocParallel::register(MulticoreParam(workers = num_cores))
  setDTthreads(num_cores)  # Set data.table threads
} else {
  message("Running in single-core mode")
}

# Create necessary directories
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("cache", showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE, recursive = TRUE)

# Memory tracking function
log_memory_usage <- function(step) {
  mem_used <- pryr::mem_used()
  gc_stats <- gc()
  message(sprintf("[%s] Memory used: %.2f GB", step, mem_used/1024^3))
  return(invisible(NULL))
}

# Function to get clinical data with optimizations
get_clinical_data <- function(cancer_type) {
  log_memory_usage("Starting clinical data retrieval")
  
  message(sprintf("Getting clinical data for %s", cancer_type))
  clinical <- GDCquery_clinic(paste0("TCGA-", cancer_type))
  
  if (nrow(clinical) == 0) {
    stop(sprintf("No clinical data found for %s", cancer_type))
  }
  
  # Process survival information
  clinical <- clinical %>%
    mutate(
      deceased = ifelse(vital_status == "Alive", FALSE, TRUE),
      overall_survival = ifelse(vital_status == "Alive",
                              days_to_last_follow_up,
                              days_to_death)
    ) %>%
    filter(!is.na(overall_survival), overall_survival > 0)
  
  message(sprintf("Processed clinical data summary:"))
  message(sprintf("- Total patients: %d", nrow(clinical)))
  message(sprintf("- Patients with death events: %d", sum(clinical$deceased)))
  message(sprintf("- Median survival time: %.1f days", median(clinical$overall_survival)))
  
  # Add debug information
  message("Sample of clinical data case IDs:")
  print(head(clinical$submitter_id))
  
  log_memory_usage("Finished clinical data processing")
  return(clinical)
}

# Function to get expression data with optimizations
get_expression_data <- function(cancer_type, gene_name) {
  log_memory_usage("Starting expression data retrieval")
  
  cache_file <- file.path("cache", paste0("expression_data_", cancer_type, "_", gene_name, ".rds"))
  
  if (file.exists(cache_file)) {
    message("Loading cached expression data...")
    cached_data <- readRDS(cache_file)
    # Verify the cached data has proper TCGA IDs
    if (!grepl("^TCGA-", cached_data$case_id[1])) {
      message("Cached data has incorrect ID format. Regenerating...")
      unlink(cache_file)
    } else {
      return(cached_data)
    }
  }
  
  message("Getting expression data from recount3...")
  projects <- available_projects()
  
  project_info <- subset(projects, 
                        project == toupper(cancer_type) &
                        file_source == "tcga" &
                        project_type == "data_sources")
  
  if (nrow(project_info) != 1) {
    message("\nAvailable TCGA projects:")
    print(subset(projects, file_source == "tcga" & project_type == "data_sources"))
    stop("Could not find unique project for ", cancer_type)
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
  
  # Get the TCGA barcodes
  tcga_barcodes <- colData(rse_gene)$tcga.tcga_barcode
  if (is.null(tcga_barcodes)) {
    message("Available colData columns:")
    print(colnames(colData(rse_gene)))
    stop("Could not find tcga.tcga_barcode in the data")
  }
  
  # Create data frame with expression data
  expression_data <- data.frame(
    case_id = tcga_barcodes,
    counts = assay(rse_gene)[gene_idx, ],
    stringsAsFactors = FALSE
  )
  
  # Clean case IDs to match TCGA format (TCGA-XX-XXXX)
  expression_data$case_id <- substr(expression_data$case_id, 1, 12)
  
  # Handle duplicates
  if (anyDuplicated(expression_data$case_id)) {
    expression_data <- expression_data %>%
      group_by(case_id) %>%
      summarise(counts = mean(counts, na.rm = TRUE)) %>%
      ungroup()
  }
  
  message("Sample of expression data case IDs:")
  print(head(expression_data$case_id))
  
  saveRDS(expression_data, cache_file)
  log_memory_usage("Finished expression data processing")
  return(expression_data)
}

# Function to perform survival analysis with optimizations
perform_survival_analysis <- function(cancer_type, gene_name, grouping_method = "median") {
  log_memory_usage("Starting survival analysis")
  
  # Get clinical data
  clinical_data <- get_clinical_data(cancer_type)
  message(sprintf("Clinical data: %d samples", nrow(clinical_data)))
  
  # Get expression data
  tryCatch({
    molecular_data <- get_expression_data(cancer_type, gene_name)
    message(sprintf("Expression data: %d samples", nrow(molecular_data)))
    
    if (nrow(molecular_data) == 0) {
      stop("No samples in expression data")
    }
    
    # Make sure case IDs are in the same format
    clinical_data$submitter_id <- substr(clinical_data$submitter_id, 1, 12)
    molecular_data$case_id <- substr(molecular_data$case_id, 1, 12)
    
    # Debug case IDs
    message("\nSample clinical IDs:")
    print(head(clinical_data$submitter_id, 10))
    message("\nSample molecular IDs:")
    print(head(molecular_data$case_id, 10))
    
    message("\nChecking for any matching IDs:")
    matching_ids <- intersect(clinical_data$submitter_id, molecular_data$case_id)
    message(sprintf("Found %d matching IDs", length(matching_ids)))
    if (length(matching_ids) > 0) {
        message("First few matching IDs:")
        print(head(matching_ids))
    }
    
    message("\nChecking ID formats:")
    message("Clinical ID format example:", clinical_data$submitter_id[1])
    message("Molecular ID format example:", molecular_data$case_id[1])
    message("\nUnique clinical ID patterns:")
    print(unique(substr(clinical_data$submitter_id, 1, 4)))
    message("\nUnique molecular ID patterns:")
    print(unique(substr(molecular_data$case_id, 1, 4)))
    
    # Merge the data
    merged_data <- merge(molecular_data, clinical_data,
                        by.x = "case_id", by.y = "submitter_id",
                        all = FALSE)  # Only keep matching cases
    
    message(sprintf("Merged data samples: %d", nrow(merged_data)))
    
    if (nrow(merged_data) == 0) {
      stop("No samples after merging clinical and expression data")
    }
    
    # Group samples based on expression values
    if (grouping_method == "median") {
      median_value <- median(merged_data$counts, na.rm = TRUE)
      message(sprintf("Median expression value: %.2f", median_value))
      merged_data$group <- ifelse(merged_data$counts >= median_value, "High", "Low")
    } else if (grouping_method == "quartile") {
      quartiles <- quantile(merged_data$counts, probs = c(0.25, 0.75), na.rm = TRUE)
      message(sprintf("Expression quartiles: Q1=%.2f, Q3=%.2f", quartiles[1], quartiles[2]))
      merged_data$group <- case_when(
        merged_data$counts <= quartiles[1] ~ "Low",
        merged_data$counts >= quartiles[2] ~ "High",
        TRUE ~ "Medium"
      )
    }
    
    # Remove NA groups and ensure all required columns are non-NA
    merged_data <- subset(merged_data,
                         !is.na(group) & 
                         !is.na(overall_survival) & 
                         !is.na(deceased) & 
                         overall_survival > 0)
    
    if (nrow(merged_data) == 0) {
      stop("No samples after grouping and NA removal")
    }
    
    # Print group sizes
    message("Sample sizes per group:")
    print(table(merged_data$group))
    
    # Calculate group statistics
    group_stats <- aggregate(counts ~ group, data = merged_data, 
                           FUN = function(x) c(mean = mean(x), 
                                             median = median(x),
                                             min = min(x),
                                             max = max(x)))
    
    message("\nGroup statistics:")
    print(group_stats)
    
    # Create subtitle with ranges
    subtitle <- paste(
      sapply(1:nrow(group_stats), function(i) {
        sprintf("%s: %.1f-%.1f (median: %.1f)",
                group_stats$group[i],
                group_stats$counts[i, "min"],
                group_stats$counts[i, "max"],
                group_stats$counts[i, "median"])
      }),
      collapse = "\n"
    )
    
    # Fit survival model
    fit <- survfit(Surv(overall_survival, deceased) ~ group, data = merged_data)
    
    # Create plot with optimized settings
    plot <- ggsurvplot(
      fit,
      data = merged_data,
      pval = TRUE,
      risk.table = TRUE,
      tables.height = 0.3,
      title = sprintf("%s survival by %s expression", cancer_type, gene_name),
      subtitle = subtitle,
      palette = "Dark2",
      xlab = "Time (days)",
      ylab = "Overall survival probability",
      font.main = c(16, "bold"),
      font.x = c(14),
      font.y = c(14),
      font.tickslab = c(12),
      risk.table.fontsize = 4,
      break.time.by = 365,  # Break x-axis by years
      legend.title = "Expression Level"
    )
    
    # Perform statistical tests
    logrank_test <- survdiff(Surv(overall_survival, deceased) ~ group, data = merged_data)
    cox_model <- coxph(Surv(overall_survival, deceased) ~ counts, data = merged_data)
    
    # Save plots with compression
    plot_base <- file.path("results", sprintf("%s_%s_expression_survival", cancer_type, gene_name))
    ggsave(paste0(plot_base, ".pdf"), plot = plot$plot, width = 10, height = 8, compress = TRUE)
    ggsave(paste0(plot_base, ".png"), plot = plot$plot, width = 10, height = 8, dpi = 300)
    
    message(sprintf("Saved plots to %s.{pdf,png}", plot_base))
    
    # Save results
    results <- list(
      plot = plot,
      fit = fit,
      data = merged_data,
      logrank = logrank_test,
      cox = cox_model,
      group_stats = group_stats
    )
    
    # Save results with compression
    saveRDS(results, file.path("results", sprintf("%s_%s_results.rds", cancer_type, gene_name)), 
            compress = TRUE)
    
    # Save summary to CSV
    summary_df <- data.frame(
      cancer_type = cancer_type,
      gene = gene_name,
      total_samples = nrow(merged_data),
      high_group = sum(merged_data$group == "High"),
      low_group = sum(merged_data$group == "Low"),
      medium_group = sum(merged_data$group == "Medium", na.rm = TRUE),
      logrank_pvalue = 1 - pchisq(logrank_test$chisq, df = length(logrank_test$n) - 1),
      cox_pvalue = summary(cox_model)$coefficients[5],
      cox_hr = exp(coef(cox_model)[1])
    )
    
    write.csv(summary_df, 
              file = file.path("results", sprintf("%s_%s_summary.csv", cancer_type, gene_name)),
              row.names = FALSE)
    
    log_memory_usage("Finished survival analysis")
    return(results)
    
  }, error = function(e) {
    message(sprintf("Error in %s analysis: %s", gene_name, e$message))
    return(NULL)
  })
}

# Main analysis with optimizations
main <- function() {
  log_memory_usage("Starting main analysis")
  
  cancer_type <- NULL
  results_list <- list()
  
  array_task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "-1"))
  if (array_task_id >= 0) {
    cancer_types <- c("ACC", "UVM", "SKCM", "LLG", "GBM")
    if (array_task_id < length(cancer_types)) {
      cancer_type <- cancer_types[array_task_id + 1]
      message(sprintf("\nProcessing cancer type: %s", cancer_type))
    }
  }
  
  if (!is.null(cancer_type)) {
    # Run analysis for both SRRM3 and SRRM4
    
    # 1. SRRM3 Expression Analysis
    tryCatch({
      message("\nRunning SRRM3 expression analysis...")
      results_list[["SRRM3"]] <- perform_survival_analysis(cancer_type, "SRRM3", "median")
    }, error = function(e) {
      message(sprintf("Error in SRRM3 analysis: %s", e$message))
    })
    
    # Clear memory between analyses
    gc()
    
    # 2. SRRM4 Expression Analysis
    tryCatch({
      message("\nRunning SRRM4 expression analysis...")
      results_list[["SRRM4"]] <- perform_survival_analysis(cancer_type, "SRRM4", "median")
    }, error = function(e) {
      message(sprintf("Error in SRRM4 analysis: %s", e$message))
    })
    
    # Save combined results with compression
    saveRDS(results_list, 
            file.path("results", sprintf("%s_combined_results.rds", cancer_type)),
            compress = TRUE)
    
    message(sprintf("\nCompleted all analyses for %s at %s", 
                   cancer_type, 
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  }
  
  log_memory_usage("Finished main analysis")
}

# Run the main function
main()
