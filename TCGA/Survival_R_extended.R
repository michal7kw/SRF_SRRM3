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

# At the top of the file, add:
if (!dir.exists("cache")) {
  dir.create("cache")
}

clear_cache <- function(cancer_type = NULL, gene = NULL) {
  if (!is.null(cancer_type) && !is.null(gene)) {
    # Remove specific cache file
    cache_file <- file.path("cache", paste0("expression_data_", cancer_type, "_", gene, ".rds"))
    if (file.exists(cache_file)) {
      unlink(cache_file)
      message("Removed cache file: ", cache_file)
    }
  } else {
    # Remove all cache files
    if (dir.exists("cache")) {
      unlink("cache/*")
      message("Cleared all cache files")
    }
  }
}

# clear_cache("BRCA", "SRRM3")

# Add this helper function at the top
find_tcga_project <- function(cancer_type) {
  message("Getting available projects...")
  projects <- available_projects()
  
  # Print all available projects for debugging
  message("All available projects with 'tcga' in name:")
  tcga_projects <- grep("tcga", projects$project, value=TRUE, ignore.case=TRUE)
  print(head(tcga_projects))
  
  # Look for project in data_sources/tcga
  project_info <- subset(
    projects,
    file_source == "tcga" &
    project_type == "data_sources"
  )
  
  message("Available TCGA projects:")
  print(head(project_info))
  
  # Find specific cancer type
  cancer_matches <- grep(tolower(cancer_type), tolower(project_info$project), value=TRUE)
  
  if (length(cancer_matches) == 0) {
    stop(paste("No matching TCGA project found for", cancer_type))
  }
  
  message("Found matching project: ", cancer_matches[1])
  return(subset(project_info, project == cancer_matches[1]))
}

# Function to get and process clinical data
get_clinical_data <- function(cancer_type) {
  clinical <- GDCquery_clinic(paste0("TCGA-", cancer_type))
  
  # Process survival information
  clinical <- clinical %>%
    mutate(
      deceased = ifelse(vital_status == "Alive", FALSE, TRUE),
      overall_survival = ifelse(vital_status == "Alive",
                                days_to_last_follow_up,
                                days_to_death)
    )
  
  # Add debug information
  message(sprintf("Clinical data has %d patients", nrow(clinical)))
  message("Sample of clinical data case IDs:")
  print(head(clinical$submitter_id))
  
  return(clinical)
}

# Function to get expression data
get_expression_data <- function(cancer_type, gene_name) {
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
  
  # Make sure we get exactly one project
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
    verbose = TRUE
  )
  
  # Get gene expression
  gene_idx <- which(rowData(rse_gene)$gene_name == gene_name)
  if (length(gene_idx) == 0) {
    stop(paste("Gene", gene_name, "not found in dataset"))
  }
  
  # Get the TCGA barcodes from tcga.tcga_barcode column
  tcga_barcodes <- colData(rse_gene)$tcga.tcga_barcode
  if (is.null(tcga_barcodes)) {
    message("Available colData columns:")
    print(colnames(colData(rse_gene)))
    stop("Could not find tcga.tcga_barcode in the data")
  }
  
  # Create data frame with expression data
  gene_data <- data.frame(
    case_id = tcga_barcodes,
    counts = assay(rse_gene)[gene_idx, ],
    stringsAsFactors = FALSE
  )
  
  # Clean case IDs - make sure they match TCGA format (TCGA-XX-XXXX)
  gene_data$case_id <- substr(gene_data$case_id, 1, 12)
  
  # Add debug information
  message("Sample of expression data case IDs after cleaning:")
  print(head(gene_data$case_id))
  
  # Cache results
  saveRDS(gene_data, cache_file)
  return(gene_data)
}

# Add this function for handling GDC queries with retries
query_gdc_with_retry <- function(cancer_type, data_type) {
  message("Querying GDC data with retry...")
  
  # Set timeout options
  old_timeout <- getOption("timeout")
  options(timeout = 300) # Set timeout to 5 minutes
  
  while(attempt <= max_attempts) {
    tryCatch({
      if(data_type == "Splice Junction Quantification") {
        # For junction data, use RNA-Seq Alignments
        query <- GDCquery(
          project = paste0("TCGA-", cancer_type),
          data.category = "Transcriptome Profiling",
          experimental.strategy = "RNA-Seq",
          workflow.type = "HTSeq - Counts",  # Changed from STAR - Counts
          data.type = "Gene Expression Quantification",  # Changed to get alignment data
          sample.type = c("Primary Tumor", "Solid Tissue Normal"),
          access = "open"
        )
      } else {
        # For other data types, use original parameters
        query <- GDCquery(
          project = paste0("TCGA-", cancer_type),
          data.category = "Transcriptome Profiling",
          experimental.strategy = "RNA-Seq",
          workflow.type = "STAR - Counts",
          data.type = data_type,
          sample.type = "Primary Tumor",
          access = "open"
        )
      }
      
      # Get sample info
      samples <- getResults(query)
      message(sprintf("Found %d samples", nrow(samples)))
      
      # For testing, limit to first 100 samples
      if(nrow(samples) > 100) {
        message("Limiting to first 100 samples for testing")
        query$results <- query$results[1:100,]
      }
      
      return(query)
    }, error = function(e) {
      message(sprintf("Attempt %d failed: %s", attempt, e$message))
      if(attempt == max_attempts) {
        stop("Max retry attempts reached")
      }
      Sys.sleep(5)
      attempt <<- attempt + 1
    })
  }
}

# Function to calculate PSI values
calculate_psi_values <- function(cancer_type) {
  cache_file <- file.path("cache", paste0("psi_values_", cancer_type, ".rds"))
  
  if (file.exists(cache_file)) {
    message("Loading cached PSI values...")
    return(readRDS(cache_file))
  }
  
  # Define SRRM3 information first
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
  
  # Make sure we get exactly one project - match the project name directly
  project_info <- subset(projects, 
                        project == toupper(cancer_type) &  # Projects are in uppercase (BRCA)
                        file_source == "tcga" &
                        project_type == "data_sources")
  
  if (nrow(project_info) != 1) {
    message("\nAvailable TCGA projects:")
    print(subset(projects, file_source == "tcga" & project_type == "data_sources"))
    stop("Could not find unique project for ", cancer_type)
  }
  
  message(sprintf("Found project: %s", project_info$project[1]))
  
  # Create RSE object with correct parameters from working implementation
  rse <- create_rse(
    project_info[1, ],  # Ensure single row
    type = "jxn",
    jxn_format = "UNIQUE",
    verbose = TRUE
  )
  
  # Get junction data
  junction_counts <- assay(rse)
  jxn_coords <- rowRanges(rse)
  
  # Find relevant junctions
  message("Finding relevant junctions...")
  chr7_jxns <- which(seqnames(jxn_coords) == "chr7")
  srrm3_region <- which(
    start(jxn_coords) >= (SRRM3_INFO$exon15$start - 10000) &
    end(jxn_coords) <= (SRRM3_INFO$exon15$end + 10000) &
    seqnames(jxn_coords) == "chr7"
  )
  
  message(sprintf("Found %d junctions in SRRM3 region", length(srrm3_region)))
  
  # Get relevant junctions
  jxn_coords <- jxn_coords[srrm3_region]
  junction_counts <- junction_counts[srrm3_region, ]
  
  # Find inclusion and exclusion junctions
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
  
  # Calculate PSI values with correct IDs
  psi_values <- data.frame(
    case_id = colData(rse)$tcga.tcga_barcode,  # Use tcga.tcga_barcode instead of external_id
    psi = sapply(seq_len(ncol(junction_counts)), function(i) {
      inclusion_reads <- sum(junction_counts[inclusion_jxns, i])
      exclusion_reads <- sum(junction_counts[exclusion_jxns, i])
      total_reads <- inclusion_reads + exclusion_reads
      
      if(total_reads >= 10) {
        return((inclusion_reads / total_reads) * 100)
      } else {
        return(NA)
      }
    })
  )
  
  # Clean case IDs to match TCGA format
  psi_values$case_id <- substr(psi_values$case_id, 1, 12)
  
  message("Caching PSI values...")
  saveRDS(psi_values, cache_file)
  return(psi_values)
}

# Main function for survival analysis
perform_survival_analysis <- function(cancer_type, analysis_type = "expression", 
                                    gene = "SRRM3", stratification = "median") {
  message(sprintf("\nStarting survival analysis for %s cancer, %s analysis of %s", 
                 cancer_type, analysis_type, gene))
  
  # Get clinical data
  message("Getting clinical data...")
  clinical_data <- get_clinical_data(cancer_type)
  message(sprintf("Retrieved clinical data for %d patients", nrow(clinical_data)))
  
  # Get analysis data based on type
  if(analysis_type == "expression") {
    message("Getting expression data...")
    analysis_data <- get_expression_data(cancer_type, gene)
    value_col <- "counts"
  } else if(analysis_type == "PSI") {
    message("Getting PSI values...")
    analysis_data <- calculate_psi_values(cancer_type)
    if(is.null(analysis_data)) {
      stop("Failed to get PSI values")
    }
    value_col <- "psi"
  } else {
    stop("Invalid analysis type. Use 'expression' or 'PSI'")
  }
  
  message(sprintf("Retrieved %s data for %d samples", 
                 analysis_type, nrow(analysis_data)))
  
  # Ensure case_id is unique before merging using dplyr operations
  message("Processing analysis data...")
  analysis_data <- analysis_data %>%
    group_by(case_id) %>%
    summarise(!!sym(value_col) := mean(!!sym(value_col), na.rm = TRUE)) %>%
    ungroup()
  
  message(sprintf("Processed data has %d unique cases", nrow(analysis_data)))
  
  # Determine stratification threshold
  message("Calculating stratification...")
  if(stratification == "median") {
    median_value <- median(analysis_data[[value_col]], na.rm = TRUE)
    message(sprintf("Median value: %.2f", median_value))
    analysis_data$strata <- ifelse(analysis_data[[value_col]] >= median_value, "HIGH", "LOW")
  } else if(stratification == "quartile") {
    quartiles <- quantile(analysis_data[[value_col]], probs = c(0.25, 0.75), na.rm = TRUE)
    message(sprintf("Quartiles: Q1=%.2f, Q3=%.2f", quartiles[1], quartiles[2]))
    analysis_data$strata <- case_when(
      analysis_data[[value_col]] <= quartiles[1] ~ "LOW",
      analysis_data[[value_col]] >= quartiles[2] ~ "HIGH",
      TRUE ~ "INTERMEDIATE"
    )
  }
  
  # Print strata distribution
  strata_counts <- table(analysis_data$strata)
  message("Strata distribution:")
  print(strata_counts)
  
  # Merge with clinical data - add debugging
  message("\nBefore merging:")
  message("Clinical data dimensions: ", nrow(clinical_data), " x ", ncol(clinical_data))
  message("Analysis data dimensions: ", nrow(analysis_data), " x ", ncol(analysis_data))
  message("\nSample clinical IDs:")
  print(head(clinical_data$submitter_id))
  message("\nSample analysis IDs:")
  print(head(analysis_data$case_id))
  
  # Make sure case IDs are in the same format
  clinical_data$submitter_id <- substr(clinical_data$submitter_id, 1, 12)
  analysis_data$case_id <- substr(analysis_data$case_id, 1, 12)
  
  # Merge with clinical data
  message("Merging with clinical data...")
  merged_data <- merge(analysis_data, clinical_data, 
                      by.x = 'case_id', by.y = 'submitter_id',
                      all = FALSE)  # Only keep matching cases
  
  message(sprintf("Final merged dataset has %d patients", nrow(merged_data)))
  
  if (nrow(merged_data) == 0) {
    stop("No matching patients found after merging clinical and expression data!")
  }
  
  # Create survival object and fit
  message("Fitting survival model...")
  fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = merged_data)
  
  # Add summary statistics calculation
  message("Calculating group statistics...")
  group_stats <- merged_data %>%
    group_by(strata) %>%
    summarise(
      mean_val = mean(!!sym(value_col), na.rm = TRUE),
      median_val = median(!!sym(value_col), na.rm = TRUE),
      sd_val = sd(!!sym(value_col), na.rm = TRUE),
      min_val = min(!!sym(value_col), na.rm = TRUE),
      max_val = max(!!sym(value_col), na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Create subtitle with ranges
  subtitle <- paste(
    sapply(1:nrow(group_stats), function(i) {
      group <- group_stats$strata[i]
      if(analysis_type == "PSI") {
        sprintf("%s: PSI %.1f-%.1f%% (median: %.1f%%)",
                group,
                group_stats$min_val[i],
                group_stats$max_val[i],
                group_stats$median_val[i])
      } else {
        sprintf("%s: Expr %.1f-%.1f (median: %.1f)",
                group,
                group_stats$min_val[i],
                group_stats$max_val[i],
                group_stats$median_val[i])
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
    title = paste(gene, analysis_type, "Survival Analysis"),
    subtitle = subtitle,  # Add the subtitle with ranges
    legend.title = paste(analysis_type, "Level"),
    palette = c("HIGH" = "#E7B800", "LOW" = "#2E9FDF", "INTERMEDIATE" = "#868686"),
    font.subtitle = c(10, "plain"),  # Adjust subtitle font size
    break.time.by = 365  # Break x-axis by years
  )
  
  # Perform statistical tests
  message("Performing statistical tests...")
  logrank_test <- survdiff(Surv(overall_survival, deceased) ~ strata, data = merged_data)
  cox_model <- coxph(Surv(overall_survival, deceased) ~ get(value_col), data = merged_data)
  
  message("Analysis complete!")
  
  return(list(
    plot = surv_plot,
    fit = fit,
    data = merged_data,
    logrank = logrank_test,
    cox = cox_model,
    group_stats = group_stats  # Add group statistics to returned results
  ))
}

# For SRRM3 expression analysis:
results_exp <- perform_survival_analysis("BRCA", "expression", "SRRM3", "median")
results_exp$plot
# View results
summary(results_exp$cox)  # View Cox model results
results_exp$logrank  # View log-rank test results

# For SRRM3 PSI analysis:
results_psi <- perform_survival_analysis("BRCA", "PSI", "SRRM3", "quartile")
results_psi$plot

# For SRRM4 expression analysis:
results_srrm4 <- perform_survival_analysis("BRCA", "expression", "SRRM4", "median")
results_srrm4$plot

# Create results directory if it doesn't exist
if (!dir.exists("results_ext")) {
  dir.create("results_ext")
}

# Function to safely save plot to PDF
save_plot_to_pdf <- function(plot, filename, width = 10, height = 8) {
  pdf(filename, width = width, height = height, onefile = TRUE)
  tryCatch({
    print(plot)
  }, finally = {
    dev.off()  # Ensure device is closed even if there's an error
  })
}

# Save SRRM3 expression plot
if (!is.null(results_exp)) {
  save_plot_to_pdf(
    results_exp$plot,
    file.path("results_ext", "BRCA_SRRM3_expression_survival.pdf")
  )
  saveRDS(results_exp, file.path("results_ext", "BRCA_SRRM3_expression_results.rds"))
}

# Save SRRM3 PSI plot
if (!is.null(results_psi)) {
  save_plot_to_pdf(
    results_psi$plot,
    file.path("results_ext", "BRCA_SRRM3_PSI_survival.pdf")
  )
  saveRDS(results_psi, file.path("results_ext", "BRCA_SRRM3_PSI_results.rds"))
}

# Save SRRM4 expression plot
if (!is.null(results_srrm4)) {
  save_plot_to_pdf(
    results_srrm4$plot,
    file.path("results_ext", "BRCA_SRRM4_expression_survival.pdf")
  )
  saveRDS(results_srrm4, file.path("results_ext", "BRCA_SRRM4_expression_results.rds"))
}

message("All plots and results have been saved to the results_ext directory")

# Function to safely save plot to PNG
save_plot_to_png <- function(plot, filename, width = 10, height = 8) {
  tryCatch({
    png(filename, width = width*100, height = height*100, res = 100)
    print(plot)
    dev.off()
    message("Successfully saved plot to: ", filename)
  }, error = function(e) {
    warning("Failed to save plot to ", filename, ": ", e$message)
  })
}

# Save SRRM3 expression results
if (!is.null(results_exp)) {
  message("\nSaving SRRM3 expression analysis results...")
  save_plot_to_png(
    results_exp$plot,
    file.path("results_ext", "BRCA_SRRM3_expression_survival.png")
  )
}

# Save SRRM3 PSI results
if (!is.null(results_psi)) {
  message("\nSaving SRRM3 PSI analysis results...")
  save_plot_to_png(
    results_psi$plot,
    file.path("results_ext", "BRCA_SRRM3_PSI_survival.png")
  )
}

# Save SRRM4 expression results
if (!is.null(results_srrm4)) {
  message("\nSaving SRRM4 expression analysis results...")
  save_plot_to_png(
    results_srrm4$plot,
    file.path("results_ext", "BRCA_SRRM4_expression_survival.png")
  )
}

message("\nAll plots have been saved to the results_ext directory")

# Function to print summary of saved files
print_saved_files_summary <- function() {
  files <- list.files("results_ext", pattern = "\\.png$", full.names = TRUE)
  if (length(files) > 0) {
    message("\nSaved files:")
    for (file in files) {
      message("- ", basename(file))
    }
  } else {
    message("\nNo files were saved in the results_ext directory")
  }
}

# Print summary of saved files
print_saved_files_summary()