# Load required libraries
library(survival)
library(survminer)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(recount3)
library(GenomicRanges)
library(progress)
library(futile.logger)
library(R.utils)
library(viridis)

# Define SRRM3 information
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
  ),
  transcripts = list(
    with_exon15 = "NM_001291831.2",    # 16 exons
    without_exon15 = "NM_001110199.3"   # 15 exons
  )
)

# Function to find relevant junctions for exon 15
find_exon15_junctions <- function(jxn_coords) {
  exon_start <- SRRM3_INFO$exon15$start
  exon_end <- SRRM3_INFO$exon15$end
  
  # Find upstream junction (ending at exon start)
  upstream_jxns <- which(abs(end(jxn_coords) - exon_start) <= 5)
  
  # Find downstream junction (starting at exon end)
  downstream_jxns <- which(abs(start(jxn_coords) - exon_end) <= 5)
  
  # Find exclusion junctions (those that skip exon 15)
  exclusion_jxns <- which(
    start(jxn_coords) < (exon_start - 5) & 
      end(jxn_coords) > (exon_end + 5)
  )
  
  # Combine inclusion junctions
  inclusion_jxns <- unique(c(upstream_jxns, downstream_jxns))
  
  return(list(
    inclusion = inclusion_jxns,
    exclusion = exclusion_jxns,
    details = list(
      upstream = upstream_jxns,
      downstream = downstream_jxns
    )
  ))
}

# Function to calculate PSI values
calculate_psi_values <- function(rse) {
  if(is.null(rse)) return(NULL)
  
  # Get junction counts
  junction_counts <- assay(rse)
  jxn_coords <- rowRanges(rse)
  
  # Find relevant junctions
  junctions <- find_exon15_junctions(jxn_coords)
  
  if(length(junctions$inclusion) == 0 || length(junctions$exclusion) == 0) {
    return(NULL)
  }
  
  # Calculate PSI for each sample
  psi_values <- sapply(seq_len(ncol(junction_counts)), function(i) {
    inclusion_reads <- sum(junction_counts[junctions$inclusion, i])
    exclusion_reads <- sum(junction_counts[junctions$exclusion, i])
    total_reads <- inclusion_reads + exclusion_reads
    
    if(total_reads >= 10) {  # Minimum coverage threshold
      psi <- (inclusion_reads / total_reads) * 100
      return(psi)
    } else {
      return(NA)
    }
  })
  
  return(psi_values)
}

# Add after library imports:
if (!dir.exists("cache")) {
  dir.create("cache")
}

# Add utility functions
get_cached_data <- function(cache_file) {
  if (file.exists(cache_file)) {
    message("Loading cached data from: ", cache_file)
    return(readRDS(cache_file))
  }
  return(NULL)
}

save_cached_data <- function(data, cache_file) {
  message("Saving data to cache: ", cache_file)
  saveRDS(data, cache_file)
}

# Main function for survival analysis
perform_survival_analysis <- function(cancer_type, analysis_type = "expression", gene = "SRRM3") {
  # Create cache filename based on parameters
  cache_file <- file.path("cache", 
                         paste0("survival_analysis_", cancer_type, "_", 
                               analysis_type, "_", gene, ".rds"))
  
  # Check cache first
  cached_results <- get_cached_data(cache_file)
  if (!is.null(cached_results)) {
    return(cached_results)
  }
  
  if(analysis_type == "PSI") {
    # Cache for recount3 data
    recount3_cache <- file.path("cache", paste0("recount3_", cancer_type, ".rds"))
    rse <- get_cached_data(recount3_cache)
    
    if (is.null(rse)) {
      projects <- available_projects()
      # Make sure we get exactly one project - note the changed format
      project_info <- subset(projects, project == paste0("TCGA/", cancer_type))
      if (nrow(project_info) != 1) {
        stop("Could not find unique project for ", cancer_type, 
             ". Available TCGA projects: ", 
             paste(grep("TCGA/", projects$project, value = TRUE), collapse = ", "))
      }
      
      message("Creating RSE object for ", project_info$project)
      rse <- create_rse(
        project_info[1, ],  # Ensure single row
        type = "jxn",
        jxn_format = "UNIQUE",
        verbose = TRUE
      )
      save_cached_data(rse, recount3_cache)
    }
    
    # Calculate PSI values
    psi_values <- calculate_psi_values(rse)
    
    if(is.null(psi_values)) {
      stop("Could not calculate PSI values")
    }
    
    # Create analysis dataset
    analysis_data <- data.frame(
      sample = colnames(rse),
      psi = psi_values
    )
    
    # Define high/low PSI groups using median
    median_psi <- median(analysis_data$psi, na.rm = TRUE)
    analysis_data$group <- ifelse(analysis_data$psi > median_psi, "High PSI", "Low PSI")
    
  } else {
    # Cache for expression data
    exp_cache <- file.path("cache", paste0("expression_", cancer_type, ".rds"))
    exp_data <- get_cached_data(exp_cache)
    
    if (is.null(exp_data)) {
      query_exp <- GDCquery(
        project = paste0("TCGA-", cancer_type),
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = "STAR - Counts"
      )
      
      GDCdownload(query_exp)
      exp_data <- GDCprepare(query_exp)
      save_cached_data(exp_data, exp_cache)
    }
    
    # Get expression data using TCGAbiolinks
    query_exp <- GDCquery(
      project = paste0("TCGA-", cancer_type),
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
    
    GDCdownload(query_exp)
    exp_data <- GDCprepare(query_exp)
    
    # Extract gene expression
    gene_exp <- assay(exp_data)[grep(gene, rowData(exp_data)$gene_name), ]
    
    if(length(gene_exp) == 0) {
      stop("Gene not found in dataset")
    }
    
    # Create analysis dataset
    analysis_data <- data.frame(
      sample = colnames(exp_data),
      expression = as.numeric(gene_exp)
    )
    
    # Define high/low expression groups using median
    median_exp <- median(analysis_data$expression)
    analysis_data$group <- ifelse(analysis_data$expression > median_exp, "High Expression", "Low Expression")
  }
  
  # Get clinical data
  clinical <- GDCquery_clinic(paste0("TCGA-", cancer_type))
  
  # Clean sample names to match clinical data
  analysis_data$sample <- substr(analysis_data$sample, 1, 12)
  
  # Merge with clinical data
  merged_data <- merge(
    analysis_data,
    clinical[, c("submitter_id", "days_to_death", "days_to_last_follow_up", "vital_status")],
    by.x = "sample",
    by.y = "submitter_id"
  )
  
  # Calculate survival time and status
  merged_data <- merged_data %>%
    mutate(
      survival_time = ifelse(vital_status == "Dead", 
                             days_to_death,
                             days_to_last_follow_up),
      status = ifelse(vital_status == "Dead", 1, 0)
    )
  
  # Create survival object
  surv_obj <- Surv(merged_data$survival_time, merged_data$status)
  
  # Fit survival curve
  fit <- survfit(surv_obj ~ group, data = merged_data)
  
  # Create plot
  plot <- ggsurvplot(
    fit,
    data = merged_data,
    pval = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    xlab = "Time (days)",
    ylab = "Overall Survival Probability",
    title = paste0(gene, " ", analysis_type, " Survival Analysis"),
    legend.title = ifelse(analysis_type == "PSI", "PSI Level", "Expression Level"),
    palette = c("#E7B800", "#2E9FDF")
  )
  
  results <- list(plot = plot, fit = fit, data = merged_data)
  save_cached_data(results, cache_file)
  return(results)
}

# Example usage:
# For PSI analysis:
# results <- perform_survival_analysis("BRCA", "PSI", "SRRM3")
# results$plot

# For expression analysis:
# results <- perform_survival_analysis("BRCA", "expression", "SRRM3")
# results$plot

# For PSI-based survival analysis:
psi_results <- perform_survival_analysis("BRCA", "PSI", "SRRM3")
psi_results$plot

# For expression-based survival analysis:
exp_results <- perform_survival_analysis("BRCA", "expression", "SRRM3")
exp_results$plot