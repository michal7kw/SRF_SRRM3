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

# At the top of the file, add:
if (!dir.exists("cache")) {
  dir.create("cache")
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
  
  return(clinical)
}

# Function to get expression data
get_expression_data <- function(cancer_type, gene_name, sample_types = c("Primary Tumor", "Solid Tissue Normal")) {
  # Create cache filename
  cache_file <- file.path("cache", paste0("expression_data_", cancer_type, ".rds"))
  
  # Check if cached data exists
  if (file.exists(cache_file)) {
    message("Loading cached expression data...")
    tcga_data <- readRDS(cache_file)
  } else {
    message("Downloading expression data...")
    # Build query
    query <- GDCquery(
      project = paste0("TCGA-", cancer_type),
      data.category = "Transcriptome Profiling",
      experimental.strategy = "RNA-Seq",
      workflow.type = "STAR - Counts",
      data.type = "Gene Expression Quantification",
      sample.type = sample_types,
      access = "open"
    )
    
    # Download and prepare data
    GDCdownload(query)
    tcga_data <- GDCprepare(query, summarizedExperiment = TRUE)
    
    # Cache the data
    saveRDS(tcga_data, cache_file)
  }
  
  # Extract count matrix and metadata
  count_matrix <- assay(tcga_data, "unstranded")
  gene_metadata <- as.data.frame(rowData(tcga_data))
  coldata <- as.data.frame(colData(tcga_data))
  
  # VST transformation
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = coldata,
                                design = ~ 1)
  
  # Filter low count genes
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  # Perform VST
  vsd <- vst(dds, blind = FALSE)
  expression_matrix_vst <- assay(vsd)
  
  # Process gene expression data
  gene_data <- expression_matrix_vst %>%
    as.data.frame() %>%
    rownames_to_column(var = 'gene_id') %>%
    gather(key = 'case_id', value = 'counts', -gene_id) %>%
    left_join(., gene_metadata, by = "gene_id") %>%
    filter(gene_name == gene_name)
  
  # Clean case IDs
  gene_data$case_id <- gsub('-01.*', '', gene_data$case_id)
  
  return(gene_data)
}

# Function to calculate PSI values (from previous code)
calculate_psi_values <- function(cancer_type) {
  # SRRM3 information
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
  
  # Get junction data using recount3
  projects <- available_projects()
  project_info <- subset(projects, project == paste0("TCGA_", cancer_type))
  
  # Create RSE object
  rse <- create_rse(
    project_info,
    type = "jxn",
    jxn_format = "UNIQUE",
    verbose = TRUE
  )
  
  # Calculate PSI values
  junction_counts <- assay(rse)
  jxn_coords <- rowRanges(rse)
  
  # Find relevant junctions
  exon_start <- SRRM3_INFO$exon15$start
  exon_end <- SRRM3_INFO$exon15$end
  
  inclusion_jxns <- which(
    (abs(end(jxn_coords) - exon_start) <= 5) |
      (abs(start(jxn_coords) - exon_end) <= 5)
  )
  
  exclusion_jxns <- which(
    start(jxn_coords) < (exon_start - 5) &
      end(jxn_coords) > (exon_end + 5)
  )
  
  # Calculate PSI for each sample
  psi_values <- data.frame(
    case_id = colnames(rse),
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
  
  # Clean case IDs
  psi_values$case_id <- substr(psi_values$case_id, 1, 12)
  
  return(psi_values)
}

# Main function for survival analysis
perform_survival_analysis <- function(cancer_type, analysis_type = "expression", 
                                      gene = "SRRM3", stratification = "median") {
  # Get clinical data
  clinical_data <- get_clinical_data(cancer_type)
  
  # Get analysis data based on type
  if(analysis_type == "expression") {
    analysis_data <- get_expression_data(cancer_type, gene)
    value_col <- "counts"
  } else if(analysis_type == "PSI") {
    analysis_data <- calculate_psi_values(cancer_type)
    value_col <- "psi"
  } else {
    stop("Invalid analysis type. Use 'expression' or 'PSI'")
  }
  
  # Determine stratification threshold
  if(stratification == "median") {
    median_value <- median(analysis_data[[value_col]], na.rm = TRUE)
    analysis_data$strata <- ifelse(analysis_data[[value_col]] >= median_value, "HIGH", "LOW")
  } else if(stratification == "quartile") {
    quartiles <- quantile(analysis_data[[value_col]], probs = c(0.25, 0.75), na.rm = TRUE)
    analysis_data$strata <- case_when(
      analysis_data[[value_col]] <= quartiles[1] ~ "LOW",
      analysis_data[[value_col]] >= quartiles[2] ~ "HIGH",
      TRUE ~ "INTERMEDIATE"
    )
  }
  
  # Merge with clinical data
  merged_data <- merge(analysis_data, clinical_data, 
                       by.x = 'case_id', by.y = 'submitter_id')
  
  # Create survival object and fit
  fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = merged_data)
  
  # Create plot
  surv_plot <- ggsurvplot(
    fit,
    data = merged_data,
    pval = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    xlab = "Time (days)",
    ylab = "Overall Survival Probability",
    title = paste(gene, analysis_type, "Survival Analysis"),
    legend.title = paste(analysis_type, "Level"),
    palette = c("HIGH" = "#E7B800", "LOW" = "#2E9FDF", "INTERMEDIATE" = "#868686")
  )
  
  # Perform statistical tests
  logrank_test <- survdiff(Surv(overall_survival, deceased) ~ strata, data = merged_data)
  
  # Cox proportional hazards model
  cox_model <- coxph(Surv(overall_survival, deceased) ~ get(value_col), data = merged_data)
  
  return(list(
    plot = surv_plot,
    fit = fit,
    data = merged_data,
    logrank = logrank_test,
    cox = cox_model
  ))
}

# Example usage:
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