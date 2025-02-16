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
  library(gridExtra)
  library(grid)
  library(cowplot)
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
  
  # Construct path to clinical data directory
  clinical_dir <- file.path("../../TCGA_clinical", 
                           sprintf("clinical.project-tcga-%s.2025-02-12", 
                                 tolower(cancer_type)))
  
  # Check if directory exists
  if (!dir.exists(clinical_dir)) {
    stop(sprintf("Clinical data directory not found for %s", cancer_type))
  }
  
  # Read clinical.tsv and follow_up.tsv files
  clinical <- read.delim(file.path(clinical_dir, "clinical.tsv"), 
                        sep="\t", stringsAsFactors=FALSE, na.strings=c("'--", "--", "[Not Available]", "[Not Applicable]"))
  follow_up <- read.delim(file.path(clinical_dir, "follow_up.tsv"), 
                         sep="\t", stringsAsFactors=FALSE, na.strings=c("'--", "--", "[Not Available]", "[Not Applicable]"))
  
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
      gender = gender,
      tumor_grade = tumor_grade,  # Include tumor grade
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
  
  # Check if tumor grade information is available
  grade_counts <- table(clinical_clean$tumor_grade)
  message("\nTumor grade distribution:")
  print(grade_counts)
  
  # Try to map different grading systems to G2/G3
  clinical_clean <- clinical_clean %>%
    mutate(
      tumor_grade = case_when(
        # Standard grading
        tumor_grade %in% c("G2", "G3") ~ tumor_grade,
        # Alternative grading systems
        tumor_grade %in% c("Grade II", "grade 2", "2") ~ "G2",
        tumor_grade %in% c("Grade III", "grade 3", "3") ~ "G3",
        # For some cancers, intermediate/high might correspond to G2/G3
        tolower(tumor_grade) == "intermediate" ~ "G2",
        tolower(tumor_grade) == "high" ~ "G3",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(tumor_grade))
  
  # Print summary statistics
  message(sprintf("\nProcessed clinical data summary:"))
  message(sprintf("- Total patients: %d", nrow(clinical_clean)))
  message(sprintf("- G2 patients: %d", sum(clinical_clean$tumor_grade == "G2", na.rm = TRUE)))
  message(sprintf("- G3 patients: %d", sum(clinical_clean$tumor_grade == "G3", na.rm = TRUE)))
  message(sprintf("- Male patients: %d", sum(tolower(clinical_clean$gender) == "male", na.rm = TRUE)))
  message(sprintf("- Female patients: %d", sum(tolower(clinical_clean$gender) == "female", na.rm = TRUE)))
  message(sprintf("- Patients with death events: %d", sum(clinical_clean$deceased)))
  message(sprintf("- Median survival time: %.1f days", 
                 median(clinical_clean$overall_survival, na.rm = TRUE)))
  
  if (nrow(clinical_clean) == 0) {
    message(sprintf("Warning: No patients with G2/G3 grading found for %s", cancer_type))
    message("Original grade values found:")
    print(grade_counts)
  }
  
  return(clinical_clean)
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
    expression = assay(rse_gene, "raw_counts")[gene_idx, ]
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
# Helper Functions for Distribution Analysis
#####################################################################
determine_optimal_grouping <- function(values) {
  # Calculate distribution characteristics
  q1 <- quantile(values, 0.25)
  q2 <- quantile(values, 0.5)
  q3 <- quantile(values, 0.75)
  iqr <- q3 - q1
  
  # Calculate density
  d <- density(values)
  peaks <- findPeaks(d$y)
  valleys <- findValleys(d$y)
  
  # Criteria for using two groups:
  # 1. Clear bimodality (2 peaks with significant valley)
  # 2. Large gaps in distribution
  # 3. Small middle group
  
  middle_group_size <- sum(values > q1 & values < q3)
  total_size <- length(values)
  middle_group_proportion <- middle_group_size / total_size
  
  use_two_groups <- length(peaks) >= 2 || 
                    middle_group_proportion < 0.3 ||
                    (q3 - q1) > 2 * mad(values)
  
  message(sprintf("\n=== Group Distribution Analysis ==="))
  message(sprintf("Number of density peaks: %d", length(peaks)))
  message(sprintf("Middle group proportion: %.2f", middle_group_proportion))
  message(sprintf("IQR/MAD ratio: %.2f", (q3 - q1) / mad(values)))
  message(sprintf("Selected grouping: %s", if(use_two_groups) "two groups" else "three groups"))
  
  return(use_two_groups)
}

# Helper functions for peak detection
findPeaks <- function(x) {
  # Find local maxima
  peaks <- which(diff(sign(diff(c(-Inf, x, -Inf)))) == -2)
  return(peaks)
}

findValleys <- function(x) {
  # Find local minima
  valleys <- which(diff(sign(diff(c(Inf, x, Inf)))) == 2)
  return(valleys)
}

#####################################################################
# High SRRM3 PSI Analysis with Grade Stratification
#####################################################################
perform_high_srrm3_psi_analysis_by_grade <- function(cancer_type, expression_threshold = 0.75,
                                                    grouping_method = "quartile") {
  tryCatch({
    # Get expression data first to include in summary
    expression_data <- get_expression_data(cancer_type, "SRRM3")
    expression_cutoff <- quantile(expression_data$expression, probs = expression_threshold, na.rm = TRUE)
    message(sprintf("Expression cutoff (%.2f percentile): %.2f", expression_threshold * 100, expression_cutoff))
    
    high_expr_samples <- expression_data %>%
      filter(expression >= expression_cutoff) %>%
      pull(case_id)
    
    message("Number of high expression samples: ", length(high_expr_samples))
    
    if (length(high_expr_samples) < 10) {
      stop("Insufficient samples with high SRRM3 expression")
    }
    
    clinical_data <- get_clinical_data(cancer_type)
    message("Number of clinical samples: ", nrow(clinical_data))
    
    psi_data <- get_psi_data(cancer_type, "SRRM3", high_expr_samples)
    message("Number of PSI samples: ", nrow(psi_data))
    
    # Merge all data (clinical, PSI, and expression)
    merged_data <- clinical_data %>% 
      mutate(merge_id = substr(case_id, 1, 12)) %>%
      inner_join(psi_data %>% mutate(merge_id = substr(case_id, 1, 12)), by = "merge_id") %>%
      inner_join(expression_data %>% mutate(merge_id = substr(case_id, 1, 12)), by = "merge_id")
    
    message("Number of samples after merging: ", nrow(merged_data))
    
    if (nrow(merged_data) == 0) {
      stop("No samples after merging clinical and molecular data")
    }
    
    # Determine optimal grouping based on PSI distribution
    message("\n=== Grouping Samples ===")
    use_two_groups <- determine_optimal_grouping(merged_data$psi)
    
    if (use_two_groups) {
      # Use median split for two groups
      median_val <- median(merged_data$psi, na.rm = TRUE)
      message(sprintf("Using two-group split at median %.2f", median_val))
      merged_data$group <- if_else(merged_data$psi > median_val, "High", "Low")
    } else {
      # Use quartiles for three groups
      quartiles <- quantile(merged_data$psi, probs = c(0.25, 0.75), na.rm = TRUE)
      message(sprintf("Using three-group split at quartiles: Q1=%.2f, Q3=%.2f", 
                     quartiles[1], quartiles[2]))
      merged_data$group <- case_when(
        merged_data$psi <= quartiles[1] ~ "Low",
        merged_data$psi >= quartiles[2] ~ "High",
        TRUE ~ "Medium"
      )
    }
    
    merged_data <- merged_data %>% 
      filter(!is.na(group))
    
    # Create results directory
    results_dir <- file.path("results", "high_srrm3_analysis_by_grade")
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Create a list to store all plots
    plot_list <- list()
    
    # Function to create survival plot
    create_survival_plot <- function(data, title, grade) {
      # Filter data for specific grade
      grade_data <- data %>% filter(tumor_grade == grade)
      
      if (nrow(grade_data) == 0) {
        message(sprintf("No data available for %s in %s", grade, title))
        return(NULL)
      }
      
      # Calculate counts
      counts <- grade_data %>%
        group_by(group) %>%
        summarize(n = n(), .groups = "drop") %>%
        mutate(label = sprintf("%s (n=%d)", group, n)) %>%
        pull(label) %>%
        paste(collapse = ", ")
      
      # Create survival fit
      fit <- survfit(Surv(time = overall_survival, event = deceased) ~ group,
                    data = grade_data)
      
      # Create plot
      plot <- ggsurvplot(
        fit, data = grade_data,
        pval = TRUE,
        risk.table = FALSE,  # Remove risk table
        title = sprintf("%s - %s\n%s", title, grade, counts),
        legend.title = "PSI Group",
        xlab = "Time (days)",
        break.time.by = 365.25,  # Break x-axis by years
        font.main = c(14, "bold"),  # Increase title font size
        font.x = c(12),            # Increase x-axis font size
        font.y = c(12),            # Increase y-axis font size
        font.legend = c(10),       # Adjust legend font size
        legend.labs = levels(as.factor(grade_data$group)),  # Clean legend labels
        palette = "npg",           # Use Nature publishing group color palette
        ggtheme = theme_bw() +     # Use black and white theme
                 theme(plot.title = element_text(hjust = 0.5),  # Center title
                       legend.position = "right")
      )
      
      return(plot)
    }
    
    # Create plots for all patients and by gender for each grade
    plots_g2 <- list()
    plots_g3 <- list()
    
    # All patients - G2
    plot_all_g2 <- create_survival_plot(
      merged_data,
      sprintf("%s - All Patients", cancer_type),
      "G2"
    )
    if (!is.null(plot_all_g2)) plots_g2$all <- plot_all_g2
    
    # All patients - G3
    plot_all_g3 <- create_survival_plot(
      merged_data,
      sprintf("%s - All Patients", cancer_type),
      "G3"
    )
    if (!is.null(plot_all_g3)) plots_g3$all <- plot_all_g3
    
    # Male patients - G2
    plot_male_g2 <- create_survival_plot(
      merged_data %>% filter(tolower(gender) == "male"),
      sprintf("%s - Male Patients", cancer_type),
      "G2"
    )
    if (!is.null(plot_male_g2)) plots_g2$male <- plot_male_g2
    
    # Male patients - G3
    plot_male_g3 <- create_survival_plot(
      merged_data %>% filter(tolower(gender) == "male"),
      sprintf("%s - Male Patients", cancer_type),
      "G3"
    )
    if (!is.null(plot_male_g3)) plots_g3$male <- plot_male_g3
    
    # Female patients - G2
    plot_female_g2 <- create_survival_plot(
      merged_data %>% filter(tolower(gender) == "female"),
      sprintf("%s - Female Patients", cancer_type),
      "G2"
    )
    if (!is.null(plot_female_g2)) plots_g2$female <- plot_female_g2
    
    # Female patients - G3
    plot_female_g3 <- create_survival_plot(
      merged_data %>% filter(tolower(gender) == "female"),
      sprintf("%s - Female Patients", cancer_type),
      "G3"
    )
    if (!is.null(plot_female_g3)) plots_g3$female <- plot_female_g3
    
    # Save G2 plots if any exist
    if (length(plots_g2) > 0) {
      combined_plot_g2 <- arrange_ggsurvplots(
        plots_g2,
        ncol = 1, nrow = length(plots_g2),
        print = FALSE
      )
      
      ggsave(
        file.path(results_dir, sprintf("%s_combined_survival_G2.pdf", cancer_type)),
        combined_plot_g2,
        width = 8,
        height = 6 * length(plots_g2)
      )
      ggsave(
        file.path(results_dir, sprintf("%s_combined_survival_G2.png", cancer_type)),
        combined_plot_g2,
        width = 8,
        height = 6 * length(plots_g2),
        dpi = 300
      )
    }
    
    # Save G3 plots if any exist
    if (length(plots_g3) > 0) {
      combined_plot_g3 <- arrange_ggsurvplots(
        plots_g3,
        ncol = 1, nrow = length(plots_g3),
        print = FALSE
      )
      
      ggsave(
        file.path(results_dir, sprintf("%s_combined_survival_G3.pdf", cancer_type)),
        combined_plot_g3,
        width = 8,
        height = 6 * length(plots_g3)
      )
      ggsave(
        file.path(results_dir, sprintf("%s_combined_survival_G3.png", cancer_type)),
        combined_plot_g3,
        width = 8,
        height = 6 * length(plots_g3),
        dpi = 300
      )
    }
    
    # Create summary statistics
    summary_stats <- data.frame(
      cancer_type = cancer_type,
      total_samples = nrow(merged_data),
      g2_samples = sum(merged_data$tumor_grade == "G2"),
      g3_samples = sum(merged_data$tumor_grade == "G3"),
      male_g2 = sum(merged_data$tumor_grade == "G2" & tolower(merged_data$gender) == "male"),
      male_g3 = sum(merged_data$tumor_grade == "G3" & tolower(merged_data$gender) == "male"),
      female_g2 = sum(merged_data$tumor_grade == "G2" & tolower(merged_data$gender) == "female"),
      female_g3 = sum(merged_data$tumor_grade == "G3" & tolower(merged_data$gender) == "female"),
      expression_threshold = expression_cutoff,
      expression_threshold_percentile = expression_threshold,
      grouping_method = if(use_two_groups) "two-level" else "three-level"
    )
    
    write.csv(summary_stats,
              file.path(results_dir, sprintf("%s_grade_summary.csv", cancer_type)),
              row.names = FALSE)
    
    return(list(
      plots_g2 = if(length(plots_g2) > 0) plots_g2 else NULL,
      plots_g3 = if(length(plots_g3) > 0) plots_g3 else NULL,
      data = merged_data,
      summary = summary_stats
    ))
    
  }, error = function(e) {
    message(sprintf("Error in high SRRM3 PSI analysis for %s: %s", cancer_type, e$message))
    return(NULL)
  })
}

#####################################################################
# Main Function
#####################################################################
main <- function() {
    # Get cancer type from environment variable
    cancer_type <- Sys.getenv("CANCER_TYPE")
    if (nchar(cancer_type) == 0) {
        stop("CANCER_TYPE environment variable not set")
    }
    
    message(sprintf("\nAnalyzing %s with grade stratification", cancer_type))
    
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
    
    results <- perform_high_srrm3_psi_analysis_by_grade(cancer_type, expression_threshold = 0.75)
    
    if (is.null(results)) {
        message("No results generated")
        quit(status = 1)
    }
    
    return(results)
}

# Execute main function if script is run directly
if (sys.nframe() == 0) {
    main()
} 