# Load additional required libraries
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
library(mclust)      # For mixture modeling
library(timeROC)     # For time-dependent ROC
library(plotly)      # For interactive plots
library(cmprsk)      # For competing risks
library(rms)         # For advanced survival modeling
library(forestplot)  # For forest plots

# Add after library imports:
if (!dir.exists("cache")) {
  dir.create("cache")
}

# Add these utility functions
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

# Add this debugging function after the utility functions
debug_available_projects <- function() {
  projects <- available_projects()
  message("\nAll available projects:")
  print(head(projects))
  
  message("\nProjects containing 'TCGA' (case-insensitive):")
  tcga_projects <- grep("TCGA", projects$project, value = TRUE, ignore.case = TRUE)
  print(tcga_projects)
  
  message("\nProjects containing 'breast' (case-insensitive):")
  breast_projects <- grep("breast", projects$project, value = TRUE, ignore.case = TRUE)
  print(breast_projects)
  
  return(projects)
}

# Simplify find_tcga_project function
find_tcga_project <- function(cancer_type) {
  projects <- available_projects()
  
  # Try just the cancer type code (e.g., "BRCA", "LUAD", etc.)
  project_info <- subset(projects, project == cancer_type)
  
  if (nrow(project_info) == 1) {
    message("Found project: ", project_info$project)
    return(project_info)
  }
  
  # If not found, show available projects for debugging
  message("\nAll available projects:")
  print(head(projects))
  
  stop("Could not find project for ", cancer_type, 
       ". Please use the cancer type code directly (e.g., 'BRCA' for breast cancer)")
}

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

# New function for advanced PSI filtering
filter_psi_data <- function(rse, psi_values, min_coverage = 10, min_quality = 20, 
                            min_samples_per_group = 30) {
  # Get junction quality scores
  quality_scores <- rowData(rse)$score  # Adjust based on actual quality metric available
  
  # Filter based on quality
  high_quality_junctions <- quality_scores >= min_quality
  
  # Calculate technical variance
  tech_var <- apply(assay(rse), 2, function(x) var(log2(x + 1)))
  high_var_samples <- tech_var > quantile(tech_var, 0.95)
  
  # Apply filters
  filtered_psi <- psi_values
  filtered_psi[high_var_samples] <- NA
  
  # Ensure minimum group sizes
  valid_samples <- sum(!is.na(filtered_psi)) >= (2 * min_samples_per_group)
  
  if (!valid_samples) {
    warning("Insufficient samples after filtering")
    return(NULL)
  }
  
  return(filtered_psi)
}

# Function for advanced grouping
determine_groups <- function(values, method = "quartile") {
  if (method == "quartile") {
    cuts <- quantile(values, probs = c(0.25, 0.75), na.rm = TRUE)
    groups <- case_when(
      values <= cuts[1] ~ "Low",
      values >= cuts[2] ~ "High",
      TRUE ~ "Medium"
    )
  } else if (method == "mixture") {
    mod <- Mclust(values, G = 2)
    groups <- ifelse(mod$classification == 1, "Low", "High")
  } else if (method == "median") {
    median_val <- median(values, na.rm = TRUE)
    groups <- ifelse(values > median_val, "High", "Low")
  }
  return(groups)
}

# Function for advanced statistical analysis
perform_advanced_statistics <- function(data, time_col, status_col, 
                                        predictor_col, clinical_vars) {
  # Prepare survival data
  surv_obj <- Surv(data[[time_col]], data[[status_col]])
  
  # Univariate Cox model
  univ_cox <- coxph(surv_obj ~ data[[predictor_col]])
  
  # Multivariate Cox model
  formula_str <- paste("surv_obj ~", 
                       paste(c(predictor_col, clinical_vars), collapse = " + "))
  multi_cox <- coxph(as.formula(formula_str), data = data)
  
  # Time-dependent coefficients test
  zph_test <- cox.zph(multi_cox)
  
  # C-index calculation
  c_index <- concordance(multi_cox)
  
  # Competing risks analysis (if death causes are available)
  if ("death_cause" %in% names(data)) {
    cr_data <- crr(data[[time_col]], 
                   data$death_cause, 
                   data[c(predictor_col, clinical_vars)])
  } else {
    cr_data <- NULL
  }
  
  return(list(
    univariate = univ_cox,
    multivariate = multi_cox,
    proportional_hazards_test = zph_test,
    concordance = c_index,
    competing_risks = cr_data
  ))
}

# Function for advanced visualizations
create_advanced_plots <- function(data, survival_fit, stats_results) {
  # Basic survival plot with enhanced features
  surv_plot <- ggsurvplot(
    survival_fit,
    data = data,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    ncensor.plot = TRUE,
    surv.median.line = "hv",
    palette = "npg",
    ggtheme = theme_minimal()
  )
  
  # Time-dependent ROC curve
  roc_obj <- timeROC(
    T = data$survival_time,
    delta = data$status,
    marker = data$predictor_value,
    cause = 1,
    times = c(365, 730, 1095)  # 1, 2, and 3 years
  )
  
  roc_plot <- plot(roc_obj)
  
  # Forest plot for multivariate analysis
  forest_data <- data.frame(
    hr = exp(coef(stats_results$multivariate)),
    lower = exp(confint(stats_results$multivariate)[,1]),
    upper = exp(confint(stats_results$multivariate)[,2])
  )
  
  forest_plot <- forestplot(
    forest_data,
    title = "Hazard Ratios from Multivariate Analysis"
  )
  
  # Waterfall plot
  waterfall_plot <- ggplot(data, aes(x = reorder(sample_id, predictor_value), 
                                     y = predictor_value)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Distribution of Values Across Samples",
         x = "Samples",
         y = "Value")
  
  return(list(
    survival = surv_plot,
    roc = roc_plot,
    forest = forest_plot,
    waterfall = waterfall_plot
  ))
}

# Enhanced main analysis function
perform_enhanced_survival_analysis <- function(cancer_type, 
                                               analysis_type = "expression",
                                               gene = "SRRM3",
                                               grouping_method = "quartile",
                                               min_coverage = 10,
                                               min_quality = 20,
                                               min_samples_per_group = 30) {
  # Create cache filename based on parameters
  cache_file <- file.path("cache", 
                         paste0("survival_analysis_", cancer_type, "_", 
                               analysis_type, "_", gene, ".rds"))
  
  # Check cache first
  cached_results <- get_cached_data(cache_file)
  if (!is.null(cached_results)) {
    return(cached_results)
  }
  
  # Get and process data (PSI or expression)
  if(analysis_type == "PSI") {
    # Cache for recount3 data
    recount3_cache <- file.path("cache", paste0("recount3_", cancer_type, ".rds"))
    rse <- get_cached_data(recount3_cache)
    
    if (is.null(rse)) {
      # Use the new function to find the project
      project_info <- find_tcga_project(cancer_type)
      
      message("Creating RSE object for ", project_info$project)
      rse <- create_rse(
        project_info[1, ],
        type = "jxn",
        jxn_format = "UNIQUE",
        verbose = TRUE
      )
      save_cached_data(rse, recount3_cache)
    }
    
    # [Previous PSI calculation code]
    psi_values <- calculate_psi_values(rse)
    values <- filter_psi_data(rse, psi_values, min_coverage, min_quality, min_samples_per_group)
    value_type <- "PSI"
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
    
    # [Previous expression calculation code]
    values <- get_expression_data(cancer_type, gene)
    value_type <- "Expression"
  }
  
  # Determine groups using advanced method
  groups <- determine_groups(values, method = grouping_method)
  
  # Get clinical data with extended variables
  clinical <- get_extended_clinical_data(cancer_type)
  
  # Combine all data
  analysis_data <- prepare_analysis_data(values, groups, clinical)
  
  # Perform advanced statistical analyses
  stats_results <- perform_advanced_statistics(
    analysis_data,
    time_col = "survival_time",
    status_col = "status",
    predictor_col = value_type,
    clinical_vars = c("age", "stage", "grade")
  )
  
  # Create survival object and fit
  surv_obj <- Surv(analysis_data$survival_time, analysis_data$status)
  surv_fit <- survfit(surv_obj ~ groups, data = analysis_data)
  
  # Generate all plots
  plots <- create_advanced_plots(analysis_data, surv_fit, stats_results)
  
  # Cross-validation
  cv_results <- perform_cross_validation(analysis_data, stats_results$multivariate)
  
  results <- list(
    data = analysis_data,
    statistics = stats_results,
    plots = plots,
    cross_validation = cv_results
  )
  save_cached_data(results, cache_file)
  return(results)
}

# Example usage:
# results <- perform_enhanced_survival_analysis(
#   cancer_type = "BRCA",
#   analysis_type = "PSI",
#   grouping_method = "quartile",
#   min_coverage = 10,
#   min_quality = 20,
#   min_samples_per_group = 30
# )
#
# # Access results
# results$plots$survival  # View survival plot
# results$statistics$multivariate  # View multivariate analysis results
# results$cross_validation  # View cross-validation results

results <- perform_enhanced_survival_analysis(
  cancer_type = "BRCA",
  analysis_type = "PSI",
  grouping_method = "quartile"
)

# View specific results
results$plots$survival  # Survival curves
results$statistics$multivariate  # Multivariate analysis
results$plots$forest  # Forest plot of hazard ratios