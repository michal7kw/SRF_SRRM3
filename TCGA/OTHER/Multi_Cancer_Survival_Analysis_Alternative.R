# Load the same libraries as in Survival_R_Alternative_extended.R
source("Survival_R_Alternative_extended.R")

# Define cancer types to analyze
CANCER_TYPES <- c(
  "ACC",  # Adrenocortical carcinoma
  "UVM",  # Uveal melanoma
  "SKCM", # Skin cutaneous melanoma
  "LGG",  # Brain lower grade glioma (note: fixed typo from LLG)
  "GBM"   # Glioblastoma
)

# Add this function at the beginning after loading libraries
setup_logging <- function() {
  log_dir <- "logs/progress"
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Setup futile.logger to write to both console and file
  log_file <- file.path(log_dir, paste0("analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
  flog.appender(appender.tee(log_file))
  flog.threshold(INFO)
}

# Add this function to check sample sizes
check_sample_sizes <- function(values, groups) {
  sample_counts <- table(groups)
  flog.info("Sample sizes per group:")
  print(sample_counts)
  
  if (any(sample_counts < 10)) {
    flog.warn("Warning: Some groups have fewer than 10 samples")
    return(FALSE)
  }
  return(TRUE)
}

#' Perform survival analysis across multiple cancer types
#' @param cancer_types Vector of cancer type codes
#' @param analysis_type Either "PSI" or "expression"
#' @param gene Gene name to analyze
#' @param grouping_method Method for grouping samples ("quartile", "mixture", or "median")
#' @param min_samples Minimum samples per group
#' @return List of results for each cancer type
perform_multi_cancer_analysis <- function(
    cancer_types = CANCER_TYPES,
    analysis_type = "PSI",
    gene = "SRRM3",
    grouping_method = "quartile",
    min_samples = 10
) {
  setup_logging()
  flog.info("Starting multi-cancer analysis")
  flog.info("Parameters: analysis_type=%s, gene=%s, grouping_method=%s, min_samples=%s", 
            analysis_type, gene, grouping_method, min_samples)
  
  # Create results directory if it doesn't exist
  results_dir <- file.path("results", "multi_cancer_analysis")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Initialize results list
  all_results <- list()
  
  # Analyze each cancer type
  for (cancer_type in cancer_types) {
    flog.info("Starting analysis for cancer type: %s", cancer_type)
    
    tryCatch({
      # Record memory usage
      flog.info("Memory usage before analysis: %s MB", 
                format(gc()["Ncells", "used"] / 1024^2, digits=2))
      
      # Perform analysis for current cancer type using wrapped function
      results <- perform_enhanced_survival_analysis_with_validation(
        cancer_type = cancer_type,
        analysis_type = analysis_type,
        gene = gene,
        grouping_method = grouping_method,
        min_samples_per_group = min_samples
      )
      
      # Check sample sizes before proceeding
      if (!check_sample_sizes(results$data$value, results$groups)) {
        flog.warn("Skipping %s due to insufficient sample sizes", cancer_type)
        next
      }
      
      # Store results
      all_results[[cancer_type]] <- results
      
      # Save individual plots
      flog.info("Saving plots for %s", cancer_type)
      save_cancer_specific_plots(results, cancer_type, results_dir)
      
      flog.info("Completed analysis for %s", cancer_type)
      
    }, error = function(e) {
      flog.error("Error analyzing %s: %s", cancer_type, e$message)
      all_results[[cancer_type]] <- list(error = e$message)
    })
    
    # Clean up memory
    gc()
  }
  
  flog.info("Generating comparative analyses")
  comparative_results <- generate_comparative_analyses(all_results)
  
  flog.info("Saving combined results")
  save_combined_results(all_results, comparative_results, results_dir)
  
  flog.info("Analysis complete")
  return(list(
    individual_results = all_results,
    comparative_results = comparative_results
  ))
}

#' Save plots for a specific cancer type
#' @param results Results from survival analysis
#' @param cancer_type Cancer type code
#' @param results_dir Directory to save results
save_cancer_specific_plots <- function(results, cancer_type, results_dir) {
  # Create cancer-specific directory
  cancer_dir <- file.path(results_dir, cancer_type)
  dir.create(cancer_dir, showWarnings = FALSE)
  
  # Save survival plot
  survival_plot <- results$plots$survival
  pdf(file.path(cancer_dir, "survival_plot.pdf"), width = 10, height = 8)
  print(survival_plot)
  dev.off()
  
  # Save forest plot
  forest_plot <- results$plots$forest
  pdf(file.path(cancer_dir, "forest_plot.pdf"), width = 10, height = 8)
  print(forest_plot)
  dev.off()
  
  # Save waterfall plot
  waterfall_plot <- results$plots$waterfall
  pdf(file.path(cancer_dir, "waterfall_plot.pdf"), width = 10, height = 8)
  print(waterfall_plot)
  dev.off()
}

#' Generate comparative analyses across cancer types
#' @param all_results List of results for each cancer type
#' @return List of comparative analyses
generate_comparative_analyses <- function(all_results) {
  # Extract hazard ratios and p-values
  hr_data <- data.frame(
    cancer_type = character(),
    hazard_ratio = numeric(),
    ci_lower = numeric(),
    ci_upper = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (cancer_type in names(all_results)) {
    results <- all_results[[cancer_type]]
    if (!is.null(results$statistics$univariate)) {
      model <- results$statistics$univariate
      hr <- exp(coef(model))
      ci <- exp(confint(model))
      p_val <- summary(model)$coefficients[,5]
      
      hr_data <- rbind(hr_data, data.frame(
        cancer_type = cancer_type,
        hazard_ratio = hr,
        ci_lower = ci[1],
        ci_upper = ci[2],
        p_value = p_val
      ))
    }
  }
  
  # Create forest plot comparing all cancer types
  forest_plot <- create_comparative_forest_plot(hr_data)
  
  # Create survival comparison plot
  survival_comparison <- create_survival_comparison_plot(all_results)
  
  return(list(
    hazard_ratio_data = hr_data,
    forest_plot = forest_plot,
    survival_comparison = survival_comparison
  ))
}

#' Create forest plot comparing hazard ratios across cancer types
#' @param hr_data Data frame with hazard ratio information
#' @return ggplot object
create_comparative_forest_plot <- function(hr_data) {
  ggplot(hr_data, aes(y = cancer_type, x = hazard_ratio)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_log10() +
    theme_minimal() +
    labs(
      title = "Hazard Ratios Across Cancer Types",
      x = "Hazard Ratio (log scale)",
      y = "Cancer Type"
    )
}

#' Create survival comparison plot across cancer types
#' @param all_results List of results for each cancer type
#' @return ggplot object
create_survival_comparison_plot <- function(all_results) {
  # Combine survival data from all cancer types
  combined_data <- data.frame()
  
  for (cancer_type in names(all_results)) {
    if (!is.null(all_results[[cancer_type]]$data)) {
      data <- all_results[[cancer_type]]$data
      data$cancer_type <- cancer_type
      combined_data <- rbind(combined_data, data)
    }
  }
  
  # Create survival plot using survminer
  fit <- survfit(Surv(survival_time, status) ~ cancer_type, data = combined_data)
  
  survplot <- ggsurvplot(
    fit,
    data = combined_data,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.col = "strata",
    linetype = "strata",
    surv.median.line = "hv",
    ggtheme = theme_minimal(),
    palette = "Dark2",
    title = "Survival Curves Across Cancer Types",
    xlab = "Time (days)",
    ylab = "Survival Probability"
  )
  
  return(survplot)
}

#' Save combined results and comparative analyses
#' @param all_results Individual cancer type results
#' @param comparative_results Comparative analyses
#' @param results_dir Directory to save results
save_combined_results <- function(all_results, comparative_results, results_dir) {
  # Save comparative forest plot
  pdf(file.path(results_dir, "comparative_forest_plot.pdf"), width = 12, height = 8)
  print(comparative_results$forest_plot)
  dev.off()
  
  # Save survival comparison plot (modified for survminer output)
  pdf(file.path(results_dir, "survival_comparison_plot.pdf"), width = 12, height = 10)
  print(comparative_results$survival_comparison$plot)
  dev.off()
  
  # Save risk table
  pdf(file.path(results_dir, "survival_risk_table.pdf"), width = 12, height = 4)
  print(comparative_results$survival_comparison$table)
  dev.off()
  
  # Save hazard ratio data
  write.csv(
    comparative_results$hazard_ratio_data,
    file.path(results_dir, "hazard_ratios.csv"),
    row.names = FALSE
  )
  
  # Save complete results as RDS
  saveRDS(
    list(
      individual_results = all_results,
      comparative_results = comparative_results
    ),
    file.path(results_dir, "complete_results.rds")
  )
}

# Add this function to prepare and validate clinical data
prepare_clinical_data <- function(cancer_type, values) {
  flog.info("Preparing clinical data for %s", cancer_type)
  
  # Print structure of input values
  flog.info("Structure of input values:")
  print(str(values))
  
  tryCatch({
    # Get clinical data
    clinical <- get_extended_clinical_data(cancer_type)
    
    if (is.null(clinical)) {
      flog.error("No clinical data found for %s", cancer_type)
      return(NULL)
    }
    
    # Ensure we have the required columns
    if (!all(c("survival_time", "status") %in% names(clinical))) {
      flog.error("Missing required columns in clinical data")
      flog.info("Available columns: %s", paste(names(clinical), collapse = ", "))
      return(NULL)
    }
    
    # Print sample IDs from both datasets
    flog.info("Sample IDs in values: %s", paste(head(names(values)), collapse = ", "))
    flog.info("Sample IDs in clinical: %s", paste(head(rownames(clinical)), collapse = ", "))
    
    # Try to standardize sample IDs
    names(values) <- gsub("\\.[0-9]+$", "", names(values))  # Remove version numbers
    rownames(clinical) <- gsub("\\.[0-9]+$", "", rownames(clinical))
    
    # Match samples
    common_samples <- intersect(rownames(clinical), names(values))
    
    if (length(common_samples) == 0) {
      flog.error("No matching samples between clinical data and values")
      return(NULL)
    }
    
    flog.info("Found %d matching samples", length(common_samples))
    
    # Create analysis data frame
    analysis_data <- data.frame(
      sample_id = common_samples,
      survival_time = clinical[common_samples, "survival_time"],
      status = clinical[common_samples, "status"],
      value = values[common_samples],
      stringsAsFactors = FALSE
    )
    
    # Remove any NA values
    complete_cases <- complete.cases(analysis_data)
    if (!all(complete_cases)) {
      flog.warn("Removing %d samples with missing values", sum(!complete_cases))
      analysis_data <- analysis_data[complete_cases, ]
    }
    
    flog.info("Final analysis data structure:")
    print(str(analysis_data))
    
    return(analysis_data)
    
  }, error = function(e) {
    flog.error("Error in prepare_clinical_data: %s", e$message)
    print(sys.calls())
    return(NULL)
  })
}

# Wrap the original function to add our clinical data validation
perform_enhanced_survival_analysis_with_validation <- function(
    cancer_type,
    analysis_type = "PSI",
    gene = "SRRM3",
    grouping_method = "quartile",
    min_samples_per_group = 10
) {
  # Call the original function
  results <- perform_enhanced_survival_analysis(
    cancer_type = cancer_type,
    analysis_type = analysis_type,
    gene = gene,
    grouping_method = grouping_method,
    min_samples_per_group = min_samples_per_group
  )
  
  # Add our additional validation
  if (!is.null(results$data)) {
    analysis_data <- prepare_clinical_data(cancer_type, results$data$values)
    
    if (!is.null(analysis_data)) {
      # Update the results with validated data
      results$data <- analysis_data
      results$groups <- determine_groups(analysis_data$value, method = grouping_method)
      
      # Update survival analysis with validated data
      surv_obj <- Surv(analysis_data$survival_time, analysis_data$status)
      results$survival_fit <- survfit(surv_obj ~ results$groups, data = analysis_data)
      
      # Update plots and statistics
      results$plots <- create_advanced_plots(analysis_data, results$survival_fit)
      results$statistics <- perform_advanced_statistics(analysis_data)
    }
  }
  
  return(results)
}

# Add this function that was missing
prepare_analysis_data <- function(values, groups, clinical_data) {
  flog.info("Preparing analysis data")
  
  # Print input data summaries
  flog.info("Number of values: %d", length(values))
  flog.info("Number of groups: %d", length(groups))
  flog.info("Number of clinical samples: %d", nrow(clinical_data))
  
  # Create the basic data frame
  analysis_data <- data.frame(
    sample_id = names(values),
    value = values,
    group = groups,
    stringsAsFactors = FALSE
  )
  
  # Print sample IDs for debugging
  flog.info("First few sample IDs in values:")
  print(head(names(values)))
  flog.info("First few sample IDs in clinical data:")
  print(head(rownames(clinical_data)))
  
  # Add clinical data
  analysis_data$survival_time <- clinical_data$survival_time[match(analysis_data$sample_id, rownames(clinical_data))]
  analysis_data$status <- clinical_data$status[match(analysis_data$sample_id, rownames(clinical_data))]
  
  # Check data types
  flog.info("Data types in analysis_data:")
  print(sapply(analysis_data, class))
  
  # Remove any NA values
  complete_cases <- complete.cases(analysis_data)
  if (!all(complete_cases)) {
    flog.warn("Removing %d samples with missing values", sum(!complete_cases))
    analysis_data <- analysis_data[complete_cases, ]
  }
  
  # Validate final data
  if (nrow(analysis_data) == 0) {
    flog.error("No valid samples after merging clinical data")
    return(NULL)
  }
  
  # Ensure survival_time is numeric
  analysis_data$survival_time <- as.numeric(analysis_data$survival_time)
  analysis_data$status <- as.numeric(analysis_data$status)
  
  flog.info("Final analysis data structure:")
  print(str(analysis_data))
  
  return(analysis_data)
}

# Modify get_extended_clinical_data function to use correct TCGA column names
get_extended_clinical_data <- function(cancer_type) {
  flog.info("Getting clinical data for %s", cancer_type)
  
  tryCatch({
    project_info <- find_tcga_project(cancer_type)
    flog.info("Creating RSE for clinical data with project info:")
    print(project_info)
    
    # Create RSE
    rse <- create_rse(
      project = project_info,
      type = "jxn"
    )
    
    flog.info("RSE object created, extracting clinical data")
    clinical <- colData(rse)
    
    flog.info("Available columns in raw clinical data:")
    print(names(clinical))
    
    clinical <- as.data.frame(clinical)
    
    # Look for TCGA-specific survival columns
    time_col <- NULL
    status_col <- NULL
    
    # Check for days_to_death and days_to_last_follow_up in different locations
    if ("tcga.gdc_cases.diagnoses.days_to_death" %in% names(clinical)) {
      time_col <- "tcga.gdc_cases.diagnoses.days_to_death"
    } else if ("tcga.xml_days_to_death" %in% names(clinical)) {
      time_col <- "tcga.xml_days_to_death"
    } else if ("tcga.cgc_case_days_to_death" %in% names(clinical)) {
      time_col <- "tcga.cgc_case_days_to_death"
    }
    
    # For patients without death event, use days_to_last_follow_up
    followup_col <- NULL
    if ("tcga.gdc_cases.diagnoses.days_to_last_follow_up" %in% names(clinical)) {
      followup_col <- "tcga.gdc_cases.diagnoses.days_to_last_follow_up"
    } else if ("tcga.xml_days_to_last_followup" %in% names(clinical)) {
      followup_col <- "tcga.xml_days_to_last_followup"
    } else if ("tcga.cgc_case_days_to_last_follow_up" %in% names(clinical)) {
      followup_col <- "tcga.cgc_case_days_to_last_follow_up"
    }
    
    # Look for vital status in different locations
    if ("tcga.gdc_cases.diagnoses.vital_status" %in% names(clinical)) {
      status_col <- "tcga.gdc_cases.diagnoses.vital_status"
    } else if ("tcga.xml_vital_status" %in% names(clinical)) {
      status_col <- "tcga.xml_vital_status"
    } else if ("tcga.cgc_case_vital_status" %in% names(clinical)) {
      status_col <- "tcga.cgc_case_vital_status"
    }
    
    if (is.null(time_col) || is.null(status_col)) {
      flog.error("Could not find survival time or status columns")
      return(NULL)
    }
    
    # Extract time values, using days_to_death if available, otherwise days_to_last_follow_up
    time_values <- clinical[[time_col]]
    if (!is.null(followup_col)) {
      followup_values <- clinical[[followup_col]]
      # Replace NA death times with followup times
      time_values[is.na(time_values)] <- followup_values[is.na(time_values)]
    }
    
    flog.info("Raw time values (first few):")
    print(head(time_values))
    
    # Convert time to numeric, handling special cases
    clinical$survival_time <- as.numeric(as.character(time_values))
    
    # Handle NA values in time
    if (all(is.na(clinical$survival_time))) {
      flog.error("All survival times are NA after conversion")
      return(NULL)
    }
    
    # Convert vital status to binary (1 = deceased, 0 = alive)
    status_values <- clinical[[status_col]]
    if (is.factor(status_values)) {
      status_values <- as.character(status_values)
    }
    
    flog.info("Unique status values:")
    print(table(status_values))
    
    clinical$status <- case_when(
      tolower(as.character(status_values)) %in% 
        c("1", "true", "yes", "dead", "deceased", "died") ~ 1,
      tolower(as.character(status_values)) %in% 
        c("0", "false", "no", "alive", "living") ~ 0,
      TRUE ~ NA_real_
    )
    
    # Print summaries
    flog.info("Summary of converted survival times:")
    print(summary(clinical$survival_time))
    flog.info("Summary of converted status:")
    print(table(clinical$status, useNA = "ifany"))
    
    # Validate
    valid_samples <- !is.na(clinical$survival_time) & !is.na(clinical$status)
    valid_count <- sum(valid_samples)
    
    flog.info("Found %d valid samples", valid_count)
    
    if (valid_count == 0) {
      flog.error("No valid survival data found for %s", cancer_type)
      return(NULL)
    }
    
    # Keep only valid samples
    clinical <- clinical[valid_samples, ]
    
    # Ensure rownames are preserved
    rownames(clinical) <- rownames(clinical[valid_samples, ])
    
    flog.info("Successfully prepared clinical data for %s with %d valid samples", 
              cancer_type, nrow(clinical))
    
    return(clinical)
    
  }, error = function(e) {
    flog.error("Error getting clinical data for %s: %s", cancer_type, e$message)
    flog.error("Error traceback:")
    print(sys.calls())
    return(NULL)
  })
}

# Example usage:
# results <- perform_multi_cancer_analysis(
#   cancer_types = CANCER_TYPES,
#   analysis_type = "PSI",
#   gene = "SRRM3",
#   grouping_method = "quartile"
# ) 