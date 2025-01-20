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

#' Perform survival analysis across multiple cancer types
#' @param cancer_types Vector of cancer type codes
#' @param analysis_type Either "PSI" or "expression"
#' @param gene Gene name to analyze
#' @param grouping_method Method for grouping samples ("quartile", "mixture", or "median")
#' @return List of results for each cancer type
perform_multi_cancer_analysis <- function(
    cancer_types = CANCER_TYPES,
    analysis_type = "PSI",
    gene = "SRRM3",
    grouping_method = "quartile"
) {
  # Create results directory if it doesn't exist
  results_dir <- file.path("results", "multi_cancer_analysis")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Initialize results list
  all_results <- list()
  
  # Create progress bar
  pb <- progress_bar$new(
    format = "Analyzing cancer type :what [:bar] :percent eta: :eta",
    total = length(cancer_types)
  )
  
  # Analyze each cancer type
  for (cancer_type in cancer_types) {
    pb$tick(tokens = list(what = cancer_type))
    
    tryCatch({
      # Perform analysis for current cancer type
      results <- perform_enhanced_survival_analysis(
        cancer_type = cancer_type,
        analysis_type = analysis_type,
        gene = gene,
        grouping_method = grouping_method
      )
      
      # Store results
      all_results[[cancer_type]] <- results
      
      # Save individual plots
      save_cancer_specific_plots(results, cancer_type, results_dir)
      
    }, error = function(e) {
      message("\nError analyzing ", cancer_type, ": ", e$message)
      all_results[[cancer_type]] <- list(error = e$message)
    })
  }
  
  # Generate comparative analyses
  comparative_results <- generate_comparative_analyses(all_results)
  
  # Save combined results
  save_combined_results(all_results, comparative_results, results_dir)
  
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
  
  # Create survival plot
  ggplot(combined_data, aes(x = survival_time, color = cancer_type)) +
    stat_survival(aes(time = survival_time, status = status), 
                 data = combined_data) +
    theme_minimal() +
    labs(
      title = "Survival Curves Across Cancer Types",
      x = "Time (days)",
      y = "Survival Probability"
    )
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
  
  # Save survival comparison plot
  pdf(file.path(results_dir, "survival_comparison_plot.pdf"), width = 12, height = 8)
  print(comparative_results$survival_comparison)
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

# Example usage:
# results <- perform_multi_cancer_analysis(
#   cancer_types = CANCER_TYPES,
#   analysis_type = "PSI",
#   gene = "SRRM3",
#   grouping_method = "quartile"
# ) 