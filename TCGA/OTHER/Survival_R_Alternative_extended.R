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
# library(timeROC)     # For time-dependent ROC
library(plotly)      # For interactive plots
library(cmprsk)      # For competing risks
# library(rms)         # For advanced survival modeling
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

# Modify find_tcga_project function to handle TCGA's specific format
find_tcga_project <- function(cancer_type) {
  message("Getting available projects...")
  projects <- available_projects()
  
  # Debug output
  message("\nSearching for cancer type: ", cancer_type)
  message("Total number of projects: ", nrow(projects))
  
  # First look for TCGA projects specifically
  tcga_projects <- subset(projects, 
                         file_source == "tcga" & 
                         project_type == "data_sources")
  
  message("\nFound ", nrow(tcga_projects), " TCGA projects")
  
  if (nrow(tcga_projects) > 0) {
    message("\nSample of TCGA projects:")
    print(head(tcga_projects))
    
    # Look for exact match first
    exact_match <- subset(tcga_projects, project == cancer_type)
    if (nrow(exact_match) == 1) {
      message("\nFound exact match:")
      print(exact_match)
      return(exact_match)
    }
    
    # If no exact match, try case-insensitive search
    matches <- tcga_projects[grep(cancer_type, tcga_projects$project, ignore.case = TRUE), ]
    if (nrow(matches) > 0) {
      message("\nFound matching project:")
      print(matches[1, ])
      return(matches[1, ])
    }
  }
  
  # If no match found, show what's available
  message("\nNo match found. Available TCGA projects:")
  print(head(tcga_projects, n = 10))
  
  stop(paste("Could not find TCGA project for", cancer_type, 
             ". Please check the project list above and verify the correct project name."))
}

# Add a helper function to explore available projects
explore_projects <- function() {
  projects <- available_projects()
  
  message("Total number of projects: ", nrow(projects))
  
  # Look for TCGA projects
  tcga_projects <- projects[grep("tcga", projects$project, ignore.case = TRUE), ]
  message("\nTCGA projects found: ", nrow(tcga_projects))
  
  if (nrow(tcga_projects) > 0) {
    message("\nSample of TCGA projects:")
    print(head(tcga_projects))
    
    # Try to identify the pattern used for TCGA projects
    message("\nUnique TCGA project patterns:")
    patterns <- unique(gsub("[A-Z]+$", "", tcga_projects$project))
    print(patterns)
  }
  
  return(list(
    all_projects = projects,
    tcga_projects = tcga_projects
  ))
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

# Add a new function to process junction data in chunks
process_junction_data <- function(rse, chunk_size = 1000) {
  message("Processing junction data in chunks...")
  
  total_features <- nrow(rse)
  total_chunks <- ceiling(total_features / chunk_size)
  
  pb <- progress_bar$new(
    format = "Processing chunk :current/:total [:bar] :percent eta: :eta",
    total = total_chunks
  )
  
  processed_data <- list()
  
  for (i in seq(1, total_features, by = chunk_size)) {
    chunk_end <- min(i + chunk_size - 1, total_features)
    chunk <- rse[i:chunk_end, ]
    
    # Process chunk here
    jxn_coords <- rowRanges(chunk)
    junctions <- find_exon15_junctions(jxn_coords)
    
    if (length(junctions$inclusion) > 0 || length(junctions$exclusion) > 0) {
      processed_data[[length(processed_data) + 1]] <- list(
        coords = jxn_coords,
        junctions = junctions
      )
    }
    
    pb$tick()
  }
  
  return(processed_data)
}

# Modify calculate_psi_values to use chunked processing
calculate_psi_values <- function(rse) {
  if(is.null(rse)) return(NULL)
  
  message("Starting PSI calculation...")
  
  # Process data in chunks
  processed_chunks <- process_junction_data(rse)
  
  if (length(processed_chunks) == 0) {
    message("No relevant junctions found")
    return(NULL)
  }
  
  # Combine results from all chunks
  all_inclusion <- unique(unlist(lapply(processed_chunks, function(x) x$junctions$inclusion)))
  all_exclusion <- unique(unlist(lapply(processed_chunks, function(x) x$junctions$exclusion)))
  
  message("Calculating PSI values for ", ncol(rse), " samples...")
  
  # Get junction counts
  junction_counts <- assay(rse)
  
  # Calculate PSI for each sample with progress bar
  pb <- progress_bar$new(
    format = "Calculating PSI :current/:total [:bar] :percent eta: :eta",
    total = ncol(junction_counts)
  )
  
  psi_values <- sapply(seq_len(ncol(junction_counts)), function(i) {
    inclusion_reads <- sum(junction_counts[all_inclusion, i])
    exclusion_reads <- sum(junction_counts[all_exclusion, i])
    total_reads <- inclusion_reads + exclusion_reads
    
    pb$tick()
    
    if(total_reads >= 10) {  # Minimum coverage threshold
      return((inclusion_reads / total_reads) * 100)
    } else {
      return(NA)
    }
  })
  
  message("PSI calculation complete")
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
  # roc_obj <- timeROC(
  #   T = data$survival_time,
  #   delta = data$status,
  #   marker = data$predictor_value,
  #   cause = 1,
  #   times = c(365, 730, 1095)  # 1, 2, and 3 years
  # )
  
  # roc_plot <- plot(roc_obj)
  
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
    # roc = roc_plot,
    forest = forest_plot,
    waterfall = waterfall_plot
  ))
}

# Add a helper function to retry downloads
retry_download <- function(url, destfile, max_tries = 3) {
  for (i in 1:max_tries) {
    tryCatch({
      download.file(url, destfile, quiet = TRUE)
      return(TRUE)
    }, error = function(e) {
      message("Download attempt ", i, " failed: ", e$message)
      if (i == max_tries) stop("Maximum retry attempts reached")
      Sys.sleep(2^i)  # Exponential backoff
    })
  }
}

# Modify perform_enhanced_survival_analysis to handle the RSE creation more carefully
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
  
  if(analysis_type == "PSI") {
    # Cache for recount3 data
    recount3_cache <- file.path("cache", paste0("recount3_", cancer_type, ".rds"))
    rse <- get_cached_data(recount3_cache)
    
    if (is.null(rse)) {
      tryCatch({
        # Use the modified find_tcga_project function
        project_info <- find_tcga_project(cancer_type)
        
        message("Creating RSE object for ", project_info$project)
        message("Project info:")
        print(project_info)
        
        # Create RSE with correct parameters
        rse <- create_rse(
          project = project_info,
          type = "jxn",
          jxn_format = "UNIQUE",
          verbose = TRUE
        )
        
        if (is.null(rse)) {
          stop("Failed to create RSE object")
        }
        
        save_cached_data(rse, recount3_cache)
      }, error = function(e) {
        message("Error creating RSE object: ", e$message)
        message("\nAttempting to get more information about the error...")
        message("Available arguments for create_rse:")
        print(args(create_rse))
        message("\nProject search details:")
        print(explore_projects())
        stop(e)
      })
    }
    
    # Continue with PSI calculation
    message("Calculating PSI values...")
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

# Comment out the example execution
# results <- perform_enhanced_survival_analysis(
#   cancer_type = "BRCA",
#   analysis_type = "PSI",
#   grouping_method = "quartile"
# )

# # View specific results
# results$plots$survival  # Survival curves
# results$statistics$multivariate  # Multivariate analysis
# results$plots$forest  # Forest plot of hazard ratios