# 1. Package Management and Environment Setup
required_packages <- c(
  "survival", "survminer", "TCGAbiolinks", "SummarizedExperiment",
  "tidyverse", "recount3", "GenomicRanges", "progress", "futile.logger",
  "R.utils", "viridis", "mclust", "timeROC", "plotly", "cmprsk", "rms",
  "forestplot", "R6", "digest"  # Added digest as it's used in the code
)

# Environment validation and package loading
validate_environment <- function() {
  missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
  if (length(missing_packages) > 0) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
  }
  
  # Load all required packages
  for(pkg in required_packages) {
    library(pkg, character.only = TRUE)
  }
  
  # Also validate required directories
  required_dirs <- c("cache", "logs", "output")
  for(dir in required_dirs) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  }
}

# Run validation first
validate_environment()

# Initialize logging after packages are loaded
flog.appender(appender.file("./logs/survival_analysis.log"))
flog.threshold(DEBUG)

# 2. Constants and Configuration
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

# 3. Cache Management System
CacheManager <- R6::R6Class(
  "CacheManager",
  public = list(
    cache_dir = "cache",
    max_age_days = 30,
    
    initialize = function(cache_dir = "cache", max_age_days = 30) {
      self$cache_dir <- cache_dir
      self$max_age_days <- max_age_days
      if (!dir.exists(self$cache_dir)) {
        dir.create(self$cache_dir, recursive = TRUE)
      }
    },
    
    get_data = function(key) {
      cache_file <- file.path(self$cache_dir, paste0(key, ".rds"))
      
      if (file.exists(cache_file)) {
        file_age <- difftime(Sys.time(), file.mtime(cache_file), units = "days")
        if (file_age < self$max_age_days) {
          futile.logger::flog.info("Loading cached data: %s", key)
          return(readRDS(cache_file))
        }
      }
      return(NULL)
    },
    
    save_data = function(data, key) {
      cache_file <- file.path(self$cache_dir, paste0(key, ".rds"))
      futile.logger::flog.info("Caching data: %s", key)
      saveRDS(data, cache_file)
    }
  )
)

# 4. Error Handling
with_error_handling <- function(expr, message = "Operation failed") {
  tryCatch(
    expr,
    error = function(e) {
      error_msg <- sprintf("%s: %s", message, e$message)
      futile.logger::flog.error(error_msg)
      stop(error_msg)
    },
    warning = function(w) {
      warning_msg <- sprintf("Warning in %s: %s", message, w$message)
      futile.logger::flog.warn(warning_msg)
    }
  )
}

# 5. Core Functions for Junction Analysis and Project Management
find_tcga_project_optimized <- function(cancer_type) {
  cache_key <- paste("project_info", cancer_type, sep = "_")
  cached_data <- cache_mgr$get_data(cache_key)
  if (!is.null(cached_data)) return(cached_data)
  
  with_error_handling({
    projects <- available_projects()
    tcga_projects <- subset(projects, 
                           file_source == "tcga" & 
                           project_type == "data_sources")
    
    # Look for exact match first
    exact_match <- subset(tcga_projects, project == cancer_type)
    if (nrow(exact_match) == 1) {
      cache_mgr$save_data(exact_match[1, ], cache_key)
      return(exact_match[1, ])
    }
    
    # Try with TCGA prefix
    tcga_name <- paste0("TCGA-", cancer_type)
    exact_match <- subset(tcga_projects, project == tcga_name)
    if (nrow(exact_match) == 1) {
      cache_mgr$save_data(exact_match[1, ], cache_key)
      return(exact_match[1, ])
    }
    
    stop(paste("Could not find TCGA project for", cancer_type))
  }, "Project search failed")
}

find_exon15_junctions <- function(jxn_coords) {
  exon_start <- SRRM3_INFO$exon15$start
  exon_end <- SRRM3_INFO$exon15$end
  
  # Debug coordinates
  futile.logger::flog.debug("Looking for junctions around exon 15: %d-%d", exon_start, exon_end)
  
  # Find upstream junction (ending at exon start)
  upstream_jxns <- which(abs(end(jxn_coords) - exon_start) <= 5)  # Allow 5bp flexibility
  
  # Find downstream junction (starting at exon end)
  downstream_jxns <- which(abs(start(jxn_coords) - exon_end) <= 5)  # Allow 5bp flexibility
  
  # Find exclusion junctions (those that skip exon 15)
  exclusion_jxns <- which(
    start(jxn_coords) < (exon_start - 5) & 
    end(jxn_coords) > (exon_end + 5)
  )
  
  # Debug junction counts
  futile.logger::flog.debug("Upstream junctions found: %d", length(upstream_jxns))
  futile.logger::flog.debug("Downstream junctions found: %d", length(downstream_jxns))
  futile.logger::flog.debug("Exclusion junctions found: %d", length(exclusion_jxns))
  
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

calculate_exon15_psi <- function(rse) {
  if(is.null(rse)) return(NULL)
  
  cache_key <- digest::digest(list(
    coords = as.data.frame(rowRanges(rse)),
    counts = assay(rse)
  ))
  cache_file <- file.path(cache_mgr$cache_dir, paste0("psi_", cache_key, ".rds"))
  
  cached_data <- cache_mgr$get_data(cache_file)
  if (!is.null(cached_data)) return(cached_data)
  
  # Get junction counts
  junction_counts <- assay(rse)
  jxn_coords <- rowRanges(rse)
  
  # Find relevant junctions
  junctions <- find_exon15_junctions(jxn_coords)
  
  if(length(junctions$inclusion) == 0 || length(junctions$exclusion) == 0) {
    futile.logger::flog.warn("Missing junctions - Inclusion: %d, Exclusion: %d", 
                            length(junctions$inclusion), length(junctions$exclusion))
    return(NULL)
  }
  
  # Calculate PSI for each sample
  psi_values <- sapply(seq_len(ncol(junction_counts)), function(i) {
    inclusion_reads <- sum(junction_counts[junctions$inclusion, i])
    exclusion_reads <- sum(junction_counts[junctions$exclusion, i])
    total_reads <- inclusion_reads + exclusion_reads
    
    if(total_reads >= 10) {  # Minimum coverage threshold
      return((inclusion_reads / total_reads) * 100)
    } else {
      return(NA)
    }
  })
  
  results <- list(
    psi = psi_values,
    junctions = list(
      inclusion = rownames(junction_counts)[junctions$inclusion],
      exclusion = rownames(junction_counts)[junctions$exclusion]
    )
  )
  
  cache_mgr$save_data(results, cache_file)
  
  futile.logger::flog.info("PSI calculation complete:")
  futile.logger::flog.info("Total samples: %d", length(psi_values))
  futile.logger::flog.info("Non-NA samples: %d", sum(!is.na(psi_values)))
  futile.logger::flog.info("Mean PSI: %.2f", mean(psi_values, na.rm = TRUE))
  
  return(results)
}

create_rse_safe <- function(project_info) {
  cache_file <- file.path(cache_mgr$cache_dir, 
                         paste0("rse_", project_info$project, ".rds"))
  
  cached_data <- cache_mgr$get_data(cache_file)
  if (!is.null(cached_data)) return(cached_data)
  
  with_error_handling({
    futile.logger::flog.info("Creating RSE for project: %s", project_info$project)
    
    # Create region around exon 15
    region <- GRanges(
      seqnames = SRRM3_INFO$gene$chr,
      ranges = IRanges(
        start = SRRM3_INFO$exon15$start - 5000,
        end = SRRM3_INFO$exon15$end + 5000
      )
    )
    
    # Use UNIQUE junction format for TCGA
    rse <- create_rse(
      project_info,
      type = "jxn",
      jxn_format = "UNIQUE",
      verbose = TRUE
    )
    
    # Filter to relevant region
    relevant_rows <- which(
      as.character(seqnames(rowRanges(rse))) == as.character(seqnames(region)) &
        start(rowRanges(rse)) >= start(region) &
        end(rowRanges(rse)) <= end(region)
    )
    
    if(length(relevant_rows) == 0) {
      futile.logger::flog.warn("No features found in the specified region")
      return(NULL)
    }
    
    rse_filtered <- rse[relevant_rows, ]
    rm(rse)
    gc()
    
    cache_mgr$save_data(rse_filtered, cache_file)
    return(rse_filtered)
  }, "RSE creation failed")
}

# 6. Data Retrieval Functions
get_psi_data_optimized <- function(cancer_type, min_coverage = 10) {
  if (!is.character(cancer_type) || length(cancer_type) != 1) {
    stop("cancer_type must be a single character string")
  }
  if (!is.numeric(min_coverage) || min_coverage < 0) {
    stop("min_coverage must be a positive number")
  }
  cache_key <- paste("psi_data", cancer_type, min_coverage, sep = "_")
  cached_data <- cache_mgr$get_data(cache_key)
  if (!is.null(cached_data)) return(cached_data)
  
  with_error_handling({
    # Get project info
    project_info <- find_tcga_project_optimized(cancer_type)
    
    # Create RSE object with region filtering
    rse <- create_rse_safe(project_info)
    
    if(is.null(rse)) {
      stop("Failed to create RSE object")
    }
    
    # Calculate PSI values
    psi_result <- calculate_exon15_psi(rse)
    
    if(is.null(psi_result)) {
      stop("Failed to calculate PSI values")
    }
    
    # Create return data frame
    result_df <- data.frame(
      sample = colnames(rse),
      value = psi_result$psi,
      stringsAsFactors = FALSE
    )
    
    cache_mgr$save_data(result_df, cache_key)
    return(result_df)
  }, "PSI data retrieval failed")
}

get_expression_data_optimized <- function(cancer_type, gene) {
  if (!is.character(cancer_type) || length(cancer_type) != 1) {
    stop("cancer_type must be a single character string")
  }
  if (!is.character(gene) || length(gene) != 1) {
    stop("gene must be a single character string")
  }
  cache_key <- paste("expression", cancer_type, gene, sep = "_")
  cached_data <- cache_mgr$get_data(cache_key)
  if (!is.null(cached_data)) return(cached_data)
  
  with_error_handling({
    query <- GDCquery(
      project = paste0("TCGA-", cancer_type),
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      genes = gene
    )
    
    GDCdownload(query, method = "api", directory = "cache")
    exp_data <- GDCprepare(query, directory = "cache")
    
    # Extract and process expression data
    expression_values <- data.frame(
      sample = colnames(exp_data),
      expression = as.numeric(assay(exp_data)[grep(gene, rowData(exp_data)$gene_name), ]),
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        sample = substr(sample, 1, 12),
        expression = log2(expression + 1)
      )
    
    cache_mgr$save_data(expression_values, cache_key)
    return(expression_values)
  }, paste("Expression data retrieval failed for", gene))
}

get_clinical_data_optimized <- function(cancer_type) {
  cache_key <- paste("clinical", cancer_type, sep = "_")
  cached_data <- cache_mgr$get_data(cache_key)
  if (!is.null(cached_data)) return(cached_data)
  
  with_error_handling({
    clinical <- GDCquery_clinic(paste0("TCGA-", cancer_type))
    
    # Process both overall survival and PFS
    clinical <- clinical %>%
      mutate(
        # Overall survival
        os_time = ifelse(vital_status == "Alive",
                        days_to_last_follow_up,
                        days_to_death),
        os_status = ifelse(vital_status == "Dead", 1, 0),
        
        # Progression-free survival
        pfs_time = ifelse(!is.na(days_to_progression),
                         days_to_progression,
                         days_to_last_follow_up),
        pfs_status = ifelse(!is.na(days_to_progression), 1, 0)
      ) %>%
      select(submitter_id, os_time, os_status, pfs_time, pfs_status,
             age_at_diagnosis, gender, ajcc_pathologic_stage,
             tumor_grade, everything())
    
    cache_mgr$save_data(clinical, cache_key)
    return(clinical)
  }, "Clinical data retrieval failed")
}

# 7. Survival Analysis Functions
prepare_survival_data <- function(values, clinical, surv_vars) {
  # Add validation
  if (!all(c("sample", "value") %in% names(values))) {
    stop("values must contain 'sample' and 'value' columns")
  }
  if (!all(c("submitter_id", surv_vars$time, surv_vars$status) %in% names(clinical))) {
    stop("Missing required clinical columns")
  }
  
  with_error_handling({
    # Clean and prepare the data
    analysis_data <- values %>%
      mutate(
        sample = substr(sample, 1, 12)  # Ensure consistent sample ID format
      ) %>%
      inner_join(clinical, by = c("sample" = "submitter_id")) %>%
      filter(!is.na(!!sym(surv_vars$time)) & !is.na(!!sym(surv_vars$status)))
    
    # Add grouping based on values
    quartiles <- quantile(analysis_data$value, probs = c(0.25, 0.75), na.rm = TRUE)
    analysis_data$group <- case_when(
      analysis_data$value <= quartiles[1] ~ "Low",
      analysis_data$value >= quartiles[2] ~ "High",
      TRUE ~ "Medium"
    )
    
    return(analysis_data)
  }, "Data preparation failed")
}

analyze_survival_data <- function(data, grouping_method, survival_type, analysis_type) {
  # Validate inputs
  if (nrow(data) < 10) {
    stop("Insufficient samples for analysis (minimum 10 required)")
  }
  
  if (!all(c("group", paste0(tolower(survival_type), "_time"), 
             paste0(tolower(survival_type), "_status")) %in% names(data))) {
    stop("Missing required columns in data")
  }
  
  # Create survival object
  surv_obj <- Surv(data[[paste0(tolower(survival_type), "_time")]],
                   data[[paste0(tolower(survival_type), "_status")]])
  
  # Fit survival curves
  fit <- survfit(surv_obj ~ group, data = data)
  
  # Cox proportional hazards model
  cox <- coxph(surv_obj ~ group + age_at_diagnosis + ajcc_pathologic_stage,
               data = data)
  
  # Create enhanced survival plot
  plot <- ggsurvplot(
    fit,
    data = data,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    ncensor.plot = TRUE,
    surv.median.line = "hv",
    palette = "npg",
    xlab = paste("Time (days)"),
    ylab = paste("Probability of", survival_type),
    title = paste(analysis_type, survival_type, "Analysis"),
    legend.title = "Group",
    risk.table.height = 0.25,
    ggtheme = theme_minimal()
  )
  
  # Additional statistics
  stats <- list(
    cox_model = cox,
    logrank = survdiff(surv_obj ~ group, data = data),
    median_survival = summary(fit)$table,
    concordance = concordance(cox)
  )
  
  # Validate results before returning
  if (is.null(fit) || is.null(cox) || is.null(plot)) {
    stop("Analysis failed to produce valid results")
  }
  
  # Add analysis metadata
  results <- list(
    plot = plot,
    statistics = stats,
    data = data,
    fit = fit,
    metadata = list(
      analysis_date = Sys.time(),
      samples = nrow(data),
      events = sum(data[[paste0(tolower(survival_type), "_status")]]),
      analysis_type = analysis_type,
      survival_type = survival_type,
      grouping_method = grouping_method
    )
  )
  
  return(results)
}

# 8. Main Analysis Function
perform_survival_analysis_enhanced <- function(cancer_type, 
                                             analysis_type = c("SRRM3_expression", 
                                                             "SRRM3_PSI", 
                                                             "SRRM4_expression"),
                                             survival_type = c("OS", "PFS"),
                                             grouping_method = c("quartile", "median", "mixture")) {
  # Add validation for cancer_type
  if (!is.character(cancer_type) || length(cancer_type) != 1) {
    stop("cancer_type must be a single character string")
  }
  
  # Validate environment first
  validate_environment()
  
  # Input validation
  analysis_type <- match.arg(analysis_type)
  survival_type <- match.arg(survival_type)
  grouping_method <- match.arg(grouping_method)
  
  # Create cache key
  cache_key <- paste("survival_analysis", cancer_type, analysis_type, 
                    survival_type, grouping_method, sep = "_")
  
  # Try to get cached results
  cached_results <- cache_mgr$get_data(cache_key)
  if (!is.null(cached_results)) return(cached_results)
  
  # Get data based on analysis type
  values <- with_error_handling({
    switch(analysis_type,
           SRRM3_expression = get_expression_data_optimized(cancer_type, "SRRM3"),
           SRRM3_PSI = get_psi_data_optimized(cancer_type),
           SRRM4_expression = get_expression_data_optimized(cancer_type, "SRRM4"))
  }, paste("Failed to get", analysis_type, "data"))
  
  if (is.null(values) || nrow(values) == 0) {
    stop("No data retrieved for ", analysis_type)
  }
  
  # Get clinical data
  clinical <- get_clinical_data_optimized(cancer_type)
  
  # Prepare survival variables
  surv_vars <- if(survival_type == "OS") {
    list(time = "os_time", status = "os_status", 
         label = "Overall Survival")
  } else {
    list(time = "pfs_time", status = "pfs_status", 
         label = "Progression-Free Survival")
  }
  
  # Prepare analysis data
  analysis_data <- prepare_survival_data(values, clinical, surv_vars)
  
  # Perform analysis
  results <- analyze_survival_data(
    analysis_data,
    grouping_method = grouping_method,
    survival_type = survival_type,
    analysis_type = analysis_type
  )
  
  # Cache results
  cache_mgr$save_data(results, cache_key)
  
  return(results)
}

# 9. Convenience Wrapper Function
run_survival_analysis <- function(cancer_type = "BRCA", 
                                retry_failed = TRUE, 
                                max_retries = 3) {
  futile.logger::flog.info("Starting survival analysis for %s", cancer_type)
  
  analysis_types <- c("SRRM3_expression", "SRRM3_PSI", "SRRM4_expression")
  survival_types <- c("OS", "PFS")
  total_analyses <- length(analysis_types) * length(survival_types)
  
  pb <- progress::progress_bar$new(
    format = "Analysis [:bar] :percent eta: :eta",
    total = total_analyses
  )
  
  results <- list()
  failed_analyses <- list()
  
  for(analysis_type in analysis_types) {
    for(survival_type in survival_types) {
      key <- paste(analysis_type, survival_type, sep = "_")
      futile.logger::flog.info("Running analysis: %s", key)
      
      success <- FALSE
      retries <- 0
      
      while (!success && retries < max_retries) {
        tryCatch({
          results[[key]] <- perform_survival_analysis_enhanced(
            cancer_type = cancer_type,
            analysis_type = analysis_type,
            survival_type = survival_type,
            grouping_method = "quartile"
          )
          success <- TRUE
        }, error = function(e) {
          retries <- retries + 1
          if (retries >= max_retries) {
            failed_analyses[[key]] <- e$message
            futile.logger::flog.error("Failed analysis %s after %d retries: %s", 
                                    key, max_retries, e$message)
          } else {
            futile.logger::flog.warn("Retry %d for %s: %s", 
                                   retries, key, e$message)
            Sys.sleep(retries * 2)  # Exponential backoff
          }
        })
      }
      
      pb$tick()
    }
  }
  
  if (length(failed_analyses) > 0) {
    results$failed_analyses <- failed_analyses
  }
  
  # Add summary for successful analyses
  successful_analyses <- results[!names(results) %in% c("failed_analyses", "summary")]
  if (length(successful_analyses) > 0) {
    results$summary <- lapply(successful_analyses, function(x) x$metadata)
  }
  
  return(results)
}

# Initialize cache manager
cache_mgr <- CacheManager$new()

# Add after Error Handling section
save_analysis_results <- function(results, cancer_type) {
  if (is.null(results)) {
    stop("No results to save")
  }
  
  output_dir <- file.path("output", cancer_type)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save plots
  for (analysis_name in names(results)) {
    if (analysis_name == "summary") next
    
    analysis <- results[[analysis_name]]
    if (!is.null(analysis$plot)) {
      filename <- file.path(output_dir, paste0(analysis_name, "_survival.pdf"))
      ggsave(filename, analysis$plot$plot, width = 10, height = 8)
    }
  }
  
  # Save statistics
  stats_file <- file.path(output_dir, "analysis_statistics.rds")
  saveRDS(results, stats_file)
  
  # Save summary as CSV
  summary_df <- do.call(rbind, lapply(results$summary, as.data.frame))
  write.csv(summary_df, 
            file.path(output_dir, "analysis_summary.csv"), 
            row.names = FALSE)
}

cleanup_cache <- function(max_age_days = 30) {
  cache_dir <- cache_mgr$cache_dir
  files <- list.files(cache_dir, full.names = TRUE)
  
  for (file in files) {
    file_age <- difftime(Sys.time(), file.mtime(file), units = "days")
    if (file_age > max_age_days) {
      unlink(file)
      futile.logger::flog.info("Removed old cache file: %s", basename(file))
    }
  }
}

# Example usage with all features:
results <- run_survival_analysis("BRCA")
save_analysis_results(results, "BRCA")
# cleanup_cache()  # Optional cleanup of old cache files 