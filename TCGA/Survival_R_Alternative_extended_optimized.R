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
library(mclust)
library(timeROC)
library(plotly)
library(cmprsk)
library(rms)
library(forestplot)
library(R6)
library(digest)

# Initialize logging and directories
if (!dir.exists("logs")) dir.create("logs")
flog.appender(appender.file("logs/survival_analysis.log"))
flog.threshold(DEBUG)

# Create necessary directories
required_dirs <- c("cache", "logs", "results_opt")
for(dir in required_dirs) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
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

# Add the find_exon15_junctions function
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

# Cache Manager class (from your working implementation)
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
    
    save_data = function(key, data) {
      cache_file <- file.path(self$cache_dir, paste0(key, ".rds"))
      saveRDS(data, cache_file)
      futile.logger::flog.info("Saved data to cache: %s", key)
    }
  )
)

# Initialize cache manager
cache_mgr <- CacheManager$new()

# Function to get expression data (based on your SRRM3_and_SRRM4_expressions.R)
get_expression_data_optimized <- function(cancer_type, gene, timeout = 300) {
  futile.logger::flog.info("Starting expression data retrieval for %s in %s", gene, cancer_type)
  
  cache_key <- paste("expression", cancer_type, gene, sep = "_")
  cached_data <- cache_mgr$get_data(cache_key)
  if (!is.null(cached_data)) {
    futile.logger::flog.info("Using cached expression data")
    return(cached_data)
  }
  
  with_error_handling({
    # Set timeout
    old_timeout <- options("timeout")
    on.exit(options(old_timeout))
    options(timeout = timeout)
    
    # First try to get project info
    project <- paste0("TCGA-", cancer_type)
    
    # Create query
    query <- GDCquery(
      project = project,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      sample.type = c("Primary Tumor")
    )
    
    # Download and prepare data
    GDCdownload(query)
    expr_data <- GDCprepare(query)
    
    # Get gene expression
    gene_idx <- which(rowData(expr_data)$gene_name == gene)
    if (length(gene_idx) == 0) {
      stop(paste("Gene", gene, "not found in dataset"))
    }
    
    expression_data <- data.frame(
      sample_id = colnames(expr_data),
      expression = assay(expr_data)[gene_idx, ],
      stringsAsFactors = FALSE
    )
    
    # Cache the results
    cache_mgr$save_data(cache_key, expression_data)
    
    return(expression_data)
  })
}

# Function to get PSI values
get_psi_data_optimized <- function(cancer_type) {
  futile.logger::flog.info("Starting PSI data retrieval for %s", cancer_type)
  
  cache_key <- paste("psi", cancer_type, sep = "_")
  cached_data <- cache_mgr$get_data(cache_key)
  if (!is.null(cached_data)) {
    futile.logger::flog.info("Using cached PSI data")
    return(cached_data)
  }
  
  # Create query for junction data
  project <- paste0("TCGA-", cancer_type)
  
  # Create region around exon 15
  region <- GRanges(
    seqnames = SRRM3_INFO$gene$chr,
    ranges = IRanges(
      start = SRRM3_INFO$exon15$start - 5000,
      end = SRRM3_INFO$exon15$end + 5000
    )
  )
  
  # Get junction data using recount3
  rse <- create_rse(
    project = project,
    type = "jxn",
    jxn_format = "UNIQUE",
    verbose = TRUE
  )
  
  # Filter to relevant region
  relevant_rows <- which(
    as.character(seqnames(rowRanges(rse))) == as.character(seqnames(region)) &
      start(rowRanges(rse)) >= (start(region)) &
      end(rowRanges(rse)) <= (end(region))
  )
  
  if(length(relevant_rows) == 0) {
    futile.logger::flog.error("No relevant junctions found for %s", cancer_type)
    return(NULL)
  }
  
  rse_filtered <- rse[relevant_rows, ]
  
  # Calculate PSI values
  psi_values <- calculate_exon15_psi(rse_filtered)
  
  if(is.null(psi_values)) {
    futile.logger::flog.error("Failed to calculate PSI values for %s", cancer_type)
    return(NULL)
  }
  
  # Create data frame with results
  cancer_data <- data.frame(
    sample_id = colnames(rse_filtered),
    psi = psi_values,
    stringsAsFactors = FALSE
  ) %>%
    # Convert TCGA sample IDs to match clinical data format
    mutate(
      sample_id = substr(sample_id, 1, 12)  # Take first 12 characters of TCGA ID
    )
  
  # Cache the results
  cache_mgr$save_data(cache_key, cancer_data)
  
  return(cancer_data)
}

# Function to save analysis results
save_analysis_results <- function(results, cancer_type, analysis_type, survival_type) {
  if (is.null(results)) {
    futile.logger::flog.warn("No results to save for %s %s %s", 
                            cancer_type, analysis_type, survival_type)
    return(FALSE)
  }
  
  tryCatch({
    # Create results directory if it doesn't exist
    results_dir <- file.path("results_opt", cancer_type)
    if (!dir.exists(results_dir)) {
      dir.create(results_dir, recursive = TRUE)
    }
    
    # Base filename
    base_filename <- sprintf("%s_%s_%s", cancer_type, analysis_type, survival_type)
    
    # Save plot if it exists
    if (!is.null(results$plot)) {
      plot_file <- file.path(results_dir, paste0(base_filename, "_survival_plot.pdf"))
      futile.logger::flog.info("Saving plot to %s", plot_file)
      
      tryCatch({
        pdf(plot_file, width = 10, height = 8)
        print(results$plot)
        dev.off()
      }, error = function(e) {
        futile.logger::flog.error("Failed to save plot: %s", e$message)
      })
    }
    
    # Save data
    data_file <- file.path(results_dir, paste0(base_filename, "_survival_data.rds"))
    saveRDS(results, data_file)
    
    # Create and save summary statistics
    if (!is.null(results$fit) && !is.null(results$data)) {
      fit_summary <- summary(results$fit)
      
      # Calculate hazard ratio properly
      coxph_fit <- coxph(Surv(time, event) ~ strata, data = results$data)
      hr_data <- summary(coxph_fit)
      
      summary_data <- data.frame(
        cancer_type = cancer_type,
        analysis_type = analysis_type,
        survival_type = survival_type,
        n_samples = nrow(results$data),
        n_events = sum(results$data$event),
        median_survival_low = fit_summary$table["Low", "median"],
        median_survival_high = fit_summary$table["High", "median"],
        pvalue = results$plot$plot.data$pval,
        hr = exp(coef(coxph_fit))[1],
        hr_conf_low = exp(confint(coxph_fit))[1],
        hr_conf_high = exp(confint(coxph_fit))[2],
        low_group_size = sum(results$data$strata == "Low"),
        high_group_size = sum(results$data$strata == "High"),
        analysis_date = format(Sys.time(), "%Y-%m-%d")
      )
      
      summary_file <- file.path(results_dir, paste0(base_filename, "_summary.csv"))
      write.csv(summary_data, summary_file, row.names = FALSE)
      
      futile.logger::flog.info("Saved results to %s", results_dir)
      return(TRUE)
    }
  }, error = function(e) {
    futile.logger::flog.error("Failed to save results: %s", e$message)
    return(FALSE)
  })
}

# Function to perform survival analysis
perform_survival_analysis_enhanced <- function(cancer_type, 
                                            analysis_type = c("SRRM3_expression", "SRRM3_PSI", "SRRM4_expression"),
                                            survival_type = c("OS", "PFS")) {
  
  futile.logger::flog.info("Starting %s analysis for %s", analysis_type, cancer_type)
  
  # Create cache key for the entire analysis
  cache_key <- paste("survival", cancer_type, analysis_type, survival_type, sep = "_")
  cached_results <- cache_mgr$get_data(cache_key)
  if (!is.null(cached_results)) {
    return(cached_results)
  }
  
  # Get appropriate data based on analysis type
  if (grepl("expression", analysis_type)) {
    gene <- sub("_expression", "", analysis_type)
    molecular_data <- get_expression_data_optimized(cancer_type, gene)
    value_col <- "expression"
  } else {
    molecular_data <- get_psi_data_optimized(cancer_type)
    value_col <- "psi"
  }
  
  # Get clinical data
  clinical <- GDCquery_clinic(paste0("TCGA-", cancer_type))
  
  # Process survival data based on type
  if (survival_type == "OS") {
    # Overall Survival
    clinical <- clinical %>%
      mutate(
        time = ifelse(vital_status == "Alive",
                     days_to_last_follow_up,
                     days_to_death),
        event = ifelse(vital_status == "Dead", 1, 0)
      ) %>%
      filter(!is.na(time), time > 0)  # Remove invalid times
      
  } else {
    # Progression-Free Survival - using disease-free survival data
    clinical <- clinical %>%
      mutate(
        time = days_to_last_follow_up,  # Default to last follow-up
        event = 0  # Default to no event
      ) %>%
      mutate(
        # Update time if progression event exists
        time = case_when(
          !is.na(days_to_progression) ~ days_to_progression,
          !is.na(days_to_recurrence) ~ days_to_recurrence,
          TRUE ~ time
        ),
        # Update event status
        event = case_when(
          !is.na(progression_status) & progression_status == "YES" ~ 1,
          !is.na(disease_free_status) & disease_free_status == "Recurred/Progressed" ~ 1,
          TRUE ~ event
        )
      ) %>%
      filter(!is.na(time), time > 0)  # Remove invalid times
  }
  
  # Check if we have enough data
  if (nrow(clinical) < 10) {
    futile.logger::flog.error("Insufficient clinical data for %s analysis", survival_type)
    return(NULL)
  }
  
  # Merge clinical and molecular data
  analysis_data <- merge(molecular_data, clinical, 
                        by.x = "sample_id", 
                        by.y = "submitter_id")
  
  if (nrow(analysis_data) < 10) {
    futile.logger::flog.error("Insufficient matched data after merging")
    return(NULL)
  }
  
  # Create stratification based on median
  analysis_data <- analysis_data %>%
    mutate(
      strata = factor(ifelse(get(value_col) > median(get(value_col)), 
                            "High", "Low"),
                     levels = c("Low", "High"))
    )
  
  # Create survival object and fit
  surv_obj <- Surv(analysis_data$time, analysis_data$event)
  fit <- survfit(surv_obj ~ strata, data = analysis_data)
  
  # Create plot
  plot <- ggsurvplot(
    fit,
    data = analysis_data,
    pval = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    xlab = paste(survival_type, "Time (days)"),
    ylab = paste(survival_type, "Probability"),
    title = paste(analysis_type, "Survival Analysis"),
    subtitle = paste("Cancer type:", cancer_type),
    palette = c("#E7B800", "#2E9FDF"),
    legend.title = "Expression Level",
    legend.labs = c("Low", "High"),
    font.main = c(14, "bold"),
    font.x = c(12),
    font.y = c(12),
    font.legend = c(10)
  )
  
  results <- list(
    plot = plot,
    fit = fit,
    data = analysis_data,
    survival_type = survival_type,
    analysis_type = analysis_type
  )
  
  # Save results
  save_analysis_results(results, cancer_type, analysis_type, survival_type)
  
  # Cache results
  cache_mgr$save_data(cache_key, results)
  
  return(results)
}

# Main execution
analysis_types <- c("SRRM3_expression", "SRRM3_PSI", "SRRM4_expression")
survival_types <- c("OS", "PFS")

# Initialize progress bar
total_analyses <- length(analysis_types) * length(survival_types)
pb <- progress::progress_bar$new(
  format = "Analysis [:bar] :percent | ETA: :eta | :current/:total | :message",
  total = total_analyses,
  clear = FALSE
)

# Before running analyses, add validation
cancer_type <- "BRCA"
futile.logger::flog.info("Validating data availability for %s", cancer_type)

# Validate clinical data
clinical_data <- tryCatch({
  GDCquery_clinic(paste0("TCGA-", cancer_type))
}, error = function(e) {
  futile.logger::flog.error("Failed to retrieve clinical data: %s", e$message)
  return(NULL)
})

if (is.null(clinical_data) || nrow(clinical_data) == 0) {
  stop("No clinical data available for ", cancer_type)
}

# Validate PSI data
psi_data <- get_psi_data_optimized(cancer_type)
if (is.null(psi_data)) {
  futile.logger::flog.warn("No PSI data available for %s", cancer_type)
}

# Run all analyses
for(analysis_type in analysis_types) {
  for(survival_type in survival_types) {
    tryCatch({
      message(sprintf("\nRunning %s_%s analysis...", analysis_type, survival_type))
      pb$tick(tokens = list(message = sprintf("Running %s_%s", analysis_type, survival_type)))
      
      results <- perform_survival_analysis_enhanced(
        cancer_type = "BRCA",
        analysis_type = analysis_type,
        survival_type = survival_type
      )
      
      if (!is.null(results)) {
        message(sprintf("Completed %s %s analysis", analysis_type, survival_type))
      } else {
        message(sprintf("Failed %s %s analysis", analysis_type, survival_type))
      }
      
    }, error = function(e) {
      message(sprintf("Error in %s %s analysis: %s", 
                     analysis_type, survival_type, e$message))
    })
    
    # Add a small delay between analyses
    Sys.sleep(1)
  }
}

# Create summary of all analyses
create_analysis_summary <- function() {
  tryCatch({
    # Get all summary files
    results_files <- list.files("results_opt", 
                              pattern = "*_summary.csv$", 
                              recursive = TRUE, 
                              full.names = TRUE)
    
    if (length(results_files) == 0) {
      futile.logger::flog.error("No summary files found in results_opt directory")
      return(NULL)
    }
    
    # Read and combine summaries
    summaries <- lapply(results_files, function(file) {
      tryCatch({
        read.csv(file)
      }, error = function(e) {
        futile.logger::flog.warn("Failed to read %s: %s", file, e$message)
        return(NULL)
      })
    })
    
    # Remove NULL entries
    summaries <- summaries[!sapply(summaries, is.null)]
    
    if (length(summaries) == 0) {
      futile.logger::flog.error("No valid summary data found")
      return(NULL)
    }
    
    # Combine all summaries
    all_summaries <- do.call(rbind, summaries)
    
    # Save complete summary
    output_file <- file.path("results_opt", "complete_survival_analysis_summary.csv")
    write.csv(all_summaries, output_file, row.names = FALSE)
    
    futile.logger::flog.info("Created complete summary at %s", output_file)
    return(all_summaries)
    
  }, error = function(e) {
    futile.logger::flog.error("Failed to create analysis summary: %s", e$message)
    return(NULL)
  })
}

# Main execution remains the same, but add error handling for the final summary
tryCatch({
  summary_results <- create_analysis_summary()
  if (!is.null(summary_results)) {
    cat("Analysis complete. Check results_opt directory for outputs.\n")
  } else {
    cat("Analysis complete but summary creation failed. Check individual results.\n")
  }
}, error = function(e) {
  cat("Error in final summary creation:", e$message, "\n")
})

calculate_exon15_psi <- function(rse) {
  if(is.null(rse)) return(NULL)
  
  # Get junction counts and coordinates
  junction_counts <- assay(rse)
  jxn_coords <- rowRanges(rse)
  
  # Find relevant junctions
  junctions <- find_exon15_junctions(jxn_coords)
  
  if(length(junctions$inclusion) == 0 || length(junctions$exclusion) == 0) {
    return(NULL)
  }
  
  # Calculate PSI for each sample
  psi_values <- sapply(seq_len(ncol(junction_counts)), function(i) {
    # Sum reads supporting exon 15 inclusion
    inclusion_reads <- sum(junction_counts[junctions$inclusion, i])
    
    # Sum reads supporting exon 15 skipping
    exclusion_reads <- sum(junction_counts[junctions$exclusion, i])
    
    # Total reads covering this splicing event
    total_reads <- inclusion_reads + exclusion_reads
    
    # Calculate PSI if we have sufficient coverage
    if(total_reads >= 10) {  # Minimum coverage threshold
      psi <- (inclusion_reads / total_reads) * 100
      return(psi)
    } else {
      return(NA)
    }
  })
  
  return(psi_values)
}