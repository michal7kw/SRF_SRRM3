### Output #################################
# "SRRM3_and_SRRM4_expressions_distributions.pdf"
# "SRRM3_and_SRRM4_expressions_correlation_matrix.pdf"
# "SRRM3_and_SRRM4_expressions_correlation_matrix_simple.pdf"
# "SRRM3_and_SRRM4_expressions_summary.pdf"
# "SRRM3_and_SRRM4_expressions_summary_data.csv"
# "SRRM3_and_SRRM4_expressions_TCGA.rds"
###########################################

library(recount3)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(SummarizedExperiment)
library(pheatmap)
library(DESeq2)
library(survival)
library(survminer)
library(corrplot)
library(reshape2)
library(futile.logger)

# Source common theme
source("common_theme.R")

# Set up logging
flog.appender(appender.file("./logs/SRRM3_and_SRRM4_expressions_TCGA.log"))
flog.threshold(DEBUG)

# Define cache directory
CACHE_DIR <- "../cache"
if (!dir.exists(CACHE_DIR)) {
  dir.create(CACHE_DIR)
}

# Function to get cached data
get_cached_data <- function(cache_file) {
  if (file.exists(cache_file)) {
    flog.info("Loading cached data from %s", cache_file)
    return(readRDS(cache_file))
  }
  return(NULL)
}

# Function to save cached data
save_cached_data <- function(data, cache_file) {
  flog.info("Saving data to cache: %s", cache_file)
  saveRDS(data, cache_file)
}

# Define gene information
SRRM3_INFO <- list(
  gene = list(
    name = "SRRM3",
    chr = "chr7",
    start = 76201896,
    end = 76287287,
    ensembl_id = "ENSG00000177679"
  )
)

SRRM4_INFO <- list(
  gene = list(
    name = "SRRM4",
    chr = "chr12",
    start = 118587761,
    end = 118708623,
    ensembl_id = "ENSG00000139767"
  )
)

# Function to get available TCGA projects
get_tcga_projects <- function() {
  message("Checking available projects...")
  projects <- available_projects()
  tcga_projects <- subset(projects, file_source == "tcga")
  return(tcga_projects)
}

# Modify get_cancer_type function to match other TCGA scripts
get_cancer_type <- function(project_info, use_full_names = FALSE) {
  if (use_full_names) {
    # Get histological diagnosis from clinical data
    cancer_type <- project_info$tcga.cgc_case_histological_diagnosis
    if (is.null(cancer_type) || is.na(cancer_type) || cancer_type == "") {
      # Fallback to project code if no histological diagnosis
      cancer_type <- sub("^TCGA-", "", project_info$project)
    }
  } else {
    # Use short TCGA code (e.g., "ACC", "BLCA")
    cancer_type <- sub("^TCGA-", "", project_info$project)
  }
  return(cancer_type)
}

# Modified get_expression_data function
get_expression_data <- function(project_info) {
  project_name <- project_info$project
  cache_file <- file.path(CACHE_DIR, paste0("srrm_expr_", project_name, ".rds"))
  
  # Check cache first
  cached_data <- get_cached_data(cache_file)
  if (!is.null(cached_data)) {
    return(cached_data)
  }
  
  flog.info("\nProcessing project: %s", project_name)
  
  tryCatch({
    # Create RSE object
    rse_gene <- create_rse(subset(available_projects(), 
                                 project == project_info$project & 
                                   file_source == "tcga"), 
                          type = "gene")
    
    # Get SRRM3 and SRRM4 expression
    srrm3_idx <- which(rowData(rse_gene)$gene_name == "SRRM3")
    srrm4_idx <- which(rowData(rse_gene)$gene_name == "SRRM4")
    
    if(length(srrm3_idx) == 0 || length(srrm4_idx) == 0) {
      flog.warn("SRRM3 or SRRM4 not found in %s", project_name)
      return(NULL)
    }
    
    # Extract expression data and clinical info
    data <- data.frame(
      sample_id = colnames(rse_gene),
      srrm3_expr = assay(rse_gene)[srrm3_idx, ],
      srrm4_expr = assay(rse_gene)[srrm4_idx, ],
      stringsAsFactors = FALSE
    )
    
    # Add cancer type information
    data$cancer_type_short <- sub("^TCGA-", "", project_info$project)
    data$cancer_type_full <- NA
    
    # Get histological diagnosis from clinical data
    clinical <- as.data.frame(colData(rse_gene))
    if ("tcga.cgc_case_histological_diagnosis" %in% colnames(clinical)) {
      hist_diag <- clinical$tcga.cgc_case_histological_diagnosis
      if (!all(is.na(hist_diag))) {
        data$cancer_type_full <- hist_diag
      }
    }
    
    # If no histological diagnosis available, use short name
    data$cancer_type_full[is.na(data$cancer_type_full)] <- data$cancer_type_short[is.na(data$cancer_type_full)]
    
    # Cache the results
    save_cached_data(data, cache_file)
    
    return(data)
  }, error = function(e) {
    flog.error("Error processing project %s: %s", project_name, e$message)
    return(NULL)
  })
}

# Function to analyze expression patterns
analyze_expression_patterns <- function(results) {
  # Create two analysis data frames
  analysis_data_short <- data.frame()
  analysis_data_full <- data.frame()
  
  for(result in results) {
    if(is.null(result)) next
    
    # Calculate statistics
    stats <- data.frame(
      srrm3_mean = mean(result$srrm3_expr),
      srrm3_sd = sd(result$srrm3_expr),
      srrm4_mean = mean(result$srrm4_expr),
      srrm4_sd = sd(result$srrm4_expr),
      correlation = cor(result$srrm3_expr, result$srrm4_expr),
      sample_size = length(result$srrm3_expr)
    )
    
    # Add to both datasets with appropriate cancer type
    analysis_data_short <- rbind(
      analysis_data_short,
      cbind(cancer_type = result$cancer_type_short, stats)
    )
    
    analysis_data_full <- rbind(
      analysis_data_full,
      cbind(cancer_type = result$cancer_type_full, stats)
    )
  }
  
  return(list(
    short_names = analysis_data_short,
    full_names = analysis_data_full
  ))
}

# Visualization functions


plot_expression_summary <- function(analysis_data, use_full_names = FALSE) {
  # Select appropriate data
  data_to_plot <- if(use_full_names) analysis_data$full_names else analysis_data$short_names
  name_type <- if(use_full_names) "full_names" else "short_names"
  
  if(is.null(data_to_plot) || nrow(data_to_plot) == 0) {
    stop(sprintf("No data available for %s", name_type))
  }
  
  # Create plot data
  plot_data <- data_to_plot %>%
    tidyr::pivot_longer(
      cols = c("srrm3_mean", "srrm4_mean"),
      names_to = "gene",
      values_to = "mean"
    ) %>%
    mutate(
      gene = ifelse(gene == "srrm3_mean", "SRRM3", "SRRM4"),
      sd = ifelse(gene == "SRRM3", srrm3_sd, srrm4_sd),
      se = sd / sqrt(sample_size)
    )
  
  # Create plot
  p <- ggplot(plot_data, 
              aes(x = reorder(cancer_type, mean), 
                  y = mean, 
                  fill = gene)) +
    geom_bar(stat = "identity", 
             position = position_dodge(0.9)) +
    geom_errorbar(aes(ymin = mean - se, 
                      ymax = mean + se),
                  position = position_dodge(0.9),
                  width = 0.25) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "SRRM3 and SRRM4 Expression by Cancer Type",
         x = "Cancer Type",
         y = "Mean Expression (normalized counts)",
         fill = "Gene")
  
  # Save plot
  ggsave(sprintf("./output/SRRM3_and_SRRM4_expressions_%s.pdf", name_type),
         p, width = 15, height = 8)
  
  return(p)
}

# Modified main_analysis function
main_analysis <- function() {
  # Check for cached final results
  final_cache_file <- file.path(CACHE_DIR, "SRRM3_and_SRRM4_expressions.rds")
  cached_results <- get_cached_data(final_cache_file)
  
  if (!is.null(cached_results)) {
    flog.info("Using cached final results")
    return(cached_results)
  }
  
  # Get TCGA projects
  projects <- get_tcga_projects()
  flog.info("Found %d TCGA projects", nrow(projects))
  
  # Process each project
  results <- lapply(seq_len(nrow(projects)), function(i) {
    flog.info("Processing project %d/%d", i, nrow(projects))
    get_expression_data(projects[i,])
  })
  
  # Combine all results
  results <- do.call(rbind, results[!sapply(results, is.null)])
  
  if(is.null(results) || nrow(results) == 0) {
    stop("No valid results to analyze")
  }
  
  # Calculate statistics for both naming schemes
  analysis_data <- list(
    short_names = results %>%
      group_by(cancer_type_short) %>%
      summarise(
        srrm3_mean = mean(srrm3_expr, na.rm = TRUE),
        srrm3_sd = sd(srrm3_expr, na.rm = TRUE),
        srrm4_mean = mean(srrm4_expr, na.rm = TRUE),
        srrm4_sd = sd(srrm4_expr, na.rm = TRUE),
        correlation = cor(srrm3_expr, srrm4_expr, use = "pairwise.complete.obs"),
        sample_size = n(),
        .groups = "drop"
      ) %>%
      rename(cancer_type = cancer_type_short),
    
    full_names = results %>%
      group_by(cancer_type_full) %>%
      summarise(
        srrm3_mean = mean(srrm3_expr, na.rm = TRUE),
        srrm3_sd = sd(srrm3_expr, na.rm = TRUE),
        srrm4_mean = mean(srrm4_expr, na.rm = TRUE),
        srrm4_sd = sd(srrm4_expr, na.rm = TRUE),
        correlation = cor(srrm3_expr, srrm4_expr, use = "pairwise.complete.obs"),
        sample_size = n(),
        .groups = "drop"
      ) %>%
      rename(cancer_type = cancer_type_full)
  )
  
  # Save raw data and analysis results
  final_results <- list(
    raw_results = results,
    analysis_data = analysis_data
  )
  
  save_cached_data(final_results, final_cache_file)
  
  return(final_results)
}

# Run the analysis with caching
results <- main_analysis()

# Create plots with both naming schemes
plot_short <- plot_expression_summary(results$analysis_data, use_full_names = FALSE)
plot_full <- plot_expression_summary(results$analysis_data, use_full_names = TRUE)

# plot_correlation_matrix(results$analysis_data)
# plot_correlation_matrix <- function(analysis_data) {
#   # Create correlation matrix
#   cor_matrix <- matrix(
#     analysis_data$correlation,
#     nrow = 1,
#     dimnames = list("SRRM3-SRRM4", analysis_data$project)
#   )
  
#   # Create custom color palette from blue to red
#   col_palette <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  
#   # Plot main correlation matrix
#   pdf("./output/SRRM3_and_SRRM4_expressions_correlation_matrix.pdf", width = 12, height = 4)  # Adjusted dimensions for horizontal layout
#   par(mar = c(4, 8, 4, 4))  # Adjusted margins
  
#   corrplot(cor_matrix,  # Removed t() to make it horizontal
#            method = "color",
#            col = col_palette,
#            addCoef.col = "black",
#            number.cex = 0.65,      # Smaller coefficient size
#            tl.col = "black",
#            tl.srt = 0,            # Horizontal text
#            tl.cex = 0.7,          # Smaller label size
#            cl.cex = 0.7,          # Smaller color scale labels
#            cl.ratio = 0.15,       # Thinner color key
#            title = "SRRM3-SRRM4 Expression Correlation by Cancer Type",
#            mar = c(0,0,2,0),
#            addgrid.col = "gray80",
#            is.corr = FALSE)
  
#   # Add significance stars with adjusted positioning
#   if("sample_size" %in% colnames(analysis_data)) {
#     p_values <- sapply(1:nrow(analysis_data), function(i) {
#       n <- analysis_data$sample_size[i]
#       r <- analysis_data$correlation[i]
#       if(!is.na(r) && !is.na(n) && n > 2) {
#         z <- 0.5 * log((1 + r)/(1 - r))
#         p <- 2 * (1 - pnorm(abs(z) * sqrt(n - 3)))
#         return(p)
#       }
#       return(NA)
#     })
    
#     stars <- sapply(p_values, function(p) {
#       if(is.na(p)) return("")
#       if(p < 0.001) return("***")
#       if(p < 0.01) return("**")
#       if(p < 0.05) return("*")
#       return("")
#     })
    
#     text(x = rep(1, ncol(cor_matrix)), 
#          y = 1:ncol(cor_matrix), 
#          labels = stars,
#          pos = 4, 
#          cex = 0.7,
#          offset = 0.5)
#   }
  
#   # Add legend with adjusted positioning
#   if(any(p_values < 0.05, na.rm = TRUE)) {
#     legend("bottom", 
#            legend = c("p < 0.05 *", "p < 0.01 **", "p < 0.001 ***"),
#            bty = "n",
#            horiz = TRUE,
#            cex = 0.6,
#            inset = c(0, -0.1),
#            xpd = TRUE)
#   }
  
#   dev.off()
  
#   # Create simplified version
#   pdf("./output/SRRM3_and_SRRM4_expressions_correlation_matrix_simple.pdf", width = 10, height = 3)
#   par(mar = c(4, 8, 4, 4))
  
#   corrplot(cor_matrix, 
#            method = "color",
#            col = col_palette,
#            addCoef.col = "black",
#            number.cex = 0.65,
#            tl.col = "black",
#            tl.srt = 0,
#            tl.cex = 0.7,
#            cl.cex = 0.7,
#            title = "SRRM3-SRRM4 Correlation",
#            addgrid.col = "gray80",
#            is.corr = FALSE)
  
#   dev.off()
# }


# plot_expression_summary(results$analysis_data)
# plot_expression_distributions <- function(results) {
#   # Prepare data for plotting
#   plot_data <- data.frame()
  
#   for(result in results) {
#     if(is.null(result)) next
    
#     project_data <- data.frame(
#       project = result$project,
#       SRRM3 = result$srrm3_expr,
#       SRRM4 = result$srrm4_expr
#     )
    
#     plot_data <- rbind(plot_data, 
#                        data.frame(
#                          project = project_data$project,
#                          gene = "SRRM3",
#                          expression = project_data$SRRM3
#                        ),
#                        data.frame(
#                          project = project_data$project,
#                          gene = "SRRM4",
#                          expression = project_data$SRRM4
#                        ))
#   }
  
#   # Create violin plot
#   ggplot(plot_data, aes(x = project, y = log2(expression + 1), fill = gene)) +
#     geom_violin(position = position_dodge(width = 0.7), alpha = 0.7) +
#     geom_boxplot(position = position_dodge(width = 0.7), width = 0.1, alpha = 0.7) +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     labs(
#       title = "SRRM3 and SRRM4 Expression Distribution Across Cancer Types",
#       x = "Cancer Type",
#       y = "log2(normalized counts + 1)",
#       fill = "Gene"
#     )
    
#   ggsave("./output/SRRM3_and_SRRM4_expressions_distributions.pdf", width = 15, height = 8)
# }