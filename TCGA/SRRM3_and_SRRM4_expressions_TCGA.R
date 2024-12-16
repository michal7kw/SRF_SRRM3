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

# Modified get_expression_data function with caching
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
    
    # Extract expression data
    expr_data <- list(
      project = project_name,
      srrm3_expr = assay(rse_gene)[srrm3_idx, ],
      srrm4_expr = assay(rse_gene)[srrm4_idx, ],
      clinical = colData(rse_gene)
    )
    
    # Cache the results
    save_cached_data(expr_data, cache_file)
    
    return(expr_data)
  }, error = function(e) {
    flog.error("Error processing project %s: %s", project_name, e$message)
    return(NULL)
  })
}

# Function to analyze expression patterns
analyze_expression_patterns <- function(results) {
  # Remove NULL results
  results <- results[!sapply(results, is.null)]
  
  if(length(results) == 0) {
    stop("No valid results to analyze")
  }
  
  # Create data frame for analysis
  analysis_data <- data.frame(
    project = character(),
    srrm3_mean = numeric(),
    srrm3_sd = numeric(),
    srrm4_mean = numeric(),
    srrm4_sd = numeric(),
    correlation = numeric(),
    sample_size = numeric()
  )
  
  # Calculate statistics for each project
  for(result in results) {
    project_stats <- data.frame(
      project = result$project,
      srrm3_mean = mean(result$srrm3_expr, na.rm = TRUE),
      srrm3_sd = sd(result$srrm3_expr, na.rm = TRUE),
      srrm4_mean = mean(result$srrm4_expr, na.rm = TRUE),
      srrm4_sd = sd(result$srrm4_expr, na.rm = TRUE),
      correlation = cor(result$srrm3_expr, result$srrm4_expr, 
                        use = "pairwise.complete.obs"),
      sample_size = length(result$srrm3_expr)
    )
    analysis_data <- rbind(analysis_data, project_stats)
  }
  
  return(analysis_data)
}

# Visualization functions


plot_expression_summary <- function(analysis_data) {
  # Calculate standard error
  analysis_data$srrm3_se <- analysis_data$srrm3_sd / sqrt(analysis_data$sample_size)
  analysis_data$srrm4_se <- analysis_data$srrm4_sd / sqrt(analysis_data$sample_size)
  
  # First, create long format for means
  means_long <- reshape2::melt(analysis_data, 
                               id.vars = "project",
                               measure.vars = c("srrm3_mean", "srrm4_mean"),
                               variable.name = "gene",
                               value.name = "mean")
  
  # Create long format for standard errors
  se_long <- reshape2::melt(analysis_data, 
                            id.vars = "project",
                            measure.vars = c("srrm3_se", "srrm4_se"),
                            variable.name = "gene",
                            value.name = "se")
  
  # Combine the data
  plot_data <- means_long
  plot_data$se <- se_long$se
  plot_data$sample_size <- rep(analysis_data$sample_size, each = 2)
  
  # Clean up gene names
  plot_data$gene <- gsub("_mean", "", plot_data$gene)
  plot_data$gene <- toupper(gsub("srrm", "SRRM", plot_data$gene))
  
  # Order projects by mean expression of SRRM3
  project_order <- analysis_data[order(-analysis_data$srrm3_mean), "project"]
  plot_data$project <- factor(plot_data$project, levels = project_order)
  
  # Create sample size data frame for labels
  sample_labels <- unique(analysis_data[, c("project", "sample_size")])
  sample_labels$project <- factor(sample_labels$project, levels = project_order)
  
  # Calculate max y value for each project to position labels
  max_values <- aggregate(mean + se ~ project, data = plot_data, max)
  sample_labels$y_pos <- max_values$`mean + se`
  
  # Calculate total samples for subtitle
  total_samples <- sum(analysis_data$sample_size)
  sample_text <- sprintf("Total samples: %d", total_samples)
  
  # Create enhanced plot
  p <- ggplot(plot_data, 
              aes(x = project, y = mean, fill = gene)) +
    # Add bars
    geom_bar(stat = "identity", 
             position = position_dodge(width = 0.8),
             width = 0.7,
             alpha = 0.8) +
    # Add error bars using standard error
    geom_errorbar(aes(ymin = mean - se, 
                      ymax = mean + se),
                  position = position_dodge(width = 0.8),
                  width = 0.25,
                  color = "black") +
    # Add sample size labels above bars
    geom_text(data = sample_labels,
              aes(x = project, 
                  y = max(y_pos),
                  label = sprintf("%d", sample_size)),
              inherit.aes = FALSE,
              vjust = -0.5,
              size = 3,
              fontface = "bold",
              color = "black") +
    # Customize theme
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, 
                                 hjust = 1, 
                                 size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, 
                                face = "bold", 
                                hjust = 0.5),
      plot.subtitle = element_text(size = 10, 
                                   hjust = 0.5),
      legend.position = "top",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", 
                                  fill = NA, 
                                  linewidth = 0.5)
    ) +
    # Customize labels
    labs(
      title = "Mean Expression Levels of SRRM3 and SRRM4 Across Cancer Types",
      subtitle = paste("Error bars represent ± standard error of the mean\n", sample_text),
      x = "Cancer Type",
      y = "Mean Expression (normalized counts)",
      fill = "Gene"
    ) +
    # Customize colors
    scale_fill_manual(values = c("SRRM3" = "#2166AC", 
                                 "SRRM4" = "#B2182B")) +
    # Adjust y-axis limits to accommodate labels
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))  # Increased top expansion for labels
  
  # Save plot with high resolution
  ggsave("./output/SRRM3_and_SRRM4_expressions_summary.pdf", 
         plot = p,
         width = 15, 
         height = 8,
         dpi = 300)
  
  # Save plot data
  write.csv(plot_data, 
            "./output/SRRM3_and_SRRM4_expressions_summary_data.csv", 
            row.names = FALSE)
  
  return(p)
}

# Modified main analysis function with caching
main_analysis <- function() {
  # Check for cached final results
  final_cache_file <- file.path(CACHE_DIR, "SRRM3_and_SRRM4_expressions_TCGA.rds")
  cached_results <- get_cached_data(final_cache_file)
  
  if (!is.null(cached_results)) {
    flog.info("Using cached final results")
    return(cached_results)
  }
  
  # Get TCGA projects
  projects <- get_tcga_projects()
  flog.info("Found %d TCGA projects", nrow(projects))
  
  # Process each project with progress tracking
  total_projects <- nrow(projects)
  results <- lapply(seq_len(total_projects), function(i) {
    flog.info("Processing project %d/%d", i, total_projects)
    get_expression_data(projects[i,])
  })
  
  # Analyze expression patterns
  analysis_data <- analyze_expression_patterns(results)
  
  # Generate visualizations
  # plot_expression_distributions(results)
  # plot_correlation_matrix(analysis_data)
  # plot_expression_summary(analysis_data)
  
  # Save final results
  final_results <- list(
    raw_results = results,
    analysis_data = analysis_data
  )
  save_cached_data(final_results, final_cache_file)
  
  return(final_results)
}

# Run the analysis with caching
results <- main_analysis()

plot_expression_summary(results$analysis_data)

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