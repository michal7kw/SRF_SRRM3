# Analysis of SRRM3 isoform expression across cancer types
library(recount3)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(SummarizedExperiment)
library(futile.logger)
library(viridis)

# Set up logging
flog.appender(appender.file("srrm3_cancer_analysis.log"))
flog.threshold(DEBUG)

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

# Simplified RSE creation function
create_rse_safe <- function(project_info) {  
  tryCatch({
    # Create region around exon 15
    region <- GRanges(
      seqnames = SRRM3_INFO$gene$chr,
      ranges = IRanges(
        start = SRRM3_INFO$exon15$start - 5000,  # Include surrounding region
        end = SRRM3_INFO$exon15$end + 5000
      )
    )
    
    rse <- create_rse(
      project_info,
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
      return(NULL)
    }
    
    rse_filtered <- rse[relevant_rows, ]
    rm(rse)
    gc()
    
    return(rse_filtered)
  }, error = function(e) {
    flog.error("Error creating RSE: %s", e$message)
    return(NULL)
  })
}

# Calculate PSI for exon 15
calculate_exon15_psi <- function(rse) {
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
    # Get read counts
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
  
  return(list(
    psi = psi_values,
    junctions = list(
      inclusion = rownames(junction_counts)[junctions$inclusion],
      exclusion = rownames(junction_counts)[junctions$exclusion]
    )
  ))
}

# Function to get cancer type from TCGA project name
get_cancer_type <- function(project) {
  gsub("TCGA-", "", project)
}

# Function to analyze SRRM3 isoforms across cancer types
analyze_cancer_types <- function() {
  # Get available TCGA projects for human
  projects <- available_projects(organism = "human")
  
  # Filter for TCGA projects only
  tcga_projects <- projects[projects$file_source == "tcga", ]
  
  # Print available projects for debugging
  flog.info("Found %d TCGA projects", nrow(tcga_projects))
  for(i in 1:nrow(tcga_projects)) {
    flog.info("Project %d: %s", i, tcga_projects$project[i])
  }
  
  # Store results
  results_list <- list()
  
  # Progress tracking
  total_projects <- nrow(tcga_projects)
  flog.info("Starting analysis of %d TCGA projects", total_projects)
  
  if(total_projects == 0) {
    flog.error("No TCGA projects found")
    return(NULL)
  }
  
  for(i in 1:nrow(tcga_projects)) {
    project_info <- tcga_projects[i,]
    project_name <- project_info$project
    cancer_type <- get_cancer_type(project_name)
    
    flog.info("==========================================")
    flog.info("Processing project %d/%d: %s", i, total_projects, project_name)
    
    tryCatch({
      rse <- create_rse_safe(project_info)
      if(is.null(rse) || ncol(rse) == 0 || nrow(rse) == 0) {
        flog.warn("Invalid RSE object for %s", project_name)
        next
      }
      
      psi_result <- calculate_exon15_psi(rse)
      if(is.null(psi_result) || all(is.na(psi_result$psi))) {
        flog.warn("No valid PSI values for %s", project_name)
        next
      }
      
      # Store results
      results_list[[project_name]] <- data.frame(
        cancer_type = cancer_type,
        sample_id = colnames(rse),
        psi = psi_result$psi,
        stringsAsFactors = FALSE
      )
      
      flog.info("Processed %s: %d samples", project_name, sum(!is.na(psi_result$psi)))
      
    }, error = function(e) {
      flog.error("Error processing %s: %s", project_name, e$message)
    })
  }
  
  if(length(results_list) == 0) {
    flog.error("No results were generated - all projects failed processing")
    return(NULL)
  }
  
  all_results <- do.call(rbind, results_list)
  all_results <- all_results[!is.na(all_results$psi), ]
  
  if(nrow(all_results) == 0) {
    flog.error("No valid results after removing NA values")
    return(NULL)
  }
  
  flog.info("Analysis complete: %d samples across %d cancer types", 
            nrow(all_results), length(unique(all_results$cancer_type)))
  
  return(all_results)
}

# Function to plot cancer type distributions
plot_cancer_distributions <- function(results) {
  # Calculate summary statistics for cancer types with at least 10 samples
  summary_stats <- results %>%
    group_by(cancer_type) %>%
    summarise(
      median_psi = median(psi, na.rm = TRUE),
      mean_psi = mean(psi, na.rm = TRUE),
      sd_psi = sd(psi, na.rm = TRUE),
      n_samples = n(),
      .groups = "drop"
    ) %>%
    filter(n_samples >= 10)
  
  # Create violin plot
  p1 <- ggplot(results %>% filter(cancer_type %in% summary_stats$cancer_type), 
               aes(x = reorder(cancer_type, psi, FUN = median), y = psi)) +
    geom_violin(aes(fill = cancer_type), alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    scale_fill_viridis(discrete = TRUE) +
    labs(
      title = "SRRM3 Isoform Distribution Across Cancer Types",
      subtitle = "Ratio of long (16 exons) to short (15 exons) isoform",
      x = "Cancer Type",
      y = "PSI Value"
    )
  
  # Create summary plot
  p2 <- ggplot(summary_stats, 
               aes(x = reorder(cancer_type, median_psi), y = median_psi)) +
    geom_point(aes(size = n_samples), alpha = 0.7) +
    geom_errorbar(aes(ymin = median_psi - sd_psi, 
                     ymax = median_psi + sd_psi), 
                 width = 0.2) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = "SRRM3 Isoform Distribution Summary by Cancer Type",
      subtitle = "Median PSI values with standard deviation",
      x = "Cancer Type",
      y = "Median PSI Value",
      size = "Number of Samples"
    )
  
  return(list(distribution_plot = p1, summary_plot = p2))
}

# Main execution
flog.info("Starting SRRM3 isoform analysis across cancer types")

results <- analyze_cancer_types()

if(!is.null(results) && nrow(results) > 0) {
  # Save results
  write.csv(results, "srrm3_isoforms_cancer_types.csv", row.names = FALSE)
  
  # Create and save plots
  plots <- plot_cancer_distributions(results)
  ggsave("srrm3_isoforms_distribution.pdf", plots$distribution_plot, width = 12, height = 8)
  ggsave("srrm3_isoforms_summary.pdf", plots$summary_plot, width = 12, height = 8)
  
  flog.info("Analysis complete. Results and plots saved.")
} else {
  flog.error("Analysis failed - no valid results to plot")
}
