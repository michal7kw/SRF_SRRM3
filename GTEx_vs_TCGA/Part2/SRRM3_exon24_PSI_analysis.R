# Load required libraries
library(recount3)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(SummarizedExperiment)

# Define SRRM3 information with genomic coordinates for exon 24
SRRM3_INFO <- list(
  gene = list(
    name = "SRRM3",
    chr = "chr7",
    start = 76201896,  # GRCh38 coordinates
    end = 76287287
  ),
  exon24 = list(
    start = 76283524,  # Genomic coordinates from GTF
    end = 76283602,
    length = 79
  )
)

# Function to calculate PSI values
calculate_PSI <- function(inclusion_counts, skipping_counts) {
    psi <- inclusion_counts / (inclusion_counts + skipping_counts)
    return(psi)
}

# Function to get junction data from recount3
get_junction_data <- function(project_info) {
  message("\nGetting junction data for project: ", project_info$project)
  
  # Get junction-level data
  rse_jxn <- create_rse(project_info, type = "jxn")
  
  # Extract junction data and coordinates
  jxn_coords <- rowRanges(rse_jxn)
  
  # Create junction info dataframe first
  junction_info <- data.frame(
    junction = paste0(seqnames(jxn_coords), ":",
                     start(jxn_coords), "-",
                     end(jxn_coords)),
    start = start(jxn_coords),
    end = end(jxn_coords),
    strand = strand(jxn_coords)
  )
  
  # Filter for relevant junctions before normalizing
  relevant_junctions <- which(
    (junction_info$end == SRRM3_INFO$exon24$start & junction_info$strand == "+") |
    (junction_info$start == SRRM3_INFO$exon24$end & junction_info$strand == "+") |
    (junction_info$junction == paste0(SRRM3_INFO$gene$chr, ":",
                                    SRRM3_INFO$exon24$start - 1, "-",
                                    SRRM3_INFO$exon24$end + 1))
  )
  
  # Extract only relevant junctions
  junction_counts <- assay(rse_jxn)[relevant_junctions, , drop = FALSE]
  
  # Calculate normalized counts (CPM) for only relevant junctions
  totals <- colSums(assay(rse_jxn))
  normalized_counts <- t(t(junction_counts) / totals) * 1e6
  
  # Convert to regular matrix, then to data frame (much smaller now)
  junction_counts_df <- as.data.frame(as.matrix(normalized_counts))
  colnames(junction_counts_df) <- colnames(junction_counts)
  
  return(list(
    junctions = junction_info[relevant_junctions, ],
    counts = junction_counts_df
  ))
}

# Main analysis function
run_analysis <- function() {
  # Get projects
  projects <- available_projects()
  
  # Get TCGA GBM data
  tcga_project <- subset(projects, project == "GBM" & file_source == "tcga")
  tcga_data <- get_junction_data(tcga_project)
  
  # Get GTEx brain data
  gtex_project <- subset(projects, project == "BRAIN" & file_source == "gtex")
  gtex_data <- get_junction_data(gtex_project)
  
  # Define inclusion and skipping junctions
  inclusion_upstream <- SRRM3_INFO$exon24$start - 1
  inclusion_downstream <- SRRM3_INFO$exon24$end + 1
  skipping_junction <- paste0(SRRM3_INFO$gene$chr, ":", 
                            inclusion_upstream, "-",
                            inclusion_downstream)
  
  # Calculate PSI for TCGA samples
  tcga_psi <- data.frame(
    sample_id = colnames(tcga_data$counts)
  ) %>%
    rowwise() %>%
    mutate(
      inclusion_counts = (
        mean(c(
          sum(tcga_data$counts[tcga_data$junctions$end == SRRM3_INFO$exon24$start & 
                              tcga_data$junctions$strand == "+", sample_id]),
          sum(tcga_data$counts[tcga_data$junctions$start == SRRM3_INFO$exon24$end & 
                              tcga_data$junctions$strand == "+", sample_id])
        ))
      ),
      skipping_counts = sum(tcga_data$counts[tcga_data$junctions$junction == skipping_junction, sample_id]),
      PSI = calculate_PSI(inclusion_counts, skipping_counts)
    )
  
  # Calculate PSI for GTEx samples
  gtex_psi <- data.frame(
    sample_id = colnames(gtex_data$counts)
  ) %>%
    rowwise() %>%
    mutate(
      inclusion_counts = (
        mean(c(
          sum(gtex_data$counts[gtex_data$junctions$end == SRRM3_INFO$exon24$start & 
                              gtex_data$junctions$strand == "+", sample_id]),
          sum(gtex_data$counts[gtex_data$junctions$start == SRRM3_INFO$exon24$end & 
                              gtex_data$junctions$strand == "+", sample_id])
        ))
      ),
      skipping_counts = sum(gtex_data$counts[gtex_data$junctions$junction == skipping_junction, sample_id]),
      PSI = calculate_PSI(inclusion_counts, skipping_counts)
    )
  
  # Add dataset information
  tcga_psi$dataset <- "TCGA_Glioma"
  gtex_psi$dataset <- "GTEx_Brain"
  
  # Combine the datasets
  combined_psi <- bind_rows(tcga_psi, gtex_psi)
  
  # Remove NA values and validate data before statistical test
  combined_psi_clean <- combined_psi %>%
    filter(!is.na(PSI)) %>%
    filter(!is.infinite(PSI))
  
  # Check if we have enough data for comparison
  if (n_distinct(combined_psi_clean$dataset) == 2 &&
      all(table(combined_psi_clean$dataset) > 0)) {
    
    # Statistical comparison
    wilcox_test <- wilcox.test(
      PSI ~ dataset, 
      data = combined_psi_clean,
      alternative = "two.sided"
    )
    
    wilcox_pval <- wilcox_test$p.value
    test_subtitle <- paste("Wilcoxon test p-value:", format(wilcox_pval, digits = 3))
  } else {
    warning("Insufficient data for statistical comparison")
    wilcox_test <- NULL
    test_subtitle <- "Insufficient data for statistical comparison"
  }
  
  # Create visualization
  plot <- ggplot(combined_psi_clean, aes(x = dataset, y = PSI, fill = dataset)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
    theme_minimal() +
    labs(
      title = "SRRM3 Exon 24 PSI Values Comparison",
      subtitle = test_subtitle,
      y = "Percent Spliced In (PSI)",
      x = "Dataset"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  # Save results
  write.csv(combined_psi, "SRRM3_exon24_PSI_results.csv", row.names = FALSE)
  ggsave("SRRM3_exon24_PSI_comparison.pdf", width = 8, height = 6, plot = plot)
  
  # Print summary statistics
  summary_stats <- combined_psi %>%
    group_by(dataset) %>%
    summarize(
      mean_PSI = mean(PSI, na.rm = TRUE),
      median_PSI = median(PSI, na.rm = TRUE),
      sd_PSI = sd(PSI, na.rm = TRUE),
      n_samples = n()
    )
  
  print(summary_stats)
  print(paste("Wilcoxon test p-value:", wilcox_test$p.value))
  
  return(list(
    data = combined_psi,
    plot = plot,
    stats = summary_stats,
    wilcox_test = wilcox_test
  ))
}

# Run analysis
results <- run_analysis()