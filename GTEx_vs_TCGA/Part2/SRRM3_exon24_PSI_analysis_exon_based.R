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
calculate_PSI <- function(exon_coverage, flanking_coverage) {
    psi <- exon_coverage / (exon_coverage + flanking_coverage)
    return(psi)
}

# Function to get exon coverage data from recount3
get_exon_coverage <- function(project_info) {
  message("\nGetting exon coverage data for project: ", project_info$project)
  
  # Get exon-level data
  rse_exon <- create_rse(project_info, type = "exon")
  
  # Extract coverage data
  exon_counts <- assay(rse_exon)
  exon_coords <- rowRanges(rse_exon)
  
  # Calculate normalized counts (RPM)
  totals <- colSums(exon_counts)
  normalized_counts <- t(t(exon_counts) / totals) * 1e6
  
  # Create exon info dataframe with more detailed information
  exon_info <- data.frame(
    exon_id = names(exon_coords),
    chr = as.character(seqnames(exon_coords)),
    start = start(exon_coords),
    end = end(exon_coords),
    width = width(exon_coords),
    strand = as.character(strand(exon_coords))
  )
  
  # Print debugging information
  message("Number of exons found: ", nrow(exon_info))
  message("Chromosome names present: ", paste(unique(exon_info$chr), collapse=", "))
  
  # Check if our target exon is present
  target_exons <- exon_info %>%
    filter(
      chr == SRRM3_INFO$gene$chr,
      start >= SRRM3_INFO$gene$start,
      end <= SRRM3_INFO$gene$end
    )
  
  message("Number of exons in SRRM3 region: ", nrow(target_exons))
  if(nrow(target_exons) > 0) {
    message("Found exons in SRRM3 region:")
    print(head(target_exons))
  }
  
  # Convert the normalized counts matrix to a data frame and set column names
  # This makes the coverage data easier to work with downstream and ensures
  # sample IDs are preserved as column names. The coverage values are in RPM
  # (Reads Per Million) to account for differences in sequencing depth
  coverage_df <- as.data.frame(normalized_counts)
  colnames(coverage_df) <- colnames(exon_counts)
  
  # Create unique row names by adding a suffix to duplicates
  exon_ids <- exon_info$exon_id
  duplicated_ids <- duplicated(exon_ids) | duplicated(exon_ids, fromLast = TRUE)
  if(any(duplicated_ids)) {
    message("Found ", sum(duplicated_ids)/2, " duplicate exon IDs")
    # Add numeric suffix to duplicates
    exon_ids[duplicated_ids] <- make.unique(exon_ids[duplicated_ids], sep = "_dup")
  }
  
  # Assign the unique row names
  rownames(coverage_df) <- exon_ids
  # Update exon_info to maintain consistency
  exon_info$exon_id <- exon_ids
  
  return(list(
    exons = exon_info,
    coverage = coverage_df
  ))
}

# Add this function before run_analysis()
create_comparison_plot <- function(combined_psi) {
  # Calculate mean PSI for each dataset
  mean_psi <- combined_psi %>%
    group_by(dataset) %>%
    summarize(mean_PSI = mean(PSI, na.rm = TRUE)) %>%
    pivot_wider(names_from = dataset, 
                values_from = mean_PSI,
                names_prefix = "mean_PSI_")
  
  # Calculate -log10(p-value) for coloring
  if(n_distinct(combined_psi$dataset) == 2) {
    wilcox_test <- wilcox.test(
      PSI ~ dataset, 
      data = combined_psi,
      alternative = "two.sided"
    )
    neg_log10_pval <- -log10(wilcox_test$p.value)
  } else {
    neg_log10_pval <- NA
  }
  
  # Get the maximum PSI value to determine if we need the full 0-1 range
  max_psi <- max(c(mean_psi$mean_PSI_GTEx_Brain, mean_psi$mean_PSI_TCGA_Glioma))
  
  # Create scatter plot
  scatter_plot <- ggplot() +
    # Add diagonal line
    geom_abline(slope = 1, intercept = 0, color = "black") +
    # Add point for SRRM3
    geom_point(data = mean_psi, 
              aes(x = mean_PSI_GTEx_Brain, 
                  y = mean_PSI_TCGA_Glioma),
              color = "black",
              size = 3) +
    # # Add gene label
    # geom_text(data = mean_psi,
    #           aes(x = mean_PSI_GTEx_Brain, 
    #               y = mean_PSI_TCGA_Glioma,
    #               label = "SRRM3"),
    #           vjust = -1,
    #           size = 4) +
    # Customize appearance with adjusted axis ranges
    scale_x_continuous(
      limits = c(0, max_psi*1.2),  # Adjust based on your data
      breaks = seq(0, max_psi*1.2, 0.005),
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      limits = c(0, max_psi*1.2),  # Adjust based on your data
      breaks = seq(0, max_psi*1.2, 0.005),
      expand = c(0.01, 0)
    ) +
    # Add more gridlines for better readability
    theme_minimal() +
    theme(
      panel.grid.minor = element_line(color = "grey90"),
      panel.grid.major = element_line(color = "grey85")
    ) +
    labs(
      x = "Mean PSI (GTEx)",
      y = "Mean PSI (TCGA)",
      title = "SRRM3 Exon 24 PSI Comparison",
      subtitle = ifelse(!is.na(neg_log10_pval),
                       paste("-log10(p-value) =", round(neg_log10_pval, 2)),
                       "Statistical test not performed")
    ) +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  return(scatter_plot)
}

# Function to calculate PSI for a dataset
calculate_dataset_PSI <- function(data, sample_ids) {
  # Get all exons that overlap with the SRRM3 gene region
  srrm3_exons <- data$exons %>%
    filter(
      chr == SRRM3_INFO$gene$chr,
      start >= SRRM3_INFO$gene$start,
      end <= SRRM3_INFO$gene$end
    ) %>%
    # Sort exons by their start position to see them in genomic order
    arrange(start) %>%
    # Remove any duplicate exon entries that might exist in the annotation
    # This ensures we don't double-count any exons in our analysis
    distinct(start, end, .keep_all = TRUE)
  
  
  message("All SRRM3 exons (duplicates removed):")
  print(srrm3_exons)
  
  # Find target exon (exon 24)
  target_exon_idx <- which(
    data$exons$chr == SRRM3_INFO$gene$chr &
      data$exons$start == SRRM3_INFO$exon24$start &
      data$exons$end == SRRM3_INFO$exon24$end &
      data$exons$strand == "+"
  )[1]  # Take first match if multiple exist
  
  if(length(target_exon_idx) == 0) {
    message("Target exon not found!")
    return(data.frame(
      sample_id = sample_ids,
      exon_coverage = NA,
      flanking_coverage = NA,
      PSI = NA
    ))
  }
  
  message("Found target exon: ", 
          data$exons$chr[target_exon_idx], ":",
          data$exons$start[target_exon_idx], "-",
          data$exons$end[target_exon_idx])
  
  # Find flanking exons (exon 23 and exon 25)
  upstream_exon_idx <- which(
    data$exons$chr == SRRM3_INFO$gene$chr &
      data$exons$start == 76282964 &  # exon 23 coordinates
      data$exons$end == 76283101 &
      data$exons$strand == "+"
  )[1]  # Take first match if multiple exist
  
  downstream_exon_idx <- which(
    data$exons$chr == SRRM3_INFO$gene$chr &
      data$exons$start == 76285615 &  # exon 25 coordinates
      data$exons$end == 76287286 &
      data$exons$strand == "+"
  )[1]  # Take first match if multiple exist
  
  message("Found upstream exon (exon 23): ", 
          if(length(upstream_exon_idx) > 0) 
            paste0(data$exons$start[upstream_exon_idx], "-", 
                   data$exons$end[upstream_exon_idx]) 
          else "NONE")
  
  message("Found downstream exon (exon 25): ", 
          if(length(downstream_exon_idx) > 0) 
            paste0(data$exons$start[downstream_exon_idx], "-", 
                   data$exons$end[downstream_exon_idx]) 
          else "NONE")
  
  # Calculate PSI for each sample
  psi_df <- data.frame(
    sample_id = sample_ids
  ) %>%
    rowwise() %>%
    mutate(
      exon_coverage = data$coverage[target_exon_idx, sample_id],
      upstream_coverage = if(length(upstream_exon_idx) == 1) 
        data$coverage[upstream_exon_idx, sample_id] 
      else NA,
      downstream_coverage = if(length(downstream_exon_idx) == 1) 
        data$coverage[downstream_exon_idx, sample_id] 
      else NA,
      flanking_coverage = mean(c(upstream_coverage, downstream_coverage), 
                               na.rm = TRUE),
      PSI = if(!is.na(flanking_coverage) && flanking_coverage > 0) 
        calculate_PSI(exon_coverage, flanking_coverage) 
      else NA
    ) %>%
    select(sample_id, exon_coverage, flanking_coverage, PSI)
  
  return(psi_df)
}


################################################################################
############################ Main analysis function ############################
################################################################################
# run_analysis <- function() {
# Get projects
projects <- available_projects()

# Get TCGA GBM data
tcga_project <- subset(projects, project == "GBM" & file_source == "tcga")
tcga_data <- get_exon_coverage(tcga_project)
saveRDS(tcga_data, "tcga_data.rds")
tcga_data <- readRDS("tcga_data.rds")

# Get GTEx brain data
gtex_project <- subset(projects, project == "BRAIN" & file_source == "gtex")
gtex_data <- get_exon_coverage(gtex_project)
saveRDS(gtex_data, "gtex_data.rds")
gtex_data <- readRDS("gtex_data.rds")

# Calculate PSI for both datasets
tcga_psi <- calculate_dataset_PSI(tcga_data, colnames(tcga_data$coverage))
gtex_psi <- calculate_dataset_PSI(gtex_data, colnames(gtex_data$coverage))

# Add dataset information
tcga_psi$dataset <- "TCGA_Glioma"
gtex_psi$dataset <- "GTEx_Brain"

# Combine the datasets
combined_psi <- bind_rows(tcga_psi, gtex_psi)

# Remove rows with NA PSI values
combined_psi <- combined_psi %>% filter(!is.na(PSI))

# Check if we have data for both datasets
datasets_present <- unique(combined_psi$dataset)

if (length(datasets_present) == 2) {
  # Perform statistical comparison
  wilcox_test <- wilcox.test(
    PSI ~ dataset, 
    data = combined_psi,
    alternative = "two.sided"
  )
  wilcox_p_value <- wilcox_test$p.value
} else {
  wilcox_test <- NULL
  wilcox_p_value <- NA
  warning("Cannot perform Wilcoxon test: data available for only one dataset or no data available.")
}

# Create visualization
plot <- ggplot(combined_psi, aes(x = dataset, y = PSI, fill = dataset)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "SRRM3 Exon 24 PSI Values Comparison (Exon-based)",
    subtitle = ifelse(!is.na(wilcox_p_value),
                      paste("Wilcoxon test p-value:", format(wilcox_p_value, digits = 3)),
                      "Wilcoxon test not performed: insufficient data"),
    y = "Percent Spliced In (PSI)",
    x = "Dataset"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Save results
write.csv(combined_psi, "SRRM3_exon24_PSI_results_exon_based.csv", row.names = FALSE)
ggsave("SRRM3_exon24_PSI_comparison_exon_based.pdf", width = 8, height = 6, plot = plot)

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
if (!is.null(wilcox_test)) {
  print(paste("Wilcoxon test p-value:", wilcox_p_value))
} else {
  print("Wilcoxon test not performed: insufficient data")
}

# Create scatter plot
scatter_plot <- create_comparison_plot(combined_psi)

# Save scatter plot
ggsave("SRRM3_exon24_PSI_scatter_comparison.pdf", 
        plot = scatter_plot, 
        width = 8, 
        height = 8)

# Add scatter_plot to the return list
# return(list(
#   data = combined_psi,
#   plot = plot,
#   scatter_plot = scatter_plot,
#   stats = summary_stats,
#   wilcox_test = wilcox_test
# ))
# }