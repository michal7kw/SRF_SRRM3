library(recount3)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(SummarizedExperiment)

# Function to calculate PSI values
calculate_psi <- function(project_info, gene_coords) {
  message("\nCalculating PSI for project: ", project_info$project)
  
  # Get junction-level data
  rse_jxn <- create_rse(project_info, type = "jxn")
  
  # Define the SRRM3 region
  region <- GRanges(
    seqnames = "chr7",
    ranges = IRanges(
      start = gene_coords$start,
      end = gene_coords$end
    )
  )
  
  # Get overlapping junctions
  overlaps <- findOverlaps(rowRanges(rse_jxn), region)
  
  if(length(overlaps) == 0) {
    message("No junctions found in the specified region")
    return(NULL)
  }
  
  # Extract junction data
  jxn_data <- rse_jxn[queryHits(overlaps),]
  junction_counts <- assay(jxn_data)
  jxn_coords <- rowRanges(jxn_data)
  
  # Create junction identifiers
  jxn_ids <- paste0(
    seqnames(jxn_coords), ":",
    start(jxn_coords), "-",
    end(jxn_coords)
  )
  
  # Calculate PSI for each junction
  psi_matrix <- matrix(NA, nrow = length(jxn_ids), ncol = ncol(junction_counts))
  rownames(psi_matrix) <- jxn_ids
  colnames(psi_matrix) <- colnames(junction_counts)
  
  # For each junction
  for(i in seq_along(jxn_ids)) {
    # Find alternative junctions sharing either 5' or 3' splice site
    alt_5prime <- which(start(jxn_coords) == start(jxn_coords[i]) & 
                          seq_along(jxn_coords) != i)
    alt_3prime <- which(end(jxn_coords) == end(jxn_coords[i]) & 
                          seq_along(jxn_coords) != i)
    
    # Calculate PSI for each sample
    for(j in seq_len(ncol(junction_counts))) {
      # Junction reads
      inclusion_reads <- junction_counts[i, j]
      
      # Alternative junction reads
      alt_5prime_reads <- if(length(alt_5prime) > 0) 
        sum(junction_counts[alt_5prime, j]) else 0
      alt_3prime_reads <- if(length(alt_3prime) > 0) 
        sum(junction_counts[alt_3prime, j]) else 0
      
      # Calculate PSI
      total_reads <- inclusion_reads + alt_5prime_reads + alt_3prime_reads
      if(total_reads > 0) {
        psi_matrix[i, j] <- inclusion_reads / total_reads * 100
      }
    }
  }
  
  return(list(
    psi_values = psi_matrix,
    junction_info = data.frame(
      junction_id = jxn_ids,
      start = start(jxn_coords),
      end = end(jxn_coords),
      stringsAsFactors = FALSE
    )
  ))
}

# SRRM3 coordinates
srrm3_coords <- list(
  start = 76201900,
  end = 76287288
)

# Get projects
projects <- available_projects()

# Calculate PSI values for both datasets
message("Calculating TCGA PSI values...")
tcga_psi <- calculate_psi(
  subset(projects, project == "GBM" & file_source == "tcga"),
  srrm3_coords
)

message("Calculating GTEx PSI values...")
gtex_psi <- calculate_psi(
  subset(projects, project == "BRAIN" & file_source == "gtex"),
  srrm3_coords
)

# Create comparison data frame
create_psi_comparison <- function(tcga_psi, gtex_psi) {
  # Find common junctions
  common_junctions <- intersect(
    rownames(tcga_psi$psi_values),
    rownames(gtex_psi$psi_values)
  )
  
  # Calculate mean PSI values
  comparison_df <- data.frame(
    junction_id = common_junctions,
    tcga_psi = rowMeans(tcga_psi$psi_values[common_junctions,], na.rm = TRUE),
    gtex_psi = rowMeans(gtex_psi$psi_values[common_junctions,], na.rm = TRUE)
  )
  
  # Calculate statistics
  comparison_df$p_value <- sapply(common_junctions, function(j) {
    wilcox.test(
      tcga_psi$psi_values[j,],
      gtex_psi$psi_values[j,],
      na.rm = TRUE
    )$p.value
  })
  
  # Add additional metrics
  comparison_df <- comparison_df %>%
    mutate(
      delta_psi = tcga_psi - gtex_psi,
      significance = -log10(p_value),
      is_significant = p_value < 0.05 & abs(delta_psi) > 10
    )
  
  return(comparison_df)
}

# Create comparison
comparison_data <- create_psi_comparison(tcga_psi, gtex_psi)

# Create scatter plot
p <- ggplot(comparison_data, aes(x = gtex_psi, y = tcga_psi)) +
  # Base points
  geom_point(aes(color = significance), alpha = 0.7, size = 3) +
  
  # Reference line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  
  # Add labels for significant changes
  geom_label_repel(
    data = subset(comparison_data, is_significant),
    aes(label = junction_id),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.5,
    max.overlaps = 20
  ) +
  
  # Styling
  scale_color_viridis_c(name = "-log10(p-value)") +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  labs(
    title = "SRRM3 PSI Comparison",
    subtitle = "TCGA Glioma vs GTEx Brain",
    x = "Mean PSI (GTEx)",
    y = "Mean PSI (TCGA)",
    caption = paste0(
      "Significant changes (p < 0.05, |ΔPSI| > 10): ",
      sum(comparison_data$is_significant)
    )
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

# Save plot
pdf("srrm3_psi_comparison.pdf", width = 12, height = 10)
print(p)
dev.off()

# Create summary table
psi_summary <- comparison_data %>%
  arrange(p_value) %>%
  mutate(
    significant = p_value < 0.05 & abs(delta_psi) > 10
  )

# Save results
write.csv(psi_summary, "srrm3_psi_summary.csv", row.names = FALSE)

# Print significant changes
message("\nSignificant PSI changes (p < 0.05, |ΔPSI| > 10):")
print(subset(psi_summary, significant))