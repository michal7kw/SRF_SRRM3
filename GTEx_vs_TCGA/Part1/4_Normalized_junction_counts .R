library(recount3)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)  # For better label placement
library(SummarizedExperiment)

# Function to calculate junction usage with better annotation
calculate_junction_usage <- function(project_info, gene_coords) {
  message("\nCalculating junction usage for project: ", project_info$project)
  
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
  
  # Get junction coordinates
  jxn_coords <- rowRanges(jxn_data)
  
  # Create detailed junction annotations
  junction_info <- data.frame(
    junction_id = paste0(
      seqnames(jxn_coords), ":",
      start(jxn_coords), "-",
      end(jxn_coords)
    ),
    start_exon = paste0("Exon_", seq_along(start(jxn_coords))),
    end_exon = paste0("Exon_", seq_along(end(jxn_coords))),
    start_pos = start(jxn_coords),
    end_pos = end(jxn_coords),
    width = width(jxn_coords)
  )
  
  # Calculate normalized junction usage
  junction_totals <- colSums(junction_counts)
  normalized_counts <- t(t(junction_counts) / junction_totals) * 1e6  # Per million scaling
  
  return(list(
    junction_info = junction_info,
    counts = junction_counts,
    normalized_counts = normalized_counts,
    coords = jxn_coords
  ))
}

# SRRM3 coordinates
srrm3_coords <- list(
  start = 76201900,
  end = 76287288
)

# Get projects
projects <- available_projects()

# Get data
tcga_data <- calculate_junction_usage(
  subset(projects, project == "GBM" & file_source == "tcga"),
  srrm3_coords
)

gtex_data <- calculate_junction_usage(
  subset(projects, project == "BRAIN" & file_source == "gtex"),
  srrm3_coords
)

# Function to create detailed comparison data
create_detailed_comparison <- function(tcga_data, gtex_data) {
  message("\nCreating detailed comparison...")
  
  # Find common junctions
  common_junctions <- intersect(
    tcga_data$junction_info$junction_id,
    gtex_data$junction_info$junction_id
  )
  
  message("Found ", length(common_junctions), " common junctions")
  
  # Get indices for common junctions
  tcga_idx <- match(common_junctions, tcga_data$junction_info$junction_id)
  gtex_idx <- match(common_junctions, gtex_data$junction_info$junction_id)
  
  # Calculate statistics
  comparison_df <- data.frame(
    junction_id = common_junctions,
    tcga_usage = rowMeans(tcga_data$normalized_counts[tcga_idx,], na.rm = TRUE),
    gtex_usage = rowMeans(gtex_data$normalized_counts[gtex_idx,], na.rm = TRUE)
  ) %>%
    mutate(
      # Calculate additional metrics
      p_value = sapply(seq_along(common_junctions), function(i) {
        tcga_vals <- tcga_data$normalized_counts[tcga_idx[i],]
        gtex_vals <- gtex_data$normalized_counts[gtex_idx[i],]
        wilcox.test(tcga_vals, gtex_vals)$p.value
      }),
      log2FC = log2((tcga_usage + 1)/(gtex_usage + 1)),
      significance = -log10(p_value),
      
      # Add junction details
      start_pos = tcga_data$junction_info$start_pos[tcga_idx],
      end_pos = tcga_data$junction_info$end_pos[tcga_idx],
      
      # Create junction labels
      junction_label = paste0(
        "J", seq_along(common_junctions), "\n",
        "(", start_pos, "-", end_pos, ")"
      ),
      
      # Flag highly significant differences
      is_significant = p_value < 0.05 & abs(log2FC) > 1
    )
  
  return(comparison_df)
}

# Create detailed comparison
comparison_data <- create_detailed_comparison(tcga_data, gtex_data)

# Create enhanced scatter plot
p <- ggplot(comparison_data, aes(x = gtex_usage, y = tcga_usage)) +
  # Base points
  geom_point(aes(color = significance), alpha = 0.7, size = 3) +
  
  # Reference line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  
  # Add labels for significant junctions
  geom_label_repel(
    data = subset(comparison_data, is_significant),
    aes(label = junction_label),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.5,
    force = 5,
    max.overlaps = 20
  ) +
  
  # Styling
  scale_color_viridis_c(name = "-log10(p-value)") +
  theme_minimal() +
  labs(
    title = "SRRM3 Junction Usage Comparison",
    subtitle = paste0(
      "TCGA Glioma (n=", ncol(tcga_data$counts), 
      ") vs GTEx Brain (n=", ncol(gtex_data$counts), ")"
    ),
    x = "Mean Junction Usage (GTEx)",
    y = "Mean Junction Usage (TCGA)",
    caption = paste0(
      "Total junctions: ", nrow(comparison_data), "\n",
      "Significant differential usage: ", 
      sum(comparison_data$is_significant)
    )
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

# Save plot
pdf("4_srrm3_labeled_junction_comparison.pdf", width = 12, height = 10)
print(p)
dev.off()

# Create detailed summary table
junction_summary <- comparison_data %>%
  arrange(p_value) %>%
  mutate(
    fold_change = 2^log2FC,
    significant = p_value < 0.05
  ) %>%
  select(
    junction_id,
    start_pos,
    end_pos,
    tcga_usage,
    gtex_usage,
    fold_change,
    p_value,
    significant
  )

# Save detailed results
write.csv(junction_summary, "4_srrm3_junction_summary.csv", row.names = FALSE)

# Print summary of significant changes
message("\nTop differential junctions:")
print(head(junction_summary, n = 10))