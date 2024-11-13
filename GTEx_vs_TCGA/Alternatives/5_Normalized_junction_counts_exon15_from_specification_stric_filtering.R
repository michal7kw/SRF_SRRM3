# Load required libraries
library(recount3)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(SummarizedExperiment)

# Define SRRM3 information with exact coordinates
SRRM3_INFO <- list(
  gene = list(
    name = "SRRM3",
    chr = "chr7",
    start = 76201896,  # GRCh38 coordinates
    end = 76287287
  ),
  transcripts = list(
    with_exon15 = "NM_001291831.2",    # 16 exons
    without_exon15 = "NM_001110199.3"   # 15 exons
  ),
  exon15 = list(
    transcript_start = 1945,   # Position in transcript
    transcript_end = 2023,     # Position in transcript
    length = 79               # Exact length from transcript comparison
  )
)

# Function to get junction data
get_junction_data <- function(project_info) {
  message("\nGetting junction data for project: ", project_info$project)
  
  # Get junction-level data
  rse_jxn <- create_rse(project_info, type = "jxn")
  
  # Define SRRM3 region
  region <- GRanges(
    seqnames = SRRM3_INFO$gene$chr,
    ranges = IRanges(
      start = SRRM3_INFO$gene$start,
      end = SRRM3_INFO$gene$end
    )
  )
  
  # Get overlapping junctions
  overlaps <- findOverlaps(rowRanges(rse_jxn), region)
  
  if(length(overlaps) == 0) {
    stop("No junctions found in SRRM3 region")
  }
  
  # Extract junction data
  jxn_data <- rse_jxn[queryHits(overlaps),]
  junction_counts <- assay(jxn_data)
  jxn_coords <- rowRanges(jxn_data)
  
  # Calculate normalized counts
  totals <- colSums(junction_counts)
  normalized_counts <- t(t(junction_counts) / totals) * 1e6
  
  # Create junction info
  junction_info <- data.frame(
    junction_id = paste0(
      seqnames(jxn_coords), ":",
      start(jxn_coords), "-",
      end(jxn_coords)
    ),
    start_pos = start(jxn_coords),
    end_pos = end(jxn_coords),
    width = width(jxn_coords)
  )
  
  return(list(
    junction_info = junction_info,
    normalized_counts = normalized_counts
  ))
}

# Function to analyze junctions
analyze_junctions <- function(tcga_data, gtex_data) {
  # Find common junctions
  common_junctions <- intersect(
    tcga_data$junction_info$junction_id,
    gtex_data$junction_info$junction_id
  )
  
  tcga_idx <- match(common_junctions, tcga_data$junction_info$junction_id)
  gtex_idx <- match(common_junctions, gtex_data$junction_info$junction_id)
  
  comparison <- data.frame(
    junction_id = common_junctions,
    start_pos = tcga_data$junction_info$start_pos[tcga_idx],
    end_pos = tcga_data$junction_info$end_pos[tcga_idx],
    width = tcga_data$junction_info$width[tcga_idx],
    tcga_usage = rowMeans(tcga_data$normalized_counts[tcga_idx,], na.rm = TRUE),
    gtex_usage = rowMeans(gtex_data$normalized_counts[gtex_idx,], na.rm = TRUE)
  ) %>%
    mutate(
      p_value = sapply(seq_along(common_junctions), function(i) {
        wilcox.test(
          tcga_data$normalized_counts[tcga_idx[i],],
          gtex_data$normalized_counts[gtex_idx[i],]
        )$p.value
      }),
      log2FC = log2((tcga_usage + 1)/(gtex_usage + 1)),
      
      # More precise junction classification
      junction_type = case_when(
        # Exact 79bp junctions within the expected region
        abs(width - 79) <= 2 & 
          start_pos >= 76283000 & 
          end_pos <= 76284000 ~ "Exon 15 (79bp)",
        
        # Skipping junctions that span the exon 15 region
        width > 200 & 
          start_pos < 76283000 & 
          end_pos > 76284000 ~ "Skips Exon 15",
        
        # Everything else
        TRUE ~ "Other Junctions"
      ),
      
      junction_label = paste0(
        "J", seq_along(common_junctions), "\n",
        format(start_pos, big.mark = ","), "-",
        format(end_pos, big.mark = ","), "\n",
        "Width: ", width, "bp"
      ),
      
      is_significant = p_value < 0.05 & abs(log2FC) > 1
    )
  
  # Print summary of junctions by type
  cat("\nJunction Summary:\n")
  print(comparison %>%
          group_by(junction_type) %>%
          summarize(
            count = n(),
            avg_width = mean(width),
            significant = sum(is_significant),
            mean_tcga = mean(tcga_usage),
            mean_gtex = mean(gtex_usage)
          ))
  
  return(comparison)
}

# Update visualization to be more informative
create_visualization <- function(comparison_data) {
  ggplot(comparison_data, 
         aes(x = gtex_usage, y = tcga_usage, color = junction_type)) +
    geom_point(aes(size = -log10(p_value)), alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    geom_label_repel(
      data = subset(comparison_data, 
                    junction_type != "Other Junctions" & is_significant),
      aes(label = junction_label),
      size = 3,
      box.padding = 0.5,
      force = 5
    ) +
    scale_color_manual(
      values = c(
        "Exon 15 (79bp)" = "blue",
        "Skips Exon 15" = "red",
        "Other Junctions" = "gray70"
      )
    ) +
    scale_size_continuous(range = c(2, 6)) +
    theme_minimal() +
    labs(
      title = "SRRM3 Junction Usage Comparison",
      subtitle = paste0(
        "Points represent individual splice junctions\n",
        "Blue: Potential exon 15 junctions (≈79bp)\n",
        "Red: Junctions that skip exon 15\n",
        "Grey: Other splice junctions in SRRM3"
      ),
      x = "Mean Junction Usage (GTEx)",
      y = "Mean Junction Usage (TCGA)",
      color = "Junction Type",
      size = "-log10(p-value)"
    )
}

# Main analysis function
run_analysis <- function() {
  # Get projects
  projects <- available_projects()
  
  # Get data
  tcga_data <- get_junction_data(
    subset(projects, project == "GBM" & file_source == "tcga")
  )
  
  gtex_data <- get_junction_data(
    subset(projects, project == "BRAIN" & file_source == "gtex")
  )
  
  # Analyze junctions
  comparison_data <- analyze_junctions(tcga_data, gtex_data)
  
  # Create and save visualization
  plot <- create_visualization(comparison_data)
  pdf("srrm3_isoform_analysis_strict_filtering.pdf", width = 12, height = 10)
  print(plot)
  dev.off()
  
  # Save results
  write.csv(comparison_data, "srrm3_junction_analysis_strict_filtering.csv", row.names = FALSE)
  
  # Print summary
  cat("\nJunction Analysis Summary:\n")
  print(comparison_data %>%
          group_by(junction_type) %>%
          summarize(
            count = n(),
            significant = sum(is_significant),
            mean_tcga = mean(tcga_usage),
            mean_gtex = mean(gtex_usage)
          ))
  
  return(list(data = comparison_data, plot = plot))
}

# Run analysis
results <- run_analysis()