library(ggplot2)
library(dplyr)
library(tidyr)
library(SummarizedExperiment)

# Function to process expression data for each dataset
get_isoform_expression <- function(structure, source) {
  # Get expression data for each exon
  expr_data <- structure$expression_data
  ranges <- structure$raw_ranges
  
  # Get transcript IDs for each exon
  transcript_ids <- mcols(ranges)$transcript_id
  
  # Calculate mean expression per transcript across all samples
  expr_by_transcript <- data.frame(
    transcript_id = transcript_ids,
    expression = expr_data
  ) %>%
    group_by(transcript_id) %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    gather(sample, expression, -transcript_id)
  
  # Add source information
  expr_by_transcript$source <- source
  
  return(expr_by_transcript)
}

# Process TCGA data
tcga_expr <- get_isoform_expression(tcga_structure, "TCGA GBM")

# Process GTEx data
gtex_expr <- get_isoform_expression(gtex_structure, "GTEx Brain")

# Combine the data
all_expr <- rbind(tcga_expr, gtex_expr)

# Calculate summary statistics
summary_stats <- all_expr %>%
  group_by(transcript_id, source) %>%
  summarise(
    mean_expr = mean(expression, na.rm = TRUE),
    median_expr = median(expression, na.rm = TRUE),
    sd_expr = sd(expression, na.rm = TRUE),
    n = n()
  )

# Create violin plot with box plot overlay
p1 <- ggplot(all_expr, aes(x = transcript_id, y = expression, fill = source)) +
  geom_violin(position = position_dodge(width = 0.7), alpha = 0.6) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.7), alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "SRRM3 Isoform Expression Comparison",
    subtitle = "TCGA GBM vs GTEx Brain",
    x = "Transcript ID",
    y = "Expression Level",
    fill = "Source"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "top"
  ) +
  scale_fill_manual(values = c("GTEx Brain" = "#2ecc71", "TCGA GBM" = "#e74c3c"))

# Create bar plot of mean expression with error bars
p2 <- ggplot(summary_stats, aes(x = transcript_id, y = mean_expr, fill = source)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr),
    position = position_dodge(width = 0.8),
    width = 0.25
  ) +
  theme_minimal() +
  labs(
    title = "Mean SRRM3 Isoform Expression",
    subtitle = "With Standard Deviation",
    x = "Transcript ID",
    y = "Mean Expression Level",
    fill = "Source"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "top"
  ) +
  scale_fill_manual(values = c("GTEx Brain" = "#2ecc71", "TCGA GBM" = "#e74c3c"))

# Perform statistical tests
stat_tests <- all_expr %>%
  group_by(transcript_id) %>%
  summarise(
    wilcox_p = wilcox.test(
      expression[source == "TCGA GBM"],
      expression[source == "GTEx Brain"]
    )$p.value,
    log2fc = log2(mean(expression[source == "TCGA GBM"]) / 
                    mean(expression[source == "GTEx Brain"]))
  )

# Print statistical results
message("\nStatistical Analysis:")
print(stat_tests)

# Save plots
pdf("3_srrm3_isoform_expression.pdf", width = 12, height = 10)
print(p1)
print(p2)
dev.off()

# Create a summary table
summary_table <- summary_stats %>%
  spread(source, mean_expr) %>%
  left_join(stat_tests, by = "transcript_id") %>%
  arrange(wilcox_p)

# Save results
save(summary_stats, stat_tests, summary_table, 
     file = "3_srrm3_expression_analysis.RData")

# Print summary table
message("\nExpression Summary Table:")
print(summary_table)