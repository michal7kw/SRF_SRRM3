#####################################################################
# Combine High SRRM3 PSI Analysis Results
#####################################################################
# This script combines and summarizes results from multiple cancer types
# analyzed in the High SRRM3 PSI Analysis pipeline

library(tidyverse)

results_dir <- "results/"

# Combine summaries
summaries <- list.files(results_dir, 
                       pattern = "*_summary.csv",
                       full.names = TRUE) %>%
    lapply(read.csv) %>%
    bind_rows()

# Add additional statistics
summaries <- summaries %>%
    mutate(
        psi_analysis_rate = analyzed_samples / high_expr_samples,
        survival_ratio = median_survival_high / median_survival_low
    )

# Save combined results
write.csv(summaries,
          file.path(results_dir, "combined_summary.csv"),
          row.names = FALSE)

# Create summary visualization
p <- ggplot(summaries, aes(x = cancer_type)) +
    geom_bar(aes(y = high_expr_samples), stat = "identity") +
    geom_bar(aes(y = analyzed_samples), stat = "identity", alpha = 0.5) +
    theme_minimal() +
    labs(title = "Sample counts by cancer type",
         y = "Number of samples",
         x = "Cancer type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(results_dir, "sample_counts.pdf"), p, width = 10, height = 6) 