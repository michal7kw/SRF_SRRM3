# Load required libraries
library(recount3)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(SummarizedExperiment)
library(progress)
library(futile.logger)
library(R.utils)
library(viridis)

# Set up logging
flog.appender(appender.file("srrm3_tcga_analysis.log"))
flog.threshold(DEBUG)

# Define SRRM3 information (keep the same as in GTEx script)
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
  
  # Debug coordinates
  flog.debug("Looking for junctions around exon 15: %d-%d", exon_start, exon_end)
  
  # Find upstream junction (ending at exon start)
  upstream_jxns <- which(abs(end(jxn_coords) - exon_start) <= 5)  # Allow 5bp flexibility
  
  # Find downstream junction (starting at exon end)
  downstream_jxns <- which(abs(start(jxn_coords) - exon_end) <= 5)  # Allow 5bp flexibility
  
  # Find exclusion junctions (those that skip exon 15)
  exclusion_jxns <- which(
    start(jxn_coords) < (exon_start - 5) & 
    end(jxn_coords) > (exon_end + 5)
  )
  
  # Debug junction coordinates
  if(length(upstream_jxns) > 0) {
    flog.debug("Upstream junction coordinates:")
    for(i in upstream_jxns) {
      flog.debug("Junction %d: %d-%d", i, start(jxn_coords[i]), end(jxn_coords[i]))
    }
  }
  
  if(length(downstream_jxns) > 0) {
    flog.debug("Downstream junction coordinates:")
    for(i in downstream_jxns) {
      flog.debug("Junction %d: %d-%d", i, start(jxn_coords[i]), end(jxn_coords[i]))
    }
  }
  
  if(length(exclusion_jxns) > 0) {
    flog.debug("Exclusion junction coordinates:")
    for(i in exclusion_jxns) {
      flog.debug("Junction %d: %d-%d", i, start(jxn_coords[i]), end(jxn_coords[i]))
    }
  }
  
  # Combine inclusion junctions
  inclusion_jxns <- unique(c(upstream_jxns, downstream_jxns))
  
  flog.debug("Upstream junctions found: %d", length(upstream_jxns))
  flog.debug("Downstream junctions found: %d", length(downstream_jxns))
  flog.debug("Total inclusion junctions: %d", length(inclusion_jxns))
  flog.debug("Exclusion junctions found: %d", length(exclusion_jxns))
  
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
    flog.info("Creating RSE for project: %s", project_info$project)
    
    # Create region around exon 15
    region <- GRanges(
      seqnames = SRRM3_INFO$gene$chr,
      ranges = IRanges(
        start = SRRM3_INFO$exon15$start - 5000,  # Include surrounding region
        end = SRRM3_INFO$exon15$end + 5000
      )
    )
    
    # Use UNIQUE junction format for GTEx and TCGA
    jxn_format <- if(project_info$file_source %in% c("gtex", "tcga")) "UNIQUE" else "ALL"
    
    rse <- create_rse(
      project_info,
      type = "jxn",
      jxn_format = jxn_format,
      verbose = TRUE
    )
    
    # Filter to relevant region
    relevant_rows <- which(
      as.character(seqnames(rowRanges(rse))) == as.character(seqnames(region)) &
        start(rowRanges(rse)) >= (start(region)) &
        end(rowRanges(rse)) <= (end(region))
    )
    
    if(length(relevant_rows) == 0) {
      flog.warn("No features found in the specified region")
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
    flog.warn("Missing junctions - Inclusion: %d, Exclusion: %d", 
              length(junctions$inclusion), length(junctions$exclusion))
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

  # psi_values <- sapply(seq_len(ncol(junction_counts)), function(i) {
  #   # Separate upstream and downstream inclusion reads
  #   upstream_reads <- sum(junction_counts[junctions$details$upstream, i])
  #   downstream_reads <- sum(junction_counts[junctions$details$downstream, i])
  #   exclusion_reads <- sum(junction_counts[junctions$exclusion, i])
    
  #   if(exclusion_reads >= 10) {  # Changed threshold condition
  #       psi <- 100 * ((upstream_reads + downstream_reads)/2) / exclusion_reads
  #       return(psi)
  #   } else {
  #       return(NA)
  #   }
  # })
  
  # Log summary statistics
  flog.info("PSI calculation complete:")
  flog.info("Total samples: %d", length(psi_values))
  flog.info("Non-NA samples: %d", sum(!is.na(psi_values)))
  flog.info("Mean PSI: %.2f", mean(psi_values, na.rm = TRUE))
  
  return(list(
    psi = psi_values,
    junctions = list(
      inclusion = rownames(junction_counts)[junctions$inclusion],
      exclusion = rownames(junction_counts)[junctions$exclusion]
    )
  ))
}

# New function to process TCGA metadata
process_tcga_metadata <- function(rse) {
  # Extract metadata
  metadata <- colData(rse) %>% 
    as.data.frame()
  
  # Extract cancer type from TCGA metadata
  metadata$cancer_type <- metadata$tcga.cgc_case_histological_diagnosis
  
  # Clean up cancer type names if needed
  metadata$cancer_type <- gsub("_", " ", metadata$cancer_type)
  
  # Log unique cancer types found
  unique_types <- unique(metadata$cancer_type)
  flog.info("Found %d unique cancer types:", length(unique_types))
  for(type in unique_types) {
    flog.info("  - %s", type)
  }
  
  return(metadata)
}

# Get available TCGA projects
projects <- available_projects()
tcga_projects <- subset(projects, file_source == "tcga")

# Initialize lists to store results
all_psi_data <- list()
all_metadata <- list()

# Process each TCGA project
for(i in 1:nrow(tcga_projects)) {
  project_name <- tcga_projects$project[i]
  
  flog.info("Processing TCGA project: %s", project_name)
  
  # Create RSE object
  rse <- create_rse_safe(tcga_projects[i,])
  
  if(is.null(rse)) {
    flog.error("Failed to create RSE for project %s", project_name)
    next
  }
  
  # Calculate PSI values
  psi_result <- calculate_exon15_psi(rse)
  
  if(is.null(psi_result)) {
    flog.error("Failed to calculate PSI values for project %s", project_name)
    next
  }
  
  # Process metadata
  metadata <- process_tcga_metadata(rse)
  
  # Store results
  all_psi_data[[project_name]] <- psi_result$psi
  all_metadata[[project_name]] <- metadata
}

# Combine all data
plot_data <- data.frame(
  sample_id = unlist(lapply(names(all_psi_data), function(proj) {
    paste0(proj, "_", seq_along(all_psi_data[[proj]]))
  })),
  PSI = unlist(all_psi_data),
  cancer_type = unlist(lapply(names(all_metadata), function(proj) {
    all_metadata[[proj]]$cancer_type
  }))
) %>%
  filter(!is.na(PSI)) %>%  # Remove NA values
  filter(cancer_type != "") # Remove empty cancer types

# Calculate summary statistics
summary_stats <- plot_data %>%
  group_by(cancer_type) %>%
  summarize(
    mean_PSI = mean(PSI),
    median_PSI = median(PSI),
    sd_PSI = sd(PSI),
    n_samples = n(),
    .groups = 'drop'
  )

# Create violin plot
p1 <- ggplot(plot_data, aes(x = reorder(cancer_type, PSI, FUN = median), 
                           y = PSI, 
                           fill = cancer_type)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "SRRM3 Exon 15 PSI Values Across Cancer Types",
    y = "Percent Spliced In (PSI)",
    x = "Cancer Type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(size = 12)
  )

# Create scatter plot with means and error bars
p2 <- ggplot(summary_stats, 
             aes(x = reorder(cancer_type, mean_PSI), 
                 y = mean_PSI, 
                 color = cancer_type)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_PSI - sd_PSI, 
                    ymax = mean_PSI + sd_PSI), 
                width = 0.2) +
  theme_minimal() +
  labs(
    title = "Mean SRRM3 Exon 15 PSI Values by Cancer Type",
    y = "Mean PSI",
    x = "Cancer Type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(size = 12)
  )

# Try to save plots
tryCatch({
  ggsave("SRRM3_exon15_PSI_cancer_types_violin.pdf", p1, width = 15, height = 8)
  ggsave("SRRM3_exon15_PSI_cancer_types_means.pdf", p2, width = 15, height = 8)
}, error = function(e) {
  flog.error("Error saving plots: %s", e$message)
})

# # Perform statistical test (Kruskal-Wallis)
# kw_test <- tryCatch({
#   kruskal.test(PSI ~ cancer_type, data = plot_data)
# }, error = function(e) {
#   flog.error("Error in Kruskal-Wallis test: %s", e$message)
#   return(NULL)
# })

# # If Kruskal-Wallis is significant, perform pairwise Wilcoxon tests
# if(!is.null(kw_test) && kw_test$p.value < 0.05) {
#   pairwise_tests <- tryCatch({
#     pairwise.wilcox.test(
#       plot_data$PSI,
#       plot_data$cancer_type,
#       p.adjust.method = "BH"
#     )
#   }, error = function(e) {
#     flog.error("Error in pairwise Wilcoxon tests: %s", e$message)
#     return(NULL)
#   })
# }

# Try to save summary statistics
tryCatch({
  write.csv(summary_stats, "SRRM3_exon15_cancer_types_summary.csv", row.names = FALSE)
}, error = function(e) {
  flog.error("Error saving summary statistics: %s", e$message)
})

# # Print results
# flog.info("Summary statistics by cancer type:")
# print(summary_stats)
# 
# if(!is.null(kw_test)) {
#   flog.info("Kruskal-Wallis test p-value: %f", kw_test$p.value)
#   
#   if(exists("pairwise_tests") && !is.null(pairwise_tests)) {
#     flog.info("Significant differences between cancer types detected")
#     print(pairwise_tests)
#   }
# } 