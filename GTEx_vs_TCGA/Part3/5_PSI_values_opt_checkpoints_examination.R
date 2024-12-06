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
flog.appender(appender.file("srrm3_analysis.log"))
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
  checkpoint_name <- paste0(project_info$project, "_rse")
  rse_filtered <- load_checkpoint(checkpoint_name)
  
  if(!is.null(rse_filtered)) {
    flog.info("Loaded RSE from checkpoint")
    return(rse_filtered)
  }
  
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
    
    save_checkpoint(rse_filtered, checkpoint_name)
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

# Main analysis
main_analysis <- function() {
  # Get projects
  projects <- available_projects()
  
  # Process TCGA data
  tcga_project <- subset(projects, project == "GBM" & file_source == "tcga")
  tcga_rse <- create_rse_safe(tcga_project)
  tcga_psi <- calculate_exon15_psi(tcga_rse)
  
  # Process GTEx data
  gtex_project <- subset(projects, project == "BRAIN" & file_source == "gtex")
  gtex_rse <- create_rse_safe(gtex_project)
  gtex_psi <- calculate_exon15_psi(gtex_rse)
  
  # Create comparison plot
  if(!is.null(tcga_psi) && !is.null(gtex_psi)) {
    # Prepare data for plotting
    plot_data <- data.frame(
      dataset = c(rep("TCGA", length(tcga_psi$psi)), 
                 rep("GTEx", length(gtex_psi$psi))),
      PSI = c(tcga_psi$psi, gtex_psi$psi)
    )
    
    # Remove NA values
    plot_data <- plot_data[!is.na(plot_data$PSI),]
    
    # Create plot
    p <- ggplot(plot_data, aes(x = dataset, y = PSI, fill = dataset)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
      theme_minimal() +
      labs(
        title = "SRRM3 Exon 15 PSI Values",
        y = "Percent Spliced In (PSI)",
        x = "Dataset"
      ) +
      theme(legend.position = "none")
    
    # Save results
    ggsave("SRRM3_exon15_PSI_comparison.pdf", p, width = 8, height = 6)
    
    # Perform statistical test
    test_result <- wilcox.test(
      PSI ~ dataset, 
      data = plot_data,
      alternative = "two.sided"
    )
    
    # Save summary statistics
    summary_stats <- plot_data %>%
      group_by(dataset) %>%
      summarize(
        mean_PSI = mean(PSI),
        median_PSI = median(PSI),
        sd_PSI = sd(PSI),
        n_samples = n()
      )
    
    write.csv(summary_stats, "SRRM3_exon15_summary.csv", row.names = FALSE)
    
    # Print results
    print(summary_stats)
    print(paste("Wilcoxon test p-value:", test_result$p.value))
  }
}

# Run analysis
main_analysis()
