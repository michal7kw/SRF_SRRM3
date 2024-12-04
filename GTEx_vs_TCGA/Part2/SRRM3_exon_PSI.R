# Libraries
library(recount3)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(SummarizedExperiment)
library(progress)
library(futile.logger)

# Set up logging
flog.appender(appender.file("srrm3_exon_analysis.log"))
flog.threshold(DEBUG)

# Read and process exon coordinates
exon_coords <- read.delim("SRRM3_exon_coordinates.tsv", stringsAsFactors = FALSE)

# Define SRRM3 information
SRRM3_INFO <- list(
  gene = list(
    name = "SRRM3",
    chr = "chr7",
    start = as.numeric(min(exon_coords$start)),  # Ensure numeric
    end = as.numeric(max(exon_coords$end)),      # Ensure numeric
    strand = "+"
  ),
  exon15 = {
    # Find exon 15 coordinates from the data
    exon15_row <- exon_coords[exon_coords$exon_number == 15, ]
    list(
      start = as.numeric(exon15_row$start),
      end = as.numeric(exon15_row$end),
      length = as.numeric(exon15_row$exon_length)
    )
  },
  # Define constitutive exons (exons present in all transcripts)
  constitutive_exons = {
    # Find exons that appear in all transcripts
    exon_counts <- table(paste(exon_coords$start, exon_coords$end))
    max_count <- max(exon_counts)
    constitutive_pairs <- names(exon_counts[exon_counts == max_count])
    
    # Convert to list format
    lapply(strsplit(constitutive_pairs, " "), function(pair) {
      list(
        start = as.numeric(pair[1]),
        end = as.numeric(pair[2])
      )
    })
  }
)

# Function to create RSE object with exon counts
create_exon_rse <- function(project_info) {
  flog.info("Creating RSE for project: %s", project_info$project)
  
  # Create RSE with exon-level counts
  rse <- create_rse(
    project_info,
    type = "exon",
    verbose = TRUE
  )
  
  # Filter to SRRM3 region
  chr <- SRRM3_INFO$gene$chr
  relevant_rows <- which(
    as.character(seqnames(rowRanges(rse))) == chr &
      start(rowRanges(rse)) >= SRRM3_INFO$gene$start &
      end(rowRanges(rse)) <= SRRM3_INFO$gene$end
  )
  
  if(length(relevant_rows) == 0) {
    flog.warn("No SRRM3 exons found in the dataset")
    return(NULL)
  }
  
  rse_filtered <- rse[relevant_rows, ]
  return(rse_filtered)
}

# Function to identify exon 15 and constitutive exons
identify_exons <- function(rse) {
  exon_ranges <- rowRanges(rse)
  
  # Create GRanges for exon 15
  exon15_range <- GRanges(
    seqnames = SRRM3_INFO$gene$chr,
    ranges = IRanges(
      start = SRRM3_INFO$exon15$start[[1]],
      end = SRRM3_INFO$exon15$end[[1]]
    ),
    strand = SRRM3_INFO$gene$strand
  )
  
  # Find exon 15 in the data
  exon15_idx <- which(overlapsAny(exon_ranges, exon15_range))
  
  if(length(exon15_idx) == 0) {
    flog.error("Exon 15 not found in the dataset")
    return(NULL)
  }
  
  # Create GRanges for constitutive exons
  constitutive_ranges <- GRangesList(
    lapply(SRRM3_INFO$constitutive_exons, function(exon) {
      GRanges(
        seqnames = SRRM3_INFO$gene$chr,
        ranges = IRanges(
          start = exon$start,
          end = exon$end
        ),
        strand = SRRM3_INFO$gene$strand
      )
    })
  )
  
  # Find constitutive exons in the data
  constitutive_idx <- which(overlapsAny(exon_ranges, unlist(constitutive_ranges)))
  
  if(length(constitutive_idx) == 0) {
    flog.error("No constitutive exons found in the dataset")
    return(NULL)
  }
  
  return(list(
    exon15_idx = exon15_idx,
    constitutive_idx = constitutive_idx
  ))
}

# Function to calculate PSI using exon counts
calculate_exon_psi <- function(rse, exon_indices) {
  flog.info("Calculating PSI values")
  
  # Get counts
  counts <- assay(rse)
  
  # Get exon 15 counts - ensure we get a numeric vector
  exon15_counts <- counts[exon_indices$exon15_idx[[1]], ]
  
  # Get mean counts of constitutive exons - ensure we get numeric values
  constitutive_idx <- unlist(exon_indices$constitutive_idx)
  constitutive_counts <- colMeans(counts[constitutive_idx, , drop = FALSE])
  
  # Calculate PSI
  psi <- (exon15_counts / constitutive_counts) * 100
  
  # Handle edge cases
  psi[is.infinite(psi)] <- NA
  psi[psi > 100] <- 100
  
  return(psi)
}

# Function to process a single project
process_project <- function(project_info) {
  flog.info("Processing project: %s", project_info$project)
  
  # Create RSE object
  rse <- create_exon_rse(project_info)
  if(is.null(rse)) {
    return(NULL)
  }
  
  # Identify relevant exons
  exon_indices <- identify_exons(rse)
  if(is.null(exon_indices)) {
    return(NULL)
  }
  
  # Calculate PSI
  psi_values <- calculate_exon_psi(rse, exon_indices)
  
  # Create results object
  results <- list(
    project = project_info$project,
    psi_values = psi_values,
    sample_info = colData(rse),
    mean_psi = mean(psi_values, na.rm = TRUE),
    median_psi = median(psi_values, na.rm = TRUE),
    sd_psi = sd(psi_values, na.rm = TRUE)
  )
  
  return(results)
}

# Function to visualize results
plot_psi_distribution <- function(results) {
  # Create data frame for plotting
  plot_data <- data.frame(
    project = results$project,
    psi = results$psi_values
  )
  
  # Create violin plot
  p <- ggplot(plot_data, aes(x = project, y = psi)) +
    geom_violin(fill = "lightblue", alpha = 0.5) +
    geom_boxplot(width = 0.1) +
    theme_minimal() +
    labs(
      title = "SRRM3 Exon 15 PSI Distribution",
      x = "Project",
      y = "PSI (%)"
    )
  
  # Save plot
  ggsave(paste0("SRRM3_exon15_PSI_", results$project, ".pdf"), 
         p, width = 8, height = 6)
  
  return(p)
}

# Main execution function
main <- function() {
  flog.info("Starting SRRM3 exon 15 PSI analysis")
  
  # Get available projects
  projects <- available_projects()
  tcga_projects <- subset(projects, file_source == "tcga")
  
  # Process each project
  results_list <- list()
  for(i in seq_len(nrow(tcga_projects))) {
    results_list[[i]] <- process_project(tcga_projects[i,])
  }
  
  # Remove NULL results
  results_list <- results_list[!sapply(results_list, is.null)]
  
  # Create summary table
  summary_df <- do.call(rbind, lapply(results_list, function(x) {
    data.frame(
      project = x$project,
      mean_psi = x$mean_psi,
      median_psi = x$median_psi,
      sd_psi = x$sd_psi,
      n_samples = length(x$psi_values)
    )
  }))
  
  # Save results
  write.csv(summary_df, "SRRM3_exon15_PSI_summary.csv", row.names = FALSE)
  saveRDS(results_list, "SRRM3_exon15_PSI_full_results.rds")
  
  # Create visualizations
  lapply(results_list, plot_psi_distribution)
  
  flog.info("Analysis completed successfully")
}

# Run the analysis with error handling
tryCatch({
  main()
}, error = function(e) {
  flog.error("Analysis failed: %s", e$message)
  stop(e)
}) 