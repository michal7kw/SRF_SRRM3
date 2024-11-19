# analysis of SRRM3 isoforms and SRRM4 across TCGA tumor types
library(recount3)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(SummarizedExperiment)
library(pheatmap)
library(DESeq2)
library(biomaRt)

# Define gene information
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
    transcript_start = 1945,
    transcript_end = 2023,
    length = 79
  )
)

SRRM4_INFO <- list(
  gene = list(
    name = "SRRM4",
    chr = "chr12",
    start = 118587761,  # GRCh38 coordinates
    end = 118708623
  )
)

# Function to get available TCGA projects
get_tcga_projects <- function() {
  message("Checking available projects...")
  projects <- available_projects()
  tcga_projects <- subset(projects, file_source == "tcga")
  return(tcga_projects)
}

# Function to get junction data
get_junction_data <- function(project_info) {
  message("\nGetting junction data for project: ", project_info$project)
  
  # Create RSE object with the proper project information
  rse_jxn <- create_rse(subset(available_projects(), 
                              project == project_info$project & 
                                file_source == "tcga"), 
                       type = "jxn")
  
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

# Function to calculate PSI for exon 15
calculate_psi <- function(junction_data) {
  # Find junctions that correspond to exon 15 inclusion/exclusion
  junctions <- junction_data$junction_info
  
  # Identify inclusion and exclusion junctions based on coordinates
  inclusion_jxn <- which(
    junctions$end_pos >= SRRM3_INFO$gene$start + SRRM3_INFO$exon15$transcript_start &
    junctions$start_pos <= SRRM3_INFO$gene$start + SRRM3_INFO$exon15$transcript_end
  )
  
  exclusion_jxn <- which(
    junctions$end_pos < SRRM3_INFO$gene$start + SRRM3_INFO$exon15$transcript_start |
    junctions$start_pos > SRRM3_INFO$gene$start + SRRM3_INFO$exon15$transcript_end
  )
  
  if(length(inclusion_jxn) == 0 || length(exclusion_jxn) == 0) {
    warning("No inclusion or exclusion junctions found for exon 15")
    return(rep(NA, ncol(junction_data$normalized_counts)))
  }
  
  # Calculate PSI
  inclusion_counts <- colSums(junction_data$normalized_counts[inclusion_jxn, , drop = FALSE])
  exclusion_counts <- colSums(junction_data$normalized_counts[exclusion_jxn, , drop = FALSE])
  
  psi <- inclusion_counts / (inclusion_counts + exclusion_counts)
  return(psi)
}

# Function to process junction data for a specific project
process_project_data <- function(project_info) {
  message("\nProcessing project: ", project_info$project)
  
  tryCatch({
    # Get SRRM3 junction data
    srrm3_data <- get_junction_data(project_info)
    
    # Get SRRM4 gene expression data
    rse_gene <- create_rse(subset(available_projects(), 
                                project == project_info$project & 
                                  file_source == "tcga"), 
                         type = "gene")
    
    srrm4_idx <- which(rowData(rse_gene)$gene_name == "SRRM4")
    if(length(srrm4_idx) == 0) {
      warning("SRRM4 gene not found in ", project_info$project)
      srrm4_counts <- matrix(NA, nrow = 1, ncol = ncol(rse_gene))
    } else {
      srrm4_counts <- assay(rse_gene)[srrm4_idx, , drop = FALSE]
    }
    
    # Calculate PSI (Percent Spliced In) for SRRM3 exon 15
    psi_values <- calculate_psi(srrm3_data)
    
    # Combine data
    result <- list(
      project = project_info$project,
      srrm3_junction_data = srrm3_data,
      srrm4_expression = srrm4_counts,
      psi_values = psi_values
    )
    
    return(result)
  }, error = function(e) {
    warning("Error processing project ", project_info$project, ": ", e$message)
    return(NULL)
  })
}

# Main analysis pipeline
main_analysis <- function() {
  # Get all TCGA projects
  projects <- get_tcga_projects()
  
  # Process each project
  results <- list()
  for(i in 1:nrow(projects)) {
    results[[i]] <- process_project_data(projects[i,])
  }
  
  # Combine results across projects
  combined_results <- combine_project_results(results)
  
  # Generate visualizations
  plot_isoform_distribution(combined_results)
  plot_srrm4_correlation(combined_results)
  plot_cancer_type_comparison(combined_results)
  
  return(combined_results)
}

# Function to combine results across projects
combine_project_results <- function(results) {
  # Remove NULL results from failed processing
  results <- results[!sapply(results, is.null)]
  
  if(length(results) == 0) {
    stop("No valid results to combine")
  }
  
  # Get project names
  project_names <- sapply(results, function(x) x$project)
  
  # Initialize matrices for combined data
  n_projects <- length(results)
  
  # Combine PSI values
  psi_matrix <- matrix(NA, nrow = n_projects, ncol = 1)
  rownames(psi_matrix) <- project_names
  
  # Combine SRRM4 expression
  srrm4_matrix <- matrix(NA, nrow = n_projects, ncol = 1)
  rownames(srrm4_matrix) <- project_names
  
  # Fill matrices with mean values for each project
  for(i in seq_len(n_projects)) {
    # PSI values
    psi_vals <- results[[i]]$psi_values
    psi_matrix[i, 1] <- mean(psi_vals, na.rm = TRUE)
    
    # SRRM4 expression
    srrm4_expr <- results[[i]]$srrm4_expression
    srrm4_matrix[i, 1] <- mean(as.numeric(srrm4_expr), na.rm = TRUE)
  }
  
  return(list(
    psi = psi_matrix,
    srrm4 = srrm4_matrix,
    projects = project_names
  ))
}

# Visualization functions
plot_isoform_distribution <- function(combined_results) {
  # Create data frame for plotting
  plot_data <- data.frame(
    cancer_type = combined_results$projects,
    psi_value = combined_results$psi[,1]
  )
  
  ggplot(plot_data, aes(x = cancer_type, y = psi_value)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "SRRM3 Exon 15 Inclusion Levels Across Cancer Types",
      x = "Cancer Type",
      y = "PSI (Percent Spliced In)"
    )
  
  ggsave("output/SRRM3_isoform_distribution.pdf", width = 12, height = 8)
}

plot_srrm4_correlation <- function(combined_results) {
  # Create data frame for correlation plot
  cor_data <- data.frame(
    cancer_type = combined_results$projects,
    srrm4_expr = combined_results$srrm4[,1],
    srrm3_psi = combined_results$psi[,1]
  )
  
  ggplot(cor_data, aes(x = srrm4_expr, y = srrm3_psi, label = cancer_type)) +
    geom_point() +
    geom_text_repel() +
    geom_smooth(method = "lm", se = TRUE) +
    theme_minimal() +
    labs(
      title = "Correlation between SRRM4 Expression and SRRM3 Exon 15 Inclusion",
      x = "SRRM4 Expression (mean normalized counts)",
      y = "SRRM3 Exon 15 PSI"
    )
  
  ggsave("output/SRRM4_SRRM3_correlation.pdf", width = 8, height = 8)
}

plot_cancer_type_comparison <- function(combined_results) {
  # Create matrix for heatmap
  plot_matrix <- rbind(
    combined_results$psi[,1],
    1 - combined_results$psi[,1],
    combined_results$srrm4[,1]
  )
  
  rownames(plot_matrix) <- c(
    "SRRM3_with_exon15",
    "SRRM3_without_exon15",
    "SRRM4"
  )
  colnames(plot_matrix) <- combined_results$projects
  
  # Scale the data for better visualization
  plot_matrix_scaled <- t(scale(t(plot_matrix)))
  
  pheatmap(
    plot_matrix_scaled,
    main = "Expression Patterns Across Cancer Types",
    filename = "output/expression_heatmap.pdf",
    width = 12,
    height = 6,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation"
  )
}

# Run the analysis
results <- main_analysis()

# Save results
saveRDS(results, "output/SRRM3_SRRM4_analysis_results.rds")

