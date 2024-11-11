library(recount3)
library(GenomicRanges)
library(dplyr)
library(SummarizedExperiment)

# First, let's see what projects are available
message("Checking available projects...")
projects <- available_projects()
message("\nTCGA projects:")
tcga_projects <- subset(projects, file_source == "tcga")
print(tcga_projects$project[1:5])  # Just print first few to check
message("\nGTEx projects:")
gtex_projects <- subset(projects, file_source == "gtex")
print(gtex_projects$project[1:5])  # Just print first few to check

get_expression_data <- function(project = c("TCGA", "GTEx"), 
                                chromosome, 
                                start_pos, 
                                end_pos) {
  project <- match.arg(project)
  
  tryCatch({
    # Create genomic ranges object
    region <- GRanges(
      seqnames = chromosome,
      ranges = IRanges(start = start_pos, end = end_pos)
    )
    
    # Get proper project name based on what's available
    if (project == "TCGA") {
      project_name <- tcga_projects$project[1]  # Use first available TCGA project
      rse <- create_rse(project_name)
    } else {
      project_name <- gtex_projects$project[1]  # Use first available GTEx project
      rse <- create_rse(project_name)
    }
    
    # Filter for our region of interest
    overlaps <- findOverlaps(rowRanges(rse), region)
    
    if (length(overlaps) == 0) {
      stop("No data found in specified region")
    }
    
    # Extract relevant data
    expression_data <- assay(rse)[queryHits(overlaps), ]
    
    # Add gene annotations
    gene_info <- as.data.frame(rowRanges(rse)[queryHits(overlaps)])
    
    # Add which project was used
    message("Successfully retrieved data from project: ", project_name)
    
    return(list(
      expression = expression_data,
      annotations = gene_info,
      project = project_name
    ))
    
  }, error = function(e) {
    message("Error occurred: ", e$message)
    message("\nPlease check:")
    message("1. Internet connection")
    message("2. Coordinates (", chromosome, ":", start_pos, "-", end_pos, ")")
    message("3. Project availability")
    if(project == "TCGA") {
      message("Available TCGA projects:")
      print(tcga_projects$project[1:5])
    } else {
      message("Available GTEx projects:")
      print(gtex_projects$project[1:5])
    }
    return(NULL)
  })
}

# Now let's try getting the data
message("\nAttempting to get TCGA data...")
tcga_data <- get_expression_data(
  project = "TCGA",
  chromosome = "chr7",
  start_pos = 76282502,
  end_pos = 76287486
)

message("\nAttempting to get GTEx data...")
gtex_data <- get_expression_data(
  project = "GTEx",
  chromosome = "chr7",
  start_pos = 76282502,
  end_pos = 76287486
)

# Save and summarize results if successful
if (!is.null(tcga_data)) {
  save(tcga_data, file = "srrm3_tcga_data3.RData")
  message("\nTCGA Data Summary:")
  message("Project used: ", tcga_data$project)
  message("Dimensions: ", paste(dim(tcga_data$expression), collapse=" x "))
  message("Number of genes: ", nrow(tcga_data$annotations))
}

if (!is.null(gtex_data)) {
  save(gtex_data, file = "srrm3_gtex_data3.RData")
  message("\nGTEx Data Summary:")
  message("Project used: ", gtex_data$project)
  message("Dimensions: ", paste(dim(gtex_data$expression), collapse=" x "))
  message("Number of genes: ", nrow(gtex_data$annotations))
}