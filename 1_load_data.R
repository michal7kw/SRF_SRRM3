library(recount3)
library(GenomicRanges)
library(dplyr)
library(SummarizedExperiment)

# First, let's get available projects
message("Checking available projects...")
projects <- available_projects()

# Filter for brain-related projects
message("\nFiltering for brain-related projects...")
# TCGA brain cancer projects (GBM: Glioblastoma, LGG: Lower Grade Glioma)
tcga_brain_projects <- subset(projects, 
                              file_source == "tcga" & 
                                project %in% c("GBM", "LGG"))
message("\nTCGA brain cancer projects:")
print(tcga_brain_projects$project)

# GTEx brain tissue projects
gtex_brain_projects <- subset(projects, 
                              file_source == "gtex" & 
                                grepl("BRAIN", project, ignore.case = TRUE))
message("\nGTEx brain tissue projects:")
print(gtex_brain_projects$project)

get_brain_expression_data <- function(project = c("TCGA", "GTEx"), 
                                      chromosome, 
                                      start_pos, 
                                      end_pos,
                                      brain_project = NULL) {
  project <- match.arg(project)
  
  tryCatch({
    # Create genomic ranges object
    region <- GRanges(
      seqnames = chromosome,
      ranges = IRanges(start = start_pos, end = end_pos)
    )
    
    # Get proper project name and create RSE object
    if (project == "TCGA") {
      if (is.null(brain_project)) {
        project_name <- tcga_brain_projects$project[1]  # Default to first brain cancer project
      } else {
        if (!brain_project %in% tcga_brain_projects$project) {
          stop("Invalid TCGA brain project. Choose from: ", 
               paste(tcga_brain_projects$project, collapse = ", "))
        }
        project_name <- brain_project
      }
      proj_info <- subset(projects, project == project_name & file_source == "tcga")
    } else {
      if (is.null(brain_project)) {
        project_name <- gtex_brain_projects$project[1]  # Default to first brain tissue project
      } else {
        if (!brain_project %in% gtex_brain_projects$project) {
          stop("Invalid GTEx brain project. Choose from: ", 
               paste(gtex_brain_projects$project, collapse = ", "))
        }
        project_name <- brain_project
      }
      proj_info <- subset(projects, project == project_name & file_source == "gtex")
    }
    
    message("Creating RSE object for project: ", project_name)
    rse <- create_rse(proj_info)
    
    # Ensure the RSE object was created successfully
    if (is.null(rse)) {
      stop("Failed to create RSE object")
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
      message("Available TCGA brain projects:")
      print(tcga_brain_projects$project)
    } else {
      message("Available GTEx brain projects:")
      print(gtex_brain_projects$project)
    }
    return(NULL)
  })
}

# Now let's try getting the data
message("\nAttempting to get TCGA brain cancer data (GBM)...")
tcga_brain_data <- get_brain_expression_data(
  project = "TCGA",
  chromosome = "chr7",
  start_pos = 76282502,
  end_pos = 76287486,
  brain_project = "GBM"  # Specifically get Glioblastoma data
)

message("\nAttempting to get GTEx brain tissue data...")
gtex_brain_data <- get_brain_expression_data(
  project = "GTEx",
  chromosome = "chr7",
  start_pos = 76282502,
  end_pos = 76287486  # Will use first available brain tissue project by default
)

# Save and summarize results if successful
if (!is.null(tcga_brain_data)) {
  save(tcga_brain_data, file = "srrm3_tcga_brain_data.RData")
  message("\nTCGA Brain Data Summary:")
  message("Project used: ", tcga_brain_data$project)
  message("Dimensions: ", paste(dim(tcga_brain_data$expression), collapse=" x "))
  message("Number of genes: ", nrow(tcga_brain_data$annotations))
}

if (!is.null(gtex_brain_data)) {
  save(gtex_brain_data, file = "srrm3_gtex_brain_data.RData")
  message("\nGTEx Brain Data Summary:")
  message("Project used: ", gtex_brain_data$project)
  message("Dimensions: ", paste(dim(gtex_brain_data$expression), collapse=" x "))
  message("Number of genes: ", nrow(gtex_brain_data$annotations))
}