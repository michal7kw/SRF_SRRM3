# Install required packages if not present
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("sparseMatrixStats")

# Load and verify key packages
required_packages <- c(
    "dplyr", "tidyverse", "survival", "survminer",
    "recount3", "biomaRt", "TCGAbiolinks",
    "SummarizedExperiment", "DESeq2", "sparseMatrixStats"
)

for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
        BiocManager::install(pkg)
    }
} 