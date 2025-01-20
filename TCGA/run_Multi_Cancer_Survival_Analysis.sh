#!/bin/bash
#SBATCH --job-name=Multi_Cancer_Survival_Analysis
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/Multi_Cancer_Survival_Analysis_%A_%a.err"
#SBATCH --output="logs/Multi_Cancer_Survival_Analysis_%A_%a.out"
#SBATCH --array=0-4  # One job per cancer type (5 cancer types)

# Load required modules and activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3/TCGA/
# Create logs directory if it doesn't exist
mkdir -p logs

# Create a temporary R script for this specific array job
TEMP_SCRIPT="temp_analysis_${SLURM_ARRAY_TASK_ID}.R"

# Define cancer types array
declare -a CANCER_TYPES=("ACC" "UVM" "SKCM" "LGG" "GBM")

# Get the cancer type for this array job
CANCER_TYPE=${CANCER_TYPES[$SLURM_ARRAY_TASK_ID]}

# Create the R script for this specific cancer type
cat << EOF > $TEMP_SCRIPT
# Source the main analysis script
source("Multi_Cancer_Survival_Analysis.R")

# Create results directory
if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

# Run analysis for single cancer type
results <- perform_multi_cancer_analysis(
  cancer_types = c("$CANCER_TYPE"),
  analysis_type = "PSI",
  gene = "SRRM3",
  grouping_method = "quartile"
)

# Save results for this cancer type
saveRDS(results, file.path("results", paste0("${CANCER_TYPE}_results.rds")))
EOF

# Run the R script
Rscript $TEMP_SCRIPT

# Clean up temporary script
rm $TEMP_SCRIPT

# If this is the last array job, combine all results
if [ $SLURM_ARRAY_TASK_ID -eq 4 ]; then
    # Create the combination R script
    cat << EOF > combine_results.R
    # Source the main analysis script
    source("Multi_Cancer_Survival_Analysis.R")

    # List all individual result files
    result_files <- list.files("results", pattern = "*_results.rds", full.names = TRUE)
    
    # Load and combine all results
    all_results <- list()
    for (file in result_files) {
        cancer_type <- gsub("_results.rds", "", basename(file))
        all_results[[cancer_type]] <- readRDS(file)
    }
    
    # Generate and save comparative analyses
    comparative_results <- generate_comparative_analyses(all_results)
    save_combined_results(all_results, comparative_results, "results/multi_cancer_analysis")
    
    # Clean up individual result files
    file.remove(result_files)
EOF

    # Run the combination script
    Rscript combine_results.R
    rm combine_results.R
fi



    