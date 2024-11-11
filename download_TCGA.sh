#!/bin/bash
#SBATCH --job-name=TCGA_download
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=INFINITE
#SBATCH --exclusive
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --error="./logs/TCGA_download.err"
#SBATCH --output="./logs/TCGA_download.out"

echo "my job start now" > ./logs/TCGA_download.log;
date >> ./logs/TCGA_download.log;

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

Rscript ./TCGA/download_tcga.R

echo "my job end now" >> ./logs/TCGA_download.log;
date >> ./logs/TCGA_download.log;
