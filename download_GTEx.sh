#!/bin/bash
#SBATCH --job-name=GTEx_download
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=INFINITE
#SBATCH --exclusive
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --error="./logs/GTEx_download.err"
#SBATCH --output="./logs/GTEx_download.out"

echo "my job start now" > ./logs/GTEx_download.log;
date >> ./logs/GTEx_download.log;

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

Rscript ./GTEx/download_gtex.R

echo "my job end now" >> ./logs/GTEx_download.log;
date >> ./logs/GTEx_download.log;
