from Survival_in_py_functions import * 
import os

os.chdir('/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3/TCGA')

data, pval = run_survival_analysis(
    'tcga_clinical.csv',
    'tcga_expression.csv',
    'SRRM3',
    method='expression'
)

data, pval = run_survival_analysis(
    'tcga_clinical.csv',
    'tcga_psi.csv',
    'SRRM3',
    method='psi'
)

data, pval = run_survival_analysis(
    'tcga_clinical.csv',
    'tcga_expression.csv',
    'SRRM4',
    method='expression'
)