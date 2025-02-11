import pandas as pd
import numpy as np
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def prepare_survival_data(clinical_data, molecular_data, gene_id, method='expression'):
  """
    Prepare survival data by merging clinical and molecular data
    
    Parameters:
    -----------
    clinical_data : pandas DataFrame
        TCGA clinical data with columns: 'patient_id', 'days_to_death'/'days_to_last_followup', 'vital_status'
    molecular_data : pandas DataFrame
        Either gene expression or PSI data with patient IDs as columns
    gene_id : str
        Gene identifier (e.g., 'SRRM3' or 'SRRM4')
    method : str
        'expression' or 'psi' to indicate the type of molecular data
    
    Returns:
    --------
    DataFrame with survival data and molecular measurements
    """
  # Merge clinical and molecular data
  if method == 'expression':
    mol_data = molecular_data.loc[gene_id]
  else:  # PSI
    mol_data = molecular_data.loc[molecular_data['gene_id'] == gene_id]

  # Create survival time and event columns
  survival_data = pd.DataFrame({
    'time': clinical_data['days_to_death'].fillna(clinical_data['days_to_last_followup']),
    'event': clinical_data['vital_status'].map({'Dead': 1, 'Alive': 0}),
    'molecular_value': mol_data
  })

  # Remove any rows with missing data
  survival_data = survival_data.dropna()

  # Stratify patients into high/low groups based on median
  median_value = survival_data['molecular_value'].median()
  survival_data['group'] = survival_data['molecular_value'].map(
    lambda x: 'High' if x > median_value else 'Low'
  )

  return survival_data

def plot_survival_curve(survival_data, gene_id, method='expression'):
  """
    Plot Kaplan-Meier survival curves and calculate log-rank test
    
    Parameters:
    -----------
    survival_data : pandas DataFrame
        Prepared survival data from prepare_survival_data function
    gene_id : str
        Gene identifier for plot title
    method : str
        'expression' or 'psi' to indicate type of analysis
    """
  # Initialize KM model
  kmf = KaplanMeierFitter()

  # Set up plot
  plt.figure(figsize=(10, 6))

  # Plot high group
  high_mask = survival_data['group'] == 'High'
  kmf.fit(
    survival_data[high_mask]['time'],
    survival_data[high_mask]['event'],
    label=f'High {method}'
  )
  kmf.plot()

  # Plot low group
  low_mask = survival_data['group'] == 'Low'
  kmf.fit(
    survival_data[low_mask]['time'],
    survival_data[low_mask]['event'],
    label=f'Low {method}'
  )
  kmf.plot()

  # Calculate log-rank test
  from lifelines.statistics import logrank_test
  log_rank = logrank_test(
    survival_data[high_mask]['time'],
    survival_data[low_mask]['time'],
    survival_data[high_mask]['event'],
    survival_data[low_mask]['event']
  )

  # Add plot details
  plt.title(f'Survival Analysis Stratified by {gene_id} {method.capitalize()}\np = {log_rank.p_value:.2e}')
  plt.xlabel('Time (days)')
  plt.ylabel('Survival Probability')
  plt.grid(True)

  return log_rank.p_value

def run_survival_analysis(clinical_file, molecular_file, gene_id, method='expression'):
  """
    Main function to run the complete survival analysis
    
    Parameters:
    -----------
    clinical_file : str
        Path to clinical data file
    molecular_file : str
        Path to molecular data file (expression or PSI)
    gene_id : str
        Gene identifier (e.g., 'SRRM3' or 'SRRM4')
    method : str
        'expression' or 'psi' to indicate type of analysis
    """
  # Load data
  clinical_data = pd.read_csv(clinical_file)
  molecular_data = pd.read_csv(molecular_file)

  # Prepare survival data
  survival_data = prepare_survival_data(
    clinical_data,
    molecular_data,
    gene_id,
    method
  )

  # Plot and get p-value
  p_value = plot_survival_curve(survival_data, gene_id, method)

  return survival_data, p_value
