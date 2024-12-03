import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
import scanpy as sc

def analyze_gene_trajectory_anndata(adata, gene_name, use_class, cell_type=None):
    """
    Analyze expression trajectory of a gene across development stages using AnnData
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix with genes as vars and cells as obs
    gene_name : str
        Name of gene to analyze
    cell_type : str, optional
        Specific cell type to analyze. If None, analyzes across all cells
        
    Returns:
    --------
    DataFrame with mean expression per developmental stage
    """
    # Create temporary dataframe with required columns
    df = pd.DataFrame({
        'stage_id': adata.obs['stage_id'],
        'cell_type': adata.obs[use_class],
        'expression': adata[:, gene_name].X.toarray().flatten() if scipy.sparse.issparse(adata.X) else adata[:, gene_name].X
    })
    
    # Filter for cell type if specified
    if cell_type:
        df = df[df['cell_type'] == cell_type]
    
    # Calculate mean expression per stage
    stage_order = ['Fetal', 'Neonatal', 'Infancy', 'Childhood', 'Adolescence', 'Adult']
    trajectory = (df.groupby('stage_id', observed=True)['expression']
                 .agg(['mean', 'std', 'count'])
                 .reset_index())
    
    # Add standard error
    trajectory['se'] = trajectory['std'] / np.sqrt(trajectory['count'])
    
    # Ensure stages are in correct order
    trajectory['stage_id'] = pd.Categorical(
        trajectory['stage_id'], 
        categories=stage_order, 
        ordered=True
    )
    
    return trajectory.sort_values('stage_id')

def plot_gene_trajectories(adata, genes, use_class, cell_types=None, use_violin=False):
    """
    Plot gene expression trajectories across development for multiple genes and cell types
    """
    if cell_types is None:
        cell_types = ['All']
    
    # First collect valid trajectories
    valid_plots = []
    for gene in genes:
        for cell_type in cell_types:
            trajectory = analyze_gene_trajectory_anndata(
                adata, gene, use_class, 
                cell_type if cell_type != 'All' else None
            )
            if len(trajectory) > 0:  # Only store if there is data
                valid_plots.append((gene, cell_type, trajectory))
    
    if not valid_plots:  # No valid data to plot
        return None
    
    # Calculate grid dimensions based on valid plots only
    n_valid = len(valid_plots)
    plots_per_row = min(6, n_valid)
    n_rows = (n_valid + plots_per_row - 1) // plots_per_row  # Ceiling division
    
    # Create figure with exact number of needed subplots
    fig = plt.figure(figsize=(4*plots_per_row, 4*n_rows))
    
    # Create only the necessary subplots
    for idx, (gene, cell_type, trajectory) in enumerate(valid_plots):
        row = idx // plots_per_row
        col = idx % plots_per_row
        ax = fig.add_subplot(n_rows, plots_per_row, idx + 1)
        
        if use_violin:
            sns.violinplot(data=trajectory, x='stage_id', y='mean', ax=ax)
        else:
            ax.plot(range(len(trajectory)), trajectory['mean'].values, 'b-', linewidth=2)
            ax.fill_between(
                range(len(trajectory)),
                (trajectory['mean'] - trajectory['se']).values,
                (trajectory['mean'] + trajectory['se']).values,
                alpha=0.2
            )
        
        # Customize plot
        ax.set_xticks(range(len(trajectory['stage_id'])))
        ax.set_xticklabels(trajectory['stage_id'], rotation=45)
        ax.set_xlabel('Developmental Stage')
        ax.set_ylabel(f'{gene} Expression')
        
        title = f'{gene}'
        if cell_type != 'All':
            title += f'\n{cell_type}'
        ax.set_title(title)

    plt.tight_layout()
    return fig

def compare_trajectories_by_cell_type(adata, gene, cell_types):
    """
    Compare trajectories of a single gene across multiple cell types
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    gene : str
        Gene to analyze
    cell_types : list
        List of cell types to compare
    """
    plt.figure(figsize=(12, 6))
    
    for cell_type in cell_types:
        trajectory = analyze_gene_trajectory_anndata(adata, gene, cell_type)
        plt.plot(range(len(trajectory)), trajectory['mean'].values,
                label=cell_type, linewidth=2, marker='o')
    
    plt.xlabel('Developmental Stage')
    plt.ylabel(f'{gene} Expression')
    plt.title(f'{gene} Expression Trajectories by Cell Type')
    plt.xticks(range(len(trajectory)), trajectory['stage_id'], rotation=45)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    return plt

def convert_age_to_days(age_str):
    """Convert age string to days for proper sorting"""
    if 'ga' in age_str:
        # Gestational weeks - convert to negative days to place before birth
        weeks = int(age_str.replace('ga', ''))
        return -280 + (weeks * 7)  # Assuming 40 weeks (280 days) is full term
    elif 'd' in age_str:
        # Days
        return int(age_str.replace('d', ''))
    elif 'yr' in age_str:
        # Years - convert to days
        return int(age_str.replace('yr', '')) * 365
    return 0

def assign_age_group(age):
    """Convert age to broader developmental periods"""
    if age.startswith('ga'):
        return 'Prenatal'
    elif age.endswith('d'):
        days = int(age[:-1])
        if days <= 100:
            return 'Early Postnatal'
        else:
            return 'Late Postnatal'
    elif age.endswith('yr'):
        years = int(age[:-2])
        if years <= 6:
            return 'Early Childhood'
        elif years <= 12:
            return 'Late Childhood'
        else:
            return 'Adolescent/Adult'
    return 'Unknown'

# Create a function to calculate mean expression for each cell type and age
def calculate_mean_expression(adata, gene, cell_types):
    results = []
    for age in age_list:
        for cell_type in cell_types:
            mask = (adata.obs.age == age) & (adata.obs.major_clust == cell_type)
            if mask.any():
                mean_expr = adata[mask, gene].X.mean()
                results.append({
                    'Age': age,
                    'Cell_Type': cell_type,
                    'Mean_Expression': mean_expr,
                })
    return pd.DataFrame(results)