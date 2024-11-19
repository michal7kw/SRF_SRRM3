# Standard library imports
import warnings
from pathlib import Path
from typing import List, Union, Tuple

# Third party imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import anndata
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
from umap import UMAP

# Project imports
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

from typing import List, Union, Tuple, Optional

import os
# Set the current working directory
os.chdir('/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3')

def load_multi_region_data(cell_metadata: pd.DataFrame,
        abc_cache: AbcProjectCache,
        datasets: List[str], 
        genes: List[str],
        dataset_regions: Optional[List[str]] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load expression data for multiple datasets
    
    Args:
        cell_metadata: DataFrame containing cell metadata
        abc_cache: AbcProjectCache instance for accessing data
        datasets: List of dataset names (e.g., ['WMB-10X'])
        genes: List of gene names
        dataset_regions: Optional list of specific regions to load based on region_of_interest_acronym.
    """
    all_expression = []
    all_metadata = []
    
    # Get all dataset versions from metadata
    dataset_versions = cell_metadata['dataset_label'].unique()
    
    for dataset in datasets:
        # Get all versions for this dataset
        versions = [v for v in dataset_versions if v.startswith(dataset)]
        
        for version in versions:
            # Get feature matrices for this version
            feature_matrices = cell_metadata[
                cell_metadata['dataset_label'] == version
            ]['feature_matrix_label'].unique()
            
            for feature_matrix in feature_matrices:
                # Get cells for this feature matrix
                matrix_cells = cell_metadata[
                    cell_metadata['feature_matrix_label'] == feature_matrix
                ]
                
                # If specific regions requested, filter cells
                if dataset_regions is not None:
                    matrix_cells = matrix_cells[
                        matrix_cells['region_of_interest_acronym'].isin(dataset_regions)
                    ]
                    
                    if matrix_cells.empty:
                        continue
                
                print(f"Processing {feature_matrix}")
                
                try:
                    # Load expression data using the exact feature matrix label
                    adata = anndata.read_h5ad(
                        abc_cache.get_data_path(
                            directory=version,  # Use version instead of dataset
                            file_name=f"{feature_matrix}/log2"
                        ),
                        backed='r'
                    )
                    
                    # Get gene indices
                    gene_mask = [g in genes for g in adata.var['gene_symbol']]
                    if not any(gene_mask):
                        print(f"No requested genes found in {feature_matrix}")
                        continue
                        
                    # Get common cells
                    common_cells = list(set(adata.obs.index) & set(matrix_cells.index))
                    if not common_cells:
                        print(f"No matching cells found in {feature_matrix}")
                        continue
                    
                    print(f"Found {len(common_cells)} cells in {feature_matrix}")
                    
                    # Extract expression data
                    gene_ids = adata.var_names[gene_mask]
                    data = pd.DataFrame(
                        adata[common_cells, gene_ids].X.toarray(),
                        index=common_cells,
                        columns=adata.var['gene_symbol'][gene_mask]
                    )
                    
                    # Add region and dataset information
                    meta = matrix_cells.loc[common_cells].copy()
                    data['region'] = meta['region_of_interest_acronym']
                    data['dataset'] = version
                    
                    all_expression.append(data)
                    all_metadata.append(meta)
                    
                    adata.file.close()
                    del adata
                    
                except Exception as e:
                    print(f"Error processing {feature_matrix}: {str(e)}")
                    continue
    
    if not all_expression:
        error_msg = "No data could be loaded for any region. "
        error_msg += "Please check that:\n"
        error_msg += "1. The specified datasets and regions exist\n"
        error_msg += "2. The requested genes are present in the data\n" 
        error_msg += "3. There are matching cells in the metadata"
        raise ValueError(error_msg)
        
    return pd.concat(all_expression), pd.concat(all_metadata)


def plot_region_gene_comparison(expression_data: pd.DataFrame, 
                              metadata: pd.DataFrame,
                              genes: List[str],
                              group_by: str,
                              figsize: Tuple[int, int] = (15, 10),
                              top_n: int = 15):
    """
    Create separate plots for each brain region with side-by-side gene comparisons
    """
    # Prepare the data
    plot_data = pd.DataFrame({
        'Expression': pd.concat([expression_data[gene] for gene in genes]),
        'Gene': np.repeat(genes, len(expression_data)),
        'Group': pd.concat([metadata[group_by] for _ in genes]),
        'Region': pd.concat([expression_data['region'] for _ in genes])
    })
    
    # Get unique regions
    regions = sorted(plot_data['Region'].unique())
    
    # Create a figure with subplots for each region
    n_regions = len(regions)
    n_cols = 2
    n_rows = (n_regions + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(figsize[0]*n_cols, figsize[1]*n_rows))
    if n_regions == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    # Plot for each region
    for idx, region in enumerate(regions):
        if idx < len(axes):
            ax = axes[idx]
            region_data = plot_data[plot_data['Region'] == region]
            
            # Calculate mean expression per group and gene for this region
            mean_expr = region_data.groupby(['Group', 'Gene'])['Expression'].mean().reset_index()
            mean_expr_pivot = mean_expr.pivot(index='Group', columns='Gene', values='Expression')
            
            # Get top N groups for each gene separately
            top_groups_per_gene = []
            for gene in genes:
                top_for_gene = mean_expr_pivot[gene].sort_values(ascending=False).head(top_n).index
                top_groups_per_gene.extend(top_for_gene)
            
            # Get unique top groups across both genes
            top_groups = pd.Index(sorted(set(top_groups_per_gene)))
            
            # Filter data for top groups
            region_data_filtered = region_data[region_data['Group'].isin(top_groups)]
            
            # Create the box plot without outliers
            sns.boxplot(data=region_data_filtered, 
                       x='Group', 
                       y='Expression',
                       hue='Gene',
                       order=top_groups,
                       palette=['#2ecc71', '#e74c3c'],  # Green for Srrm3, Red for Srrm4
                       showfliers=False,  # Don't show outliers
                       ax=ax)
            
            # Customize the plot
            ax.set_title(f'Region: {region}', fontsize=12, pad=20)
            ax.set_xlabel('Cell Type', fontsize=10)
            ax.set_ylabel('Expression Level', fontsize=10)
            
            # Rotate x-axis labels
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
            
            # Add mean expression values above each box
            def add_mean_labels(data, ax):
                for i, group in enumerate(top_groups):
                    for j, gene in enumerate(genes):
                        mean_val = data[
                            (data['Group'] == group) & 
                            (data['Gene'] == gene)
                        ]['Expression'].mean()
                        ax.text(i, mean_val, f'{mean_val:.1f}', 
                               ha='center', va='bottom', fontsize=8)
            
            add_mean_labels(region_data_filtered, ax)
            
            # Adjust legend
            ax.legend(title='Gene', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Remove empty subplots if any
    for idx in range(len(regions), len(axes)):
        fig.delaxes(axes[idx])
    
    plt.tight_layout()
    
    # Add statistical summary
    stats_text = "Statistical Summary:\n\n"
    for region in regions:
        stats_text += f"\nRegion: {region}\n"
        region_data = plot_data[plot_data['Region'] == region]
        for group in top_groups:
            group_data = region_data[region_data['Group'] == group]
            stats_text += f"\n{group}:\n"
            for gene in genes:
                gene_data = group_data[group_data['Gene'] == gene]
                stats_text += (f"  {gene}: mean={gene_data['Expression'].mean():.2f}, "
                             f"median={gene_data['Expression'].median():.2f}, "
                             f"n={len(gene_data)}\n")
    
    print(stats_text)
    return fig, axes


def create_dimension_reduction_plot(expression_data: pd.DataFrame,
                                  metadata: pd.DataFrame,
                                  method: str = 'umap',
                                  color_by: Union[str, List[str]] = 'class_label',
                                  genes: List[str] = None,
                                  figsize: Tuple[int, int] = (20, 15)):
    """Create UMAP or t-SNE plots colored by different features"""
    # Prepare expression data
    X = StandardScaler().fit_transform(expression_data[genes])
    
    # Perform dimension reduction
    if method.lower() == 'umap':
        reducer = UMAP(random_state=42, n_neighbors=30, min_dist=0.3)
    else:
        reducer = TSNE(random_state=42, perplexity=30)
    
    print(f"Performing {method.upper()} dimension reduction...")
    embedding = reducer.fit_transform(X)
    
    # Create plots
    if isinstance(color_by, str):
        color_by = [color_by]
    
    n_plots = len(color_by) + len(genes)
    n_cols = min(3, n_plots)
    n_rows = (n_plots + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    if n_plots == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    plot_idx = 0
    
    # Plot categorical variables
    for feature in color_by:
        categories = metadata[feature].astype('category')
        scatter = axes[plot_idx].scatter(
            embedding[:, 0], 
            embedding[:, 1],
            c=categories.cat.codes,
            cmap='tab20',
            alpha=0.6,
            s=10
        )
        axes[plot_idx].set_title(f'Colored by {feature}')
        
        # Create legend with actual category names
        legend_elements = [plt.Line2D([0], [0], marker='o', color='w', 
                                    markerfacecolor=plt.cm.tab20(i/len(categories.cat.categories)), 
                                    label=cat, markersize=10)
                         for i, cat in enumerate(categories.cat.categories)]
        
        axes[plot_idx].legend(handles=legend_elements,
                            title=feature,
                            bbox_to_anchor=(1.05, 1),
                            loc='upper left',
                            ncol=1)
        plot_idx += 1
    
    # Plot gene expression
    for gene in genes:
        scatter = axes[plot_idx].scatter(
            embedding[:, 0],
            embedding[:, 1],
            c=expression_data[gene],
            cmap='viridis',
            alpha=0.6,
            s=10
        )
        axes[plot_idx].set_title(f'{gene} Expression')
        plt.colorbar(scatter, ax=axes[plot_idx])
        plot_idx += 1
    
    # Remove empty subplots
    for idx in range(plot_idx, len(axes)):
        fig.delaxes(axes[idx])
    
    plt.tight_layout()
    return fig, axes

    

def plot_region_gene_comparison_cleaned(expression_data: pd.DataFrame, 
                              metadata: pd.DataFrame,
                              genes: List[str],
                              group_by: str,
                              figsize: Tuple[int, int] = (15, 10),
                              top_n: int = 15,
                              min_expr_threshold: float = 0.1):  # Added threshold parameter
    """
    Create separate plots for each brain region with side-by-side gene comparisons
    using violin plots and expression statistics
    """
    # Prepare the data (same as before)
    plot_data = pd.DataFrame({
        'Expression': pd.concat([expression_data[gene] for gene in genes]),
        'Gene': np.repeat(genes, len(expression_data)),
        'Group': pd.concat([metadata[group_by] for _ in genes]),
        'Region': pd.concat([expression_data['region'] for _ in genes])
    })
    
    # Get unique regions and setup subplots (same as before)
    regions = sorted(plot_data['Region'].unique())
    n_regions = len(regions)
    n_cols = 2
    n_rows = (n_regions + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(figsize[0]*n_cols, figsize[1]*n_rows))
    if n_regions == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    # Plot for each region
    for idx, region in enumerate(regions):
        if idx < len(axes):
            ax = axes[idx]
            region_data = plot_data[plot_data['Region'] == region]
            
            # Calculate mean expression and % expressing cells per group and gene
            stats = []
            for group in region_data['Group'].unique():
                for gene in genes:
                    group_gene_data = region_data[
                        (region_data['Group'] == group) & 
                        (region_data['Gene'] == gene)
                    ]
                    
                    # Calculate statistics
                    total_cells = len(group_gene_data)
                    expressing_cells = len(group_gene_data[group_gene_data['Expression'] > min_expr_threshold])
                    mean_expr = group_gene_data['Expression'].mean()
                    
                    stats.append({
                        'Group': group,
                        'Gene': gene,
                        'Mean': mean_expr,
                        'Pct_expressing': (expressing_cells / total_cells) * 100
                    })
            
            stats_df = pd.DataFrame(stats)
            
            # Get top N groups based on either mean expression or % expressing
            top_groups_per_gene = []
            for gene in genes:
                gene_stats = stats_df[stats_df['Gene'] == gene]
                # You can change this to sort by 'Pct_expressing' instead
                top_for_gene = gene_stats.nlargest(top_n, 'Mean')['Group']
                top_groups_per_gene.extend(top_for_gene)
            
            top_groups = pd.Index(sorted(set(top_groups_per_gene)))
            region_data_filtered = region_data[region_data['Group'].isin(top_groups)]
            
            # Create violin plot
            sns.violinplot(data=region_data_filtered,
                          x='Group',
                          y='Expression',
                          hue='Gene',
                          order=top_groups,
                          palette=['#2ecc71', '#e74c3c'],
                          cut=0,  # Show the full range
                          scale='width',  # Scale each violin to same width
                          inner='box',  # Show box plot inside violin
                          ax=ax)
            
            # Customize the plot
            ax.set_title(f'Region: {region}', fontsize=12, pad=20)
            ax.set_xlabel('Cell Type', fontsize=10)
            ax.set_ylabel('Expression Level', fontsize=10)
            
            # Rotate x-axis labels
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
            
            # Add expression statistics above each violin
            for i, group in enumerate(top_groups):
                for j, gene in enumerate(genes):
                    stats_row = stats_df[
                        (stats_df['Group'] == group) & 
                        (stats_df['Gene'] == gene)
                    ].iloc[0]
                    
                    # Position text above each violin
                    y_pos = region_data_filtered[
                        (region_data_filtered['Group'] == group) & 
                        (region_data_filtered['Gene'] == gene)
                    ]['Expression'].max()
                    
                    # Add mean and % expressing cells
                    ax.text(i, y_pos,
                           f'μ={stats_row["Mean"]:.1f}\n{stats_row["Pct_expressing"]:.0f}%',
                           ha='center', va='bottom', fontsize=8)
            
            # Adjust legend
            ax.legend(title='Gene', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Remove empty subplots if any
    for idx in range(len(regions), len(axes)):
        fig.delaxes(axes[idx])
    
    plt.tight_layout()
    
    return fig, axes


