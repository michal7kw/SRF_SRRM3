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
import scanpy as sc
import h5py
from scipy import sparse

# Project imports
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

from typing import List, Union, Tuple, Optional

import os
# Set the current working directory
os.chdir('/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3')

from tqdm.auto import tqdm

def load_multi_region_data(cell_metadata: pd.DataFrame,
        abc_cache: AbcProjectCache,
        datasets: List[str], 
        genes: List[str],
        dataset_regions: Optional[List[str]] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load expression data for multiple datasets
    """
    print(f"\nLoading data for {len(genes)} genes across {len(datasets)} datasets")
    
    all_expression = []
    all_metadata = []
    
    # Get all dataset versions from metadata
    dataset_versions = cell_metadata['dataset_label'].unique()
    
    for dataset in datasets:
        # Get all versions for this dataset
        versions = [v for v in dataset_versions if v.startswith(dataset)]
        print(f"\nProcessing dataset: {dataset}")
        print(f"Found versions: {versions}")
        
        for version in versions:
            # Get feature matrices for this version
            feature_matrices = cell_metadata[
                cell_metadata['dataset_label'] == version
            ]['feature_matrix_label'].unique()
            
            for feature_matrix in feature_matrices:
                print(f"\nProcessing matrix: {feature_matrix}")
                
                # Get cells for this feature matrix
                matrix_cells = cell_metadata[
                    cell_metadata['feature_matrix_label'] == feature_matrix
                ]
                print(f"Initial cell count: {len(matrix_cells)}")
                
                # If specific regions requested, filter cells
                if dataset_regions is not None:
                    matrix_cells = matrix_cells[
                        matrix_cells['region_of_interest_acronym'].isin(dataset_regions)
                    ]
                    print(f"Cells after region filter: {len(matrix_cells)}")
                    
                    if matrix_cells.empty:
                        print(f"No cells found for specified regions in {feature_matrix}")
                        continue
                
                try:
                    file_path = abc_cache.get_data_path(
                        directory=version,
                        file_name=f"{feature_matrix}/log2"
                    )
                    print(f"Loading file: {file_path}")
                    
                    with h5py.File(file_path, 'r') as f:
                        # Load and process genes
                        gene_symbols = f['var']['gene_symbol'][:]
                        if isinstance(gene_symbols[0], bytes):
                            gene_symbols = [g.decode('utf-8') for g in gene_symbols]
                        
                        # Find matching genes
                        gene_indices = [i for i, g in enumerate(gene_symbols) if g in genes]
                        if not gene_indices:
                            print(f"No requested genes found in {feature_matrix}")
                            continue
                        
                        found_genes = [gene_symbols[i] for i in gene_indices]
                        print(f"Found genes: {found_genes}")
                        
                        # Get cell barcodes from metadata
                        matrix_cells_set = set(matrix_cells['cell_barcode'].values)
                        
                        # Load cell barcodes from file
                        cell_barcode_cats = f['obs']['cell_barcode']['categories'][:]
                        if isinstance(cell_barcode_cats[0], bytes):
                            cell_barcode_cats = [b.decode('utf-8') for b in cell_barcode_cats]
                        cell_barcode_codes = f['obs']['cell_barcode']['codes'][:]
                        
                        # Match cells
                        cell_indices = []
                        valid_barcodes = []
                        for idx, code in enumerate(cell_barcode_codes):
                            barcode = cell_barcode_cats[code]
                            if barcode in matrix_cells_set:
                                cell_indices.append(idx)
                                valid_barcodes.append(barcode)
                        
                        if not cell_indices:
                            print(f"No matching cells found in {feature_matrix}")
                            continue
                        
                        print(f"Found {len(cell_indices)} matching cells")
                        
                        # Load expression data
                        data = f['X/data'][:]
                        indices = f['X/indices'][:]
                        indptr = f['X/indptr'][:]
                        
                        # Create sparse matrix and extract relevant data
                        X = sparse.csr_matrix((data, indices, indptr))
                        X = X[cell_indices][:, gene_indices]
                        
                        # Convert to dense array for selected genes only
                        data = pd.DataFrame(
                            X.toarray(),
                            index=valid_barcodes,
                            columns=found_genes
                        )
                        
                        # Get corresponding metadata
                        meta = matrix_cells[
                            matrix_cells['cell_barcode'].isin(valid_barcodes)
                        ].copy()
                        
                        # Add region and dataset information
                        data['region'] = meta['region_of_interest_acronym']
                        data['dataset'] = version
                        
                        all_expression.append(data)
                        all_metadata.append(meta)
                        
                        print(f"Successfully processed {len(cell_indices)} cells")
                    
                except Exception as e:
                    print(f"Error processing {feature_matrix}: {str(e)}")
                    import traceback
                    traceback.print_exc()
                    continue
    
    if not all_expression:
        raise ValueError("No data could be loaded. Please check your input parameters.")
    
    # Combine all data
    final_expression = pd.concat(all_expression)
    final_metadata = pd.concat(all_metadata)
    
    # Ensure all requested genes are present
    for gene in genes:
        if gene not in final_expression.columns:
            final_expression[gene] = np.nan
    
    # Final cleanup
    final_expression = final_expression[genes]  # Reorder columns to match requested genes
    
    print(f"\nSuccessfully loaded data for {len(final_metadata)} cells")
    return final_expression, final_metadata

def load_multi_region_data_2(cell_metadata: pd.DataFrame,
                          abc_cache: AbcProjectCache,
                          datasets: List[str],
                          genes: List[str],
                          dataset_regions: Optional[List[str]] = None,
                          barcode_column: str = 'cell_barcode') -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load expression data for multiple datasets and merge with cell metadata.
    """
    print(f"\nLoading data for {len(genes)} genes across {len(datasets)} datasets")
    
    # Initialize empty lists to store data
    all_expression_data = []
    processed_metadata = []
    
    # Create progress bars
    dataset_pbar = tqdm(datasets, desc="Processing datasets", position=0)
    
    for dataset in dataset_pbar:
        matrices = cell_metadata[cell_metadata['dataset_label'] == dataset]['feature_matrix_label'].unique()
        dataset_pbar.write(f"\nFound {len(matrices)} matrices in {dataset}:")
        dataset_pbar.write(str(matrices))
        
        matrix_pbar = tqdm(matrices, desc=f"Processing {dataset} matrices", 
                          position=1, leave=False)
        
        for matrix in matrix_pbar:
            if dataset_regions and not any(r in matrix for r in dataset_regions):
                matrix_pbar.write(f"Skipping {matrix} - no matching regions")
                continue
                
            try:
                matrix_pbar.set_description(f"{matrix}: Loading file")
                file_path = abc_cache.get_data_path(directory=dataset, file_name=f"{matrix}/log2")
                
                # Pre-filter metadata for this matrix
                matrix_cells = cell_metadata[
                    cell_metadata['feature_matrix_label'] == matrix
                ][barcode_column].values
                
                if len(matrix_cells) == 0:
                    matrix_pbar.write(f"No cells found for {matrix}")
                    continue
                
                # Convert matrix_cells to set for faster lookup
                matrix_cells_set = set(matrix_cells)
                
                with h5py.File(file_path, 'r') as f:
                    # Load and process genes
                    gene_symbols = f['var']['gene_symbol'][:]
                    if isinstance(gene_symbols[0], bytes):
                        gene_symbols = [g.decode('utf-8') for g in gene_symbols]
                    
                    # Find matching genes
                    gene_indices = [i for i, g in enumerate(gene_symbols) if g in genes]
                    if not gene_indices:
                        matrix_pbar.write(f"No requested genes found in {matrix}")
                        continue
                    
                    found_genes = [gene_symbols[i] for i in gene_indices]
                    matrix_pbar.write(f"Found genes: {found_genes}")
                    
                    # Load cell barcodes more efficiently
                    cell_barcode_cats = f['obs']['cell_barcode']['categories'][:]
                    if isinstance(cell_barcode_cats[0], bytes):
                        cell_barcode_cats = [b.decode('utf-8') for b in cell_barcode_cats]
                    cell_barcode_codes = f['obs']['cell_barcode']['codes'][:]
                    
                    # Create mapping dictionary for faster lookup
                    matrix_pbar.set_description(f"{matrix}: Creating barcode mapping")
                    cell_indices = []
                    valid_barcodes = []
                    for idx, code in enumerate(cell_barcode_codes):
                        barcode = cell_barcode_cats[code]
                        if barcode in matrix_cells_set:
                            cell_indices.append(idx)
                            valid_barcodes.append(barcode)
                    
                    if not cell_indices:
                        matrix_pbar.write("No matching cells found!")
                        continue
                    
                    # Load expression data
                    matrix_pbar.set_description(f"{matrix}: Loading expression data")
                    data = f['X/data'][:]
                    indices = f['X/indices'][:]
                    indptr = f['X/indptr'][:]
                    
                    # Create sparse matrix and extract relevant data
                    X = sparse.csr_matrix((data, indices, indptr))
                    X = X[cell_indices][:, gene_indices]
                    
                    # Convert to dense array for selected genes only
                    matrix_exp = pd.DataFrame(
                        X.toarray(),
                        index=valid_barcodes,
                        columns=found_genes
                    )
                    
                    # Get corresponding metadata
                    matrix_meta = cell_metadata[
                        cell_metadata[barcode_column].isin(valid_barcodes)
                    ].copy()
                    
                    # Append to lists
                    all_expression_data.append(matrix_exp)
                    processed_metadata.append(matrix_meta)
                    
                    matrix_pbar.write(f"Successfully processed {len(cell_indices)} cells")
                
            except Exception as e:
                matrix_pbar.write(f"Error processing {matrix}: {str(e)}")
                import traceback
                traceback.print_exc()
                continue
            
        matrix_pbar.close()
    
    dataset_pbar.close()
    
    if not processed_metadata:
        raise ValueError("No data could be loaded. Please check your input parameters.")
    
    # Combine all data
    final_expression = pd.concat(all_expression_data)
    final_metadata = pd.concat(processed_metadata)
    
    # Ensure all requested genes are present
    for gene in genes:
        if gene not in final_expression.columns:
            final_expression[gene] = np.nan
    
    # Final cleanup
    final_expression = final_expression[genes]  # Reorder columns to match requested genes
    
    print(f"\nSuccessfully loaded data for {len(final_metadata)} cells")
    return final_expression, final_metadata


def plot_region_gene_comparison(expression_data: pd.DataFrame, 
                              metadata: pd.DataFrame,
                              genes: List[str],
                              group_by: str,
                              region_column: str = 'region',
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
        'Region': pd.concat([expression_data[region_column] for _ in genes])
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


def plot_gene_expression_by_group(expression_data, metadata, group_by='class', min_expr_threshold=0.1):
    """Plot gene expression analysis by group"""
    
    # Prepare data
    plot_data = pd.DataFrame({
        'Expression': pd.concat([expression_data[gene] for gene in genes]),
        'Gene': np.repeat(genes, len(expression_data)),
        'Group': pd.concat([metadata[group_by] for _ in genes]),
        'Region': pd.concat([expression_data['region'] for _ in genes])
    })
    
    # Calculate statistics by group
    stats = []
    for region in plot_data['Region'].unique():
        region_data = plot_data[plot_data['Region'] == region]
        for group in region_data['Group'].unique():
            for gene in genes:
                group_gene_data = region_data[
                    (region_data['Group'] == group) & 
                    (region_data['Gene'] == gene)
                ]
                
                if len(group_gene_data) > 0:
                    stats.append({
                        'Region': region,
                        'Group': group,
                        'Gene': gene,
                        'Mean': group_gene_data['Expression'].mean(),
                        'Median': group_gene_data['Expression'].median(),
                        'Std': group_gene_data['Expression'].std(),
                        'Pct_expressing': (group_gene_data['Expression'] > min_expr_threshold).mean() * 100,
                        'N_cells': len(group_gene_data)
                    })
    
    stats_df = pd.DataFrame(stats)
    
    # Plot violin plots by region
    regions = sorted(plot_data['Region'].unique())
    n_regions = len(regions)
    fig, axes = plt.subplots(2, (n_regions+1)//2, figsize=(20, 12))
    axes = axes.flatten()
    
    for idx, region in enumerate(regions):
        region_stats = stats_df[stats_df['Region'] == region]
        
        # Get top groups by mean expression
        top_groups = []
        for gene in genes:
            gene_stats = region_stats[region_stats['Gene'] == gene]
            top = gene_stats.nlargest(5, 'Mean')['Group'].tolist()
            top_groups.extend(top)
        top_groups = list(set(top_groups))
        
        region_data = plot_data[plot_data['Region'] == region]
        region_data = region_data[region_data['Group'].isin(top_groups)]
        
        sns.violinplot(data=region_data,
                      x='Group',
                      y='Expression',
                      hue='Gene',
                      ax=axes[idx],
                      cut=0,
                      density_norm='width',
                      inner='box')
        
        axes[idx].set_title(f'Region: {region}')
        # First set the ticks
        axes[idx].set_xticks(range(len(top_groups)))
        # Then set the labels
        axes[idx].set_xticklabels(top_groups, rotation=45, ha='right')
        
        # Add statistics
        for i, group in enumerate(top_groups):
            for j, gene in enumerate(genes):
                stats_row = region_stats[
                    (region_stats['Group'] == group) & 
                    (region_stats['Gene'] == gene)
                ]
                if not stats_row.empty:
                    stats_row = stats_row.iloc[0]
                    y_pos = region_data[
                        (region_data['Group'] == group) & 
                        (region_data['Gene'] == gene)
                    ]['Expression'].max()
                    if not pd.isna(y_pos):
                        axes[idx].text(i, y_pos,
                                     f'μ={stats_row["Mean"]:.1f}\n{stats_row["Pct_expressing"]:.0f}%',
                                     ha='center', va='bottom', fontsize=8)
    
    # Remove empty subplots
    for idx in range(n_regions, len(axes)):
        fig.delaxes(axes[idx])
    
    plt.suptitle(f'Gene Expression by Region and {group_by}', y=1.02)
    plt.tight_layout()
    return fig, axes, stats_df


def plot_combined_gene_expression(combined_data: pd.DataFrame, 
                                  group_by: str = 'class', 
                                  genes: List[str] = None, 
                                  min_expr_threshold: float = 0.1, 
                                  figsize: Tuple[int, int] = (20, 12)):
    """Plot gene expression analysis for combined data by group"""
    
    if genes is None:
        genes = combined_data.columns.intersection(['Srrm3', 'Srrm4']).tolist()
    
    # Prepare data
    plot_data = pd.DataFrame({
        'Expression': pd.concat([combined_data[gene] for gene in genes]),
        'Gene': np.repeat(genes, len(combined_data)),
        'Group': pd.concat([combined_data[group_by] for _ in genes]),
        'Region': pd.concat([combined_data['region'] for _ in genes])
    })
    
    # Calculate statistics by group
    stats = []
    for region in plot_data['Region'].unique():
        region_data = plot_data[plot_data['Region'] == region]
        for group in region_data['Group'].unique():
            for gene in genes:
                group_gene_data = region_data[
                    (region_data['Group'] == group) & 
                    (region_data['Gene'] == gene)
                ]
                
                if len(group_gene_data) > 0:
                    stats.append({
                        'Region': region,
                        'Group': group,
                        'Gene': gene,
                        'Mean': group_gene_data['Expression'].mean(),
                        'Median': group_gene_data['Expression'].median(),
                        'Std': group_gene_data['Expression'].std(),
                        'Pct_expressing': (group_gene_data['Expression'] > min_expr_threshold).mean() * 100,
                        'N_cells': len(group_gene_data)
                    })
    
    stats_df = pd.DataFrame(stats)
    
    # Plot violin plots by region
    regions = sorted(plot_data['Region'].unique())
    n_regions = len(regions)
    fig, axes = plt.subplots(2, (n_regions+1)//2, figsize=figsize)
    axes = axes.flatten()
    
    for idx, region in enumerate(regions):
        region_stats = stats_df[stats_df['Region'] == region]
        
        # Filter out groups with no data
        valid_groups = []
        for group in region_data['Group'].unique():
            if all(len(region_data[(region_data['Group'] == group) & 
                                 (region_data['Gene'] == gene)]) > 0 
                   for gene in genes):
                valid_groups.append(group)
        
        # Filter data to only include groups with valid data
        region_data = region_data[region_data['Group'].isin(valid_groups)]
        
        # Create violin plot only for valid groups
        violin = sns.violinplot(data=region_data,
                              x='Group',
                              y='Expression',
                              hue='Gene',
                              ax=axes[idx],
                              cut=0,
                              density_norm='width',
                              inner='box')
        
        # Set consistent y-axis limits
        axes[idx].set_ylim(plot_data['Expression'].min(), plot_data['Expression'].max() * 1.2)
        
        axes[idx].set_title(f'Region: {region}', fontsize=12, pad=20)
        axes[idx].set_xlabel('Cell Type', fontsize=10)
        axes[idx].set_ylabel('Expression Level (log2(CPM + 1))', fontsize=10)
        
        # Set x-tick labels only for valid groups
        axes[idx].set_xticks(range(len(valid_groups)))
        axes[idx].set_xticklabels(valid_groups, rotation=45, ha='right')
        
        # Prepare legend labels only for valid groups
        legend_labels = []
        legend_handles = []
        
        for group in valid_groups:
            # Create a header for each group
            legend_labels.append(f"\n{group}")
            legend_handles.append(plt.Line2D([0], [0], color='none'))
            
            for gene in genes:
                stats_row = region_stats[
                    (region_stats['Group'] == group) & 
                    (region_stats['Gene'] == gene)
                ]
                if not stats_row.empty:
                    stats_row = stats_row.iloc[0]
                    if stats_row['N_cells'] > 0:  # Only add to legend if there are cells
                        legend_labels.append(
                            f"{gene}: μ={stats_row['Mean']:.1f} ({stats_row['Pct_expressing']:.0f}%, n={stats_row['N_cells']})"
                        )
                        color = sns.color_palette()[genes.index(gene)]
                        legend_handles.append(plt.Line2D([0], [0], color=color, lw=2))
        
        # Remove existing legend
        if axes[idx].get_legend() is not None:
            axes[idx].get_legend().remove()
        
        # Add custom formatted legend
        axes[idx].legend(legend_handles, 
                        legend_labels,
                        title='Cell Type Statistics',
                        bbox_to_anchor=(1.05, 1),
                        loc='upper left',
                        frameon=True,
                        fontsize=9,
                        title_fontsize=10,
                        handlelength=1.5)
        
    # Remove empty subplots
    for idx in range(n_regions, len(axes)):
        fig.delaxes(axes[idx])
    
    plt.suptitle(f'Gene Expression by Region and {group_by}\n'
                 f'Expression values in log2(CPM + 1) scale, showing mean (μ), % expressing cells, and cell count (n)',
                 y=1.02)
    plt.tight_layout()
    return fig, axes, stats_df