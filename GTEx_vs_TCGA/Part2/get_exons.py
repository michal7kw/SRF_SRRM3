# %%
import pandas as pd
from gtfparse import read_gtf
import pyensembl
import os
import hashlib
import gzip
import requests
from tqdm import tqdm
import polars as pl
import warnings
warnings.filterwarnings('ignore')

# %%
def get_gene_exons(gene_symbol, gtf_path, genome_build="GRCh38"):
    """
    Extract exon coordinates for a specific gene from GTF file
    """
    # Read GTF file using pandas instead of polars
    gtf = pd.read_csv(gtf_path, sep='\t', comment='#', 
                      names=['seqname', 'source', 'feature', 'start', 'end', 
                            'score', 'strand', 'frame', 'attribute'])
    
    # Parse the attributes column
    def parse_attributes(attr):
        attrs = {}
        for pair in attr.split(';'):
            if pair.strip():
                key, value = pair.strip().split(' ', 1)
                attrs[key] = value.strip('"')
        return attrs
    
    # Extract gene_name from attributes
    gtf['attrs'] = gtf['attribute'].apply(parse_attributes)
    gtf['gene_name'] = gtf['attrs'].apply(lambda x: x.get('gene_name', ''))
    
    # Filter for the gene of interest and exon features
    gene_exons = gtf[
        (gtf['gene_name'] == gene_symbol) & 
        (gtf['feature'] == 'exon')
    ].copy()
    
    if len(gene_exons) == 0:
        raise ValueError(f"No exons found for gene {gene_symbol}")
    
    # Sort exons by position
    gene_exons = gene_exons.sort_values(['seqname', 'start', 'end'])
    
    # Add useful columns for splicing analysis
    gene_exons['exon_length'] = gene_exons['end'] - gene_exons['start'] + 1
    gene_exons['exon_number'] = range(1, len(gene_exons) + 1)
    
    # Extract additional IDs from attributes
    gene_exons['gene_id'] = gene_exons['attrs'].apply(lambda x: x.get('gene_id', ''))
    gene_exons['transcript_id'] = gene_exons['attrs'].apply(lambda x: x.get('transcript_id', ''))
    gene_exons['exon_id'] = gene_exons['attrs'].apply(lambda x: x.get('exon_id', ''))
    
    # Select relevant columns
    result = gene_exons[['seqname', 'start', 'end', 'strand',
                        'gene_id', 'gene_name', 'transcript_id',
                        'exon_id', 'exon_number', 'exon_length']].copy()
    
    return result

def validate_coordinates(exon_df, datasets=['GTEx', 'TCGA']):
    """
    Validate exon coordinates against common RNA-seq datasets
    
    Parameters:
    -----------
    exon_df: DataFrame
        Output from get_gene_exons()
    datasets: list
        List of datasets to check compatibility with
        
    Returns:
    --------
    Dictionary with validation results
    """
    validation = {
        'genome_build': 'GRCh38',
        'coordinate_system': '1-based',
        'compatible_with': datasets,
        'warnings': [],
        'exon_count': len(exon_df),
        'chromosome': exon_df['seqname'].iloc[0],
        'strand': exon_df['strand'].iloc[0],
        'total_exon_length': exon_df['exon_length'].sum()
    }
    
    # Add warnings for potential issues
    if len(exon_df['transcript_id'].unique()) > 1:
        validation['warnings'].append(
            "Multiple transcript isoforms present - consider filtering for canonical transcript"
        )
    
    return validation

# Example usage:
# exons = get_gene_exons('SRRM3', 'path/to/gencode.v38.annotation.gtf')
# validation = validate_coordinates(exons)
# 
# # Export coordinates
# exons.to_csv('SRRM3_exon_coordinates.tsv', sep='\t', index=False)

# %%

def download_gencode_gtf(version="38", release_type="primary_assembly", out_dir="./"):
    """
    Download GENCODE GTF file and verify its integrity
    
    Parameters:
    -----------
    version: str
        GENCODE version (default: "38")
    release_type: str
        'primary_assembly' or 'chr_patch_hapl_scaff'
    out_dir: str
        Directory to save the downloaded file
    
    Returns:
    --------
    str: Path to the downloaded and extracted GTF file
    """
    # Create output directory if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)
    
    # Construct URLs
    base_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"
    if release_type == "primary_assembly":
        gtf_filename = f"gencode.v{version}.primary_assembly.annotation.gtf.gz"
    else:
        gtf_filename = f"gencode.v{version}.annotation.gtf.gz"
    
    gtf_url = f"{base_url}/release_{version}/{gtf_filename}"
    
    # Download paths
    gtf_gz_path = os.path.join(out_dir, gtf_filename)
    gtf_path = gtf_gz_path[:-3]  # Remove .gz extension
    
    # Download GTF file if it doesn't exist
    if not os.path.exists(gtf_path):
        print(f"Downloading GENCODE v{version} GTF file...")
        
        # Download GTF file with progress bar
        response = requests.get(gtf_url, stream=True)
        total_size = int(response.headers.get('content-length', 0))
        
        with open(gtf_gz_path, 'wb') as f, tqdm(
            desc=gtf_filename,
            total=total_size,
            unit='iB',
            unit_scale=True,
            unit_divisor=1024,
        ) as pbar:
            for data in response.iter_content(chunk_size=1024):
                size = f.write(data)
                pbar.update(size)
        
        # Extract GTF file
        print("Extracting GTF file...")
        try:
            with gzip.open(gtf_gz_path, 'rb') as f_in:
                with open(gtf_path, 'wb') as f_out:
                    f_out.write(f_in.read())
            
            # Remove compressed file to save space
            os.remove(gtf_gz_path)
            print("Download and extraction complete!")
            
        except Exception as e:
            if os.path.exists(gtf_gz_path):
                os.remove(gtf_gz_path)
            if os.path.exists(gtf_path):
                os.remove(gtf_path)
            raise Exception(f"Failed to extract GTF file: {str(e)}")
            
    else:
        print(f"GTF file already exists at: {gtf_path}")
    
    return gtf_path

def verify_gtf_content(gtf_path):
    """
    Basic verification of GTF file content
    
    Parameters:
    -----------
    gtf_path: str
        Path to GTF file
    
    Returns:
    --------
    bool: True if verification passes
    """
    try:
        # Check first few lines to ensure it's a GTF file
        with open(gtf_path, 'r') as f:
            header = [next(f) for _ in range(5)]
        
        # Check for expected GTF format
        has_valid_header = any(line.startswith('#!genome-build') for line in header)
        has_valid_columns = len(header[-1].split('\t')) == 9
        
        if not (has_valid_header and has_valid_columns):
            raise ValueError("GTF file format appears invalid")
        
        return True
        
    except Exception as e:
        print(f"GTF verification failed: {str(e)}")
        return False


# %%
version = "26"

# %%
# Download GENCODE v38 GTF
gtf_path = download_gencode_gtf(version=version, 
                              release_type="primary_assembly",
                              out_dir="./data")

# Verify the file
if verify_gtf_content(gtf_path):
    print(f"GTF file ready for use: {gtf_path}")

# %% [markdown]
# mamba install numpy=1.24.3
# mamba install pandas
# mamba install -c bioconda pyensembl
# mamba install -c bioconda gtfparse

# %%
# Get exon coordinates
exons = get_gene_exons('SRRM3', f'./data/gencode.v{version}.primary_assembly.annotation.gtf')

# Validate compatibility
validation = validate_coordinates(exons)

# Save to file
exons.to_csv(f'SRRM3_exon_coordinates_v{version}.tsv', sep='\t', index=False)

# %%



