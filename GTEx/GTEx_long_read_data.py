# %% [markdown]
# # Environment

# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple
import gzip
import re

# %%
# Set the current working directory
os.chdir('/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3')

# Print the current working directory to confirm the change
print(f"Current working directory: {os.getcwd()}")


# %% [markdown]
# # Load and explore GTExdata

# %%
gencode_tpm_file = "./DATA/GTEx/quantification_gencode.tpm.txt.gz"
flair_tpm_file = "./DATA/GTEx/quantification_flair_filter.tpm.txt.gz"

# %%
# Read TPM data
print("Reading TPM data...")
gencode_tpm_data = pd.read_csv(gencode_tpm_file, sep='\t', compression='gzip')
flair_tpm_data = pd.read_csv(flair_tpm_file, sep='\t', compression='gzip')

# %%
gencode_tpm_data.head()

# %%
transcripts_of_interest = ['ENST00000611745.1', 'ENST00000612155.1']
gencode_tpm_data.loc[gencode_tpm_data['transcript'].isin(transcripts_of_interest)]

# %%
flair_tpm_data.head()

# %%
# Create histograms of gene expression levels
plt.figure(figsize=(12, 6))

# Plot histogram for GENCODE TPM data
plt.subplot(1, 2, 1)
gencode_expr = gencode_tpm_data.iloc[:, 1:].values.flatten()  # Exclude transcript column
gencode_expr = gencode_expr[gencode_expr > 0]  # Filter out zeros
plt.hist(np.log10(gencode_expr), bins=50, alpha=0.7)
plt.xlabel('log10(TPM)')
plt.ylabel('Frequency')
plt.title('GENCODE Expression Distribution')
plt.grid(True, alpha=0.3)

# Plot histogram for FLAIR TPM data  
plt.subplot(1, 2, 2)
flair_expr = flair_tpm_data.iloc[:, 1:].values.flatten()  # Exclude transcript column
flair_expr = flair_expr[flair_expr > 0]  # Filter out zeros
plt.hist(np.log10(flair_expr), bins=50, alpha=0.7)
plt.xlabel('log10(TPM)')
plt.ylabel('Frequency') 
plt.title('FLAIR Expression Distribution')
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# Print summary statistics
print("\nSummary Statistics:")
print("\nGENCODE TPM:")
print(f"Mean: {np.mean(gencode_expr):.2f}")
print(f"Median: {np.median(gencode_expr):.2f}")
print(f"Std: {np.std(gencode_expr):.2f}")

print("\nFLAIR TPM:") 
print(f"Mean: {np.mean(flair_expr):.2f}")
print(f"Median: {np.median(flair_expr):.2f}")
print(f"Std: {np.std(flair_expr):.2f}")


# %%
print(gencode_tpm_data.shape)
print(flair_tpm_data.shape)

# %%
flair_filter_transcripts_file = "./DATA/GTEx/flair_filter_transcripts.gtf.gz"
flair_filter_transcripts_data = pd.read_csv(flair_filter_transcripts_file, sep='\t', compression='gzip')
flair_filter_transcripts_data.head()

# %% [markdown]
# # Search for SRRM3 transcripts in FLAIR data

# %%
class FlairSRRM3Analyzer:
    def __init__(self, flair_gtf: str, flair_tpm: str):
        """
        Initialize analyzer with FLAIR GTF and TPM file paths
        """
        self.flair_gtf = flair_gtf
        self.flair_tpm = flair_tpm
        self.srrm3_region = "chr7"
        # Exon 15 characteristics
        self.exon15_length = 79  # c.1783-1785 suggests a 3bp exon
        self.exon15_phase = (2, 0)  # start phase 2, end phase 0
        
    def get_exon_phases(self, exons: List[dict]) -> List[Tuple[int, int]]:
        """Calculate phase for each exon"""
        phases = []
        current_phase = 0
        
        for exon in exons:
            length = exon['length']
            end_phase = (current_phase + length) % 3
            phases.append((current_phase, end_phase))
            current_phase = end_phase
            
        return phases
    
    def has_exon15_characteristics(self, exon: dict, prev_phase: int) -> bool:
        """
        Check if an exon matches exon 15 characteristics
        - Should be relatively small (around 79 bp)
        - Should have specific phase pattern
        """
        length = exon['length']
        # Allow some flexibility in length
        if abs(length - self.exon15_length) > 10:
            return False
            
        # Calculate phases
        end_phase = (prev_phase + length) % 3
        return (prev_phase, end_phase) == self.exon15_phase
    
    def parse_flair_gtf(self) -> Dict[str, dict]:
        """
        Parse FLAIR GTF file to find SRRM3 transcripts
        """
        print("Parsing FLAIR GTF file...")
        current_transcript = None
        current_exons = []
        transcript_info = {}
        
        opener = gzip.open if self.flair_gtf.endswith('.gz') else open
        
        with opener(self.flair_gtf, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                chrom, source, feature_type, start, end, score, strand, frame, attributes = fields
                start, end = int(start), int(end)
                
                # Parse attributes
                attr_dict = {}
                for attr in attributes.strip(';').split('; '):
                    key, value = attr.split(' ', 1)
                    attr_dict[key] = value.strip('"')
                
                if chrom != self.srrm3_region:
                    continue
                
                if feature_type == 'transcript':
                    # Process previous transcript
                    if current_transcript and current_exons:
                        # Calculate phases for all exons
                        phases = self.get_exon_phases(current_exons)
                        
                        # Check for exon 15-like characteristics
                        has_exon15 = False
                        exon15_index = None
                        for i in range(1, len(current_exons)-1):  # Skip first and last exon
                            if self.has_exon15_characteristics(current_exons[i], phases[i][0]):
                                has_exon15 = True
                                exon15_index = i
                                break
                        
                        transcript_info[current_transcript] = {
                            'exons': current_exons.copy(),
                            'length': sum(e['length'] for e in current_exons),
                            'strand': strand,
                            'gene_id': attr_dict['gene_id'],
                            'has_exon15': has_exon15,
                            'exon15_index': exon15_index
                        }
                    
                    # Start new transcript
                    current_transcript = attr_dict['transcript_id']
                    current_exons = []
                    
                elif feature_type == 'exon' and current_transcript:
                    exon_length = end - start + 1
                    current_exons.append({
                        'start': start,
                        'end': end,
                        'length': exon_length
                    })
        
        return transcript_info

# %%
def classify_transcripts(transcript_info: Dict[str, dict], 
                        min_length: int = 3600, 
                        max_length: int = 3800) -> Tuple[List[str], List[str]]:
    """
    Classify transcripts with strict criteria:
    - Length between 3600-3800 bp
    - Without exon 15: must have exactly 15 exons
    - With exon 15: must have exactly 16 exons
    
    Parameters:
    transcript_info: Dictionary containing transcript information
    min_length: Minimum acceptable transcript length
    max_length: Maximum acceptable transcript length
    
    Returns:
    Tuple of (transcripts_with_exon15, transcripts_without_exon15)
    """
    with_exon15 = []
    without_exon15 = []
    
    # Tracking for detailed filtering statistics
    stats = {
        'total': len(transcript_info),
        'length_filtered': 0,
        'with_exon15': {
            'total': 0,
            'correct_exon_count': 0
        },
        'without_exon15': {
            'total': 0,
            'correct_exon_count': 0
        }
    }
    
    for transcript_id, info in transcript_info.items():
        # First check length
        if not (min_length <= info['length'] <= max_length):
            continue
        
        stats['length_filtered'] += 1
        exon_count = len(info['exons'])
        
        if info['has_exon15']:
            stats['with_exon15']['total'] += 1
            if exon_count == 16:  # Must have exactly 16 exons
                stats['with_exon15']['correct_exon_count'] += 1
                with_exon15.append({
                    'transcript_id': transcript_id,
                    'length': info['length'],
                    'exon_count': exon_count,
                    'exon15_index': info['exon15_index']
                })
        else:
            stats['without_exon15']['total'] += 1
            if exon_count == 15:  # Must have exactly 15 exons
                stats['without_exon15']['correct_exon_count'] += 1
                without_exon15.append({
                    'transcript_id': transcript_id,
                    'length': info['length'],
                    'exon_count': exon_count
                })
    
    # Print detailed filtering report
    print("\nTranscript Filtering Report")
    print("=" * 50)
    print(f"Total transcripts analyzed: {stats['total']}")
    print(f"Transcripts within length range ({min_length}-{max_length} bp): {stats['length_filtered']}")
    print("\nIsoform without exon 15 (should have 15 exons):")
    print(f"- Total candidates: {stats['without_exon15']['total']}")
    print(f"- Matching all criteria: {stats['without_exon15']['correct_exon_count']}")
    print("\nIsoform with exon 15 (should have 16 exons):")
    print(f"- Total candidates: {stats['with_exon15']['total']}")
    print(f"- Matching all criteria: {stats['with_exon15']['correct_exon_count']}")
    
    # Print detailed information about matching transcripts
    print("\nMatching Transcripts Details")
    print("=" * 50)
    
    print("\nTranscripts WITH exon 15 (16 exons):")
    print("-" * 40)
    for transcript in with_exon15:
        print(f"Transcript ID: {transcript['transcript_id']}")
        print(f"Length: {transcript['length']} bp")
        print(f"Exon count: {transcript['exon_count']}")
        print(f"Exon 15 index: {transcript['exon15_index']}")
        print("-" * 40)
    
    print("\nTranscripts WITHOUT exon 15 (15 exons):")
    print("-" * 40)
    for transcript in without_exon15:
        print(f"Transcript ID: {transcript['transcript_id']}")
        print(f"Length: {transcript['length']} bp")
        print(f"Exon count: {transcript['exon_count']}")
        print("-" * 40)
    
    # Return just the transcript IDs in the final tuple
    return (
        [t['transcript_id'] for t in with_exon15],
        [t['transcript_id'] for t in without_exon15]
    )

# %%
# File paths
flair_gtf = "./DATA/GTEx/flair_filter_transcripts.gtf.gz"
flair_tpm = "./DATA/GTEx/quantification_flair_filter.tpm.txt.gz"

# Initialize analyzer
analyzer = FlairSRRM3Analyzer(flair_gtf, flair_tpm)

# Find and classify transcripts
transcript_info = analyzer.parse_flair_gtf()

# Classify transcripts with length filtering
transcript_sets = classify_transcripts(
    transcript_info,
    min_length=3600,
    max_length=3800
)

# %% [markdown]
# # Transcripts WITH exon 15 (16 exons)
# 
# ----------------------------------------
# - Transcript ID: ENST00000393651.7
# - Length: 3649 bp
# - Exon count: 16
# - Exon 15 index: 2
# ----------------------------------------
# 
# # Transcripts WITHOUT exon 15 (15 exons):
# 
# ----------------------------------------
# - Transcript ID: ENST00000462753.5
# - Length: 3719 bp
# - Exon count: 15
# ----------------------------------------

# %% [markdown]
# # Select transcripts of SRRM3 isoforms

# %%
transcript_sets_with_exon15 = transcript_sets[0]
transcript_sets_without_exon15 = transcript_sets[1]

# %%
print(transcript_sets_with_exon15)
print(transcript_sets_without_exon15)


# %%
with_exon15 = transcript_sets_with_exon15[0]
without_exon15 = transcript_sets_without_exon15[-1]
print(with_exon15)
print(without_exon15)

# %%
gencode_tpm_data.loc[gencode_tpm_data['transcript'].isin([with_exon15])]
gencode_tpm_data.loc[gencode_tpm_data['transcript'].isin([without_exon15])]

# %%
gencode_tpm_data.columns[1:10]

# %%
# Create a figure and axis with a larger size for better visibility
plt.figure(figsize=(10, 8))

# Get data for both transcripts
with_exon15_data = gencode_tpm_data.loc[gencode_tpm_data['transcript'].isin([with_exon15])].iloc[0]
without_exon15_data = gencode_tpm_data.loc[gencode_tpm_data['transcript'].isin([without_exon15])].iloc[0]

# Get expression values (excluding transcript column) and prepare data for boxplot
with_exon15_expr = with_exon15_data[1:].values
without_exon15_expr = without_exon15_data[1:].values
data = [with_exon15_expr, without_exon15_expr]

# Create boxplot
bp = plt.boxplot(data, patch_artist=True, labels=['With exon 15\n(16 exons)', 'Without exon 15\n(15 exons)'])

# Customize boxplot colors
colors = ['#3498db', '#e74c3c']
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

# Customize plot appearance
plt.ylabel('Expression (TPM)', fontsize=12)
plt.title('Expression Distribution of SRRM3 Isoforms\nAcross GTEx Samples', fontsize=14, pad=20)
plt.grid(True, linestyle='--', alpha=0.3)

# Add individual points with jitter for better visualization
for i, d in enumerate(data, 1):
    x = np.random.normal(i, 0.04, size=len(d))
    plt.scatter(x, d, alpha=0.4, s=20, color='black')

plt.tight_layout()
plt.show()

# Print summary statistics
print("\nSummary Statistics:")
print("\nWith exon 15:")
print(f"Mean: {np.mean(with_exon15_expr):.2f}")
print(f"Median: {np.median(with_exon15_expr):.2f}")
print(f"Std: {np.std(with_exon15_expr):.2f}")

print("\nWithout exon 15:")
print(f"Mean: {np.mean(without_exon15_expr):.2f}")
print(f"Median: {np.median(without_exon15_expr):.2f}")
print(f"Std: {np.std(without_exon15_expr):.2f}")

# %%
transcripts_of_interest

# %%
# Create a figure and axis with a larger size for better visibility
plt.figure(figsize=(10, 8))

# Get data for both transcripts
with_exon15_data = gencode_tpm_data.loc[gencode_tpm_data['transcript'] == "ENST00000611745.1"].iloc[0]
without_exon15_data = gencode_tpm_data.loc[gencode_tpm_data['transcript'] == "ENST00000612155.1"].iloc[0]

# Get expression values (excluding transcript column) and prepare data for boxplot
with_exon15_expr = with_exon15_data[1:].values
without_exon15_expr = without_exon15_data[1:].values
data = [with_exon15_expr, without_exon15_expr]

# Create boxplot
bp = plt.boxplot(data, patch_artist=True, labels=['With exon 15\n(16 exons)', 'Without exon 15\n(15 exons)'])

# Customize boxplot colors
colors = ['#3498db', '#e74c3c']
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

# Customize plot appearance
plt.ylabel('Expression (TPM)', fontsize=12)
plt.title('Expression Distribution of SRRM3 Isoforms\nAcross GTEx Samples', fontsize=14, pad=20)
plt.grid(True, linestyle='--', alpha=0.3)

# Add individual points with jitter for better visualization
for i, d in enumerate(data, 1):
    x = np.random.normal(i, 0.04, size=len(d))
    plt.scatter(x, d, alpha=0.4, s=20, color='black')

plt.tight_layout()
plt.show()

# Print summary statistics
print("\nSummary Statistics:")
print("\nWith exon 15:")
print(f"Mean: {np.mean(with_exon15_expr):.2f}")
print(f"Median: {np.median(with_exon15_expr):.2f}")
print(f"Std: {np.std(with_exon15_expr):.2f}")

print("\nWithout exon 15:")
print(f"Mean: {np.mean(without_exon15_expr):.2f}")
print(f"Median: {np.median(without_exon15_expr):.2f}")
print(f"Std: {np.std(without_exon15_expr):.2f}")

# %%
gencode_tpm_data.loc[gencode_tpm_data['transcript'].isin(transcripts_of_interest)]

# %% [markdown]
# # Annotate GTEx tissue types
# 

# %%
import pandas as pd

def load_gtex_sample_attributes(file_path):
    """
    Load GTEx sample attributes from a tab-separated file.
    
    Parameters:
    file_path (str): Path to the GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt file
    
    Returns:
    pandas.DataFrame: DataFrame containing the GTEx sample attributes
    """
    # Read the tab-separated file
    # The file has tab-separator and no quoted text
    df = pd.read_csv(file_path, sep='\t', quoting=3)
    
    # Clean up column names (optional)
    # Strip any whitespace from column names
    df.columns = df.columns.str.strip()
    
    # Convert numeric columns to appropriate types
    # List of columns that should be numeric
    numeric_columns = ['SMRIN', 'SMATSSCR', 'SMTSPAX', 'SMNABTCHD', 'SMGEBTCHD',
                      'SME2MPRT', 'SMCHMPRS', 'SMNUMGPS', 'SMMAPRT', 'SMEXNCRT',
                      'SM550NRM', 'SMGNSDTC', 'SMUNMPRT', 'SM350NRM', 'SMRDLGTH',
                      'SMMNCPB', 'SME1MMRT', 'SMSFLGTH', 'SMESTLBS', 'SMMPPD',
                      'SMNTERRT', 'SMRDTTL', 'SMMNCV', 'SMMPPDPR', 'SMCGLGTH',
                      'SMGAPPCT', 'SMUNPDRD', 'SMNTRNRT', 'SMMPUNRT', 'SMEXPEFF',
                      'SMMPPDUN', 'SME2MMRT', 'SME2ANTI', 'SMALTALG', 'SME2SNSE',
                      'SMMFLGTH', 'SME1ANTI', 'SMSPLTRD', 'SMBSMMRT', 'SME1SNSE',
                      'SME1PCTS', 'SMRRNART', 'SME1MPRT', 'SMNUM5CD', 'SMDPMPRT',
                      'SME2PCTS']
    
    # Convert to numeric, errors='coerce' will set invalid parsing to NaN
    for col in numeric_columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Basic data validation
    print(f"Loaded {len(df)} samples with {len(df.columns)} attributes")
    
    return df

def summarize_gtex_data(df):
    """
    Provide a basic summary of the GTEx sample attributes data.
    
    Parameters:
    df (pandas.DataFrame): DataFrame containing GTEx sample attributes
    
    Returns:
    dict: Summary statistics about the dataset
    """
    summary = {
        'total_samples': len(df),
        'unique_tissue_types': df['SMTS'].nunique(),
        'unique_detailed_tissues': df['SMTSD'].nunique(),
        'tissue_counts': df['SMTS'].value_counts().to_dict(),
        'missing_values_by_column': df.isnull().sum().to_dict()
    }
    
    return summary

# %%
file_path = "./DATA/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"

# Load the data
gtex_df = load_gtex_sample_attributes(file_path)

# %%
gtex_df.head()

# %%
gtex_df.SAMPID[:10]

# %%
# Get summary statistics
summary = summarize_gtex_data(gtex_df)

# Print some basic information
print("\nDataset Summary:")
print(f"Total samples: {summary['total_samples']}")
print(f"Unique tissue types: {summary['unique_tissue_types']}")
print(f"Unique detailed tissue types: {summary['unique_detailed_tissues']}")

print("\nSample counts by tissue type:")
for tissue, count in summary['tissue_counts'].items():
    print(f"{tissue}: {count}")

# %%
import pandas as pd
import re

def annotate_rnaseq_samples(rnaseq_sample_names, gtex_df):
    """
    Annotate RNA-seq sample names with tissue types from GTEx annotations using exact matching.
    
    Parameters:
    rnaseq_sample_names (list-like): List of RNA-seq sample names
    gtex_df (pandas.DataFrame): DataFrame containing GTEx sample annotations
    
    Returns:
    pandas.DataFrame: DataFrame with RNA-seq sample names and their tissue annotations
    """
    # Convert sample names to a Series for easier manipulation
    rnaseq_samples = pd.Series(rnaseq_sample_names)
    
    # Create a list to store annotations
    annotations = []
    
    def clean_sample_name(name):
        """Remove potential suffixes while preserving the core sample identifier."""
        return re.sub(r'_rep\d*|_ctrl|_exp\d*$', '', name)
    
    # Process each RNA-seq sample
    for sample in rnaseq_samples:
        # Clean the sample name
        cleaned_sample = clean_sample_name(sample)
        
        # Try to find an exact match (excluding potential suffixes)
        matches = gtex_df[gtex_df['SAMPID'].str.contains(f"^{cleaned_sample}$", regex=True, na=False)]
        
        if len(matches) == 1:
            # Unique exact match found
            match = matches.iloc[0]
            annotations.append({
                'rnaseq_sample': sample,
                'cleaned_sample': cleaned_sample,
                'tissue_type': match['SMTSD'],
                'broad_tissue_type': match['SMTS'],
                'matched_gtex_id': match['SAMPID'],
                'match_type': 'exact'
            })
        else:
            # Try matching without the SM part but including the R-number
            base_pattern = re.match(r'(GTEX-\w+-\d+-R\d+[ab])', cleaned_sample)
            if base_pattern:
                base_matches = gtex_df[gtex_df['SAMPID'].str.contains(base_pattern.group(1), regex=False, na=False)]
                if len(base_matches) == 1:
                    match = base_matches.iloc[0]
                    annotations.append({
                        'rnaseq_sample': sample,
                        'cleaned_sample': cleaned_sample,
                        'tissue_type': match['SMTSD'],
                        'broad_tissue_type': match['SMTS'],
                        'matched_gtex_id': match['SAMPID'],
                        'match_type': 'partial'
                    })
                elif len(base_matches) > 1:
                    annotations.append({
                        'rnaseq_sample': sample,
                        'cleaned_sample': cleaned_sample,
                        'tissue_type': None,
                        'broad_tissue_type': None,
                        'matched_gtex_id': ', '.join(base_matches['SAMPID'].tolist()),
                        'match_type': 'multiple_matches'
                    })
                else:
                    annotations.append({
                        'rnaseq_sample': sample,
                        'cleaned_sample': cleaned_sample,
                        'tissue_type': None,
                        'broad_tissue_type': None,
                        'matched_gtex_id': None,
                        'match_type': 'no_match'
                    })
            else:
                annotations.append({
                    'rnaseq_sample': sample,
                    'cleaned_sample': cleaned_sample,
                    'tissue_type': None,
                    'broad_tissue_type': None,
                    'matched_gtex_id': None,
                    'match_type': 'invalid_format'
                })
    
    # Convert annotations to DataFrame
    result_df = pd.DataFrame(annotations)
    
    # Print summary statistics
    print("\nAnnotation Summary:")
    print(f"Total samples: {len(result_df)}")
    print("\nMatch types:")
    print(result_df['match_type'].value_counts())
    
    if len(result_df[result_df['tissue_type'].notna()]) > 0:
        print("\nTissue types found:")
        print(result_df['tissue_type'].value_counts().dropna())
    
    # Print details about multiple matches
    multiple_matches = result_df[result_df['match_type'] == 'multiple_matches']
    if len(multiple_matches) > 0:
        print("\nSamples with multiple matches:")
        for _, row in multiple_matches.iterrows():
            print(f"\n{row['rnaseq_sample']}:")
            print(f"Matching GTEx IDs: {row['matched_gtex_id']}")
    
    return result_df

# %%
sample_names = gencode_tpm_data.columns

# Run annotation
annotated_samples = annotate_rnaseq_samples(sample_names, gtex_df)

# %%
# Display results
print("\nDetailed annotations:")
pd.set_option('display.max_columns', None)
print(annotated_samples[annotated_samples['match_type'].isin(['exact', 'partial'])][['rnaseq_sample', 'tissue_type', 'match_type', 'matched_gtex_id']])

# %% [markdown]
# **Tentative matches**
# 
# GTEX-1192X-0011-R10a-SM-4RXXZ:
# Matching GTEx IDs: `GTEX-1192X-0011-R10a-SM-CYKSL`, `GTEX-1192X-0011-R10a-SM-DO941`     
# 
# GTEX-14BIL-0011-R10a-SM-5EQV4:
# Matching GTEx IDs: `GTEX-14BIL-0011-R10a-SM-5SI75`, `GTEX-14BIL-0011-R10a-SM-AHZ7B`, `GTEX-14BIL-0011-R10a-SM-CYKOB`
# 
# GTEX-15DCD-0011-R10b-SM-5S51M:
# Matching GTEx IDs: `GTEX-15DCD-0011-R10b-SM-6LPII`, `GTEX-15DCD-0011-R10b-SM-CYKO7`

# %%
# Get rows for GTEX-1192X samples
print("GTEX-1192X matches:")
print(gtex_df[gtex_df['SAMPID'].isin(['GTEX-1192X-0011-R10a-SM-CYKSL', 'GTEX-1192X-0011-R10a-SM-DO941'])][['SAMPID', 'SMTSD', 'SMTS']])
print("\n")

# Get rows for GTEX-14BIL samples  
print("GTEX-14BIL matches:")
print(gtex_df[gtex_df['SAMPID'].isin(['GTEX-14BIL-0011-R10a-SM-5SI75', 'GTEX-14BIL-0011-R10a-SM-AHZ7B', 'GTEX-14BIL-0011-R10a-SM-CYKOB'])][['SAMPID', 'SMTSD', 'SMTS']])
print("\n")

# Get rows for GTEX-15DCD samples
print("GTEX-15DCD matches:")  
print(gtex_df[gtex_df['SAMPID'].isin(['GTEX-15DCD-0011-R10b-SM-6LPII', 'GTEX-15DCD-0011-R10b-SM-CYKO7'])][['SAMPID', 'SMTSD', 'SMTS']])

# %%
multiple_matches_mask = annotated_samples['match_type'] == "multiple_matches"
annotated_samples[multiple_matches_mask].head()

# %%
annotated_samples.loc[multiple_matches_mask, 'broad_tissue_type'] = "Brain - Frontal Cortex (BA9)"
annotated_samples.loc[multiple_matches_mask, 'tissue_type'] = "Brain"

annotated_samples[multiple_matches_mask].head()


# %%
annotated_samples.head()

# %%
print(annotated_samples[~annotated_samples.tissue_type.isna()].shape)
print(annotated_samples[~annotated_samples.broad_tissue_type.isna()].shape)


# %%
annotated_samples = annotated_samples.dropna()

# %%
annotated_samples.shape

# %%
annotated_samples.tissue_type.value_counts()

# %%
annotated_samples.broad_tissue_type.value_counts()

# %%
gencode_tpm_data_brain = gencode_tpm_data.loc[:, ['transcript'] + list(annotated_samples['rnaseq_sample'])]

# %%
gencode_tpm_data_brain.head()

# %% [markdown]
# # Visualize solely brain samples

# %%
# Create AnnData object from expression data
import anndata as ad
import pandas as pd

# Transpose data to have samples as rows and genes as columns
expression_matrix = gencode_tpm_data_brain.set_index('transcript').T

# Create AnnData object
adata = ad.AnnData(X=expression_matrix)

# Add tissue type information from annotated_samples as observation annotations
adata.obs['tissue_type'] = annotated_samples['tissue_type'].values
adata.obs['broad_tissue_type'] = annotated_samples['broad_tissue_type'].values

adata

# %%
print(adata[:, with_exon15].var_names)
print(adata[:, without_exon15].var_names)

# %%
# Create a figure and axis with a larger size for better visibility
plt.figure(figsize=(10, 8))

# Get data for both transcripts from adata object
with_exon15_data = adata[:, with_exon15].X.toarray().flatten()
without_exon15_data = adata[:, without_exon15].X.toarray().flatten()

# Get tissue types for coloring
tissue_types = adata.obs['tissue_type'].unique()
tissue_colors = plt.cm.Set3(np.linspace(0, 1, len(tissue_types)))
tissue_color_dict = dict(zip(tissue_types, tissue_colors))

# Create lists to store data by tissue type
with_exon15_by_tissue = []
without_exon15_by_tissue = []
labels = []

# Split data by tissue type
for tissue in tissue_types:
    tissue_mask = adata.obs['tissue_type'] == tissue
    with_exon15_by_tissue.append(adata[tissue_mask, with_exon15].X.toarray().flatten())
    without_exon15_by_tissue.append(adata[tissue_mask, without_exon15].X.toarray().flatten())

# Create labels and data in the desired order (all with exon15 first, then all without)
data = []
colors = []
labels = []

# Add all with_exon15 data first
for i, tissue in enumerate(tissue_types):
    data.append(with_exon15_by_tissue[i])
    colors.append(tissue_color_dict[tissue])
    labels.append(f"{tissue}\nWith exon 15")

# Then add all without_exon15 data
for i, tissue in enumerate(tissue_types):
    data.append(without_exon15_by_tissue[i])
    colors.append(tissue_color_dict[tissue])
    labels.append(f"{tissue}\nWithout exon 15")

# Create boxplot
positions = np.arange(1, len(data) + 1)
bp = plt.boxplot(data, patch_artist=True, tick_labels=labels, positions=positions)

# Customize boxplot colors
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

# Add individual points with jitter
for i, d in enumerate(data, 1):
    x = np.random.normal(i, 0.04, size=len(d))
    plt.scatter(x, d, alpha=0.4, s=20, color='black')

# Customize plot appearance
plt.ylabel('Expression (TPM)', fontsize=12)
plt.title('Expression Distribution of SRRM3 Isoforms by Tissue Type\nAcross GTEx Brain Samples', fontsize=14, pad=20)
plt.grid(True, linestyle='--', alpha=0.3)
plt.xticks(rotation=45, ha='right')

plt.tight_layout()
plt.show()

# Print summary statistics by tissue type
# print("\nSummary Statistics by Tissue Type:")
# for i, tissue in enumerate(tissue_types):
#     print(f"\n{tissue}:")
#     print("With exon 15:")
#     print(f"Mean: {np.mean(with_exon15_by_tissue[i]):.2f}")
#     print(f"Median: {np.median(with_exon15_by_tissue[i]):.2f}")
#     print(f"Std: {np.std(with_exon15_by_tissue[i]):.2f}")
    
#     print("\nWithout exon 15:")
#     print(f"Mean: {np.mean(without_exon15_by_tissue[i]):.2f}")
#     print(f"Median: {np.median(without_exon15_by_tissue[i]):.2f}")
#     print(f"Std: {np.std(without_exon15_by_tissue[i]):.2f}")

# %%
# Create a figure and axis with a larger size for better visibility
plt.figure(figsize=(10, 8))

# Get data for both transcripts from adata object
with_exon15_data = adata[:, "ENST00000611745.1"].X.toarray().flatten()
without_exon15_data = adata[:, "ENST00000612155.1"].X.toarray().flatten()

# Get tissue types for coloring
tissue_types = adata.obs['tissue_type'].unique()
tissue_colors = plt.cm.Set3(np.linspace(0, 1, len(tissue_types)))
tissue_color_dict = dict(zip(tissue_types, tissue_colors))

# Create lists to store data by tissue type
with_exon15_by_tissue = []
without_exon15_by_tissue = []
labels = []

# Split data by tissue type
for tissue in tissue_types:
    tissue_mask = adata.obs['tissue_type'] == tissue
    with_exon15_by_tissue.append(adata[tissue_mask, "ENST00000611745.1"].X.toarray().flatten())
    without_exon15_by_tissue.append(adata[tissue_mask, "ENST00000612155.1"].X.toarray().flatten())

# Create labels and data in the desired order (all with exon15 first, then all without)
data = []
colors = []
labels = []

# Add all with_exon15 data first
for i, tissue in enumerate(tissue_types):
    data.append(with_exon15_by_tissue[i])
    colors.append(tissue_color_dict[tissue])
    labels.append(f"{tissue}\nWith exon 15")

# Then add all without_exon15 data
for i, tissue in enumerate(tissue_types):
    data.append(without_exon15_by_tissue[i])
    colors.append(tissue_color_dict[tissue])
    labels.append(f"{tissue}\nWithout exon 15")

# Create boxplot
positions = np.arange(1, len(data) + 1)
bp = plt.boxplot(data, patch_artist=True, tick_labels=labels, positions=positions)

# Customize boxplot colors
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

# Add individual points with jitter
for i, d in enumerate(data, 1):
    x = np.random.normal(i, 0.04, size=len(d))
    plt.scatter(x, d, alpha=0.4, s=20, color='black')

# Customize plot appearance
plt.ylabel('Expression (TPM)', fontsize=12)
plt.title('Expression Distribution of SRRM3 Isoforms by Tissue Type\nAcross GTEx Brain Samples', fontsize=14, pad=20)
plt.grid(True, linestyle='--', alpha=0.3)
plt.xticks(rotation=45, ha='right')

plt.tight_layout()
plt.show()

# Print summary statistics by tissue type
# print("\nSummary Statistics by Tissue Type:")
# for i, tissue in enumerate(tissue_types):
#     print(f"\n{tissue}:")
#     print("With exon 15:")
#     print(f"Mean: {np.mean(with_exon15_by_tissue[i]):.2f}")
#     print(f"Median: {np.median(with_exon15_by_tissue[i]):.2f}")
#     print(f"Std: {np.std(with_exon15_by_tissue[i]):.2f}")
    
#     print("\nWithout exon 15:")
#     print(f"Mean: {np.mean(without_exon15_by_tissue[i]):.2f}")
#     print(f"Median: {np.median(without_exon15_by_tissue[i]):.2f}")
#     print(f"Std: {np.std(without_exon15_by_tissue[i]):.2f}")


