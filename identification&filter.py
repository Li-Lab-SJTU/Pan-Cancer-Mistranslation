# ============================================================================
# Proteomics Mistranslation Identification and Quality Filtering Pipeline
# ============================================================================
# Purpose: Identify mistranslation from proteomics data,
#          map peptides to proteins, and perform quality filtering including
#          SNP filtering and expression normalization.
# ============================================================================

import pandas as pd
import numpy as np
from Bio import SeqIO
from itertools import groupby
from operator import itemgetter
from collections import Counter


def suffix_array(text, _step=16):
    tx = text
    size = len(tx)
    step = min(max(_step, 1), len(tx))
    sa = list(range(len(tx)))
    sa.sort(key=lambda i: tx[i:i + step])
    grpstart = size * [False] + [True]  
    rsa = size * [None]
    stgrp, igrp = '', 0
    for i, pos in enumerate(sa):
        st = tx[pos:pos + step]
        if st != stgrp:
            grpstart[igrp] = (igrp < i - 1)
            stgrp = st
            igrp = i
        rsa[pos] = igrp
        sa[i] = pos
    grpstart[igrp] = (igrp < size - 1 or size == 0)
    while grpstart.index(True) < size:
        nextgr = grpstart.index(True)
        while nextgr < size:
            igrp = nextgr
            nextgr = grpstart.index(True, igrp + 1)
            glist = []
            for ig in range(igrp, nextgr):
                pos = sa[ig]
                if rsa[pos] != igrp:
                    break
                newgr = rsa[pos + step] if pos + step < size else -1
                glist.append((newgr, pos))
            glist.sort()
            for ig, g in groupby(glist, key=itemgetter(0)):
                g = [x[1] for x in g]
                sa[igrp:igrp + len(g)] = g
                grpstart[igrp] = (len(g) > 1)
                for pos in g:
                    rsa[pos] = igrp
                igrp += len(g)
        step *= 2
    del grpstart
    del rsa
    return sa

def SA_search(P, W, sa):
    lp = len(P)
    n = len(sa)
    l = 0; r = n
    while l < r:
        mid = (l + r) // 2
        a = sa[mid]
        if P > W[a : a + lp]:
            l = mid + 1
        else:
            r = mid
    s = l; r = n
    while l < r:
        mid = (l + r) // 2
        a = sa[mid]
        if P < W[a : a + lp]:
            r = mid
        else:
            l = mid + 1
    return [sa[i] for i in range(s, r)]

def find_proteins(base_seq):
    proteins = []
    uniprotIDs = []
    indices = np.searchsorted(boundaries_aa-1, SA_search(base_seq, W_aa, sa))
    for i in indices:
        proteins.append(names_list[i])
        uniprotIDs.append(uniprotID_list[i])
    proteins = " ".join(proteins)
    uniprotIDs = " ".join(uniprotIDs)
    indices = " ".join(map(str, indices))
    if proteins.strip(" ") == '':
        proteins = ''
    if uniprotIDs.strip(" ") == '':
        uniprotIDs = ''
    if indices.strip(" ") == '':
        indices = ''
    return proteins, uniprotIDs, indices

def One2three(one):
    one2three = {'A': 'Ala',
                 'C': 'Cys',
                 'D': 'Asp',
                 'E': 'Glu',
                 'F': 'Phe',
                 'G': 'Gly',
                 'H': 'His',
                 'I': 'Ile',
                 'K': 'Lys',
                 'L': 'Leu',
                 'M': 'Met',
                 'N': 'Asn',
                 'P': 'Pro',
                 'Q': 'Gln',
                 'R': 'Arg',
                 'S': 'Ser',
                 'T': 'Thr',
                 'V': 'Val',
                 'W': 'Trp',
                 'Y': 'Tyr',
                 '*': 'Ter'}
    return one2three[one]

def find_substitution_position_local(modified_seq, base_seq, index):
    """Find the position of a substitution within a protein sequence.
    
    Args:
        modified_seq: Peptide sequence with modification
        base_seq: Base sequence to search in
        index: Protein index in reference database
        
    Returns:
        1-based position in protein, or None if not found.
    """
    if len(modified_seq) > len(base_seq):
        return None
    
    # Find alignment position
    for i in range(len(base_seq) - len(modified_seq) + 1):
        diffs = [(base_seq[j], modified_seq[j], j) for j in range(len(modified_seq))
                 if base_seq[j + i] != modified_seq[j]]
        if len(diffs) == 1:
            base_pos = i + diffs[0][2]
            break
    else:
        return None
    
    # Find protein sequence start position
    protein_seq = seq_list[int(index)]
    seq_start = protein_seq.find(base_seq)
    return seq_start + base_pos + 1 if seq_start != -1 else None


def find_positions_local(modified_seq, base_seq, proteins, indices):
    """Find substitution positions for all protein matches.
    
    Args:
        modified_seq: Peptide sequence
        base_seq: Base sequence
        proteins: Protein names (unused but kept for API compatibility)
        indices: Space-separated protein indices
        
    Returns:
        Space-separated positions as strings.
    """
    positions = [str(find_substitution_position_local(modified_seq, base_seq, idx))
                 for idx in indices.split() if idx]
    return " ".join(positions)


def find_substitution(modified_seq, base_seq):
    """Identify amino acid substitution type (e.g., 'K to R').
    
    Args:
        modified_seq: Modified peptide sequence
        base_seq: Reference base sequence
        
    Returns:
        Substitution string (e.g., 'K to R') or empty string if not found.
    """
    if len(modified_seq) > len(base_seq):
        return ""
    
    for i in range(len(base_seq) - len(modified_seq) + 1):
        diffs = [(base_seq[j + i], modified_seq[j]) for j in range(len(modified_seq))
                 if base_seq[j + i] != modified_seq[j]]
        if len(diffs) == 1:
            return f"{diffs[0][0]} to {diffs[0][1]}"
    
    return ""

def union_sets(row):
    return ','.join(set(row.dropna()))

def normalize_expression(df, tmt, reference_channel):
    """Normalize TMT expression data using specified reference channel.
    
    Args:
        df: DataFrame with expression data
        tmt: TMT plex number (number of channels per group)
        reference_channel: Reference channel name for normalization
        
    Returns:
        DataFrame with normalized expression values.
    """
    # Calculate number of expression columns (excluding metadata)
    num_cols = df.shape[1] - 12 if 'error' in df.columns else df.shape[1]
    
    # Determine reference column index (0 for Ion_126.128, -1 otherwise)
    ref_idx = 0 if reference_channel == "Ion_126.128" else -1
    
    # Normalize expression by TMT groups
    for start_col in range(0, num_cols, tmt):
        end_col = min(start_col + tmt, num_cols)
        if end_col <= num_cols:
            group = df.iloc[:, start_col:end_col]
            divisor = group.iloc[:, ref_idx]
            
            # Handle zero or NaN divisors
            zero_or_nan_rows = (divisor == 0) | (divisor.isna())
            df.iloc[:, start_col:end_col] = group.div(divisor.replace(0, np.nan), axis=0)
            df.iloc[np.where(zero_or_nan_rows)[0], start_col:end_col] = 0
    
    return df

# ============================================================================
# Configuration and Data Loading
# ============================================================================
TMT_PLEX = 10
REFERENCE_CHANNEL = "Ion_131.138"
FASTA_DB = 'uniprotkb_reviewed_true_AND_model_organ.fasta'
SNP_FILE = 'snp_uniprot.txt'
INPUT_FILE = 'expression/all_raw_expression.tsv'
OUTPUT_INITIAL = 'subs_twostep.csv'
OUTPUT_FINAL = 'subs_normexpr_filt.csv'

# Load expression data
subs = pd.read_csv(INPUT_FILE, sep='\t')

# Consolidate protein columns
protein_cols = [col for col in subs.columns if 'Protein' in col]
subs['DP base sequences'] = subs[protein_cols].apply(union_sets, axis=1)
subs.drop(protein_cols, axis=1, inplace=True)
subs['DP base sequence'] = subs['DP base sequences'].apply(
    lambda x: x.split('_')[1] if '_' in x else ''
)

# ============================================================================
# Load Reference Protein Database
# ============================================================================
record_list = []
translated_record_list = []
names_list = []
seq_list = []
uniprotID_list = []
dict_uniprot = {}
boundaries_aa = [0]

with open(FASTA_DB, 'r') as fasta_handle:
    for record in SeqIO.parse(fasta_handle, 'fasta'):
        # Extract UniProt ID
        uniprotID_list.append(record.name.split('|')[1])
        
        # Extract gene name from description
        for field in record.description.split():
            if 'GN=' in field:
                record.name = field.split('=')[-1]
                break
        
        # Store sequence data
        record_list.append(record)
        seq_str = str(record.seq)
        translated_record_list.append(seq_str)
        names_list.append(record.name)
        seq_list.append(record.seq)
        dict_uniprot[record.name] = record.seq
        boundaries_aa.append(boundaries_aa[-1] + len(seq_str))

# ============================================================================
# Build Suffix Array for Sequence Search
# ============================================================================
boundaries_aa = np.array(boundaries_aa[1:])
W_aa = ''.join(translated_record_list)
sa = suffix_array(W_aa)  # Suffix array for efficient protein matching

# ============================================================================
# Protein Mapping and Substitution Identification
# ============================================================================

# Find matching proteins for each base sequence
subs[['proteins', 'uniprotIDs', 'indices']] = subs['DP base sequence'].apply(
    lambda x: pd.Series(find_proteins(x))
)

# Extract primary protein information
subs['protein'] = subs['proteins'].map(
    lambda x: x.split()[0] if x else float('nan')
)
subs['uniprotID'] = subs['uniprotIDs'].map(
    lambda x: x.split()[0] if x else float('nan')
)

# Filter to rows with valid protein mappings
subs = subs[pd.notnull(subs['protein'])]

# Identify amino acid substitutions
subs['substitution'] = subs.apply(
    lambda row: find_substitution(row['Peptide'], row['DP base sequence']),
    axis=1
)

# Extract substitution components
subs['destination'] = subs['substitution'].map(
    lambda x: x[-1] if x else False
)
subs['origin'] = subs['substitution'].map(
    lambda x: x[0] if x else False
)

# Convert to 3-letter amino acid codes
subs['destination3'] = subs['destination'].map(One2three)
subs['origin3'] = subs['origin'].map(One2three)

# Find substitution positions in protein sequences
subs['positions'] = subs.apply(
    lambda row: find_positions_local(
        row['Peptide'], row['DP base sequence'],
        row['proteins'], row['indices']
    ),
    axis=1
)

# Extract primary position
subs['position'] = subs['positions'].map(
    lambda x: int(x.split()[0]) if x and x.split() else float('nan')
)

# Create mistranslation identifier (protein_position_substitution)
subs['error'] = subs.apply(
    lambda row: f"{row['protein']}_{row['position']}_{row['substitution']}",
    axis=1
)
# ============================================================================
# Data Standardization and I/L Normalization
# ============================================================================
# Rename peptide column
subs.rename(columns={'Peptide': 'modified_sequence'}, inplace=True)

# Normalize I/L (Isoleucine/Leucine are indistinguishable in MS)
subs['error'] = subs['error'].str.replace(r'[IL]$', 'I/L', regex=True)
subs['substitution'] = subs['substitution'].str.replace(r'[IL]$', 'I/L', regex=True)
subs['destination'] = subs['destination'].replace(['I', 'L'], 'I/L')
subs['destination3'] = subs['destination3'].replace(['Ile', 'Leu'], 'Ile/Leu')

# ============================================================================
# Aggregation by Error and Expression Normalization
# ============================================================================

# Define columns for aggregation
non_numeric_cols = {
    'modified_sequence', 'DP base sequences', 'DP base sequence',
    'proteins', 'protein', 'uniprotIDs', 'uniprotID',
    'substitution', 'destination', 'origin', 'destination3', 'origin3',
    'positions', 'position', 'error', 'indices'
}
expression_cols = [col for col in subs.columns if col not in non_numeric_cols]
metadata_cols = [
    'modified_sequence', 'DP base sequences', 'DP base sequence',
    'protein', 'uniprotID', 'substitution',
    'destination', 'origin', 'destination3', 'origin3', 'position'
]

# Create aggregation dictionary
agg_dict = {col: 'first' for col in metadata_cols}
agg_dict.update({col: 'sum' for col in expression_cols})
agg_dict['modified_sequence'] = lambda x: ';'.join(set(str(v) for v in x if pd.notna(v)))
agg_dict['DP base sequences'] = lambda x: ';'.join(set(str(v) for v in x if pd.notna(v)))

# Group by error identifier and aggregate
subs = subs.groupby('error', as_index=False).agg(agg_dict)
subs = subs[expression_cols + ['error'] + metadata_cols]

# Normalize by TMT reference channel
subs = normalize_expression(subs, TMT_PLEX, REFERENCE_CHANNEL)

# Save aggregated data
subs.to_csv(OUTPUT_INITIAL, index=False)
print(f"Aggregated data saved: {OUTPUT_INITIAL}")

# ============================================================================
# Expression Normalization and Quality Filtering
# ============================================================================

# Load aggregated data
subs = pd.read_csv(OUTPUT_INITIAL, sep=',')
all_cols = subs.columns.tolist()
metadata_cols = all_cols[-12:]
expression_cols = all_cols[:-12]

# Reorder columns: metadata first, then expression
subs = subs[metadata_cols + expression_cols]

# Remove low-quality samples
subs.drop(
    [col for col in subs.columns if 'Not' in col or 'Disqualified' in col],
    axis=1,
    inplace=True
)

# Normalize expression within sample groups
sample_cols = [col for col in subs.columns if 'Tumor' in col or 'Normal' in col]
subs[sample_cols] = subs[sample_cols].replace(0, np.nan)
medians = subs[sample_cols].median()
subs[sample_cols] = subs[sample_cols] / medians

# Merge technical replicates (same base name, different sub-indices)
col_bases = [col.split('.')[0] if '.' in col else col for col in subs.columns]
duplicate_bases = [name for name, cnt in Counter(col_bases).items() if cnt > 1]

for base_name in duplicate_bases:
    matching_cols = subs.columns[subs.columns.str.startswith(base_name)]
    subs[base_name] = subs[matching_cols].mean(axis=1)
    subs = subs.drop(matching_cols[1:], axis=1)

# ============================================================================
# SNP Filtering and Final Quality Control
# ============================================================================

# Create SNP filter key (UniProt format)
subs['SNPfilter'] = subs.apply(
    lambda row: '{}:p.{}{}{}'.format(row['uniprotID'], row['origin3'], row['position'], row['destination3']),
    axis=1
)

# Load known SNPs
with open(SNP_FILE, 'r') as f:
    snp_set = set(line.strip() for line in f if line.strip())

# Normalize I/L in SNP set
snp_set = set(
    x[:-3] + 'Ile/Leu' if x[-3:] in ['Ile', 'Leu'] else x
    for x in snp_set
)

# Filter out known SNPs
original_count = len(subs)
subs_filt = subs[~subs['SNPfilter'].isin(snp_set)]
print(f"Removed {original_count - len(subs_filt)} known SNPs")

# Filter out rows with zero total expression
subs_filt = subs_filt[subs_filt.iloc[:, 12:-1].sum(axis=1) != 0]
print(f"Final dataset: {len(subs_filt)} errors")

# Save final filtered dataset
subs_filt.to_csv(OUTPUT_FINAL, index=False)
print(f"Filtered data saved: {OUTPUT_FINAL}")
