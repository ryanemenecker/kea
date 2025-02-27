'''
various utility functions.
'''
from kea.data.aa_codon_conversions import codons_to_aa, aa_to_codons
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import defaultdict

# codon to aa
def dna_to_aa(dna_sequence):
    """
    Convert a DNA sequence to amino acid sequence using the codon table
    
    Args:
        dna_sequence (str): DNA sequence string
        
    Returns:
        str: Amino acid sequence
    """
    aa_sequence = ""
    # Process DNA sequence in chunks of 3 (codons)
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        # Skip if we don't have a complete codon
        if len(codon) != 3:
            break
        # Convert codon to amino acid using the lookup table
        aa = codons_to_aa.get(codon.upper(), 'X')
        aa_sequence += aa
    return aa_sequence


def calc_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count / len(sequence)


def read_fasta(path_to_file):
    '''
    reads a fasta file and returns a list of sequences.
    '''
    with open(path_to_file, 'r') as f:
        lines = f.readlines()
    
    sequences = []
    current_sequence = ""
    
    for line in lines:
        if line.startswith('>'):
            if current_sequence:
                sequences.append(current_sequence)
                current_sequence = ""
        else:
            current_sequence += line.strip()
    
    if current_sequence:
        sequences.append(current_sequence)
    
    return sequences

def make_codon_table(sequences, require_start_codon=False):
    '''
    function that makes a codon table from a sequence or sequences
    '''
    # Initialize codon counts
    codon_counts = {}
    for aa in aa_to_codons:
        codon_counts[aa] = {codon: 0 for codon in aa_to_codons[aa]}
    if isinstance(sequences, str):
        sequences = [sequences]
    # Count codons
    for sequence in sequences:
        if require_start_codon and sequence[:3]!='ATG':
            continue
        if len(sequence) % 3 != 0:
            continue
        if not all(c in 'ACGT' for c in sequence):
            continue
        if len(sequence) < 3:
            continue
        # Process the sequence codon by codon
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i+3].upper()
            # Use direct lookup rather than nested loops
            if codon in codons_to_aa:
                aa = codons_to_aa[codon]
                if aa in codon_counts and codon in codon_counts[aa]:
                    codon_counts[aa][codon] += 1
    
    # Calculate frequencies
    codon_usage = {}
    for aa, counts in codon_counts.items():
        total = sum(counts.values())
        if total > 0:  # Avoid division by zero
            codon_usage[aa] = {codon: count/total for codon, count in counts.items()}
        else:
            codon_usage[aa] = {codon: 0.0 for codon in counts}
    
    # round values
    for aa, counts in codon_usage.items():
        for codon, freq in counts.items():
            codon_usage[aa][codon] = round(freq, 5)
    return codon_usage    

def calculate_codon_usage_table(path_to_file, require_start_codon=True):
    '''
    Calculate codon usage frequencies from coding sequences in a FASTA file.
    
    Args:
        path_to_file (str): Path to FASTA file containing coding sequences
        require_start_codon (bool): If True, sequences must start with 'ATG'
        
    Returns:
        dict: Nested dictionary of codon usage frequencies
              Format: {amino_acid: {codon: frequency}}
    '''
    try:
        # Read sequences from FASTA file
        sequences = read_fasta(path_to_file)
        
        # check sequences..
        final_sequences=[]
        # Count codons
        for sequence in sequences:
            if sequence[:3]!='ATG':
                continue
            if len(sequence) % 3 != 0:
                continue
            if not all(c in 'ACGT' for c in sequence):
                continue
            if len(sequence) < 3:
                continue
            final_sequences.append(sequence)
        
        codon_table = make_codon_table(final_sequences, require_start_codon=require_start_codon)
        return codon_table

        
    except FileNotFoundError:
        print(f"Error: File not found: {path_to_file}")
        return None
    except Exception as e:
        print(f"Error calculating codon usage: {str(e)}")
        return None


def graph_codon_usage(sequence, codon_table):
    '''
    Creates a single plot that compares codon usage between a sequence (or list of 
    sequences) and a reference codon table.
    
    Parameters
    ----------
    sequence : str or list
        DNA sequence or list of DNA sequences to analyze.
    codon_table : dict
        Reference codon usage table {aa: {codon: frequency}}
    
    Returns
    -------
    matplotlib.figure.Figure
        The generated figure object
    '''
    # Convert single sequence to list for consistent handling
    if isinstance(sequence, str):
        sequences = [sequence]
    else:
        sequences = sequence
    dna_sequences = [seq for seq in sequences if isinstance(seq, str)]

    # Calculate codon usage for the provided sequences
    sequence_codon_freqs = make_codon_table(dna_sequences)
    
    # Create a single plot with all codons
    fig, ax = plt.subplots(figsize=(16, 8))
    
    # Prepare data for plotting
    all_codons = []
    ref_values = []
    seq_values = []
    aa_boundaries = [0]  # Track where each amino acid group ends
    aa_centers = []      # Track center position for amino acid labels
    
    # Sort amino acids for consistent display
    amino_acids = sorted(codon_table.keys())
    
    # Collect data for each amino acid
    current_position = 0
    for aa in amino_acids:
        codons = sorted(codon_table[aa].keys())
        start_pos = len(all_codons)
        
        for codon in codons:
            all_codons.append(codon)
            ref_values.append(codon_table[aa].get(codon, 0))
            if aa in sequence_codon_freqs and codon in sequence_codon_freqs[aa]:
                seq_values.append(sequence_codon_freqs[aa][codon])
            else:
                seq_values.append(0)
        
        # Calculate center position for this amino acid group
        if codons:  # Only if there are codons for this amino acid
            aa_centers.append((start_pos + len(all_codons) - 1) / 2)
            # Update the boundaries after processing each amino acid
            aa_boundaries.append(len(all_codons))
    
    # Plot bars
    x = np.arange(len(all_codons))
    width = 0.35
    
    ref_bars = ax.bar(x - width/2, ref_values, width, label='Reference', color='steelblue', alpha=0.8)
    seq_bars = ax.bar(x + width/2, seq_values, width, label='Sequence', color='darkorange', alpha=0.8)
    
    # Add vertical lines to separate amino acid groups
    for boundary in aa_boundaries[1:-1]:  # Skip first and last boundaries
        ax.axvline(x=boundary - 0.5, color='gray', linestyle='--', alpha=0.3)
    
    # Add amino acid labels
    for i, aa in enumerate(amino_acids):
        if i < len(aa_centers):
            ax.text(aa_centers[i], -0.05, aa, transform=ax.get_xaxis_transform(),
                    ha='center', fontsize=12, fontweight='bold')
    
    # Customize plot
    ax.set_xticks(x)
    ax.set_xticklabels(all_codons, rotation=90, fontsize=8)  # Smaller font size for readability
    
    # Fix overlapping x-axis labels by limiting the number of visible ticks if too many
    if len(all_codons) > 40:
        # Show only a subset of ticks to avoid overcrowding
        for i, label in enumerate(ax.get_xticklabels()):
            if i % 2 != 0:  # Show every other label
                label.set_visible(False)
    
    ax.set_xlabel('Codons', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title('Codon Usage Comparison', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    
    # Handle empty values gracefully
    max_value = max(max(ref_values) if ref_values else 0, 
                    max(seq_values) if seq_values else 0)
    ax.set_ylim(0, max_value * 1.1 or 1.0)  # Fallback to 1.0 if all values are 0
    
    # Add value labels for bars with significant usage
    def add_value_labels(bars):
        threshold = 0.05  # Show more values by lowering threshold
        for bar in bars:
            height = bar.get_height()
            if height >= threshold:
                ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                        f'{height:.2f}', ha='center', va='bottom', fontsize=8)
    
    add_value_labels(ref_bars)
    add_value_labels(seq_bars)
    
    plt.tight_layout()
    
    return fig
