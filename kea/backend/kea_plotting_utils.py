'''
Code for plotting stuff.
It's nice to visuzlie GC content and codon usage.
'''

import matplotlib.pyplot as plt
from kea.backend.kea_utils import calc_gc_content, make_codon_table
import numpy as np


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

def graph_gc_content(sequences):
    '''
    plots the GC content of a list of sequences as a histogram

    Parameters
    ----------
    sequences : list
        list of sequences to plot
    
    Returns
    -------
    matplotlib.figure.Figure
        The generated figure object
    '''
    if isinstance(sequences, str):
        sequences = [sequences]
    # Calculate GC content for each sequence
    gc_contents = [calc_gc_content(seq) for seq in sequences]
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot histogram
    n, bins, patches = ax.hist(gc_contents, bins=30, color='skyblue', 
                              edgecolor='black', alpha=0.7)
    
    # Customize plot
    ax.set_xlabel('GC Content (%)', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title('Distribution of GC Content', fontsize=14, fontweight='bold')
    
    # Add mean line
    mean_gc = np.mean(gc_contents)
    ax.axvline(mean_gc, color='red', linestyle='--', alpha=0.8,
               label=f'Mean GC: {mean_gc:.1f}%')
    
    # Add legend
    ax.legend()
    
    # Set x-axis limits to 0-100 since GC content is a percentage
    ax.set_xlim(0, 100)
    
    plt.tight_layout()
    
    return fig
