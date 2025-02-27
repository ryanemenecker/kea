'''
various utility functions for library generation. 
This includes functions for adding 5' and 3' adapters, 
addinging in additional padded sequences to get to a specific
length, and more. 
'''
import random
from kea.data.aa_codon_conversions import codons_to_aa, aa_to_codons

def add_5_prime_adapter(sequence, adapter_sequence):
    """
    Add a 5' adapter sequence to the DNA sequence.
    
    Args:
        sequence (str): The original DNA sequence.
        adapter_sequence (str): The adapter sequence to add.
        
    Returns:
        str: The DNA sequence with the 5' adapter added.
    """
    return adapter_sequence + sequence

def add_3_prime_adapter(sequence, adapter_sequence):
    """
    Add a 3' adapter sequence to the DNA sequence.
    
    Args:
        sequence (str): The original DNA sequence.
        adapter_sequence (str): The adapter sequence to add.
        
    Returns:
        str: The DNA sequence with the 3' adapter added.
    """
    return sequence + adapter_sequence

def generate_barcode_sequences(num_barcodes, length, 
                               GC_content_range=(0.35, 0.45),
                               max_barcode_similarity=0.5):
    """
    Generate a list of unique barcode sequences.

    Parameters
    ----------
    num_barcodes : int
        The number of unique barcodes to generate.
    length : int
        The length of each barcode.
    GC_content_range : tuple of float
        The range of GC content for the barcodes.
    max_barcode_similarity : float
        The maximum similarity allowed between barcodes.
    """
    pass

def add_padding(sequence, target_length, pad_location=3, 
                final_GC_content_range = (0.35, 0.45),
                avoid_added_start_codons=True,
                avoid_added_stop_codons=True):
    """
    Add padding to a DNA sequence to reach a target length.
    The padding can be added to the 5' or 3' end of the sequence.
    The padding is designed to avoid creating new start or stop codons,
    and to keep the GC content within a specified range.

    Parameters
    ----------
    sequence : str
        The original DNA sequence.
    target_length : int
        The desired length of the final sequence.
    pad_location : int
        The location to add padding: 5' (5) or 3' (3).
    final_GC_content_range : tuple of float
        The range of GC content for the final sequence.
    avoid_added_start_codons : bool
        If True, avoid adding padding that creates a new start codon.
    avoid_added_stop_codons : bool
        If True, avoid adding padding that creates a new stop codon.
    """
    
