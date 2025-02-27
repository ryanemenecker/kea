'''
various utility functions for library generation. 
This includes functions for adding 5' and 3' adapters, 
addinging in additional padded sequences to get to a specific
length, and more. 
'''
import random


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

def check_for_restriction_site(sequence, restriction_site):
    """
    Check if a restriction site is present in the sequence.

    Parameters
    ----------
    sequence : str
        The DNA sequence to check.
    restriction_site : str
        The restriction site to look for.
    """
    return restriction_site in sequence

def add_padding(sequence, target_length, pad_location=3, 
                final_GC_content_range = (0.35, 0.45),
                avoid_added_start_codons=True,
                avoid_added_stop_codons=True,
                num_attempts=1000):
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
    num_attempts : int
        The number of attempts to find suitable padding.
    """
    # Calculate the number of bases to add
    num_bases_to_add = target_length - len(sequence)
    if num_bases_to_add <= 0:
        return sequence
        
    # Define start and stop codons to avoid
    start_codons = ["ATG"]
    stop_codons = ["TAA", "TAG", "TGA"]
    
    min_gc, max_gc = final_GC_content_range
    
    for _ in range(num_attempts):
        # Generate random padding
        bases = ['A', 'T', 'G', 'C']
        padding = ''.join(random.choice(bases) for _ in range(num_bases_to_add))
        
        # Create the candidate sequence with padding
        if pad_location == 5:
            candidate = padding + sequence
            junction = padding[-2:] + sequence[:1] if len(padding) >= 2 else padding + sequence[:1]
        else:  # pad_location == 3
            candidate = sequence + padding
            junction = sequence[-2:] + padding[:1] if len(sequence) >= 2 else sequence + padding[:1]
        
        # Check GC content
        gc_count = sum(1 for base in candidate if base in ['G', 'C'])
        gc_content = gc_count / len(candidate)
        
        if not (min_gc <= gc_content <= max_gc):
            continue
        
        # Check for unwanted start/stop codons at the junction
        valid_padding = True
        
        if avoid_added_start_codons:
            for i in range(len(junction) - 2):
                if junction[i:i+3] in start_codons:
                    valid_padding = False
                    break
                    
        if valid_padding and avoid_added_stop_codons:
            for i in range(len(junction) - 2):
                if junction[i:i+3] in stop_codons:
                    valid_padding = False
                    break
        
        # If we've found valid padding, return the result
        if valid_padding:
            return candidate
    
    # if no valid padding was found, raise exception
    raise ValueError("No valid padding found after maximum attempts.")


def add_padding_to_sequences(sequences, target_length, pad_location=3,
                                final_GC_content_range=(0.35, 0.45),
                                avoid_added_start_codons=True,
                                avoid_added_stop_codons=True,
                                num_attempts=1000):
        """
        Add padding to a list of DNA sequences to reach a target length.
        The padding can be added to the 5' or 3' end of the sequence.
        The padding is designed to avoid creating new start or stop codons,
        and to keep the GC content within a specified range.
    
        Parameters
        ----------
        sequences : list of str
            The original DNA sequences.
        target_length : int
            The desired length of the final sequences.
        pad_location : int
            The location to add padding: 5' (5) or 3' (3).
        final_GC_content_range : tuple of float
            The range of GC content for the final sequences.
        avoid_added_start_codons : bool
            If True, avoid adding padding that creates a new start codon.
        avoid_added_stop_codons : bool
            If True, avoid adding padding that creates a new stop codon.
        num_attempts : int
            The number of attempts to find suitable padding.
        """
        return [add_padding(seq, target_length, pad_location, 
                            final_GC_content_range,
                            avoid_added_start_codons,
                            avoid_added_stop_codons,
                            num_attempts) for seq in sequences]


