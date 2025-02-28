'''
various utility functions for library generation. 
This includes functions for adding 5' and 3' adapters, 
addinging in additional padded sequences to get to a specific
length, and more. 
'''
import random
from kea.backend.optimize_codon_usage import _calculate_gc_content as calc_gc_content
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

def create_padding_sequence(length,
                            GC_content_range, # changed from GC_content_target to GC_content_range
                            tolerance,
                            avoid_start_codon=True,
                            avoid_stop_codon=True,
                            avoid_adding=None, 
                            max_attempts=2000):
    """
    Create the padding sequence.
    This is useful for libraries that need to be a fixed size
    but have different lengths of coding sequences.

    Parameters
    ----------
    length : int
        Length of the padding sequence.
    GC_content_range : tuple of float
        Target GC content range for the padding sequence.
    tolerance : float
        Tolerance for the GC content.
    avoid_start_codon : bool
        If True, avoid creating a start codon.
    avoid_stop_codon : bool
        If True, avoid creating a stop codon.
    avoid_adding : list of str
        List of sequences to avoid adding.
        Will check for presence of these sequences in the padding.
        Will not return any padding with these sequences.
    max_attempts : int
        Maximum number of attempts to generate a valid padding sequence.
    """
    # Define start and stop codons to avoid
    start_codons = ["ATG"]
    stop_codons = ["TAA", "TAG", "TGA"]

    # Set avoid_adding to empty list if None
    if avoid_adding is None:
        avoid_adding = []
    if isinstance(avoid_adding, str):
        avoid_adding = [avoid_adding]

    # Define nucleotide weights to achieve target GC content
    # Use the midpoint of the GC content range for nucleotide weights
    GC_content_target = (GC_content_range[0] + GC_content_range[1]) / 2
    gc_weight = GC_content_target / 2  # Split between G and C
    at_weight = (1 - GC_content_target) / 2  # Split between A and T
    nucleotide_weights = {'G': gc_weight, 'C': gc_weight, 
                        'A': at_weight, 'T': at_weight}
    
    # Define alternative bases to replace problematic codons
    # These maintain similar GC content when possible
    replacements = {
        'A': ['T', 'C', 'G'],
        'T': ['A', 'C', 'G'],
        'G': ['C', 'A', 'T'],
        'C': ['G', 'A', 'T']
    }
    
    for attempt in range(max_attempts):
        # Generate a random padding sequence
        padding = ''.join(random.choices(
            list(nucleotide_weights.keys()),
            weights=list(nucleotide_weights.values()),
            k=length
        ))
        
        # Check GC content
        gc_content = calc_gc_content(padding)
        if not (GC_content_range[0] - tolerance <= gc_content <= GC_content_range[1] + tolerance):
            continue
        
        # Convert to list for easier modification
        padding_list = list(padding)
        modified = False
        
        # Fix start codons by replacing the first base of each codon
        if avoid_start_codon:
            for i in range(0, len(padding_list) - 2):
                codon = ''.join(padding_list[i:i+3])
                if codon in start_codons:
                    # Choose a replacement that isn't 'A' for the first position of ATG
                    replacement = random.choice([b for b in ['T', 'C', 'G'] if b != padding_list[i]])
                    padding_list[i] = replacement
                    modified = True
        
        # Fix stop codons by replacing the middle base of each codon
        if avoid_stop_codon:
            for i in range(0, len(padding_list) - 2):
                codon = ''.join(padding_list[i:i+3])
                if codon in stop_codons:
                    # Replace middle base (most effective for breaking stop codons)
                    replacement = random.choice([b for b in ['A', 'T', 'C', 'G'] if b != padding_list[i+1]])
                    padding_list[i+1] = replacement
                    modified = True
        
        # If modifications were made, reassemble the sequence and recheck GC content
        if modified:
            padding = ''.join(padding_list)
            gc_content = calc_gc_content(padding)
            if not (GC_content_range[0] - tolerance <= gc_content <= GC_content_range[1] + tolerance):
                continue
        
        # Check for sequences to avoid - can't easily modify these in place
        # so we'll still reject sequences containing them
        if any(seq in padding for seq in avoid_adding):
            continue
        
        # Double check no start/stop codons were missed or created during modifications
        if avoid_start_codon and any(codon in padding for codon in start_codons):
            continue
            
        if avoid_stop_codon and any(codon in padding for codon in stop_codons):
            continue
        
        # If all checks pass, return the padding sequence
        return padding

    # If no valid padding sequence is found after max attempts
    raise ValueError(f"Could not generate a suitable padding sequence after {max_attempts} attempts. "
                    f"Consider relaxing constraints or increasing the tolerance.")


def add_padding(sequence, target_length, pad_location=3,
                final_GC_content_range=(0.35, 0.45),
                avoid_added_start_codons=True,
                avoid_added_stop_codons=True,
                num_attempts=1000,
                tolerance=0.02):
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
    tolerance : float
        The tolerance for the GC content.
    """
    # calculate current sequence GC content
    current_gc_content = calc_gc_content(sequence)
    # calculate possible target_GC content values
    # based on the length of the padding, the length
    # of the sequence, and the current GC content
    
    padding_length = target_length - len(sequence)

    min_effective_padding = target_length*final_GC_content_range[0]
    max_effective_padding = target_length*final_GC_content_range[1]
    min_effective_padding -= (current_gc_content * len(sequence))
    max_effective_padding -= (current_gc_content * len(sequence))
    min_effective_padding /= padding_length
    max_effective_padding /= padding_length
    padding_gc_content_range = (min_effective_padding, max_effective_padding)
    
    # check if any target GC content values are valid
    if padding_gc_content_range[0] > 1 or padding_gc_content_range[1] < 0:
        raise ValueError("The target GC content range is not achievable with the given sequence length.")

    # create the padding sequence
    padding = create_padding_sequence(
        length=target_length - len(sequence),
        GC_content_range=padding_gc_content_range, # pass the range to create_padding_sequence
        tolerance=tolerance,
        avoid_start_codon=avoid_added_start_codons,
        avoid_stop_codon=avoid_added_stop_codons,
        max_attempts=num_attempts
    )
    # add the padding to the sequence
    if pad_location == 5:
        return padding + sequence
    elif pad_location == 3:
        return sequence + padding
    else:
        raise ValueError("pad_location must be 5 or 3.")
    

def add_padding_to_sequences(sequences, target_length, pad_location=3,
                                final_GC_content_range=(0.35, 0.45),
                                avoid_added_start_codons=True,
                                avoid_added_stop_codons=True,
                                num_attempts=1000,
                                tolerance=0.02):
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
        tolerance : float
            The tolerance for the GC content.
        """
        return [add_padding(seq, target_length, pad_location, 
                            final_GC_content_range,
                            avoid_added_start_codons,
                            avoid_added_stop_codons,
                            num_attempts, tolerance=tolerance) for seq in sequences]


