'''
various utility functions for library generation. 
This includes functions for adding 5' and 3' adapters, 
addinging in additional padded sequences to get to a specific
length, and more. 
'''
import random
from kea.backend.optimize_codon_usage import _calculate_gc_content as calc_gc_content
from kea.data.aa_codon_conversions import codons_to_aa, aa_to_codons


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


def add_padding(current_gc_content,
                padding_length,
                target_length,
                current_length,
                final_gc_range,
                avoid_adding_start_codons,
                avoid_adding_stop_codons,
                num_attempts=5000,
                tolerance=0.02,
                return_best=True):
    """
    Add padding to a DNA sequence to reach a target length.
    The padding can be added to the 5' or 3' end of the sequence.
    The padding is designed to avoid creating new start or stop codons,
    and to keep the GC content within a specified range.

    Parameters
    ----------
    current_gc_content : float
        The current GC content of the DNA sequence.
    padding_length : int
        The length of the padding to add.
    target_length : int
        The target length of the final DNA sequence.
    final_gc_range : tuple of float
        The target GC content range for the final sequence.
    num_attempts : int
        The number of attempts to find suitable padding.
    tolerance : float
        The tolerance for the GC content.
    return_best : bool
        If True, return the best sequence even if it's outside the GC range.

    """
    # calculate possible target_GC content values
    # based on the length of the padding, the length
    # of the sequence, and the current GC content
    total_length = current_length + padding_length
    
    # Calculate required GC content in padding to achieve final GC range
    min_gc_in_padding = (final_gc_range[0] * total_length - current_gc_content * current_length) / padding_length
    max_gc_in_padding = (final_gc_range[1] * total_length - current_gc_content * current_length) / padding_length
    padding_gc_content_range = (min_gc_in_padding, max_gc_in_padding)

    
    # check if any target GC content values are valid
    if return_best==False:
        if padding_gc_content_range[0] > 1 or padding_gc_content_range[1] < 0:
            raise ValueError("The target GC content range is not achievable with the given sequence length.")

    if padding_gc_content_range[0] < 0:
        padding_gc_content_range = (0, padding_gc_content_range[1])
        if padding_gc_content_range[1] < 0:
            padding_gc_content_range = (0, 0+(tolerance*3))
    if padding_gc_content_range[1] > 1:
        padding_gc_content_range = (padding_gc_content_range[0], 1)
        if padding_gc_content_range[0] > 1:
            padding_gc_content_range = (1-(tolerance*3), 1)
    


    # create the padding sequence
    padding = create_padding_sequence(
        length=padding_length,
        GC_content_range=padding_gc_content_range, # pass the range to create_padding_sequence
        tolerance=tolerance,
        avoid_start_codon=avoid_adding_start_codons,
        avoid_stop_codon=avoid_adding_stop_codons,
        max_attempts=num_attempts
    )
    return padding


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