'''
various utility functions.
'''
from kea.data.aa_codon_conversions import codons_to_aa, aa_to_codons
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import defaultdict
import random

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


def codon_from_fasta(path_to_file, require_start_codon=True):
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

def get_restriction_enzymes(path_to_file):
    '''
    Get a list of restriction enzymes from a file.
    
    Args:
        path_to_file (str): Path to the file containing restriction enzyme information
        
    Returns:
        list: List of restriction enzymes
    '''
    try:
        restriction_enzyme_dict={}
        with open(path_to_file, 'r') as f:
            lines = f.read().split('\n')
        
        for line in lines:
            site=line.split()[0]
            enzymes=line.split()[1:]
            for e in enzymes:
                restriction_enzyme_dict[e]=site
        return restriction_enzyme_dict
    except FileNotFoundError:
        print(f"Error: File not found: {path_to_file}")
        return None
    except Exception as e:
        print(f"Error reading restriction enzymes: {str(e)}")
        return None
    
def generate_random_protein_ids(length_of_id, 
                               number_of_ids,
                               alphanumeric=True,
                               cap_sensitive=False):
    """
    Generate random protein IDs.

    Parameters
    ----------
    length_of_id (int): Length of each protein ID.
    number_of_ids (int): Number of protein IDs to generate.
    alphanumeric (bool): If True, include both letters and digits. 
                         If False, include only numbers.
    cap_sensitive (bool): If True, include both uppercase and lowercase letters.
                          If False, include only uppercase letters.
    Returns
    -------
    list: A list of unique random protein IDs.
    
    Raises
    ------
    ValueError: If the number of requested IDs exceeds the possible unique combinations.
    """
    # Define the character sets
    digits = '0123456789'
    uppercase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    lowercase = 'abcdefghijklmnopqrstuvwxyz'
    
    # Select characters based on parameters
    if alphanumeric:
        if cap_sensitive:
            # Include digits, uppercase, and lowercase
            characters = digits + uppercase + lowercase
        else:
            # Include digits and uppercase only
            characters = digits + uppercase
    else:
        # Only digits, regardless of cap_sensitive
        characters = digits
    
    # Check if it's possible to generate the requested number of unique IDs
    max_possible_ids = len(characters) ** length_of_id
    if number_of_ids > max_possible_ids:
        raise ValueError(f"Cannot generate {number_of_ids} unique IDs with length {length_of_id}. "
                         f"Maximum possible is {max_possible_ids}.")
    
    unique_ids = set()
    while len(unique_ids) < number_of_ids:
        new_id = ''.join(random.choice(characters) for _ in range(length_of_id))
        unique_ids.add(new_id)
    
    return list(unique_ids)

