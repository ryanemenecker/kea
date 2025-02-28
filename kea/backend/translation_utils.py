'''
stuff for translating DNA sequences. 
'''
from kea.data.aa_codon_conversions import codons_to_aa, aa_to_codons


def find_first_orf(sequence):
    '''
    Find the first open reading frame (ORF) in a DNA sequence.
    parameters
    ----------
    sequence (str): DNA sequence string

    Returns
    --------
    int: Start index of the first ORF, or -1 if not found.
    '''
    start_codon = "ATG"
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if codon == start_codon:
            return i
    return -1

def translate_sequence(dna_sequence, return_stop_codon=False):
    '''
    Function to translate a DNA sequence. 
    The function will start translation at the first start codon
    and stop translation at the first stop codon.
    
    Parameters
    ----------
    dna_sequence (str): DNA sequence string
    
    Returns
    --------
    str: Translated amino acid sequence.
    '''
    start_index = find_first_orf(dna_sequence)
    if start_index == -1:
        return ""
    has_stop_codon = False
    protein_sequence = ""
    for i in range(start_index, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i+3]
        if codon in ['TAA', 'TAG', 'TGA']:
            has_stop_codon = True
            break
        protein_sequence += codons_to_aa[codon]
    if return_stop_codon:
        if has_stop_codon:
            return protein_sequence + "*"
    return protein_sequence


