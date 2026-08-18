'''
stuff for translating DNA sequences.
'''
from kea.data.aa_codon_conversions import codons_to_aa, aa_to_codons

# kept as a frozenset so membership tests are O(1) in the translation loop.
STOP_CODONS = frozenset(('TAA', 'TAG', 'TGA'))


def find_first_orf(sequence, frame=None):
    '''
    Find the first start codon (ATG) in a DNA sequence.

    parameters
    ----------
    sequence (str): DNA sequence string
    frame (int or None): If None, return the first ATG at any position, which is
        the scanning-ribosome model used for whole-construct checks. If 0, 1 or 2,
        only positions in that reading frame are considered, so the returned index
        is guaranteed to be a codon boundary for that frame.

    Returns
    --------
    int: Start index of the first matching ORF, or -1 if not found.
    '''
    sequence = sequence.upper()
    if frame is None:
        return sequence.find('ATG')
    if frame not in (0, 1, 2):
        raise ValueError("frame must be 0, 1, 2, or None")
    for i in range(frame, len(sequence) - 2, 3):
        if sequence[i:i + 3] == 'ATG':
            return i
    return -1


def translate_sequence(dna_sequence,
                       return_stop_codon=True,
                       return_nucleotide_sequence=False,
                       start_index=0,
                       stop_at_stop_codon=False):
    '''
    Translate a DNA sequence in the reading frame that begins at start_index.

    This translates every complete codon from start_index onwards and represents
    stop codons as "*", so a protein containing internal stop codons round-trips
    correctly. Callers that know where their coding sequence begins should pass
    start_index rather than searching for an ATG: a text search can land out of
    frame and silently shift the entire translation. Use
    translate_open_reading_frame() when you do want the scanning-ribosome model.

    Parameters
    ----------
    dna_sequence (str): DNA sequence string
    return_stop_codon (bool): If True, keep a trailing "*". If False, strip it.
    return_nucleotide_sequence (bool): If True, return the nucleotides that were
        translated instead of the amino acid sequence.
    start_index (int): Index of the first base of the first codon.
    stop_at_stop_codon (bool): If True, stop translating at the first stop codon.

    Returns
    --------
    str: Translated amino acid sequence, or the translated nucleotide sequence.

    Raises
    ------
    ValueError: If start_index is negative or the sequence contains a codon that
        is not made up of A, C, G and T.
    '''
    if start_index < 0:
        raise ValueError("start_index must be non-negative")

    dna_sequence = dna_sequence.upper()

    amino_acids = []
    codons = []
    for i in range(start_index, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i + 3]
        try:
            amino_acids.append(codons_to_aa[codon])
        except KeyError:
            raise ValueError(f"Invalid codon '{codon}' at position {i}. "
                             "Sequences must contain only A, C, G and T.") from None
        codons.append(codon)
        if stop_at_stop_codon and codon in STOP_CODONS:
            break

    if not return_stop_codon and amino_acids and amino_acids[-1] == '*':
        amino_acids.pop()
        codons.pop()

    if return_nucleotide_sequence:
        return ''.join(codons)
    return ''.join(amino_acids)


def translate_open_reading_frame(dna_sequence,
                                 return_stop_codon=True,
                                 return_nucleotide_sequence=False):
    '''
    Translate a DNA sequence the way a scanning ribosome would read it: start at
    the first ATG anywhere in the sequence and stop at the first in-frame stop
    codon.

    This is the right model for checking a whole construct (adapters + padding +
    coding sequence), because an ATG introduced by an adapter really would
    initiate translation there. It is the wrong model for verifying an isolated
    coding sequence -- use translate_sequence() for that.

    Parameters
    ----------
    dna_sequence (str): DNA sequence string
    return_stop_codon (bool): If True, include the terminating stop codon as "*".
    return_nucleotide_sequence (bool): If True, return the nucleotide sequence
        that was translated instead of the amino acid sequence.

    Returns
    --------
    str: Translated amino acid sequence, or "" if the sequence has no ATG.
    '''
    start_index = find_first_orf(dna_sequence)
    if start_index == -1:
        return ""
    return translate_sequence(dna_sequence,
                              return_stop_codon=return_stop_codon,
                              return_nucleotide_sequence=return_nucleotide_sequence,
                              start_index=start_index,
                              stop_at_stop_codon=True)
