'''
various utility functions for library generation. 
This includes functions for adding 5' and 3' adapters, 
addinging in additional padded sequences to get to a specific
length, and more. 
'''
import random
import numpy as np
from kea.backend.optimize_codon_usage import _calculate_gc_content as calc_gc_content
from kea.data.aa_codon_conversions import codons_to_aa, aa_to_codons


# A base swapped for its same-GC-class partner leaves GC content unchanged, so try
# that first when repairing a forbidden codon.
_SAME_GC_PARTNER = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def _repair_base(padding_list, pos, forbidden):
    """
    Replace padding_list[pos] with a base that leaves every 3-mer overlapping pos
    free of forbidden codons. Prefers a same-GC-class swap so GC content is
    preserved. Returns True if such a base was found (and applied), else restores
    the original base and returns False.
    """
    original = padding_list[pos]
    partner = _SAME_GC_PARTNER[original]
    others = [b for b in 'ACGT' if b != original and b != partner]
    random.shuffle(others)

    first_window = max(0, pos - 2)
    last_window = min(pos, len(padding_list) - 3)

    for base in [partner] + others:
        padding_list[pos] = base
        if all(''.join(padding_list[s:s + 3]) not in forbidden
               for s in range(first_window, last_window + 1)):
            return True
    padding_list[pos] = original
    return False


def _repair_forbidden_codons(padding_list, forbidden, max_passes=20):
    """
    Rewrite bases in place until no forbidden 3-mer remains.

    The original implementation scanned forward once and patched each hit, but a
    patch can create a new forbidden codon in a window the scan has already passed
    (e.g. fixing TAA by rewriting the middle base can yield TGA). Those escapes hit
    the final recheck and caused the whole attempt to be thrown away, which made
    padding longer than ~800 nt fail every time. Here we re-scan from just before
    each edit and repeat until a pass changes nothing.

    Returns True if the sequence is clean.
    """
    if not forbidden:
        return True

    for _ in range(max_passes):
        changed = False
        i = 0
        while i <= len(padding_list) - 3:
            if ''.join(padding_list[i:i + 3]) in forbidden:
                # Try the middle base first: it belongs to every window containing
                # this codon and cannot re-form the same codon one position over.
                if not any(_repair_base(padding_list, i + offset, forbidden)
                           for offset in (1, 0, 2)):
                    return False
                changed = True
                # An edit can invalidate the two windows starting before i.
                i = max(0, i - 2)
            else:
                i += 1
        if not changed:
            return True
    return False


def _forbidden_codon_free(bases, position, forbidden, mutable):
    """
    True if no forbidden 3-mer lying wholly inside the mutable (padding) region
    covers `position`.

    Windows that straddle the coding sequence are ignored on purpose: the CDS
    legitimately contains ATG and stop codons, and avoid_adding_start_codons only
    promises not to *add* them to the padding.
    """
    if not forbidden:
        return True
    for start in range(max(0, position - 2), min(position + 1, len(bases) - 2)):
        window = range(start, start + 3)
        if all(index in mutable for index in window):
            if ''.join(bases[start:start + 3]) in forbidden:
                return False
    return True


def _mutate_padding_base(bases, positions, mutable, forbidden):
    """
    Rewrite one of `positions` to a different base, preferring a same-GC-class
    swap so GC content is preserved. Returns True if a base was changed.
    """
    candidates = list(positions)
    random.shuffle(candidates)
    for position in candidates:
        original = bases[position]
        partner = _SAME_GC_PARTNER[original]
        others = [b for b in 'ACGT' if b != original and b != partner]
        random.shuffle(others)
        for base in [partner] + others:
            bases[position] = base
            if _forbidden_codon_free(bases, position, forbidden, mutable):
                return True
        bases[position] = original
    return False


def repair_padding_constraints(region, mutable_ranges, constraint_set,
                               forbidden=frozenset(), max_passes=200):
    """
    Rewrite padding bases until no constraint violation touches the padding.

    Padding is generated after the coding sequence has already been repaired, and
    it is otherwise arbitrary DNA -- so it happily re-introduces the very splice
    sites and polyadenylation signals the caller asked to avoid, including ones
    created across the padding/coding junction. Padding carries no translation, so
    any base in it can be rewritten freely.

    Parameters
    ----------
    region : str
        padding_5_prime + coding_sequence + padding_3_prime. This is what the
        ConstraintSet treats as the coding sequence, so its transcribed 5'/3'
        contexts are applied around it as usual.
    mutable_ranges : list of (start, stop)
        Half-open spans of `region` that may be rewritten, in region coordinates.
        Anything outside them (the coding sequence) is left untouched.
    constraint_set : ConstraintSet
        Constraints to satisfy.
    forbidden : set of str
        Codons the padding must not contain (start/stop, per the caller's flags).
    max_passes : int
        Bound on repair rounds before giving up.

    Returns
    -------
    (str, list)
        The repaired region and any violations still touching the padding.
    """
    bases = list(region)
    offset = constraint_set.coding_offset

    mutable = set()
    for start, stop in mutable_ranges:
        mutable.update(range(max(0, start), min(len(bases), stop)))

    def touching_padding(sequence):
        """Violations that overlap at least one rewritable base."""
        found = []
        for violation in constraint_set.find_violations(sequence):
            positions = [index
                         for index in range(max(0, violation.start - offset),
                                            min(len(bases), violation.end - offset))
                         if index in mutable]
            if positions:
                found.append((violation, positions))
        return found

    if not mutable:
        return region, [violation for violation, _ in touching_padding(region)]

    for _ in range(max_passes):
        current = ''.join(bases)
        actionable = touching_padding(current)
        if not actionable:
            return current, []
        # Fix the first offender; re-scanning each pass keeps the repair honest
        # about violations a rewrite creates somewhere else.
        _, positions = actionable[0]
        if not _mutate_padding_base(bases, positions, mutable, forbidden):
            break

    current = ''.join(bases)
    return current, [violation for violation, _ in touching_padding(current)]


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
    
    # Codons the padding is not allowed to contain, as a single set so the repair
    # pass can fix start and stop codons together instead of in two forward scans
    # that can each undo the other's work.
    forbidden = set()
    if avoid_start_codon:
        forbidden.update(start_codons)
    if avoid_stop_codon:
        forbidden.update(stop_codons)

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
        
        # Rewrite any forbidden codons in place, re-scanning until stable.
        padding_list = list(padding)
        if not _repair_forbidden_codons(padding_list, forbidden):
            continue
        padding = ''.join(padding_list)

        # Repairs can shift GC content, so recheck it.
        gc_content = calc_gc_content(padding)
        if not (GC_content_range[0] - tolerance <= gc_content <= GC_content_range[1] + tolerance):
            continue

        # Check for sequences to avoid - can't easily modify these in place
        # so we'll still reject sequences containing them
        if any(seq in padding for seq in avoid_adding):
            continue

        # Belt and braces: the repair pass guarantees this, but a silent regression
        # here would ship start/stop codons into every construct.
        if any(codon in padding for codon in forbidden):
            continue

        # If all checks pass, return the padding sequence
        return padding

    # If no valid padding sequence is found after max attempts
    raise ValueError(f"Could not generate a suitable padding sequence after {max_attempts} attempts. "
                    f"Consider relaxing constraints or increasing the tolerance.")


def add_padding(current_gc_content,
                padding_length,
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
    current_length : int
        The length of the sequence the padding is being added to.
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



def check_nucleotide_percent_similarity(sequences, return_max=True):
    '''
    A function that takes in a list of nucleotide sequences and does an all by all 
    comparison to find the most similar sequences.

    Parameters
    ----------
    sequences : list
        List of nucleotide sequences to compare
    return_max : bool
        If True, returns a tuple containing:
            - Array of maximum similarities for each sequence
            - Array of indices of the most similar sequences
        If False, returns the full similarity matrix

    Returns
    -------
    tuple or numpy.ndarray
        If return_max=True: (max_similarities, most_similar_indices)
        If return_max=False: Full similarity matrix
    '''
    # Validate input sequences
    if not sequences:
        raise ValueError("Empty sequence list provided")

    if return_max and len(sequences) < 2:
        raise ValueError("At least two sequences are required to find the most similar pair")

    # Convert sequences to uppercase
    sequences = [seq.upper() for seq in sequences]

    # Check if all sequences have the same length
    lengths = [len(seq) for seq in sequences]
    if len(set(lengths)) > 1:
        raise ValueError("All sequences must have the same length")
    if lengths[0] == 0:
        raise ValueError("Sequences must not be empty")

    # Check for valid nucleotides. N is permitted and scored as a mismatch, which
    # is what the encoding below implements -- rejecting it here made that branch
    # unreachable.
    valid_nucleotides = set('ATCGN')
    for seq in sequences:
        invalid = set(seq) - valid_nucleotides
        if invalid:
            raise ValueError(f"Invalid nucleotide(s) {', '.join(sorted(invalid))} "
                             f"found in sequence: {seq}")

    # Convert sequences to numpy arrays
    seq_arrays = np.array([list(seq) for seq in sequences])
    n_sequences, len_sequences = seq_arrays.shape
    encoded_array = np.zeros((n_sequences, len_sequences), dtype=np.int8)
    
    # Encode nucleotides (A=0, C=1, G=2, T=3, N=4)
    encoded_array[seq_arrays == 'C'] = 1
    encoded_array[seq_arrays == 'G'] = 2
    encoded_array[seq_arrays == 'T'] = 3
    encoded_array[seq_arrays == 'N'] = 4

    # Calculate similarity matrix
    similarity_matrix = np.zeros((n_sequences, n_sequences))
    for i in range(n_sequences):
        # Count matching positions, treating N as a mismatch
        matches = (encoded_array == encoded_array[i]) & (encoded_array != 4)
        similarity_matrix[i, :] = np.sum(matches, axis=1) / len_sequences

    if return_max:
        # Create mask to ignore self-comparisons
        mask = ~np.eye(n_sequences, dtype=bool)
        # Get maximum similarities and corresponding indices
        max_similarities = np.max(similarity_matrix * mask, axis=1)
        most_similar_indices = np.argmax(similarity_matrix * mask, axis=1)
        return max_similarities, most_similar_indices
    else:
        return similarity_matrix
    

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
    # Not implemented. This used to silently return None, which is indistinguishable
    # from "generated zero barcodes" at the call site.
    raise NotImplementedError(
        "generate_barcode_sequences is not implemented yet. Build barcodes with "
        "create_padding_sequence() and filter them with "
        "check_nucleotide_percent_similarity() in the meantime.")