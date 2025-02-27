'''
Code to optimize codon usage. The optimization problem takes two factors into account:
1) The codon usage of the host organism.
2) A range of acceptable GC content.
'''

import random
import numpy as np
from tqdm import tqdm
from kea.data.aa_codon_conversions import aa_to_codons, codons_to_aa
from kea.data.codon_tables import s288c_codon_usage



def _build_sequence(amino_acid_sequence, aa_weights):
    """Helper function to build a DNA sequence based on weighted probabilities."""
    dna_sequence = ''
    for aa in amino_acid_sequence:
        codons, weights = aa_weights[aa]
        dna_sequence += np.random.choice(codons, p=weights)
    return dna_sequence

def _calculate_gc_content(sequence):
    """Calculate GC content of a DNA sequence."""
    if not sequence:
        return 0
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count / len(sequence)

def _score_sequence(sequence, gc_range, codon_frequency_table, gc_priority=1.0):
    """Score a sequence based on codon usage and GC content."""
    gc_content = _calculate_gc_content(sequence)
    
    # Stronger GC content penalty with exponential scaling
    gc_penalty = 0
    if gc_content < gc_range[0]:
        gc_penalty = (gc_range[0] - gc_content) * 100 * gc_priority
    elif gc_content > gc_range[1]:
        gc_penalty = (gc_content - gc_range[1]) * 100 * gc_priority
    
    # Codon usage score
    usage_score = 0
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        aa = codons_to_aa[codon]
        usage_score += codon_frequency_table[aa].get(codon, 0.01)
    
    return usage_score - gc_penalty

def _adjust_gc_content(sequence, target_gc, codon_frequency_table, iterations=100):
    """Fine-tune sequence GC content by swapping synonymous codons."""
    current_gc = _calculate_gc_content(sequence)
    best_sequence = sequence
    best_gc_diff = abs(current_gc - target_gc)
    
    # Break sequence into codons
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    
    for _ in range(iterations):
        # Choose random position to modify
        pos = random.randint(0, len(codons) - 1)
        codon = codons[pos]
        aa = codons_to_aa[codon]
        
        # Get alternative codons for this amino acid
        alt_codons = [c for c in aa_to_codons[aa] if c != codon]
        if not alt_codons:
            continue
        
        # Choose replacement codon based on GC direction needed
        if current_gc < target_gc:
            # Need more GC: sort by descending GC content
            alt_codons.sort(key=lambda c: (c.count('G') + c.count('C')), reverse=True)
        else:
            # Need less GC: sort by ascending GC content
            alt_codons.sort(key=lambda c: (c.count('G') + c.count('C')))
        
        # Consider codon usage too (take top 2 candidates)
        if len(alt_codons) > 2:
            alt_codons = sorted(alt_codons[:2], 
                              key=lambda c: codon_frequency_table[aa].get(c, 0.01),
                              reverse=True)
        
        # Try the best alternative
        new_codons = codons.copy()
        new_codons[pos] = alt_codons[0]
        new_sequence = ''.join(new_codons)
        new_gc = _calculate_gc_content(new_sequence)
        
        # Keep if it's closer to target
        if abs(new_gc - target_gc) < best_gc_diff:
            best_sequence = new_sequence
            best_gc_diff = abs(new_gc - target_gc)
            codons = new_codons
            current_gc = new_gc
    
    return best_sequence

def optimize_codon_usage(amino_acid_sequence, 
                         codon_frequency_table,
                         gc_range=(0.35, 0.45), 
                         n_iter=10000,
                         usage_weight=1.0,
                         gc_weight=1.0,
                         return_best=True,
                         early_stop_threshold=0.95,
                         fine_tuning_iterations=500,
                         minimum_codon_probability=0.06,
                         show_progress_bar=True):
    '''
    Optimize codon usage for a given amino acid sequence.
    
    Parameters
    ----------
    amino_acid_sequence : str
        Amino acid sequence to optimize.
    codon_frequency_table : dict
        Codon frequency table.
        Format is {amino_acid: {codon: frequency}}.
    gc_range : tuple
        Range of acceptable GC content.
    n_iter : int
        Number of iterations for optimization.
    usage_weight : float
        Weight for codon usage optimization (higher values prioritize frequent codons).
    gc_weight : float
        Weight for GC content optimization (higher values strictly enforce GC range).
    local_optimization : bool
        If True, perform local optimization after the global search.
    local_iterations : int
        Number of iterations for local optimization.
    return_best : bool
        If True, return the best sequence even if it's outside the GC range.
    early_stop_threshold : float
        Stop optimization if a sequence with this proportion of the theoretical
        maximum score is found.
    fine_tuning_iterations : int
        Number of iterations for fine-tuning GC content after main optimization.
    minimum_codon_probability : float
        Minimum probability for a codon to be considered in the optimization.
        Values below this will not be used. This is to prevent
        very low-frequency codons from being used.
    show_progress_bar : bool
        If True, show a progress bar for the optimization process.
    Returns
    -------
    str
        Optimized DNA sequence.
    dict
        Statistics about the optimization process.
    '''
    # Input validation
    if not isinstance(amino_acid_sequence, str):
        raise TypeError("amino_acid_sequence must be a string")
    if not all(aa in aa_to_codons for aa in amino_acid_sequence):
        raise ValueError("Invalid amino acid in sequence")
    if not 0 <= gc_range[0] <= gc_range[1] <= 1:
        raise ValueError("Invalid GC range")
    
    # see if there are any codons that are below the minimum probability
    del_vals={}
    for aa in codon_frequency_table:
        for codon in codon_frequency_table[aa]:
            if codon_frequency_table[aa][codon] < minimum_codon_probability:
                if aa not in del_vals:
                    del_vals[aa]=[]
                del_vals[aa].append(codon)
    # remove those codons from the table
    for aa in del_vals:
        for codon in del_vals[aa]:
            del codon_frequency_table[aa][codon]
            # remove the codon from aa_to_codons
            aa_to_codons[aa].remove(codon)

    # Pre-calculate codon GC content
    codon_gc_content = {codon: (codon.count('G') + codon.count('C'))/3 
                       for aa in aa_to_codons for codon in aa_to_codons[aa]}
    
    # if we removed codons from the table, we need to recalculate the probabilities
    if del_vals != {}:
        for aa in codon_frequency_table:
            # check that there is at least one value
            if len(codon_frequency_table[aa]) == 0:
                raise ValueError(f"No codons left for amino acid {aa}. Need to increase minimum_codon_probability.")
            total = sum(codon_frequency_table[aa].values())
            for codon in codon_frequency_table[aa]:
                codon_frequency_table[aa][codon] = round(codon_frequency_table[aa][codon]/total, 6)


    # Pre-calculate target GC
    target_gc = sum(gc_range) / 2
    
    # Pre-calculate weights for each amino acid
    aa_weights = {}
    for aa in aa_to_codons:
        possible_codons = aa_to_codons[aa]
        usage_weights = [codon_frequency_table[aa].get(codon, 0.01) for codon in possible_codons]
        
        # Stronger GC weighting function
        gc_weights = []
        for codon in possible_codons:
            codon_gc = codon_gc_content[codon]
            # Sharper penalty for being outside target range
            if codon_gc < gc_range[0]:
                gc_factor = 1 - 2 * (gc_range[0] - codon_gc)
            elif codon_gc > gc_range[1]:
                gc_factor = 1 - 2 * (codon_gc - gc_range[1])
            else:
                # Bonus for being in range with peak at target
                gc_factor = 1 + (1 - abs(codon_gc - target_gc) * 2)
            
            gc_weights.append(max(0.01, gc_factor))
        
        # Combine weights with stronger GC emphasis
        combined_weights = np.array([uw**usage_weight * gw**(gc_weight*2) for uw, gw in zip(usage_weights, gc_weights)])
        aa_weights[aa] = (possible_codons, combined_weights/combined_weights.sum())
    
    # Generate candidates - we'll keep track of several promising sequences
    best_score = float('-inf')
    best_sequence = None
    
    # Track the best sequence that has GC content in range
    best_in_range_score = float('-inf')
    best_in_range_sequence = None
    
    # Track optimization stats
    stats = {
        'iterations': 0,
        'attempts_in_gc_range': 0,
        'final_gc': None,
        'target_gc': target_gc
    }
    if show_progress_bar:
        pbar_inner = tqdm(total=n_iter, position=1, desc='Current sequence optimization', leave=False)
    for i in range(n_iter):
        stats['iterations'] += 1
        sequence = _build_sequence(amino_acid_sequence, aa_weights)
        
        # Calculate GC content
        gc_content = _calculate_gc_content(sequence)
            
        # get score
        score = _score_sequence(sequence, gc_range, codon_frequency_table, gc_weight)
        
        # Track the overall best sequence
        if score > best_score:
            best_score = score
            best_sequence = sequence
        
        # Also track the best sequence within GC range
        if gc_range[0] <= gc_content <= gc_range[1]:
            stats['attempts_in_gc_range'] += 1
            if score > best_in_range_score:
                best_in_range_score = score
                best_in_range_sequence = sequence
                
                # Early stopping only if we found a good sequence in GC range
                if score >= len(amino_acid_sequence) * early_stop_threshold:
                    break
        if show_progress_bar:
            pbar_inner.update(1)
    
    # Phase 2: Fine-tune the best sequence for GC content
    final_sequence = best_sequence
    
    # If we found any sequence in range, prefer that
    if best_in_range_sequence and best_in_range_score > best_score * 0.85:
        final_sequence = best_in_range_sequence
    
    # Always apply fine-tuning to improve GC content
    final_sequence = _adjust_gc_content(
        final_sequence, 
        target_gc, 
        codon_frequency_table,
        iterations=fine_tuning_iterations
    )
    
    # Check if final GC is within range
    final_gc = _calculate_gc_content(final_sequence)
    stats['final_gc'] = final_gc
    
    # If GC content not within range and return_best is False, raise error
    if final_gc < gc_range[0] or final_gc > gc_range[1]:
        if not return_best:
            raise ValueError(f"No sequence found within GC range. Best achieved: {final_gc:.3f}")
    
    if show_progress_bar:
        pbar_inner.close()

    return final_sequence, stats

def optimize_codon_usage_for_library(
        list_of_sequences,
        codon_frequency_table,
        gc_range=(0.35, 0.45),
        n_iter=10000,
        usage_weight=1.0,
        gc_weight=1.0,
        return_best=True,
        early_stop_threshold=0.95,
        fine_tuning_iterations=500,
        minimum_codon_probability=0.06,
        return_statistics=False,
        show_progress_bar=True): 
    '''
    Optimize codon usage for a list of amino acid sequences.
    Parameters
    ----------
    list_of_sequences : list
        List of amino acid sequences to optimize.
    codon_frequency_table : dict
        Codon frequency table.
        Format is {amino_acid: {codon: frequency}}.
    gc_range : tuple
        Range of acceptable GC content.
    n_iter : int
        Number of iterations for optimization.
    usage_weight : float
        Weight for codon usage optimization (higher values prioritize frequent codons).
    gc_weight : float
        Weight for GC content optimization (higher values strictly enforce GC range).
    return_best : bool
        If True, return the best sequence even if it's outside the GC range.
    early_stop_threshold : float
        Stop optimization if a sequence with this proportion of the theoretical
        maximum score is found.
    fine_tuning_iterations : int
        Number of iterations for fine-tuning GC content after main optimization.
    minimum_codon_probability : float
        Minimum probability for a codon to be considered in the optimization.
        Values below this will not be used. This is to prevent
        very low-frequency codons from being used.
    return_statistics : bool
        If True, return statistics about the optimization process.
        If False, return only the optimized sequences.
    show_progress_bar : bool
        If True, show a progress bar for the optimization process.

    Returns
    -------
    list
        List of optimized DNA sequences.
    dict
        Statistics about the optimization process.
    '''
    optimized_sequences = []
    stats_list = []
    # Shoutout to Alex Holehouse for the progress bar magic. 
    if show_progress_bar:
        pbar=tqdm(total=len(list_of_sequences), position=0, desc='Progress through sequences', leave=True)
    for sequence in tqdm(list_of_sequences):
        optimized_sequence, stats = optimize_codon_usage(
            sequence,
            codon_frequency_table,
            gc_range=gc_range,
            n_iter=n_iter,
            usage_weight=usage_weight,
            gc_weight=gc_weight,
            return_best=return_best,
            early_stop_threshold=early_stop_threshold,
            fine_tuning_iterations=fine_tuning_iterations,
            minimum_codon_probability=minimum_codon_probability,
            show_progress_bar=show_progress_bar
        )
        optimized_sequences.append(optimized_sequence)
        stats_list.append(stats)
        # update if we are using a progress bar
        if show_progress_bar:
            pbar.update(1)
    # make sure to close progress bar
    if show_progress_bar:
        pbar.close()
    if return_statistics:
        return optimized_sequences, stats_list
    else:
        return optimized_sequences

