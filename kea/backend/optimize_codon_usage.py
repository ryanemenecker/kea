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
                         codon_table_obj,
                         n_iter=10000,
                         fine_tuning_iterations=500,
                         return_best=True,
                         early_stop_threshold=0.95,
                         show_progress_bar=True):
    '''
    Optimize codon usage for a given amino acid sequence.
    
    Parameters
    ----------
    amino_acid_sequence : str
        Amino acid sequence to optimize.
    codon_table_obj : CodonTable
        CodonTable object with codon usage frequencies.
    n_iter : int
        Number of iterations for optimization.
    fine_tuning_iterations : int
        Number of iterations for fine-tuning GC content after main optimization.
    return_best : bool
        If True, return the best sequence even if it's outside the GC range.
    early_stop_threshold : float
        Stop optimization if a sequence with this proportion of the theoretical
        maximum score is found.
    show_progress_bar : bool
        If True, show a progress bar for the optimization process.
    Returns
    -------
    dict
        dictionary with the amino acid sequence as the key and the optimized DNA sequence as the value.
    '''    
    # Generate candidates - we'll keep track of several promising sequences
    best_score = float('-inf')
    best_sequence = None
    
    # Track the best sequence that has GC content in range
    best_in_range_score = float('-inf')
    best_in_range_sequence = None
    # If we are tracking progress, initialize the progress bar
    if show_progress_bar:
        pbar_inner = tqdm(total=n_iter, position=1, desc='Current sequence optimization', leave=True)
    
    # main loop.
    for _ in range(n_iter):
        sequence = _build_sequence(amino_acid_sequence, codon_table_obj.aa_weights)
        
        # Calculate GC content
        gc_content = _calculate_gc_content(sequence)
            
        # get score
        score = _score_sequence(sequence, 
                                codon_table_obj.gc_range, 
                                codon_table_obj.codon_table, 
                                codon_table_obj.gc_weight)
        
        # Track the overall best sequence
        if score > best_score:
            best_score = score
            best_sequence = sequence
        
        # Also track the best sequence within GC range
        if codon_table_obj.gc_range[0] <= gc_content <= codon_table_obj.gc_range[1]:
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
        codon_table_obj.target_gc, 
        codon_table_obj.codon_table,
        iterations=fine_tuning_iterations
    )
    
    # Check if final GC is within range
    final_gc = _calculate_gc_content(final_sequence)
    
    # If GC content not within range and return_best is False, raise error
    if final_gc < codon_table_obj.gc_range[0] or final_gc > codon_table_obj.gc_range[1]:
        if not return_best:
            raise ValueError(f"No sequence found within GC range. Best achieved: {final_gc:.3f}")
    
    # Close progress bar
    if show_progress_bar:
        pbar_inner.close()
    
    # return final sequence
    return final_sequence

