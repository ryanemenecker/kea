'''
Code to optimize codon usage. The optimization problem takes two factors into account:
1) The codon usage of the host organism.
2) A range of acceptable GC content.
'''

import random
import numpy as np
import functools
from tqdm import tqdm
from kea.data.aa_codon_conversions import aa_to_codons, codons_to_aa
from kea.data.codon_tables import s288c_codon_usage


# Add a cache decorator for expensive functions
def memoize(func):
    """Simple memoization decorator for pure functions"""
    cache = {}
    @functools.wraps(func)
    def wrapper(*args):
        if args not in cache:
            cache[args] = func(*args)
        return cache[args]
    return wrapper


def _build_sequence(amino_acid_sequence, aa_weights):
    """
    Helper function to build a DNA sequence based on weighted probabilities.
    Optimized version with better performance for long sequences.
    """
    # Use a list to store codons - much more efficient than string concatenation
    dna_codons = []
    
    # Group amino acids by type to minimize lookups
    aa_groups = {}
    
    # Count amino acid frequencies to determine if prebatching is worthwhile
    aa_counts = {}
    for aa in amino_acid_sequence:
        aa_counts[aa] = aa_counts.get(aa, 0) + 1
    
    # Pre-generate codons for common amino acids (appearing at least 5 times)
    # This reduces the number of random choice operations
    batch_size = 10  # Generate codons in batches of 10
    for aa, count in aa_counts.items():
        if count >= 5:
            codons, weights = aa_weights[aa]
            # Generate multiple codons at once in batches
            prebatch_size = ((count // batch_size) + 1) * batch_size
            prebatched_codons = np.random.choice(
                codons, 
                size=prebatch_size, 
                p=weights,
                replace=True
            )
            aa_groups[aa] = (prebatched_codons, 0)  # Store with an index counter
    
    # Build the sequence
    for aa in amino_acid_sequence:
        if aa in aa_groups:
            # Use a pre-generated codon from our batch
            prebatched_codons, index = aa_groups[aa]
            dna_codons.append(prebatched_codons[index])
            # Update the index for next time
            aa_groups[aa] = (prebatched_codons, index + 1)
        else:
            # Generate a codon on-demand for less common amino acids
            codons, weights = aa_weights[aa]
            dna_codons.append(np.random.choice(codons, p=weights))
    
    # Join all codons at once - much more efficient than string concatenation in a loop
    return ''.join(dna_codons)


def _calculate_gc_content(sequence):
    """Calculate GC content of a DNA sequence more efficiently."""
    if not sequence:
        return 0
    # Count both G and C in a single pass
    gc_count = sum(1 for base in sequence if base in 'GCgc')
    return gc_count / len(sequence)


def _score_sequence(sequence, gc_range, codon_frequency_table, gc_priority=1.0, gc_tolerance=0.025):
    """Score a sequence based on codon usage and GC content with balanced weighting.
    
    Parameters:
    - sequence: DNA sequence to score
    - gc_range: Tuple of (min_gc, max_gc) acceptable range
    - codon_frequency_table: Dictionary of codon frequencies
    - gc_priority: Relative weight of GC content vs codon usage (0.0-2.0, capped if outside range)
                  1.0 = equal weight, <1 favors codon usage, >1 favors GC content
    - gc_tolerance: Allowable deviation from GC range (default 0.025 or 2.5%)
    
    Returns:
    - Combined normalized score between 0-1
    """
    if not sequence:
        return 0
    
    # Cap gc_priority to valid range (0-2)
    gc_priority = min(2.0, max(0.0, gc_priority))

    # Number of codons
    n_codons = len(sequence) // 3
    if n_codons == 0:
        return 0
        
    # Calculate GC content
    gc_content = _calculate_gc_content(sequence)
    
    # Expanded GC range with tolerance
    tolerant_min = max(0, gc_range[0] - gc_tolerance)
    tolerant_max = min(1, gc_range[1] + gc_tolerance)
    
    # GC content score calculation with tolerance
    if tolerant_min <= gc_content <= tolerant_max:
        # If within expanded range but outside original range, apply penalty
        if gc_content < gc_range[0]:
            gc_score = 1.0 - ((gc_range[0] - gc_content) / gc_tolerance)
        elif gc_content > gc_range[1]:
            gc_score = 1.0 - ((gc_content - gc_range[1]) / gc_tolerance)
        else:
            gc_score = 1.0
    else:
        if gc_content < tolerant_min:
            distance = (tolerant_min - gc_content) / tolerant_min
        else:
            distance = (gc_content - tolerant_max) / (1 - tolerant_max)
        gc_score = max(0, 1 - min(1, distance * 2))

    # Codon usage score with normalization against best possible codon
    usage_score = 0
    for i in range(0, len(sequence), 3):
        if i+3 <= len(sequence):  # Ensure complete codon
            codon = sequence[i:i+3]
            aa = codons_to_aa[codon]
            
            # Get this codon's frequency
            codon_freq = codon_frequency_table[aa].get(codon, 0.01)
            
            # Find max frequency for this amino acid (best possible codon)
            max_freq = max(codon_frequency_table[aa].values())
            
            # Normalize against the best possible codon
            if max_freq > 0:
                normalized_score = codon_freq / max_freq
            else:
                normalized_score = 0
                
            usage_score += normalized_score
    
    # Normalize usage score to 0-1 scale
    usage_score = usage_score / n_codons
    
    # Final weighted score
    return ((usage_score * (2 - gc_priority) + gc_score * gc_priority)) / 2


def _adjust_gc_content(sequence, target_gc, codon_frequency_table, iterations=100):
    """Fine-tune sequence GC content by swapping synonymous codons.
    
    Prioritizes codon changes that have minimal impact on codon usage frequency
    while moving toward the target GC content.
    """
    current_gc = _calculate_gc_content(sequence)
    best_sequence = sequence
    best_gc_diff = abs(current_gc - target_gc)
    
    # Break sequence into codons (use list comprehension instead of full list creation)
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    
    for _ in range(iterations):
        # Optimization: Instead of evaluating all positions each time, 
        # sample a subset for very long sequences
        positions_to_evaluate = range(len(codons))
        if len(codons) > 1000:
            positions_to_evaluate = random.sample(range(len(codons)), 1000)
            
        best_pos = -1
        best_replacement = None
        best_score = float('-inf')
        
        for pos in positions_to_evaluate:
            codon = codons[pos]
            aa = codons_to_aa[codon]
            current_freq = codon_frequency_table[aa].get(codon, 0.01)
            
            # Get alternative codons for this amino acid
            alt_codons = [c for c in aa_to_codons[aa] if c != codon]
            if not alt_codons:
                continue
                
            # Calculate improvement metrics for each alternative
            for alt_codon in alt_codons:
                # Calculate GC change
                codon_gc = (codon.count('G') + codon.count('C')) / 3
                alt_gc = (alt_codon.count('G') + alt_codon.count('C')) / 3
                gc_delta = alt_gc - codon_gc
                
                # Will this change move us in the right direction?
                if (current_gc < target_gc and gc_delta <= 0) or (current_gc > target_gc and gc_delta >= 0):
                    continue
                
                # Calculate frequency change (cost)
                alt_freq = codon_frequency_table[aa].get(alt_codon, 0.01)
                freq_cost = max(0, current_freq - alt_freq)  # Only consider frequency decreases as cost
                
                # Avoid zero division
                if freq_cost < 0.001:  # Almost no cost
                    efficiency = abs(gc_delta) * 1000  # Very high efficiency
                else:
                    efficiency = abs(gc_delta) / freq_cost  # GC benefit per unit of frequency cost
                
                # Calculate impact on overall GC
                new_codons = codons.copy()
                new_codons[pos] = alt_codon
                new_sequence = ''.join(new_codons)
                new_gc = _calculate_gc_content(new_sequence)
                gc_improvement = abs(current_gc - target_gc) - abs(new_gc - target_gc)
                
                # Skip if moving away from target
                if gc_improvement <= 0:
                    continue
                
                # Score combines GC improvement and codon efficiency
                # Higher score = better change
                score = gc_improvement * (1 + efficiency)
                
                if score > best_score:
                    best_score = score
                    best_pos = pos
                    best_replacement = alt_codon
        
        # Apply the best change if one was found
        if best_pos >= 0:
            codons[best_pos] = best_replacement
            new_sequence = ''.join(codons)
            new_gc = _calculate_gc_content(new_sequence)
            
            # Update tracking variables
            if abs(new_gc - target_gc) < best_gc_diff:
                best_sequence = new_sequence
                best_gc_diff = abs(new_gc - target_gc)
            
            current_gc = new_gc
        else:
            # No beneficial changes found, break the loop
            break
    
    return best_sequence


def _improve_sequence(sequence, codon_table_obj, mutation_rate=0.3):
    """
    Improve an existing sequence by targeted mutations of suboptimal codons.
    Optimized version with better performance characteristics.
    
    Parameters:
    -----------
    sequence : str
        DNA sequence to improve
    codon_table_obj : CodonTable
        Object containing codon usage frequencies and GC preferences
    mutation_rate : float
        Probability of attempting to mutate each codon (0-1)
        
    Returns:
    --------
    str
        Improved DNA sequence if improvements were found, otherwise original sequence
    """
    # Break sequence into codons - use array for faster manipulation
    codons = np.array([sequence[i:i+3] for i in range(0, len(sequence), 3)])
    n_codons = len(codons)
    
    if n_codons == 0:
        return sequence
    
    # Track if any changes were made
    improved = False
    
    # Current sequence metrics - calculate only once
    current_gc = _calculate_gc_content(sequence)
    current_score = _score_sequence(sequence, codon_table_obj.gc_range, 
                                  codon_table_obj.codon_table, codon_table_obj.gc_weight,
                                  getattr(codon_table_obj, 'gc_tolerance', 0.025))
    
    # Precompute GC content of each codon position
    codon_gc_contents = np.array([(c.count('G') + c.count('C')) / 3 for c in codons])
    
    # Direct sampling instead of shuffling - more efficient for large sequences
    # Sample ~mutation_rate fraction of positions for evaluation
    positions_to_check = np.random.choice(
        n_codons, 
        size=max(1, int(n_codons * mutation_rate)), 
        replace=False
    )
    
    # Pre-compute amino acids and their alternative codons for positions we'll check
    aa_for_positions = {}
    alt_codons_for_aa = {}
    
    for pos in positions_to_check:
        codon = codons[pos]
        try:
            aa = codons_to_aa[codon]
            aa_for_positions[pos] = aa
            
            # Cache alternative codons for this amino acid
            if aa not in alt_codons_for_aa:
                alt_codons_for_aa[aa] = [c for c in aa_to_codons[aa] if c != codon]
        except KeyError:
            # Skip invalid codons
            continue
    
    # Process positions in batches to improve memory locality
    batch_size = min(100, len(positions_to_check))
    
    for batch_start in range(0, len(positions_to_check), batch_size):
        batch_end = min(batch_start + batch_size, len(positions_to_check))
        batch_positions = positions_to_check[batch_start:batch_end]
        
        for pos in batch_positions:
            if pos not in aa_for_positions:
                continue
                
            aa = aa_for_positions[pos]
            codon = codons[pos]
            alternatives = alt_codons_for_aa.get(aa, [])
            
            if not alternatives:
                continue
            
            # Calculate the GC impact of replacing this codon without generating new sequences
            codon_gc = codon_gc_contents[pos]
            
            # Score this position's alternatives more efficiently
            best_alt_codon = None
            best_alt_score = current_score
            
            for alt_codon in alternatives:
                # Quickly calculate GC delta
                alt_gc = (alt_codon.count('G') + alt_codon.count('C')) / 3
                gc_delta = (alt_gc - codon_gc) / n_codons  # Impact on overall GC
                
                # Estimate new GC content directly
                new_gc_content = current_gc + gc_delta
                
                # Get codon frequency info
                alt_freq = codon_table_obj.codon_table[aa].get(alt_codon, 0.01)
                current_freq = codon_table_obj.codon_table[aa].get(codon, 0.01)
                
                # Quick estimate of score change
                # If it's unlikely to improve, skip detailed evaluation
                freq_improvement = alt_freq - current_freq
                
                gc_range = codon_table_obj.gc_range
                gc_tolerance = getattr(codon_table_obj, 'gc_tolerance', 0.025)
                
                # See if this change would improve GC content fit
                gc_improved = False
                if new_gc_content >= gc_range[0] - gc_tolerance and new_gc_content <= gc_range[1] + gc_tolerance:
                    if abs(new_gc_content - codon_table_obj.target_gc) < abs(current_gc - codon_table_obj.target_gc):
                        gc_improved = True
                
                # Skip detailed evaluation if neither frequency nor GC is likely to improve
                if freq_improvement <= 0 and not gc_improved:
                    continue
                
                # Generate and score the new sequence only for promising candidates
                new_codons = codons.copy()
                new_codons[pos] = alt_codon
                new_sequence = ''.join(new_codons)
                
                new_score = _score_sequence(
                    new_sequence, 
                    gc_range,
                    codon_table_obj.codon_table, 
                    codon_table_obj.gc_weight,
                    gc_tolerance
                )
                
                if new_score > best_alt_score:
                    best_alt_score = new_score
                    best_alt_codon = alt_codon
            
            # Apply change if better alternative found
            if best_alt_codon is not None and best_alt_score > current_score:
                codons[pos] = best_alt_codon
                improved = True
                # Update current score and GC content for next evaluations
                current_score = best_alt_score
                current_gc = _calculate_gc_content(''.join(codons))
    
    # Return improved sequence if changes were made
    if improved:
        return ''.join(codons)
    else:
        return sequence


def optimize_codon_usage(amino_acid_sequence, 
                         codon_table_obj,
                         n_iter=10000,
                         fine_tuning_iterations=500,
                         return_best=True,
                         early_stop_threshold=0.95,
                         show_progress_bar=True,
                         gc_tolerance=0.025):
    '''
    Optimize codon usage for a given amino acid sequence using a hybrid approach.
    
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
    gc_tolerance : float, default=0.025
        Allowable deviation from GC range (as a fraction, e.g., 0.025 = ±2.5%)
    Returns
    -------
    str
        Optimized DNA sequence.
    '''    
    # Input validation
    if not amino_acid_sequence:
        raise ValueError("Amino acid sequence cannot be empty")
    
    if not hasattr(codon_table_obj, 'codon_table') or not hasattr(codon_table_obj, 'gc_range'):
        raise ValueError("Invalid codon table object provided")
        
    if not (0 <= early_stop_threshold <= 1):
        raise ValueError("early_stop_threshold must be between 0 and 1")
    
    # Validate amino acid characters
    valid_aas = set(aa_to_codons.keys())
    invalid_aas = [aa for aa in amino_acid_sequence if aa not in valid_aas]
    if invalid_aas:
        raise ValueError(f"Invalid amino acid(s) in sequence: {', '.join(invalid_aas)}")
    
    # Phase parameters
    exploration_fraction = 0.3    # 30% of iterations for exploration
    refinement_fraction = 0.7     # 70% of iterations for refinement
    top_candidates = 10           # Number of candidates to refine
    
    # Initialize tracking variables
    best_score = float('-inf')
    best_sequence = None
    best_in_range_score = float('-inf')
    best_in_range_sequence = None
    
    # Setup progress tracking
    if show_progress_bar:
        from tqdm import tqdm
        pbar = tqdm(total=n_iter, desc='Codon optimization', leave=True)
    
    # Phase 1: Generate diverse initial candidates through random sampling
    exploration_iterations = int(n_iter * exploration_fraction)
    candidates = []
    
    for _ in range(exploration_iterations):
        # Generate a random sequence
        sequence = _build_sequence(amino_acid_sequence, codon_table_obj.aa_weights)
        
        # Calculate metrics
        gc_content = _calculate_gc_content(sequence)
        score = _score_sequence(sequence, 
                               codon_table_obj.gc_range, 
                               codon_table_obj.codon_table, 
                               codon_table_obj.gc_weight,
                               gc_tolerance)
        
        # Track candidates
        candidates.append((sequence, score, gc_content))
        
        # Track best overall and in-range sequences
        if score > best_score:
            best_score = score
            best_sequence = sequence
        
        if codon_table_obj.gc_range[0] <= gc_content <= codon_table_obj.gc_range[1]:
            if score > best_in_range_score:
                best_in_range_score = score
                best_in_range_sequence = sequence
                # Early stopping if we found an excellent sequence
                if score >= early_stop_threshold:
                    if show_progress_bar:
                        pbar.update(exploration_iterations)
                    break
        
        if show_progress_bar:
            pbar.update(1)
    
    # Select top candidates for refinement, ensuring diversity
    # Sort by score but also include some sequences with good GC content even if score is lower
    score_candidates = sorted(candidates, key=lambda x: x[1], reverse=True)[:int(top_candidates * 0.7)]
    
    # Get additional candidates with good GC content
    gc_candidates = [c for c in candidates if abs(c[2] - codon_table_obj.target_gc) < 0.05]
    gc_candidates = sorted(gc_candidates, key=lambda x: x[1], reverse=True)
    gc_candidates = [c for c in gc_candidates if c not in score_candidates][:int(top_candidates * 0.3)]
    
    # Combine candidates
    refinement_candidates = score_candidates + gc_candidates
    refinement_candidates = [c[0] for c in refinement_candidates][:top_candidates]
    
    # Make sure we have at least one candidate
    if not refinement_candidates and best_sequence:
        refinement_candidates = [best_sequence]
    
    # Phase 2: Iterative improvement of top candidates
    refinement_iterations = n_iter - exploration_iterations
    iterations_per_candidate = refinement_iterations // len(refinement_candidates) if refinement_candidates else 0
    remaining_iters = refinement_iterations
    
    for base_sequence in refinement_candidates:
        # Limit iterations for this candidate
        candidate_iters = min(iterations_per_candidate, remaining_iters)
        remaining_iters -= candidate_iters
        
        # Starting point
        current = base_sequence
        
        # Evolve this candidate through multiple iterations
        for _ in range(candidate_iters):
            # Try improving the sequence
            improved = _improve_sequence(current, codon_table_obj)
            
            # Update if improved
            if improved != current:
                current = improved
                
                # Recalculate metrics
                gc_content = _calculate_gc_content(current)
                score = _score_sequence(current, 
                                      codon_table_obj.gc_range, 
                                      codon_table_obj.codon_table, 
                                      codon_table_obj.gc_weight,
                                      gc_tolerance)
                
                # Update tracking variables
                if score > best_score:
                    best_score = score
                    best_sequence = current
                
                if codon_table_obj.gc_range[0] <= gc_content <= codon_table_obj.gc_range[1]:
                    if score > best_in_range_score:
                        best_in_range_score = score
                        best_in_range_sequence = current
                        # Early stopping if we found an excellent sequence
                        if score >= early_stop_threshold:
                            remaining_iters = 0
                            break
            
            if show_progress_bar:
                pbar.update(1)
                
        # Break if we've hit our target
        if remaining_iters <= 0:
            break
    
    # Phase 3: Fine-tune the best sequence for GC content
    final_sequence = best_sequence
    
    # If we found any sequence in range, prefer that
    if best_in_range_sequence and best_in_range_score > best_score * 0.85:
        final_sequence = best_in_range_sequence
    
    # Apply GC fine-tuning
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
        pbar.close()
    
    # Return final sequence
    return final_sequence


def calculate_codon_adaptation_score(sequence, codon_frequency_table):
    """
    Calculate a codon adaptation score for a nucleotide sequence.
    
    The score represents how well the sequence's codon usage matches the optimal
    codon usage of the organism. A score of 1.0 means the sequence uses the most 
    common codons for each amino acid, while a score of 0.0 would use the rarest codons.
    
    Parameters:
    -----------
    sequence : str
        Nucleotide sequence to analyze (must be divisible by 3)
    codon_frequency_table : dict
        Dictionary with amino acids as keys and dictionaries of {codon: frequency} as values
    
    Returns:
    --------
    float
        A normalized score between 0 and 1 representing codon optimization
    """
    # Validate sequence length
    if len(sequence) % 3 != 0:
        raise ValueError("Sequence length must be divisible by 3")
    
    if not sequence:
        return 0
    
    total_score = 0
    n_codons = len(sequence) // 3
    
    # For each codon in the sequence
    for i in range(0, len(sequence), 3):
        if i+3 <= len(sequence):
            codon = sequence[i:i+3].upper()
            
            # Skip invalid codons
            if 'N' in codon or any(n not in 'ATGC' for n in codon):
                continue
                
            # Get amino acid for this codon
            try:
                aa = codons_to_aa[codon]
                
                # Get this codon's frequency relative to other codons for same AA
                codon_freq = codon_frequency_table[aa].get(codon, 0.01)
                
                # Find max frequency for this amino acid (best possible codon)
                max_freq = max(codon_frequency_table[aa].values())
                
                # Normalize score against the best possible codon
                if max_freq > 0:
                    normalized_score = codon_freq / max_freq
                else:
                    normalized_score = 0
                    
                total_score += normalized_score
                
            except KeyError:
                # Skip invalid codons (like stop codons if not in table)
                continue
    
    # Return average score normalized by sequence length
    if n_codons > 0:
        return total_score / n_codons
    else:
        return 0