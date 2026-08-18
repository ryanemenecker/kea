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


# Per-table caches of derived codon scores. Codon tables are plain dicts and so
# unhashable, which forces keying by id() -- but id() is only unique among LIVE
# objects. Short-lived CodonTables reuse heap addresses aggressively (measured: 60
# successive tables occupied just 6 distinct addresses), so a bare id() key hands
# one table's scores to another. Each entry therefore keeps a strong reference to
# the table it was computed from, which both pins the address for the entry's
# lifetime and lets the lookup re-verify identity.
_MAX_FREQ_CACHE = {}
_NORM_SCORE_CACHE = {}

# Bound on cached tables, so a process building many libraries does not retain
# every codon table it has ever seen.
_TABLE_CACHE_LIMIT = 32


def _cache_lookup(cache, table):
    """Return the cached value computed from exactly this table, else None."""
    entry = cache.get(id(table))
    if entry is not None and entry[0] is table:
        return entry[1]
    return None


def _cache_store(cache, table, value):
    if len(cache) >= _TABLE_CACHE_LIMIT:
        cache.clear()
    cache[id(table)] = (table, value)
    return value


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


def _build_sequences_batch(amino_acid_sequence, aa_weights, num_sequences):
    """
    Build multiple DNA sequences in parallel for the same amino acid sequence.
    
    Parameters:
    -----------
    amino_acid_sequence : str
        Amino acid sequence to encode
    aa_weights : dict
        Dictionary with amino acids as keys and (codons, weights) tuples as values
    num_sequences : int
        Number of sequences to generate simultaneously
        
    Returns:
    --------
    list of str
        List of DNA sequences encoding the amino acid sequence
    """
    # Special case for single sequence
    if num_sequences == 1:
        return [_build_sequence(amino_acid_sequence, aa_weights)]
    
    # For each position in the sequence, we'll select num_sequences codons
    seq_length = len(amino_acid_sequence)
    
    # Map of amino acid to its positions in the sequence
    aa_positions = {}
    for i, aa in enumerate(amino_acid_sequence):
        if aa not in aa_positions:
            aa_positions[aa] = []
        aa_positions[aa].append(i)
    
    # Initialize array to hold all codons (num_sequences × sequence_length)
    all_codons = np.empty((num_sequences, seq_length), dtype=object)
    
    # For each amino acid, generate codons for all its positions at once
    for aa, positions in aa_positions.items():
        codons, weights = aa_weights[aa]
        # Generate codons for all occurrences of this amino acid across all sequences
        selected_codons = np.random.choice(
            codons, 
            size=(num_sequences, len(positions)), 
            p=weights,
            replace=True
        )
        
        # Place the selected codons in their corresponding positions
        for idx, pos in enumerate(positions):
            all_codons[:, pos] = selected_codons[:, idx]
    
    # Convert arrays of codons to DNA sequences
    return [''.join(seq) for seq in all_codons]


def _calculate_gc_content(sequence):
    """Calculate GC content of a DNA sequence more efficiently."""
    if not sequence:
        return 0
    # str.count runs in C. The previous generator expression walked the sequence
    # one Python-level step per base and accounted for 171 million inner
    # iterations in a 10-sequence profile.
    return (sequence.count('G') + sequence.count('C')
            + sequence.count('g') + sequence.count('c')) / len(sequence)


# G+C base count for every codon. Integer counts let callers track GC exactly
# while swapping codons, with no floating-point drift.
_CODON_GC_COUNT = {codon: codon.count('G') + codon.count('C') for codon in codons_to_aa}


def _gc_score(gc_content, gc_range, gc_tolerance, gc_centering=0.0):
    """GC term of the objective: best inside the range, decaying outside it.

    Continuous and monotonically decreasing in distance from the centre of the
    range, so a sequence closer to the target never scores worse than one further
    away.

    With gc_centering at 0 the term is flat at 1.0 everywhere inside the range,
    which leaves the optimizer indifferent about *where* in the range it lands --
    codon usage then drags it hard against whichever edge is closest to its own
    optimum (human usage plus a range of (0.30, 0.42) reliably lands on 0.419).
    Constraint repair is left with no room to shift GC. A small positive
    gc_centering tilts the plateau toward the middle so repair has headroom by
    construction, without meaningfully competing with codon usage.
    """
    low, high = gc_range
    half_width = (high - low) / 2.0
    # Score at either edge of the range. For a zero-width range the centre IS the
    # edge, so there is no tilt and the edge value stays 1.0 -- computing it here
    # keeps the outside branch continuous with the inside one in that case too.
    edge_value = 1.0 - gc_centering if half_width > 0.0 else 1.0

    if low <= gc_content <= high:
        if gc_centering <= 0.0 or half_width <= 0.0:
            return 1.0
        # 0 at the centre of the range, 1 at either edge.
        offset = abs(gc_content - (low + high) / 2.0) / half_width
        return 1.0 - gc_centering * offset

    deviation = (low - gc_content) if gc_content < low else (gc_content - high)
    if gc_tolerance <= 0:
        return 0.0
    # Scaled to meet the in-range curve exactly at the boundary, so the score stays
    # continuous and never rewards stepping outside the range.
    return edge_value / (1.0 + deviation / gc_tolerance)


def _combine_scores(usage_score, gc_score_value, gc_weight, codon_weight):
    """Weighted average of the two objective terms. Always in 0-1."""
    total_weight = gc_weight + codon_weight
    if total_weight <= 0:
        return 0.0
    return (usage_score * codon_weight + gc_score_value * gc_weight) / total_weight


def _normalized_codon_scores(codon_frequency_table):
    """
    Per-codon usage score (frequency / best frequency for its amino acid), cached
    per table.

    Precomputing this collapses the inner loop of _score_sequence from two dict
    lookups plus a division per codon down to one dict lookup, and lets
    _improve_sequence update a running total instead of rescoring whole sequences.
    """
    cached = _cache_lookup(_NORM_SCORE_CACHE, codon_frequency_table)
    if cached is None:
        max_freqs = _max_codon_frequencies(codon_frequency_table)
        cached = {}
        for codon, aa in codons_to_aa.items():
            max_freq = max_freqs.get(aa, 0.0)
            if max_freq > 0:
                # 0.01 matches the fallback _score_sequence used for codons the
                # minimum_codon_probability filter removed from the table.
                cached[codon] = codon_frequency_table.get(aa, {}).get(codon, 0.01) / max_freq
            else:
                cached[codon] = 0.0
        _cache_store(_NORM_SCORE_CACHE, codon_frequency_table, cached)
    return cached


def _max_codon_frequencies(codon_frequency_table):
    """
    Best (highest) codon frequency per amino acid, cached per table.

    _score_sequence used to recompute max(...values()) for every codon of every
    candidate, which is the same handful of numbers over and over.
    """
    cached = _cache_lookup(_MAX_FREQ_CACHE, codon_frequency_table)
    if cached is None:
        cached = {aa: (max(codons.values()) if codons else 0.0)
                  for aa, codons in codon_frequency_table.items()}
        _cache_store(_MAX_FREQ_CACHE, codon_frequency_table, cached)
    return cached


def _score_sequence(sequence, gc_range, codon_frequency_table, gc_weight=1.0,
                    gc_tolerance=0.025, codon_weight=1.0, gc_centering=0.0):
    """Score a sequence on codon usage and GC content as a weighted average.

    Parameters:
    - sequence: DNA sequence to score
    - gc_range: Tuple of (min_gc, max_gc) acceptable range
    - codon_frequency_table: Dictionary of codon frequencies
    - gc_weight: Relative weight of the GC term (>= 0)
    - gc_tolerance: GC deviation at which the GC term falls to half credit
    - codon_weight: Relative weight of the codon usage term (>= 0)

    The two terms are combined as a weighted average, so the result is always in
    0-1 and both weights are honoured at any magnitude. The previous formula,
    (usage*(2 - gc_priority) + gc*gc_priority)/2 with gc_priority capped at 2,
    silently multiplied the usage term by zero once gc_weight reached 2.

    Returns:
    - Combined normalized score between 0-1
    """
    if not sequence:
        return 0

    # Number of codons
    n_codons = len(sequence) // 3
    if n_codons == 0:
        return 0

    # GC term. Continuous and monotonically decreasing in distance outside the
    # range: full credit inside, half credit one tolerance outside, decaying from
    # there. The old piecewise version was discontinuous at the tolerance edge and
    # inverted the ranking -- GC 0.4233 scored 0.9922 while 0.4267, which is closer
    # to the range, scored 0.0667.
    gc_score = _gc_score(_calculate_gc_content(sequence), gc_range, gc_tolerance,
                         gc_centering)

    # Codon usage term, normalized against the best codon available for each
    # amino acid. The per-codon scores are precomputed and cached per table.
    norm = _normalized_codon_scores(codon_frequency_table)
    usage_score = 0.0
    for i in range(0, n_codons * 3, 3):
        usage_score += norm[sequence[i:i + 3]]

    return _combine_scores(usage_score / n_codons, gc_score, gc_weight, codon_weight)


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

    # GC content is tracked as an integer count over a fixed total length, and each
    # candidate swap is scored by an exact O(1) update. The previous version copied
    # the codon list, re-joined it and rescanned the whole string for every
    # candidate, which dominated runtime whenever a GC range was set. Using integer
    # counts (rather than adding float deltas) keeps the result bit-identical.
    total_length = len(sequence)
    if total_length == 0:
        return sequence
    codon_gc_counts = {codon: codon.count('G') + codon.count('C')
                       for aa_codons in codon_frequency_table.values()
                       for codon in aa_codons}
    current_gc_count = sum(codon_gc_counts.get(codon, codon.count('G') + codon.count('C'))
                           for codon in codons)

    for _ in range(iterations):
        # Optimization: Instead of evaluating all positions each time, 
        # sample a subset for very long sequences
        positions_to_evaluate = range(len(codons))
        if len(codons) > 1000:
            positions_to_evaluate = random.sample(range(len(codons)), 1000)
            
        best_pos = -1
        best_replacement = None
        best_score = float('-inf')
        best_gc_count = current_gc_count
        current_diff = abs(current_gc - target_gc)

        for pos in positions_to_evaluate:
            codon = codons[pos]
            aa = codons_to_aa[codon]
            current_freq = codon_frequency_table[aa].get(codon, 0.01)
            codon_gc_count = codon_gc_counts[codon]
            codon_gc = codon_gc_count / 3

            # Get alternative codons for this amino acid. These come from
            # codon_frequency_table, which has already had minimum_codon_probability
            # applied; the module-level aa_to_codons is unfiltered and would let GC
            # fine-tuning reintroduce codons the caller explicitly banned.
            alt_codons = [c for c in codon_frequency_table[aa] if c != codon]
            if not alt_codons:
                continue
                
            # Calculate improvement metrics for each alternative
            for alt_codon in alt_codons:
                # Calculate GC change
                alt_gc_count = codon_gc_counts[alt_codon]
                gc_delta = (alt_gc_count - codon_gc_count) / 3

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

                # Exact O(1) impact on overall GC, no sequence rebuild required.
                new_gc_count = current_gc_count - codon_gc_count + alt_gc_count
                new_gc = new_gc_count / total_length
                gc_improvement = current_diff - abs(new_gc - target_gc)

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
                    best_gc_count = new_gc_count

        # Apply the best change if one was found
        if best_pos >= 0:
            codons[best_pos] = best_replacement
            current_gc_count = best_gc_count
            current_gc = current_gc_count / total_length

            # Update tracking variables. Only materialize the string when this is
            # actually the best sequence seen so far.
            if abs(current_gc - target_gc) < best_gc_diff:
                best_sequence = ''.join(codons)
                best_gc_diff = abs(current_gc - target_gc)
        else:
            # No beneficial changes found, break the loop
            break

    return best_sequence


def _gc_delta_index(codons, norm, table):
    """
    Every synonymous swap in the sequence, grouped by the GC-count change it makes.

    Used to escape the GC-range boundary. Entries are ordered by usage gain, so the
    first usable one is the cheapest way to pay back a given amount of GC.
    """
    index = {}
    for position, codon in enumerate(codons):
        amino_acid = codons_to_aa.get(codon)
        if amino_acid is None:
            continue
        base_norm = norm.get(codon, 0.0)
        base_gc = _CODON_GC_COUNT[codon]
        for alternative in table[amino_acid]:
            if alternative == codon:
                continue
            delta = _CODON_GC_COUNT[alternative] - base_gc
            if delta:
                index.setdefault(delta, []).append(
                    (norm[alternative] - base_norm, position, codon, alternative))
    for entries in index.values():
        entries.sort(reverse=True)
    return index


def _improve_sequence(sequence, codon_table_obj, mutation_rate=0.3, gc_tolerance=None):
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
    n_codons = len(sequence) // 3
    if n_codons == 0:
        return sequence

    codons = [sequence[i:i + 3] for i in range(0, n_codons * 3, 3)]

    table = codon_table_obj.codon_table
    norm = _normalized_codon_scores(table)
    gc_range = codon_table_obj.gc_range
    # Fall back to the table's value only when the caller did not supply one.
    # Reading it independently here let the refinement loop score against a
    # different tolerance than the loop that ranks its output, so the two were
    # optimizing subtly different objectives.
    if gc_tolerance is None:
        gc_tolerance = getattr(codon_table_obj, 'gc_tolerance', 0.025)
    gc_weight = codon_table_obj.gc_weight
    codon_weight = codon_table_obj.usage_weight
    gc_centering = getattr(codon_table_obj, 'gc_centering', 0.0)
    if gc_weight + codon_weight <= 0:
        return sequence

    # The objective is a function of exactly two running totals: the integer G+C
    # base count and the sum of per-codon usage scores. Both update in O(1) when a
    # codon is swapped, so a candidate can be scored without rebuilding the
    # sequence. The previous version called ''.join() and rescanned every base and
    # every codon for each candidate, which was 78% of total runtime.
    total_length = n_codons * 3
    gc_count = 0
    usage_sum = 0.0
    for codon in codons:
        gc_count += _CODON_GC_COUNT[codon]
        usage_sum += norm.get(codon, 0.0)

    def score_of(usage_total, gc_total):
        return _combine_scores(usage_total / n_codons,
                               _gc_score(gc_total / total_length, gc_range, gc_tolerance,
                                         gc_centering),
                               gc_weight, codon_weight)

    current_score = score_of(usage_sum, gc_count)
    improved = False

    # Escaping the GC boundary.
    #
    # When the codon-usage optimum lies outside gc_range (human usage is ~0.59 GC
    # against a requested (0.30, 0.42), say) this hill-climb converges hard onto
    # the range edge: every remaining usage-improving single-codon swap pushes GC
    # out of range and is rejected. Measured on converged sequences, ZERO
    # single-codon improvements remained while ~3000 GC-neutral PAIR moves per
    # sequence were still available -- improve usage at one position, pay the GC
    # back at another. Those are invisible to a one-codon-at-a-time search.
    #
    # So when a swap would leave the range, look for a partner swap that brings it
    # back. The index is built once per call; entries carry the codon they were
    # computed from so a stale one is skipped rather than misapplied.
    low_gc, high_gc = gc_range
    lowest_count = low_gc * total_length
    highest_count = high_gc * total_length
    boundary_escape = (getattr(codon_table_obj, 'escape_gc_boundary', True)
                       and gc_weight > 0 and (low_gc > 0.0 or high_gc < 1.0))
    delta_index = _gc_delta_index(codons, norm, table) if boundary_escape else None

    def paired_escape(primary_position, primary_usage_sum, primary_gc_count):
        """Best (score, usage_sum, gc_count, position, codon) partner swap, or None."""
        best = None
        first = int(np.ceil(lowest_count - primary_gc_count - 1e-9))
        last = int(np.floor(highest_count - primary_gc_count + 1e-9))
        for delta in range(first, last + 1):
            if delta == 0:
                continue
            for usage_delta, position, expected, replacement in delta_index.get(delta, ()):
                if position == primary_position or codons[position] != expected:
                    continue
                candidate_usage = primary_usage_sum + usage_delta
                candidate_gc = primary_gc_count + delta
                candidate_score = score_of(candidate_usage, candidate_gc)
                if best is None or candidate_score > best[0]:
                    best = (candidate_score, candidate_usage, candidate_gc,
                            position, replacement)
                # Entries are sorted by usage gain, so the first usable one is the
                # best available for this delta.
                break
        return best

    # Direct sampling instead of shuffling - more efficient for large sequences
    positions_to_check = np.random.choice(
        n_codons,
        size=max(1, int(n_codons * mutation_rate)),
        replace=False
    )

    # Allowed codons are cached per amino acid, taken from codon_table_obj.codon_table
    # which has already had minimum_codon_probability applied. The module-level
    # aa_to_codons is unfiltered and would reintroduce banned rare codons. The cache
    # is NOT pre-filtered against the codon at one particular position, because the
    # same amino acid appears at positions holding different codons.
    codons_for_aa = {}

    for pos in positions_to_check:
        codon = codons[pos]
        aa = codons_to_aa.get(codon)
        if aa is None:
            # Skip invalid codons
            continue

        allowed_codons = codons_for_aa.get(aa)
        if allowed_codons is None:
            allowed_codons = codons_for_aa[aa] = tuple(table[aa])

        # Nothing to swap to if this amino acid only has one allowed codon.
        if len(allowed_codons) < 2:
            continue

        codon_gc = _CODON_GC_COUNT[codon]
        codon_norm = norm.get(codon, 0.0)

        best_alt_codon = None
        best_alt_score = current_score
        best_usage_sum = usage_sum
        best_gc_count = gc_count
        best_partner = None

        # Every alternative is now cheap enough to score exactly, so there is no
        # longer a heuristic pre-filter deciding which ones are "promising". That
        # filter used to skip candidates that improved the combined objective
        # without improving either term on its own.
        for alt_codon in allowed_codons:
            if alt_codon == codon:
                continue

            new_usage_sum = usage_sum - codon_norm + norm[alt_codon]
            new_gc_count = gc_count - codon_gc + _CODON_GC_COUNT[alt_codon]
            new_score = score_of(new_usage_sum, new_gc_count)

            if new_score > best_alt_score:
                best_alt_score = new_score
                best_alt_codon = alt_codon
                best_usage_sum = new_usage_sum
                best_gc_count = new_gc_count
                best_partner = None
            elif (boundary_escape
                  and norm[alt_codon] > codon_norm
                  and not (lowest_count <= new_gc_count <= highest_count)):
                # Rejected only because it left the GC range, and it does improve
                # codon usage -- worth paying the GC back somewhere else.
                partner = paired_escape(pos, new_usage_sum, new_gc_count)
                if partner is not None and partner[0] > best_alt_score:
                    best_alt_score = partner[0]
                    best_alt_codon = alt_codon
                    best_usage_sum = partner[1]
                    best_gc_count = partner[2]
                    best_partner = (partner[3], partner[4])

        # Apply change if better alternative found
        if best_alt_codon is not None:
            codons[pos] = best_alt_codon
            if best_partner is not None:
                partner_position, partner_codon = best_partner
                codons[partner_position] = partner_codon
            usage_sum = best_usage_sum
            gc_count = best_gc_count
            current_score = best_alt_score
            improved = True

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
                         early_stop_threshold=None,
                         show_progress_bar=True,
                         gc_tolerance=0.025,
                         no_progress_patience=25):
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
    early_stop_threshold : float or None, default=None
        If set, stop as soon as an in-range sequence scores at or above this value.
        Scores are the blended codon-usage/GC score from _score_sequence, which is
        capped at 1.0. Note that a threshold of ~0.95 is reached almost immediately
        for most proteins, which makes n_iter meaningless -- hence the None default.
        Convergence is handled by no_progress_patience instead.
    show_progress_bar : bool
        If True, show a progress bar for the optimization process.
    gc_tolerance : float, default=0.025
        Allowable deviation from GC range (as a fraction, e.g., 0.025 = ±2.5%)
    no_progress_patience : int, default=25
        Stop refining a candidate after this many consecutive iterations that
        produce no improvement. The unused iterations are returned to the budget
        for the remaining candidates.
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
        
    if early_stop_threshold is not None and not (0 <= early_stop_threshold <= 1):
        raise ValueError("early_stop_threshold must be between 0 and 1, or None to disable")

    if n_iter < 1:
        raise ValueError("n_iter must be at least 1")

    if no_progress_patience < 1:
        raise ValueError("no_progress_patience must be at least 1")
    
    # Validate amino acid characters against the codon table actually in use, not
    # the module-level default -- a custom table missing an amino acid would
    # otherwise pass this check and then fail with a bare KeyError in _build_sequence.
    valid_aas = set(codon_table_obj.aa_weights)
    invalid_aas = sorted({aa for aa in amino_acid_sequence if aa not in valid_aas})
    if invalid_aas:
        raise ValueError(f"Invalid amino acid(s) in sequence: {', '.join(invalid_aas)}. "
                         "The codon frequency table in use does not cover them.")
    
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
    
    # Phase 1: Generate diverse initial candidates through random sampling.
    # Floor at 1: int(n_iter * 0.3) is 0 for n_iter <= 3, which used to leave
    # best_sequence as None and crash later with an opaque TypeError.
    exploration_iterations = max(1, int(n_iter * exploration_fraction))
    candidates = []
    
    # Determine batch size based on exploration iterations
    # Smaller batches for shorter sequences, larger for longer ones
    batch_size = min(100, exploration_iterations)
    sequences_generated = 0
    
    while sequences_generated < exploration_iterations:
        current_batch_size = min(batch_size, exploration_iterations - sequences_generated)
        
        # Generate batch of sequences
        sequences_batch = _build_sequences_batch(
            amino_acid_sequence, 
            codon_table_obj.aa_weights, 
            current_batch_size
        )
        
        # Process each sequence in the batch
        for sequence in sequences_batch:
            # Calculate metrics
            gc_content = _calculate_gc_content(sequence)
            score = _score_sequence(sequence,
                                  codon_table_obj.gc_range,
                                  codon_table_obj.codon_table,
                                  codon_table_obj.gc_weight,
                                  gc_tolerance,
                                  codon_table_obj.usage_weight,
                                  getattr(codon_table_obj, 'gc_centering', 0.0))
            
            # Track candidates
            candidates.append((sequence, score, gc_content))
            
            # Update best sequences
            if score > best_score:
                best_score = score
                best_sequence = sequence
            
            if codon_table_obj.gc_range[0] <= gc_content <= codon_table_obj.gc_range[1]:
                if score > best_in_range_score:
                    best_in_range_score = score
                    best_in_range_sequence = sequence
                    # Early stopping if we found an excellent sequence
                    if early_stop_threshold is not None and score >= early_stop_threshold:
                        if show_progress_bar:
                            pbar.update(exploration_iterations - sequences_generated)
                        sequences_generated = exploration_iterations  # Force exit the loop
                        break
            
            if show_progress_bar:
                pbar.update(1)
                
        sequences_generated += len(sequences_batch)
    
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
        # _improve_sequence samples a random subset of positions each call, so a
        # single empty pass is not proof of convergence -- but a long run of them
        # is. Without this counter the loop spins on no-ops: measured 98% of 3500
        # calls returning the sequence unchanged.
        stagnant = 0

        # Evolve this candidate through multiple iterations
        for iteration in range(candidate_iters):
            # Try improving the sequence
            improved = _improve_sequence(current, codon_table_obj,
                                         gc_tolerance=gc_tolerance)

            # Update if improved
            if improved != current:
                stagnant = 0
                current = improved

                # Recalculate metrics
                gc_content = _calculate_gc_content(current)
                score = _score_sequence(current,
                                      codon_table_obj.gc_range,
                                      codon_table_obj.codon_table,
                                      codon_table_obj.gc_weight,
                                      gc_tolerance,
                                      codon_table_obj.usage_weight,
                                  getattr(codon_table_obj, 'gc_centering', 0.0))

                # Update tracking variables
                if score > best_score:
                    best_score = score
                    best_sequence = current

                if codon_table_obj.gc_range[0] <= gc_content <= codon_table_obj.gc_range[1]:
                    if score > best_in_range_score:
                        best_in_range_score = score
                        best_in_range_sequence = current
                        # Early stopping if we found an excellent sequence
                        if early_stop_threshold is not None and score >= early_stop_threshold:
                            remaining_iters = 0
                            if show_progress_bar:
                                pbar.update(candidate_iters - iteration)
                            break
            else:
                stagnant += 1
                if stagnant >= no_progress_patience:
                    # This candidate has converged. Give its unused budget back so
                    # another candidate can use it instead of burning it on no-ops.
                    remaining_iters += candidate_iters - iteration - 1
                    if show_progress_bar:
                        pbar.update(candidate_iters - iteration)
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

    # Apply GC fine-tuning, but only when the sequence is actually outside the
    # requested range. gc_range is a range, not a point: fine-tuning every sequence
    # to its midpoint collapses the whole library onto one GC value and trades away
    # codon adaptation for a change the caller never asked for.
    gc_min, gc_max = codon_table_obj.gc_range
    current_gc = _calculate_gc_content(final_sequence)
    if fine_tuning_iterations > 0 and not (gc_min <= current_gc <= gc_max):
        # Aim just inside the nearest edge, so we make the smallest change that
        # satisfies the constraint rather than overshooting to the middle.
        margin = min(gc_tolerance, (gc_max - gc_min) / 2)
        gc_target = (gc_min + margin) if current_gc < gc_min else (gc_max - margin)
        final_sequence = _adjust_gc_content(
            final_sequence,
            gc_target,
            codon_table_obj.codon_table,
            iterations=fine_tuning_iterations
        )

    # Close the progress bar before the range check below, which may raise -- a
    # tqdm bar left open here corrupts all subsequent terminal output.
    if show_progress_bar:
        pbar.close()

    # Check if final GC is within range
    final_gc = _calculate_gc_content(final_sequence)

    # If GC content not within range and return_best is False, raise error
    if final_gc < gc_min or final_gc > gc_max:
        if not return_best:
            raise ValueError(f"No sequence found within GC range. Best achieved: {final_gc:.3f}")

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