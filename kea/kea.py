"""
Functionality for building the library of nucleotide sequences 
from DNA sequences
"""
import csv
import os
import random

import numpy as np
from tqdm import tqdm

from .backend.codon_table import CodonTable
from .backend.codon_optimizer import CodonOptimizer  # Import the new class
from .backend.optimize_codon_usage import optimize_codon_usage
from .backend.translation_utils import translate_sequence
from .backend.kea_utils import generate_random_protein_ids
from .backend.sequence import Sequence
from .backend.sequence_constraints import (ConstraintRepairError, is_violation_unsatisfiable,
                                          make_sequence_constraints, repair_coding_sequence)
from .backend.library_results import (FailedSequence, LibraryResult, QualityCheckError,
                                     evaluate_sequence)
from .backend.sequence_diversity import (
    SequenceDiversityError,
    diversify_coding_sequence as _diversify_coding_sequence,
    minimum_hamming_distance as _distance_to_closest_sequence,
)
from .data.codon_tables import all_codon_tables
from .backend.sequence_report import (BASIC_COLUMNS, REPORT_COLUMNS,
                                      annotate_sequence)


from .backend.optimize_codon_usage import calculate_codon_adaptation_score



def build_library(protein_sequences,
                  codon_frequency_table,
                  adapter_5_prime=None,
                  adapter_3_prime=None,
                  total_length=None,
                  force_start_codon=True,
                  force_stop_codon=True,
                  target_gc_range=None,
                  pad_location=None,
                  avoid_adding_start_codons=True,
                  avoid_adding_stop_codons=True,
                  optimization_attempts=5000,
                  gc_finetuning_iterations=2000,
                  padding_attempts=10000,
                  gc_weight=1,
                  codon_weight=1,
                  verify_coding_sequence=True,
                  minimum_codon_probability=None,
                  show_progress=True,
                  gc_tolerance=0.025,
                  verify_unique_protein_sequences=True,
                  verify_start_stop_codons_in_adapters=False,
                  show_optimization_progress=True,
                  add_protein_identifiers=False,
                  protein_identifier_length=8,
                  return_best=False,
                  early_stop_threshold=None,
                  seed=None,
                  sequences_per_protein=1,
                  minimum_hamming_distance=0,
                  hamming_distance_attempts=100,
                  avoid_human_splice_sites=False,
                  maxent_donor_threshold=6.0,
                  maxent_acceptor_threshold=7.0,
                  avoid_premature_polyadenylation=False,
                  polyadenylation_requires_downstream_element=False,
                  max_homopolymer_length=None,
                  max_tandem_repeat_copies=None,
                  max_tandem_repeat_unit_length=6,
                  max_direct_repeat_length=None,
                  transcribed_5_prime_context=None,
                  transcribed_3_prime_context=None,
                  constraint_repair_attempts=100,
                  constraint_retry_attempts=2,
                  minimum_codon_adaptation=None,
                  gc_centering=0.0,
                  escape_gc_boundary=None,
                  on_error="skip"):
    '''
    Build a library of optimized DNA sequences that encode given protein sequences.
    
    This function transforms protein sequences into DNA sequences with optimized codon usage,
    controlled GC content, and optional features like adapters and padding to achieve
    consistent sequence lengths. The process includes validation of translation fidelity.
    
    Parameters
    ----------
    protein_sequences : str, list, or dict
        Protein sequences to encode. Can be:
        - A single sequence string
        - A list of sequence strings
        - A dictionary {name: sequence}
    
    codon_frequency_table : str or dict
        Codon usage frequencies to guide optimization:
        - If string: One of ['yeast', 's288c', 's288c_unweighted',
          'human', 'human_unweighted']
        - If dict: {amino_acid: {codon: frequency}}
    
    adapter_5_prime : str, optional
        5' adapter sequence to add to the beginning of each DNA sequence.
    
    adapter_3_prime : str, optional
        3' adapter sequence to add to the end of each DNA sequence.
    
    total_length : int, optional
        Target length for all final DNA sequences. If specified, padding will be added
        to reach this length. Must be long enough to accommodate the longest sequence.
    
    force_start_codon : bool, default=True
        If True, ensure each protein sequence begins with methionine (M).
    
    force_stop_codon : bool, default=True
        If True, ensure each protein sequence ends with a stop codon (*).
    
    target_gc_range : tuple(float, float), optional
        Desired GC content range as (min, max), both between 0 and 1. Leave unset
        (None, the default) to skip GC targeting entirely: the GC term is dropped
        from the objective and GC fine-tuning is disabled. Note that passing
        (0, 1) explicitly is NOT the same thing -- it keeps the GC machinery active
        with a range that everything satisfies, and forces gc_weight to stay
        positive.
    
    pad_location : int or None, optional
        Where to add padding to reach total_length:
        - 5: Add all padding at 5' end
        - 3: Add all padding at 3' end
        - None: Distribute padding equally at both ends
    
    avoid_adding_start_codons : bool, default=True
        If True, avoid adding start codons (ATG) in padding sequences.
    
    avoid_adding_stop_codons : bool, default=True
        If True, avoid adding stop codons (TAA, TAG, TGA) in padding sequences.
    
    optimization_attempts : int, default=5000
        Number of iterations for codon optimization process.
    
    gc_finetuning_iterations : int, default=2000
        Number of iterations for GC content fine-tuning after main optimization.
    
    padding_attempts : int, default=10000
        Maximum attempts to generate padding sequences that meet constraints.
    
    gc_weight : float, default=1
        Weight factor for GC content optimization (higher values prioritize GC content).
    
    codon_weight : float, default=1
        Weight factor for codon usage optimization (higher values prioritize natural codon usage).
    
    verify_coding_sequence : bool, default=True
        If True, verify that coding sequences translate to expected protein sequences.
    
    minimum_codon_probability : float, optional
        Minimum probability threshold for codon selection. Codons below this threshold
        won't be used. Default is 0.0 (all codons allowed).
    
    show_progress : bool, default=True
        If True, display overall progress bar.
    
    gc_tolerance : float, default=0.025
        Maximum acceptable deviation from target GC content.
    
    verify_unique_protein_sequences : bool, default=True
        If True, check that input protein sequences are unique.
    
    verify_start_stop_codons_in_adapters : bool, default=False
        If True, verify that 5' adapter and 3' adapters do not have start or stop codons. 
    
    show_optimization_progress : bool, default=True
        If True, display progress bars during sequence optimization.
    
    add_protein_identifiers : bool, default=False
        If True, generate random identifiers for protein sequences.
    
    protein_identifier_length : int, default=8
        Length of random protein identifiers if add_protein_identifiers is True.
    
    return_best : bool, default=False
        If True, return the best optimization result found, even if it falls outside
        target_gc_range. If False (the default) a sequence that cannot be brought
        inside target_gc_range raises a ValueError rather than being returned.

    seed : int, optional
        Seed for reproducible libraries. Codon sampling uses numpy's global RNG and
        padding uses the stdlib random module, so both are seeded from this value;
        seeding only one of them yourself is not enough.

    sequences_per_protein : int, default=1
        Number of independently encoded DNA sequences requested for each input
        protein. When greater than one, outputs are named ``<name>_variant_1``,
        ``<name>_variant_2``, and so on.

    minimum_hamming_distance : int, default=0
        Minimum absolute number of differing nucleotide positions required between
        every pair of coding sequences generated for the same protein. Start/stop
        codons are included when they are forced; adapters and padding are not.

    hamming_distance_attempts : int, default=100
        Maximum randomized synonymous-diversification attempts for each additional
        encoding. Each attempt preserves allowed codons and an already satisfied GC
        range, then reruns sequence-constraint repair. If no candidate reaches the
        requested all-pairs distance, the variant is reported as a
        ``sequence_diversity`` failure.

    early_stop_threshold : float or None, default=None
        If set, stop optimizing a sequence as soon as an in-range candidate reaches
        this blended codon-usage/GC score (0-1). Left unset by default because a
        threshold around 0.95 is reached almost immediately for most proteins, which
        makes optimization_attempts have no effect. Convergence is detected
        automatically instead.

    avoid_human_splice_sites : bool, default=False
        If True, synonymously remove coding-sequence windows whose human
        MaxEntScan donor score is at least ``maxent_donor_threshold`` or whose
        acceptor score is at least ``maxent_acceptor_threshold``.

    maxent_donor_threshold : float, default=6.0
        Hard MaxEntScan score cutoff for predicted human 5' splice donors.

    maxent_acceptor_threshold : float, default=7.0
        Hard MaxEntScan score cutoff for predicted human 3' splice acceptors.

    avoid_premature_polyadenylation : bool, default=False
        If True, synonymously remove AATAAA, ATTAAA, and 11 common human PAS
        variants. Hits are annotated for downstream U-rich and GU-rich context.

    polyadenylation_requires_downstream_element : bool, default=False
        If True, only treat a PAS hexamer as a hard constraint when a U-rich or
        GU-rich downstream element is also detected. The default is conservative
        and removes every supported PAS hexamer.

    max_homopolymer_length : int or None, default=None
        Maximum allowed run of one nucleotide. Leave unset to disable.

    max_tandem_repeat_copies : int or None, default=None
        Maximum consecutive copies of a 2- to
        ``max_tandem_repeat_unit_length``-nt repeat unit. Leave unset to disable.

    max_tandem_repeat_unit_length : int, default=6
        Longest repeat unit considered by the tandem-repeat constraint.

    max_direct_repeat_length : int or None, default=None
        Maximum exact substring length allowed to occur twice. Leave unset to
        disable separated direct-repeat detection.

    transcribed_5_prime_context, transcribed_3_prime_context : str, optional
        Fixed transcribed sequence flanking the CDS. These bases are included in
        every constraint scan so sites created across UTR/CDS or vector/CDS
        junctions are detected, but only the CDS is synonymously changed.

    constraint_repair_attempts : int, default=100
        Maximum rounds of synonymous single- and double-codon repair.
    
    constraint_retry_attempts : int, default=2
        Extra rounds of re-optimization to try when constraint repair fails.
        Codon optimization is stochastic and repair is a greedy local search, so a
        sequence that cannot be repaired from one starting point is often fine from
        another. Violations that no synonymous encoding can clear are detected and
        fail immediately, so retries are not spent on impossible sequences. Set to
        0 for the previous single-shot behaviour.

    escape_gc_boundary : bool or None, default=None
        Let the refinement search make paired swaps -- improve codon usage at one
        position, pay the GC back at another -- so it can keep improving once GC
        is pinned against an edge of target_gc_range.

        Without it the search stalls: on converged sequences at GC (0.30, 0.42),
        zero single-codon improvements remained while ~3000 GC-neutral pair moves
        per sequence were still available. Enabling it drives those to zero and
        lifts codon adaptation from 0.859 to 0.910 with GC still fully in range.

        None (the default) enables it only when no sequence constraints are
        requested, because maximal codon optimization is more repetitive and so
        harder to make constraint-compliant: raw optimizer output goes from 2.3 to
        7.7 violations per sequence, and on tight ranges that costs yield
        (50.7 -> 19.0 built of 60 on (0.28, 0.36)). Pass True to force it on
        anyway when codon adaptation matters more than yield, or False to disable.

    gc_centering : float, default=0.0
        How strongly the optimizer prefers the middle of target_gc_range over its
        edges, 0-1. At 0 the GC term is flat across the whole range, so codon
        usage pins GC against whichever edge is closest to its own optimum.

        Off by default, and NOT a simple dial -- the response is non-monotonic.
        Measured on 150-250 aa human sequences with constraints enabled:

          gc_centering   range (0.28, 0.36)      range (0.45, 0.60)
          0.00           75.0/80 built           CAI 0.9950  (baseline)
          0.03           55.0/80 built           --
          0.05           62.7/80 built           CAI 0.9719
          0.10           77.7/80 built           CAI 0.9740
          0.20           77.7/80 built           --

        Values around 0.03-0.05 are a hazard band: the tilt is strong enough to
        pull GC off the edge but too weak to settle it at the centre, and yield
        drops sharply (75 -> 55 of 80 above). Use 0 or >= 0.10, never in between.

        At 0.10 it helps on narrow ranges that sit below the codon table's
        preferred GC (about 0.59 for "human", 0.35 for "yeast"), both in yield
        (77.7 vs 75.0 of 80) and, on tight ranges, in codon adaptation (+0.038 on
        (0.30, 0.42)). On ranges that already contain the preferred GC it only
        costs codon adaptation (-0.020 on (0.45, 0.60)) with no yield benefit.

        Note the in-range gradient is gc_centering / half_width, so the same value
        is a mild tilt on a wide range and an overwhelming one on a narrow range.

    minimum_codon_adaptation : float, optional
        Reject a finished sequence whose codon adaptation score falls below this
        (0-1). Left unset by default, in which case the score is still measured
        and reported on every sequence, just not enforced. Useful when a tight
        target_gc_range forces rare codons and you want a floor on how much codon
        optimality you are willing to trade away.

    on_error : {"skip", "raise"}, default="skip"
        What to do when a single sequence cannot be built (no synonymous encoding
        satisfies the requested constraints, or the GC range is unreachable).

        - "skip": record the failure and carry on with the remaining sequences.
          The returned LibraryResult holds the successes and a .failures list
          naming each sequence that failed and why. A one-line warning is printed
          when show_progress is True.
        - "raise": abort the whole build on the first failure.

        "skip" is the default so that one impossible sequence near the end of a
        large library does not discard everything already built.


    Returns
    -------
    list of Sequence
        A list of Sequence objects containing optimized DNA sequences and associated
        information for each input protein sequence.
        Sequence objects contain the following attributes:
            protein_sequence : str
                The final protein sequence (including any forced start/stop)
            full_dna_sequence : str
                Complete DNA sequence including adapters and padding
            coding_sequence : str
                DNA sequence encoding the protein
            gc_content_full_sequence : float
                GC content of the complete sequence
            gc_content_coding_sequence : float
                GC content of the coding sequence only
            correct_full_translation : bool
                Whether full sequence translates correctly
            correct_coding_translation : bool
                Whether coding sequence translates correctly
        
    Examples
    --------
    >>> from kea import build_library
    >>> sequences = ["MKKFLVLLFCWAVLCEHN", "MVLSEGEWQLVLHVWAKV"]
    >>> results = build_library(sequences, "yeast", 
    ...                         target_gc_range=(0.4, 0.6),
    ...                         total_length=120)
    '''
    if on_error not in ("skip", "raise"):
        raise ValueError('on_error must be "skip" or "raise"')
    failures = []
    try:
        # Randomness is split across two global RNGs: numpy drives codon sampling
        # and the stdlib random module drives padding. Both have to be seeded for a
        # library to be reproducible.
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)

        # Validate the iteration budgets up front, before target_gc_range
        # normalization can rewrite them, and before the per-sequence loop. These
        # are configuration errors -- identical for every input -- so they must
        # escape regardless of on_error rather than being recorded as thousands of
        # identical per-sequence "failures".
        if (not isinstance(optimization_attempts, int)
                or isinstance(optimization_attempts, bool)
                or optimization_attempts < 1):
            raise ValueError("optimization_attempts must be at least 1")
        if (not isinstance(gc_finetuning_iterations, int)
                or isinstance(gc_finetuning_iterations, bool)
                or gc_finetuning_iterations < 0):
            raise ValueError("gc_finetuning_iterations must be a non-negative integer")
        if (not isinstance(padding_attempts, int)
                or isinstance(padding_attempts, bool)
                or padding_attempts < 1):
            raise ValueError("padding_attempts must be at least 1")
        if (not isinstance(constraint_repair_attempts, int)
                or isinstance(constraint_repair_attempts, bool)
                or constraint_repair_attempts < 1):
            raise ValueError("constraint_repair_attempts must be at least 1")
        if (not isinstance(constraint_retry_attempts, int)
                or isinstance(constraint_retry_attempts, bool)
                or constraint_retry_attempts < 0):
            raise ValueError("constraint_retry_attempts must be a non-negative integer")
        if (not isinstance(sequences_per_protein, int)
                or isinstance(sequences_per_protein, bool)
                or sequences_per_protein < 1):
            raise ValueError("sequences_per_protein must be at least 1")
        if (not isinstance(minimum_hamming_distance, int)
                or isinstance(minimum_hamming_distance, bool)
                or minimum_hamming_distance < 0):
            raise ValueError("minimum_hamming_distance must be a non-negative integer")
        if (not isinstance(hamming_distance_attempts, int)
                or isinstance(hamming_distance_attempts, bool)
                or hamming_distance_attempts < 1):
            raise ValueError("hamming_distance_attempts must be at least 1")

        # Validate and process input sequences first
        if isinstance(protein_sequences, str):
            sequences = [protein_sequences]
            # Honour add_protein_identifiers here too, so a single string is named
            # the same way a one-element list would be.
            if add_protein_identifiers:
                names = generate_random_protein_ids(protein_identifier_length, 1)
            else:
                names = ["Protein_0"]
        elif isinstance(protein_sequences, list):
            if not protein_sequences:
                raise ValueError("Empty sequence list provided")
            if add_protein_identifiers:
                names = generate_random_protein_ids(protein_identifier_length, len(protein_sequences))
            else:
                names = [f"Protein_{i}" for i in range(len(protein_sequences))]
            sequences = protein_sequences
        elif isinstance(protein_sequences, dict):
            if not protein_sequences:
                raise ValueError("Empty sequence dictionary provided")
            names = list(protein_sequences.keys())
            sequences = list(protein_sequences.values())
        else:
            raise ValueError("protein_sequences must be a string, list, or dict")

        # Early validation of sequence content
        if any(not isinstance(seq, str) for seq in sequences):
            raise ValueError("All sequences must be strings")
        if any(not seq for seq in sequences):
            raise ValueError("Empty sequences are not allowed")

        # Validate unique sequences if required
        if verify_unique_protein_sequences and len(sequences) != len(set(sequences)):
            raise ValueError("Protein sequences are not unique")

        # Process codon frequency table first as it's critical
        if isinstance(codon_frequency_table, str):
            codon_frequency_table = codon_frequency_table.lower()
            if codon_frequency_table not in all_codon_tables:
                raise ValueError(f"Codon frequency table '{codon_frequency_table}' not found")
            codon_frequency_table = all_codon_tables[codon_frequency_table]
        elif not isinstance(codon_frequency_table, dict):
            raise ValueError("Codon frequency table must be a string or dictionary")

        # Validate and normalize GC parameters
        if codon_weight < 0:
            raise ValueError("codon_weight must be non-negative")
        if gc_weight < 0:
            raise ValueError("gc_weight must be non-negative")

        # What the caller actually asked for, before the None -> (0, 1) rewrite
        # below. The final quality check reports against this, so "no GC target"
        # is not silently reported as "target (0, 1), satisfied".
        requested_gc_range = target_gc_range

        if target_gc_range is None:
            # No GC target: drop the GC term entirely and skip fine-tuning, so the
            # optimizer is scored purely on codon usage.
            target_gc_range = (0, 1)
            gc_finetuning_iterations = 0
            gc_weight = 0
            if codon_weight == 0:
                raise ValueError("codon_weight must be positive when target_gc_range is None, "
                                 "otherwise there is nothing to optimize for")
        else:
            if not isinstance(target_gc_range, tuple) or len(target_gc_range) != 2:
                raise ValueError("target_gc_range must be a tuple of (min, max)")
            if not (0 <= target_gc_range[0] <= target_gc_range[1] <= 1):
                raise ValueError("Invalid GC range values")

            # Ensure gc_weight is positive and non-zero when GC optimization is needed
            if gc_weight <= 0:
                raise ValueError("gc_weight must be positive when GC optimization is enabled")

        # Ensure minimum_codon_probability is a valid float
        if minimum_codon_probability is None:
            minimum_codon_probability = 0.0
        else:
            try:
                minimum_codon_probability = float(minimum_codon_probability)
                if not (0 <= minimum_codon_probability <= 1):
                    raise ValueError("minimum_codon_probability must be between 0 and 1")
            except (TypeError, ValueError):
                raise ValueError("minimum_codon_probability must be a float between 0 and 1 or None")

        # Constraints are built first: whether any are enabled decides how the
        # optimizer is told to aim within the GC range.
        constraint_set = make_sequence_constraints(
            avoid_human_splice_sites=avoid_human_splice_sites,
            maxent_donor_threshold=maxent_donor_threshold,
            maxent_acceptor_threshold=maxent_acceptor_threshold,
            avoid_premature_polyadenylation=avoid_premature_polyadenylation,
            polyadenylation_requires_downstream_element=polyadenylation_requires_downstream_element,
            max_homopolymer_length=max_homopolymer_length,
            max_tandem_repeat_copies=max_tandem_repeat_copies,
            max_tandem_repeat_unit_length=max_tandem_repeat_unit_length,
            max_direct_repeat_length=max_direct_repeat_length,
            transcribed_5_prime_context=(transcribed_5_prime_context or ''),
            transcribed_3_prime_context=(transcribed_3_prime_context or ''),
        )

        # Initialize optimizer with validated parameters
        if not 0.0 <= gc_centering <= 1.0:
            raise ValueError("gc_centering must be between 0 and 1")
        # Only meaningful when a GC range was actually requested.
        if requested_gc_range is None:
            gc_centering = 0.0

        # The paired boundary-escape move maximizes codon adaptation, which makes
        # the sequence more repetitive and so harder to make constraint-compliant
        # (measured: 2.3 -> 7.7 violations per sequence). With no constraints there
        # is nothing to lose, so enable it; with constraints, leave it to the caller.
        if escape_gc_boundary is None:
            use_boundary_escape = not constraint_set.constraints
        else:
            use_boundary_escape = bool(escape_gc_boundary)

        codon_table_obj = CodonTable(codon_frequency_table,
                                    codon_weight,
                                    gc_weight,
                                    target_gc_range,
                                    minimum_codon_probability,
                                    gc_tolerance,
                                    gc_centering,
                                    use_boundary_escape)
        
        optimizer = CodonOptimizer(codon_table_obj)

        # Normalize the adapters BEFORE validating them. Sequence() upper-cases
        # adapters on the way into the construct, so a case-sensitive check here
        # would let a lowercase "atg" through and then splice an ATG into the
        # final sequence anyway.
        adapter_5_prime = '' if adapter_5_prime is None else adapter_5_prime.upper()
        adapter_3_prime = '' if adapter_3_prime is None else adapter_3_prime.upper()

        for label, adapter in (("5'", adapter_5_prime), ("3'", adapter_3_prime)):
            invalid = set(adapter) - set('ACGT')
            if invalid:
                raise ValueError(f"{label} adapter contains non-DNA characters: "
                                 f"{', '.join(sorted(invalid))}. Only A, C, G and T are allowed.")

        if verify_start_stop_codons_in_adapters:
            for label, adapter in (("5'", adapter_5_prime), ("3'", adapter_3_prime)):
                if 'ATG' in adapter:
                    raise ValueError(f"{label} adapter cannot contain start codon (ATG) when "
                                     "verify_start_stop_codons_in_adapters is True")
                if any(stop in adapter for stop in ('TAA', 'TAG', 'TGA')):
                    raise ValueError(f"{label} adapter cannot contain stop codon (TAA, TAG, TGA) when "
                                     "verify_start_stop_codons_in_adapters is True")

        if pad_location not in (None, 3, 5):
            raise ValueError("pad_location must be 3 (3' end), 5 (5' end), or None (both ends)")


        # Validate total_length up front, against the longest sequence in the library.
        # Without this an undersized total_length produces negative padding lengths that
        # both padding branches skip, silently returning over-length sequences.
        if total_length is not None:
            if not isinstance(total_length, int) or isinstance(total_length, bool) or total_length <= 0:
                raise ValueError("total_length must be a positive integer")
            adapter_length = len(adapter_5_prime) + len(adapter_3_prime)
            longest_required = 0
            for seq in sequences:
                seq = seq.upper()
                protein_length = len(seq)
                if force_start_codon and not seq.startswith('M'):
                    protein_length += 1
                if force_stop_codon and not seq.endswith('*'):
                    protein_length += 1
                longest_required = max(longest_required, protein_length * 3 + adapter_length)
            if total_length < longest_required:
                raise ValueError(
                    f"total_length={total_length} is too short for this library. The longest "
                    f"sequence needs at least {longest_required} nt "
                    f"(coding sequence + adapters).")

        # Build each requested encoding independently. Variants of one protein share
        # an accepted-CDS list so every new coding sequence can be checked against
        # every sibling, not just against the first one generated.
        sequence_objects = []
        total_requested = len(sequences) * sequences_per_protein
        pbar = None
        try:
            if show_progress:
                pbar = tqdm(total=total_requested, position=0,
                            desc='Progress through sequences', leave=True)

            for seq, base_name in zip(sequences, names):
                accepted_coding_sequences = []
                for variant_index in range(sequences_per_protein):
                    name = (base_name if sequences_per_protein == 1
                            else f"{base_name}_variant_{variant_index + 1}")
                    stage = "sequence_setup"
                    try:
                        seq_obj = Sequence(seq, codon_table_obj,
                                         name, adapter_3_prime,
                                         adapter_5_prime,
                                         avoid_adding_start_codons,
                                         avoid_adding_stop_codons,
                                         total_length,
                                         pad_location,
                                         force_start_codon=force_start_codon,
                                         force_stop_codon=force_stop_codon,
                                         verify_coding_sequence=verify_coding_sequence,
                                         return_best=return_best,
                                         padding_attempts=padding_attempts,
                                         constraint_set=constraint_set)

                        # Codon optimization is stochastic, and constraint repair is a
                        # greedy local search from wherever it lands. A sequence that
                        # cannot be repaired from one starting point is often fine from
                        # another, so re-optimize and try again rather than giving up on
                        # the first attempt. Violations that no synonymous encoding can
                        # clear are detected and fail immediately.
                        for attempt in range(constraint_retry_attempts + 1):
                            stage = "optimization"
                            optimized_coding_seq = optimizer.optimize(
                                seq_obj.protein_sequence,
                                n_iter=optimization_attempts,
                                fine_tuning_iterations=gc_finetuning_iterations,
                                return_best=return_best,
                                show_progress_bar=show_optimization_progress,
                                early_stop_threshold=early_stop_threshold,
                                gc_tolerance=gc_tolerance
                            )

                            stage = "constraint_repair"
                            try:
                                optimized_coding_seq, constraint_violations = repair_coding_sequence(
                                    optimized_coding_seq,
                                    codon_table_obj,
                                    constraint_set,
                                    max_iterations=constraint_repair_attempts,
                                )
                            except ConstraintRepairError as repair_error:
                                impossible = [
                                    violation for violation in repair_error.remaining_violations
                                    if is_violation_unsatisfiable(
                                        repair_error.sequence, violation,
                                        codon_table_obj, constraint_set)
                                ]
                                if impossible:
                                    preview = "; ".join(
                                        f"{item.kind} at {item.start}:{item.end} ({item.sequence}) "
                                        f"score={item.score:.2f} vs threshold {item.threshold:.2f}"
                                        for item in impossible[:3]
                                    )
                                    raise ConstraintRepairError(
                                        f"'{name}' cannot satisfy the requested constraints under any "
                                        f"synonymous encoding, so retrying will not help. Relax the "
                                        f"relevant threshold or change the protein. "
                                        f"Unsatisfiable: {preview}",
                                        remaining_violations=repair_error.remaining_violations,
                                        sequence=repair_error.sequence,
                                    ) from None
                                if attempt == constraint_retry_attempts:
                                    raise
                            else:
                                break

                        if accepted_coding_sequences and minimum_hamming_distance > 0:
                            stage = "sequence_diversity"
                            best_distance = _distance_to_closest_sequence(
                                optimized_coding_seq, accepted_coding_sequences
                            )
                            if minimum_hamming_distance > len(optimized_coding_seq):
                                raise SequenceDiversityError(
                                    f"'{name}' requested minimum_hamming_distance="
                                    f"{minimum_hamming_distance}, but its coding sequence is only "
                                    f"{len(optimized_coding_seq)} nt long.",
                                    minimum_hamming_distance,
                                    hamming_distance_attempts,
                                    best_distance,
                                )

                            for _ in range(hamming_distance_attempts):
                                candidate = _diversify_coding_sequence(
                                    optimized_coding_seq,
                                    codon_table_obj,
                                    accepted_coding_sequences,
                                    minimum_hamming_distance,
                                )
                                if candidate is None:
                                    continue
                                try:
                                    candidate, candidate_violations = repair_coding_sequence(
                                        candidate,
                                        codon_table_obj,
                                        constraint_set,
                                        max_iterations=constraint_repair_attempts,
                                    )
                                except ConstraintRepairError:
                                    # This particular diversified encoding created a
                                    # hard motif; another randomized encoding may not.
                                    continue

                                candidate_distance = _distance_to_closest_sequence(
                                    candidate, accepted_coding_sequences
                                )
                                best_distance = max(best_distance, candidate_distance)
                                if candidate_distance < minimum_hamming_distance:
                                    continue
                                if minimum_codon_adaptation is not None:
                                    adaptation = calculate_codon_adaptation_score(
                                        candidate, codon_table_obj.codon_usage)
                                    if adaptation < minimum_codon_adaptation:
                                        continue
                                optimized_coding_seq = candidate
                                constraint_violations = candidate_violations
                                break
                            else:
                                raise SequenceDiversityError(
                                    f"Could not generate '{name}' at least "
                                    f"{minimum_hamming_distance} nt from every accepted encoding "
                                    f"after {hamming_distance_attempts} attempts. Best minimum "
                                    f"distance found: {best_distance}.",
                                    minimum_hamming_distance,
                                    hamming_distance_attempts,
                                    best_distance,
                                )

                        sibling_distance = _distance_to_closest_sequence(
                            optimized_coding_seq, accepted_coding_sequences
                        )
                        seq_obj.source_protein_name = base_name
                        seq_obj.variant_index = variant_index + 1
                        seq_obj.minimum_sibling_hamming_distance = sibling_distance
                        seq_obj.sequence_constraint_violations = constraint_violations
                        seq_obj.add_coding_sequence(optimized_coding_seq)

                        if total_length is not None:
                            stage = "padding"
                            seq_obj.generate_padding()

                        # Single acceptance gate. Everything is re-measured here on
                        # the final construct before the CDS joins the sibling set.
                        stage = "quality_check"
                        report = evaluate_sequence(
                            seq_obj,
                            target_gc_range=requested_gc_range,
                            total_length=total_length,
                            constraint_set=constraint_set,
                            codon_usage=codon_table_obj.codon_usage,
                            minimum_codon_adaptation=minimum_codon_adaptation,
                            gc_is_required=not return_best,
                        )
                        seq_obj.quality_report = report
                        if not report.passed:
                            raise QualityCheckError(
                                f"'{name}' failed its final quality check: "
                                + "; ".join(f"{check.name}={check.value} "
                                            f"(target {check.target})"
                                            for check in report.errors),
                                report=report,
                            )

                    except Exception as error:
                        # Deliberately not a bare except: KeyboardInterrupt and
                        # SystemExit must still stop a long library build.
                        if on_error == "raise":
                            raise
                        # Once one slot cannot be filled, later slots cannot complete
                        # the requested all-pairs set either. Record every unfilled
                        # slot so requested == succeeded + failed remains meaningful.
                        unfilled = sequences_per_protein - variant_index
                        for failed_index in range(variant_index, sequences_per_protein):
                            failed_name = (base_name if sequences_per_protein == 1
                                           else f"{base_name}_variant_{failed_index + 1}")
                            failures.append(FailedSequence(
                                name=failed_name,
                                protein_sequence=seq,
                                error=error,
                                stage=stage,
                                remaining_violations=getattr(
                                    error, "remaining_violations", []),
                            ))
                        if pbar:
                            pbar.update(unfilled)
                        break
                    else:
                        accepted_coding_sequences.append(optimized_coding_seq)
                        sequence_objects.append(seq_obj)
                        if pbar:
                            pbar.update(1)

        finally:
            if pbar:
                pbar.close()

        library = LibraryResult(sequence_objects, failures=failures, requested=total_requested)

        if failures and show_progress:
            print(f"\n{library.summary()}")

        return library

    except Exception as e:
        # Ensure progress bar is closed even if an error occurs
        if 'pbar' in locals() and pbar:
            pbar.close()
        raise


def save_library(list_of_sequence_objects, save_path, include_diagnostics=True):
    '''
    Save a library to CSV, with per-sequence design diagnostics.

    Parameters
    ----------
    list_of_sequence_objects : list of Sequence
        Sequences to write. A LibraryResult works directly.
    save_path : str
        Path to write to. A bare filename writes to the current directory.
    include_diagnostics : bool, default=True
        Write the full diagnostic report: strongest predicted splice donor and
        acceptor with their MaxEntScan scores and positions, polyadenylation
        hexamer counts, sequence-complexity measures, and the final quality-check
        outcome. Set False for the short form (sequences, GC, translation checks
        and codon adaptation only).

        Splice and polyadenylation screens run on the complete construct
        (adapters + padding + coding sequence), and are reported whether or not
        the matching constraint was enabled during the build -- so a library built
        without splice avoidance can still be screened afterwards.

        These are scores, not probabilities. See "Reading the diagnostics" in the
        README for how to interpret them.

    Returns
    -------
    None
    '''
    # An empty dirname means a bare filename like "library.csv", which is the cwd
    # and perfectly valid.
    save_dir = os.path.dirname(save_path)
    if save_dir and not os.path.isdir(save_dir):
        raise ValueError(f"Path {save_dir} is not a directory.")

    columns = REPORT_COLUMNS if include_diagnostics else BASIC_COLUMNS

    # Written through the csv module so that names containing commas or quotes are
    # escaped rather than silently shifting every subsequent column.
    with open(save_path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        for seq_obj in list_of_sequence_objects:
            record = annotate_sequence(seq_obj)
            if record.get("codon_adaptation_index") in ('', None):
                # No quality report (Sequence built outside build_library): score
                # against the organism's full codon usage, as the report would.
                record["codon_adaptation_index"] = round(calculate_codon_adaptation_score(
                    seq_obj.coding_sequence, seq_obj.codon_table.codon_usage), 4)
            writer.writerow(record)


def save_failures(library_result, save_path):
    '''
    Save the sequences that could not be built, and why.

    Parameters
    ----------
    library_result : LibraryResult
        The result of build_library.
    save_path : str
        Path to write to.

    Returns
    -------
    None
    '''
    save_dir = os.path.dirname(save_path)
    if save_dir and not os.path.isdir(save_dir):
        raise ValueError(f"Path {save_dir} is not a directory.")

    columns = ["protein_name", "protein_sequence", "stage", "reason",
               "unsatisfiable_constraints"]
    with open(save_path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for failure in getattr(library_result, 'failures', []):
            writer.writerow({
                "protein_name": failure.name,
                "protein_sequence": failure.protein_sequence,
                "stage": failure.stage,
                "reason": failure.reason,
                "unsatisfiable_constraints": "; ".join(
                    f"{item.kind}@{item.start} score={item.score:.2f}"
                    for item in failure.remaining_violations[:5]),
            })
