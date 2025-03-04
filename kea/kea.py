"""
Functionality for building the library of nucleotide sequences 
from DNA sequences
"""
import os
from tqdm import tqdm

from .backend.codon_table import CodonTable
from .backend.codon_optimizer import CodonOptimizer  # Import the new class
from .backend.optimize_codon_usage import optimize_codon_usage
from .backend.translation_utils import translate_sequence
from .backend.kea_utils import generate_random_protein_ids
from .backend.sequence import Sequence
from .data.codon_tables import all_codon_tables



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
                  return_best=False):
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
        - If string: One of ['yeast', 's288c', 's288c_unweighted']
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
        Desired GC content range as (min, max) where both values are between 0 and 1.
        Default is (0, 1) which means no specific GC content is targeted.
    
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
        If True, verify that 5' adapter ends with start codon and 3' adapter begins with stop codon.
    
    show_optimization_progress : bool, default=True
        If True, display progress bars during sequence optimization.
    
    add_protein_identifiers : bool, default=False
        If True, generate random identifiers for protein sequences.
    
    protein_identifier_length : int, default=8
        Length of random protein identifiers if add_protein_identifiers is True.
    
    return_best : bool, default=False
        If True, return the best optimization result based on constraints rather than
        the final result after all iterations.
    
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
    try:
        # Validate and process input sequences first
        if isinstance(protein_sequences, str):
            sequences = [protein_sequences]
            names = generate_random_protein_ids(protein_identifier_length, 1)
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
        if target_gc_range is None:
            target_gc_range = (0, 1)
            gc_finetuning_iterations = 0
            gc_weight = 0
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

        # Initialize optimizer with validated parameters
        codon_table_obj = CodonTable(codon_frequency_table,
                                    codon_weight,
                                    gc_weight,
                                    target_gc_range,
                                    minimum_codon_probability)
        
        optimizer = CodonOptimizer(codon_table_obj)

        if adapter_3_prime==None:
            adapter_3_prime = ''
        if adapter_5_prime==None:
            adapter_5_prime = ''

        # Create sequence objects with proper error handling
        sequence_objects = []
        pbar = None
        try:
            if show_progress:
                pbar = tqdm(total=len(sequences), position=0, desc='Progress through sequences', leave=True)

            for i, (seq, name) in enumerate(zip(sequences, names)):
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
                                 padding_attempts=padding_attempts)
                sequence_objects.append(seq_obj)

                # Optimize sequence
                optimized_coding_seq = optimizer.optimize(
                    seq_obj.protein_sequence,
                    n_iter=optimization_attempts,
                    fine_tuning_iterations=gc_finetuning_iterations,
                    return_best=return_best,
                    show_progress_bar=show_optimization_progress
                )
                
                seq_obj.add_coding_sequence(optimized_coding_seq)

                if total_length is not None:
                    seq_obj.generate_padding()

                if pbar:
                    pbar.update(1)

        finally:
            if pbar:
                pbar.close()

        return sequence_objects

    except Exception as e:
        # Ensure progress bar is closed even if an error occurs
        if 'pbar' in locals() and pbar:
            pbar.close()
        raise


def save_library(list_of_sequence_objects, save_path):
    '''
    Function to save a library of nucleotide sequences to a file.
    
    Parameters
    ----------
    list_of_sequence_objects : list of Sequence
        A list of Sequence objects containing optimized DNA sequences and associated
        information for each input protein sequence.
    save_path : str
        Path to save the library file
    
    Returns
    -------
    None
    '''
    # make sure the path excluding the file name is a dir
    if not os.path.isdir(os.path.dirname(save_path)):
        raise ValueError(f"Path {os.path.dirname(save_path)} is not a directory.")
    
    # headers
    headers = ['Protein Name', 'Protein Sequence', 'Optimized Sequence', 'Coding Sequence', 'GC content Optimized Sequence', 'GC content Coding Sequence', 'Correct Translation Optimized Sequence?', 'Correct Translation Coding Sequence?']
    # open file
    with open(save_path, 'w') as f:
        # write headers
        f.write(','.join(headers) + '\n')
        # write data
        for seq_obj in list_of_sequence_objects:
            # get data
            data = [seq_obj.name, seq_obj.protein_sequence, seq_obj.full_dna_sequence, seq_obj.coding_sequence, seq_obj.gc_content_full_sequence, seq_obj.gc_content_coding_sequence, seq_obj.correct_full_translation, seq_obj.correct_coding_translation]
            # write data
            f.write(','.join([str(d) for d in data]) + '\n')
        
    f.close()





