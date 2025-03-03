"""
Functionality for building the library of nucleotide sequences 
from DNA sequences
"""
import os
from tqdm import tqdm

from .backend.codon_table import CodonTable
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
    # get sequences as a list and names as a list
    # Check if protein_sequences is a string
    if isinstance(protein_sequences, str):
        sequences = [protein_sequences]
        names = generate_random_protein_ids(protein_identifier_length, 1)
    elif isinstance(protein_sequences, list):
        if add_protein_identifiers:
            names = generate_random_protein_ids(protein_identifier_length, len(protein_sequences))
        else:
            names = [f"Protein_{i}" for i in range(len(protein_sequences))]
        sequences=protein_sequences
    elif isinstance(protein_sequences, dict):
        names = list(protein_sequences.keys())
        sequences = list(protein_sequences.values())
    else:
        raise ValueError("protein_sequences must be a string, list, or dict.")

    # make sure names and sequences are the same length
    if len(names) != len(sequences):
        raise ValueError("Number of protein sequences does not equal to number of portein names..")
    
    # verify unique protein sequences
    if verify_unique_protein_sequences:
        if len(sequences) != len(set(sequences)):
            raise ValueError("Protein sequences are not unique.")
    

    # check codon_frequency_table
    if isinstance(codon_frequency_table, str):
        # make sure is lowercase
        codon_frequency_table = codon_frequency_table.lower()
        if codon_frequency_table not in all_codon_tables:
            raise ValueError(f"Codon frequency table {codon_frequency_table} not found.")
        else:
            codon_frequency_table = all_codon_tables[codon_frequency_table]
    elif not isinstance(codon_frequency_table, dict):
        raise ValueError("Codon frequency table must be a string or a dictionary.")
    

    # verify start and stop codons in adapters
    if verify_start_stop_codons_in_adapters:
        if adapter_5_prime != None:
            if adapter_5_prime[-3:] != 'ATG':
                raise ValueError("5' adapter does not have start codon.")
            # check for no nucleotides in adapters.
            if len(set(adapter_5_prime) - set(['A', 'T', 'G', 'C'])) > 0:
                raise ValueError("5' adapter has non-nucleotide characters.")
        else:
            adapter_5_prime=""
        if adapter_3_prime != None:
            if adapter_3_prime[:3] not in ['TAA', 'TAG', 'TGA']:
                raise ValueError("3' adapter does not have stop codon.")
            # check for no nucleotides in adapters.
            if len(set(adapter_3_prime) - set(['A', 'T', 'G', 'C'])) > 0:
                raise ValueError("3' adapter has non-nucleotide characters.")
        else:
            adapter_3_prime=""
    if adapter_3_prime==None:
        adapter_3_prime=""
    if adapter_5_prime==None:
        adapter_5_prime=""

    # make sure total_length is an int
    if total_length != None:
        # make sure specified total length is possible given the sequences and adapters
        if not isinstance(total_length, int):
            raise ValueError("total_length must be an integer.")
        if total_length <= 0:
            raise ValueError("total_length must be greater than 0.")        
        # get the longest sequence in sequences
        longest_sequence = max(sequences, key=len)
        # see if total length is possible given the longest sequence, three nucleotides
        # per amino acid, and the adapters.
        min_total_length = len(longest_sequence) * 3
        min_total_length += len(adapter_5_prime)
        min_total_length += len(adapter_3_prime)
        if total_length < min_total_length:
            raise ValueError(f"Total length {total_length} is less than minimum length {min_total_length} given the sequences and adapters.")

    # if target_gc_range is None, set to (0,1)
    if target_gc_range is None:
        target_gc_range = (0,1)    

    # make sure target_gc_range is a tuple and between 0 and 1
    if not isinstance(target_gc_range, tuple):
        raise ValueError("target_gc_range must be a tuple.")
    if target_gc_range[0] < 0 or target_gc_range[0] > 1:
        raise ValueError("target_gc_range must be between 0 and 1.")
    if target_gc_range[1] < 0 or target_gc_range[1] > 1:
        raise ValueError("target_gc_range must be between 0 and 1.")
    
    # make sure pad_location is None or an int
    if pad_location != None:
        if pad_location not in [3, 5]:
            raise ValueError("pad_location must be 3 or 5.")

    # make sure avoid_added_start codons is a bool
    if not isinstance(avoid_adding_start_codons, bool):
        raise ValueError("avoid_added_start_codons must be a boolean.")

    # make sure avoid_added_stop codons is a bool
    if not isinstance(avoid_adding_stop_codons, bool):
        raise ValueError("avoid_added_stop_codons must be a boolean.")

    # make sure optimization_attempts is an int
    if not isinstance(optimization_attempts, int):
        raise ValueError("optimization_attempts must be an integer.")
    
    # make sure optimization_attempts is greater than 0
    if optimization_attempts <= 0:
        raise ValueError("optimization_attempts must be greater than 0.")

    # make sure gc_finetuning_iterations is an int
    if not isinstance(gc_finetuning_iterations, int):
        raise ValueError("gc_finetuning_iterations must be an integer.")
    # make sure gc_finetuning_iterations is greater than 0
    if gc_finetuning_iterations <= 0:
        raise ValueError("gc_finetuning_iterations must be greater than 0.")
    
    # make sure gc_weight is a float
    if not isinstance(gc_weight, (int, float)):
        raise ValueError("gc_weight must be a float.")
    # make sure gc_weight is greater than 0
    if gc_weight <= 0:
        raise ValueError("gc_weight must be greater than 0.")
    
    # make sure codon_weight is a float
    if not isinstance(codon_weight, (int, float)):
        raise ValueError("codon_weight must be a float.")
    # make sure codon_weight is greater than 0
    if codon_weight <= 0:
        raise ValueError("codon_weight must be greater than 0.")
    
    # make sure verify_coding_sequence is a bool
    if not isinstance(verify_coding_sequence, bool):
        raise ValueError("verify_coding_sequence must be a boolean.")

    # make sure minimum_codon_probability is a float
    if minimum_codon_probability != None:
        if not isinstance(minimum_codon_probability, (int, float)):
            raise ValueError("minimum_codon_probability must be a float.")
        # make sure minimum_codon_probability is between 0 and 1
        if minimum_codon_probability < 0 or minimum_codon_probability > 1:
            raise ValueError("minimum_codon_probability must be between 0 and 1.")
    else:
        minimum_codon_probability = 0.0
    
    # make sure show_progress is a bool
    if not isinstance(show_progress, bool):
        raise ValueError("show_progress must be a boolean.")
    
    # make sure gc_tolerance is a float
    if not isinstance(gc_tolerance, (int, float)):
        raise ValueError("gc_tolerance must be a float.")
    # make sure gc_tolerance is between 0 and 1
    if gc_tolerance < 0 or gc_tolerance > 1:
        raise ValueError("gc_tolerance must be between 0 and 1.")
    
    # make sure verify_unique_protein_sequences is a bool
    if not isinstance(verify_unique_protein_sequences, bool):
        raise ValueError("verify_unique_protein_sequences must be a boolean.")
    
    # make sure verify_start_stop_codons_in_adapters is a bool
    if not isinstance(verify_start_stop_codons_in_adapters, bool):
        raise ValueError("verify_start_stop_codons_in_adapters must be a boolean.")
    
    # make sure add_protein_identifiers is a bool
    if not isinstance(add_protein_identifiers, bool):
        raise ValueError("add_protein_identifiers must be a boolean.")
    
    # make sure protein_identifier_length is an int
    if not isinstance(protein_identifier_length, int):
        raise ValueError("protein_identifier_length must be an integer.")
    
    # make sure protein_identifier_length is greater than 0
    if protein_identifier_length <= 0:
        raise ValueError("protein_identifier_length must be greater than 0.")
    

    # make sure codon_frequency_table is a dictionary
    if not isinstance(codon_frequency_table, dict):
        raise ValueError("codon_frequency_table must be a dictionary.")
    
    # make sure show_optimization_progress is a bool
    if not isinstance(show_optimization_progress, bool):
        raise ValueError("show_optimization_progress must be a boolean.")

    # initialize codon table. 
    codon_table_obj = CodonTable(codon_frequency_table,
                                    codon_weight,
                                    gc_weight,
                                    target_gc_range,
                                    minimum_codon_probability)
    
    # make everything a sequence object
    sequence_objects=[]
    for i in range(len(sequences)):
        sequence_objects.append(Sequence(sequences[i], codon_table_obj,
                                         names[i], adapter_3_prime,
                                         adapter_5_prime, 
                                         avoid_adding_start_codons,
                                         avoid_adding_stop_codons,
                                         total_length,
                                         pad_location,
                                         force_start_codon=force_start_codon,
                                         force_stop_codon=force_stop_codon,
                                         verify_coding_sequence=verify_coding_sequence,
                                         return_best=return_best,
                                         padding_attempts=padding_attempts))

    # Shoutout to Alex Holehouse for the progress bar magic. 
    if show_progress:
        pbar=tqdm(total=len(sequence_objects), position=0, desc='Progress through sequences', leave=True)


    # iterate over protein sequences in sequence_objects
    for seq_obj in sequence_objects:
        optimized_coding_seq = optimize_codon_usage(seq_obj.protein_sequence,
                                            codon_table_obj,
                                            n_iter=optimization_attempts,
                                            fine_tuning_iterations=gc_finetuning_iterations,
                                            return_best=return_best,
                                            show_progress_bar=show_optimization_progress)
        seq_obj.add_coding_sequence(optimized_coding_seq)

        # now take care of padding.
        if total_length != None:
            seq_obj.generate_padding()
        
        # update progress bar
        if show_progress:
            pbar.update(1)
    
    # close progress bar
    if show_progress:
        pbar.close()

    return sequence_objects


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





