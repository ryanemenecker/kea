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
    Function to build a library of nucleotide sequences from a list of protein sequences.
    The function will add adapters, padding, and optimize for GC content and codon usage.
    The function will also verify that the resulting sequences hold the expected coding sequences.

    Parameters
    ----------
    protein_sequences (list or dict): List or dict of protein sequences.
        if dict, expects {name: sequence}
    codon_frequency_table (str or dict): Codon frequency table to use for optimization.
        If str, must be one of ['yeast', 's288c', 's288c_unweighted']
        If dict, must be a dictionary of codon frequencies.
        Structured: {amino_acid: {codon: frequency}}
    adapter_5_prime (str): 5' adapter sequence.
    adapter_3_prime (str): 3' adapter sequence.
    total_length (int): Target length of the final nucleotide sequence.
    force_start_codon (bool): Whether to force a start codon.
    force_stop_codon (bool): Whether to force a stop codon.
    target_gc_range (tuple): Target GC content range.
    pad_location (int): Location to add padding.
    avoid_added_start_codons (bool): Whether to avoid added start codons.
    avoid_added_stop_codons (bool): Whether to avoid added stop codons.
    optimization_attempts (int): Number of attempts for optimization.
    gc_finetuning_iterations (int): Number of iterations for fine-tuning GC content after main optimization.
    gc_weight (float): Weight for GC content optimization.
    codon_weight (float): Weight for codon usage optimization.
    minimum_codon_probability (float): Minimum codon probability for optimization.
    show_progress (bool): Whether to show progress. Default is True.
    show_optimization_progress (bool): Whether to show optimization progress. Default is True.
    gc_tolerance (float): Tolerance for GC content optimization. Default is 0.025.
    verify_unique_protein_sequences (bool): Whether to verify unique protein sequences. Default is True.
    verify_start_stop_codons_in_adapters (bool): Whether to verify start and stop codons in adapters. Default is False.
    add_protein_identifiers (bool): Whether to add protein identifiers to the output. Default is False.
    protein_identifier_length (int): Length of protein identifiers. Default is 8.
    verify_coding_sequence (bool): Whether to verify the coding sequence. Default is True.
    Returns
    -------
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









def save_library(name_to_protein_and_nucleotide_dict,
                 save_path):
    '''
    Function to save a library of nucleotide sequences to a file.
    
    Parameters
    ----------
    name_to_protein_and_nucleotide_dict (dict): Dictionary of protein names to nucleotide sequences.
    save_path (str): Path to save the library.
    
    Returns
    -------
    None
    '''
    # make sure the path exists
    if not os.path.exists(save_path):
        raise ValueError(f"Path {save_path} does not exist.")
    # make sure the path excluding the file name is a dir
    if not os.path.isdir(os.path.dirname(save_path)):
        raise ValueError(f"Path {os.path.dirname(save_path)} is not a directory.")
    
    # headers
    headers = ['Protein Name', 'Protein Sequence', 'Nucleotide Sequence', 'Translated Sequence']
    # open file
    with open(save_path, 'w') as f:
        # write headers
        f.write(','.join(headers) + '\n')
        # write data
        for protein_name, protein_to_nucleotide in name_to_protein_and_nucleotide_dict.items():
            for protein_sequence, nucleotide_sequence in protein_to_nucleotide.items():
                translated_sequence = translate_sequence(nucleotide_sequence, return_stop_codon=True, return_nucleotide_sequence=True)
                f.write(f"{protein_name},{protein_sequence},{nucleotide_sequence},{translated_sequence}\n")
    f.close()


    
    

