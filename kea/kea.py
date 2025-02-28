"""
Functionality for building the library of nucleotide sequences 
from DNA sequences
"""
from .backend.optimize_codon_usage import optimize_codon_usage_for_library
from .backend.library_generation_utils import add_padding, add_3_prime_adapter, add_5_prime_adapter
from .backend.translation_utils import translate_sequence
from .backend.kea_utils import generate_random_protein_ids
from .data.codon_tables import all_codon_tables

def build_library(protein_sequences,
                  codon_frequency_table,
                  adapter_5_prime=None,
                  adapter_3_prime=None,
                  total_length=None,
                  force_start_codon=False,
                  force_stop_codon=False,
                  target_gc_range=None,
                  pad_location=None,
                  avoid_added_start_codons=True,
                  avoid_added_stop_codons=True,
                  optimization_attempts=1000,
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
                  protein_identifier_length=8):
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
    # Initialize return_dict flag
    return_dict = False
    
    # Check if protein_sequences is a string
    if isinstance(protein_sequences, str):
        protein_sequences = [protein_sequences]
        if add_protein_identifiers:
            return_dict=True
            protein_names = generate_random_protein_ids(len(protein_sequences), protein_identifier_length)
            protein_to_name_dict={}
            for i, seq in enumerate(protein_sequences):
                protein_to_name_dict[seq] = protein_names[i]

    if isinstance(protein_sequences, dict):
        # Store original dict for later reference
        original_dict = protein_sequences.copy()
        protein_names = list(original_dict.keys())
        protein_sequences = list(original_dict.values())
        protein_to_name_dict={}
        return_dict=True
        for i, seq in enumerate(protein_sequences):
            protein_to_name_dict[seq] = protein_names[i]

    if isinstance(protein_sequences, list) and add_protein_identifiers:
        protein_names = generate_random_protein_ids(len(protein_sequences), protein_identifier_length)
        protein_to_name_dict={}
        return_dict=True
        for i, seq in enumerate(protein_sequences):
            protein_to_name_dict[seq] = protein_names[i]

    # verify unique protein sequences
    if verify_unique_protein_sequences:
        if len(protein_sequences) != len(set(protein_sequences)):
            raise ValueError("Protein sequences are not unique.")
    
    # verify start and stop codons in adapters
    if verify_start_stop_codons_in_adapters:
        if adapter_5_prime != None:
            if adapter_5_prime[-3:] != 'ATG':
                raise ValueError("5' adapter does not have start codon.")
        if adapter_3_prime != None:
            if adapter_3_prime[:3] not in ['TAA', 'TAG', 'TGA']:
                raise ValueError("3' adapter does not have stop codon.")
    
    # check codon_frequency_table
    if isinstance(codon_frequency_table, str):
        # make sure is lowercase
        codon_frequency_table= codon_frequency_table.lower()
        if codon_frequency_table not in all_codon_tables:
            raise ValueError(f"Codon frequency table {codon_frequency_table} not found.")
        else:
            codon_frequency_table = all_codon_tables[codon_frequency_table]
    elif not isinstance(codon_frequency_table, dict):
        raise ValueError("Codon frequency table must be a string or a dictionary.")

    
    # if target_gc_range is None, set to (0,1)
    if target_gc_range is None:
        target_gc_range = (0,1)
    
    # see if we need to add start codon
    if force_start_codon:
        temp=[]
        for seq in protein_sequences:
            if seq[0] != 'M':
                seq = 'M' + seq
            temp.append(seq)
        protein_sequences = temp
    # see if we need to add stop codon
    if force_stop_codon:
        temp=[]
        for seq in protein_sequences:
            if seq[-1] != '*':
                seq = seq + '*'
            temp.append(seq)
        protein_sequences = temp

    # optimize the sequences
    seq_to_nt_dict = optimize_codon_usage_for_library(protein_sequences,
                                                        codon_frequency_table,
                                                        gc_range=target_gc_range,
                                                        n_iter=optimization_attempts,
                                                        gc_weight=gc_weight,
                                                        usage_weight=codon_weight, 
                                                        return_best=True,
                                                        minimum_codon_probability=minimum_codon_probability,
                                                        show_progress_bar=show_progress,
                                                        return_protien_nucleotide_dict=True,
                                                        show_optimization_progress=show_optimization_progress)
    
    # see if we need to verify the coding sequence
    if verify_coding_sequence:
        for protein, nt_seq in seq_to_nt_dict.items():
            if protein != translate_sequence(nt_seq, return_stop_codon=True):
                raise ValueError(f"Protein sequence {protein} does not match nucleotide sequence {nt_seq}")
    
    # now see if we need to add padding. 
    if total_length != None:
        # make new dict
        new_dict = {}
        # iterate over seq_to_nt_dict
        for protein, nt_seq in seq_to_nt_dict.items():
            # figure out how much padding we need to add. This will be length
            # of the nt_seq - length of 3' adapter if there is one - length of 5' adapter
            # if there is one.
            if adapter_3_prime!=None:
                adapter_3_prime_length = len(adapter_3_prime)
            else:
                adapter_3_prime_length = 0
            if adapter_5_prime!=None:
                adapter_5_prime_length = len(adapter_5_prime)
            else:
                adapter_5_prime_length = 0
            padding_length = total_length - len(nt_seq) - adapter_3_prime_length - adapter_5_prime_length
            # if padding length is less than 0, raise error
            if padding_length < 0:
                raise ValueError(f"Total length {total_length} is less than length of nucleotide sequence {nt_seq}")
            # if padding length is greater than 0...
            if padding_length > 0:
                # if pad_location is None, add padding to each end. 
                if pad_location == None:
                    # figure out length of 5' padding
                    five_prime_padding_length = padding_length // 2
                    three_prime_padding_length = padding_length - five_prime_padding_length
                    # add padding to each end
                    nt_seq = add_padding(nt_seq, len(nt_seq)+five_prime_padding_length, 
                                         pad_location=5,
                                         avoid_added_start_codons=avoid_added_start_codons,
                                         avoid_added_stop_codons=avoid_added_stop_codons,
                                         tolerance=gc_tolerance)
                    nt_seq = add_padding(nt_seq, len(nt_seq)+three_prime_padding_length,    
                                         pad_location=3,
                                         avoid_added_start_codons=avoid_added_start_codons,
                                         avoid_added_stop_codons=avoid_added_stop_codons,
                                         tolerance=gc_tolerance)
                else:
                    # add padding to specified location
                    nt_seq = add_padding(nt_seq, len(nt_seq)+padding_length, 
                                         pad_location=pad_location,
                                         avoid_added_start_codons=avoid_added_start_codons,
                                         avoid_added_stop_codons=avoid_added_stop_codons,
                                         tolerance=gc_tolerance)
                    
            
            # see if we need to verify. Do this at each sequence to fail early. 
            if verify_coding_sequence:
                # grab protein.
                if protein != translate_sequence(nt_seq, return_stop_codon=True):
                    raise ValueError(f"Protein sequence {protein} does not match nucleotide sequence {nt_seq}")
                        
            # add nt_seq to new dict
            new_dict[protein] = nt_seq

        # set seq_to_nt_dict to new dict
        seq_to_nt_dict = new_dict

    # make new dict
    new_dict = {}
    # iterate over nt_to_seq_dict
    for protein, nt_seq in seq_to_nt_dict.items():
        # now see if we need to add adapters
        if adapter_5_prime != None:
            nt_seq = add_5_prime_adapter(nt_seq, adapter_5_prime)
        if adapter_3_prime != None:
            nt_seq = add_3_prime_adapter(nt_seq, adapter_3_prime)
        new_dict[protein] = nt_seq
    
    # if we are returning a dict, format it as name:{protein_sequence:nt_sequence}
    final_dict = {}  # Initialize final_dict
    if return_dict:
        for protein, nt_seq in new_dict.items():  # Fixed: iterate over new_dict instead of final_dict
            final_dict[protein_to_name_dict[protein]] = {protein: nt_seq}
    else:
        final_dict = new_dict

    # make sure len new_dict is same as protein_sequences
    if len(final_dict) != len(protein_sequences):
        raise ValueError("Number of input protein sequences does not match number of generated nucleotide sequences")
    
    return final_dict  # Added explicit return statement

