# Kea

## Overview

Kea is a Python library for codon optimization that also accounts for GC content variability. 
It helps you design DNA sequences that encode your target proteins while controlling for:
- Codon usage frequencies in your organism of interest
- GC content within a desired range
- Consistent total sequence length
- Proper integration of adapters
- Avoidance of unwanted sequence features

## Important Note
KEA was originally made as an internal mini project to make going from protein sequence --> library easier. I do plan on improving KEA and adding additional functionality as time allows. However, it is currently focused on specific things that I need it for.

## Installation

```bash
pip install git+https://git@github.com/ryanemenecker/kea.git
```



## Basic Usage

The main functionality of Kea is provided through the `build_library()` function. Here are some examples of how to use it:

### Simple Example

```python
from kea import build_library

# Single protein sequence
sequence = "MKKFLVLLFCWAVLCEHN"
results = build_library(sequence, "yeast")
```

### Multiple Sequences with Options

```python
sequences = [
    "MKKFLVLLFCWAVLCEHN",
    "MVLSEGEWQLVLHVWAKV"
]

results = build_library(
    sequences,
    "yeast",
    target_gc_range=(0.4, 0.6),
    total_length=120,
    adapter_5_prime="GGTCTC",
    adapter_3_prime="GAGACC"
)
```

### Named Sequences

```python
sequences = {
    "protein1": "MKKFLVLLFCWAVLCEHN",
    "protein2": "MVLSEGEWQLVLHVWAKV"
}

results = build_library(
    sequences,
    "yeast",
    target_gc_range=(0.45, 0.55)
)
```

### Saving Results

```python
from kea import build_library, save_library

results = build_library(sequences, "yeast")
save_library(results, "optimized_sequences.csv")
```

## Key Parameters

- `protein_sequences`: Input protein sequences (string, list, or dict)
- `codon_frequency_table`: Codon usage table ("yeast", "s288c", "s288c_unweighted", or custom dict)
- `adapter_5_prime`: 5' adapter sequence to add to each DNA sequence (optional)
- `adapter_3_prime`: 3' adapter sequence to add to each DNA sequence (optional)
- `total_length`: Target length for all final DNA sequences (optional)
- `force_start_codon`: Ensure sequences begin with start codon (default: True)
- `force_stop_codon`: Ensure sequences end with stop codon (default: True)
- `target_gc_range`: Desired GC content range as tuple (min, max) (default: (0,1))
- `pad_location`: Where to add padding (3, 5, or None for both ends)
- `avoid_adding_start_codons`: Avoid start codons in padding (default: True)
- `avoid_adding_stop_codons`: Avoid stop codons in padding (default: True)
- `optimization_attempts`: Number of iterations for codon optimization (default: 5000)
- `gc_finetuning_iterations`: Number of iterations for GC content fine-tuning (default: 2000)
- `padding_attempts`: Maximum attempts to generate padding sequences (default: 10000)
- `gc_weight`: Weight factor for GC content optimization (default: 1)
- `codon_weight`: Weight factor for codon usage optimization (default: 1)
- `verify_coding_sequence`: Verify translation of coding sequences (default: True)
- `minimum_codon_probability`: Minimum probability threshold for codon selection (default: None)
- `show_progress`: Display overall progress bar (default: True)
- `gc_tolerance`: Maximum acceptable deviation from target GC content (default: 0.025)
- `verify_unique_protein_sequences`: Check that input protein sequences are unique (default: True)
- `verify_start_stop_codons_in_adapters`: Verify adapter sequence requirements (default: False)
- `show_optimization_progress`: Display progress during sequence optimization (default: True)
- `add_protein_identifiers`: Generate random identifiers for protein sequences (default: False)
- `protein_identifier_length`: Length of random protein identifiers (default: 8)
- `return_best`: Return best optimization result instead of final result (default: False)

## Output

The function returns a list of Sequence objects containing:

- `protein_sequence`: Final protein sequence
- `full_dna_sequence`: Complete DNA sequence with adapters and padding
- `coding_sequence`: DNA sequence encoding the protein
- `gc_content_full_sequence`: GC content of complete sequence
- `gc_content_coding_sequence`: GC content of coding sequence
- `correct_full_translation`: Translation verification of full sequence
- `correct_coding_translation`: Translation verification of coding sequence

#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.10.

### Copyright

Copyright (c) 2025, Ryan Emenecker WUSM Holehouse Lab