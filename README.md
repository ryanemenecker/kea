<div align="center">

# Kea

**Codon optimization for protein libraries — with control over GC content.**

[![CI](https://github.com/ryanemenecker/kea/actions/workflows/CI.yaml/badge.svg)](https://github.com/ryanemenecker/kea/actions/workflows/CI.yaml)
[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

</div>

---

Kea is a Python library that turns protein sequences into codon-optimized DNA
sequences. Unlike most codon optimizers, it jointly optimizes for **organism
codon usage** *and* a **target GC-content range**, making it well suited to
designing balanced oligo libraries.

## Features

- 🧬 **Codon optimization** against organism-specific usage tables (yeast / *S. cerevisiae* S288C built in, or supply your own).
- ⚖️ **GC-content control** — optimize sequences to fall within a target GC range, with tunable tolerance.
- 📏 **Fixed-length libraries** — automatically pad sequences to a common length while avoiding spurious start/stop codons.
- 🔗 **Adapter integration** — add 5′ and 3′ adapters, with optional validation against unwanted start/stop codons.
- ✅ **Built-in verification** — every sequence is translated back and checked against the target protein.
- 💾 **One-line export** to CSV for downstream ordering and analysis.

## Installation

```bash
pip install git+https://github.com/ryanemenecker/kea.git
```

Kea supports Python 3.8+ and depends only on `numpy`, `tqdm`, and `matplotlib`.

## Quick Start

```python
from kea import build_library

results = build_library("MKKFLVLLFCWAVLCEHN", "yeast")

seq = results[0]
print(seq.coding_sequence)               # codon-optimized DNA
print(seq.gc_content_coding_sequence)    # GC content of the coding sequence
print(seq.correct_coding_translation)    # True — verified round-trip
```

`build_library` accepts a single sequence, a list, or a `{name: sequence}`
dictionary, and returns one `Sequence` object per input protein.

## Usage

### Multiple sequences with options

```python
sequences = [
    "MKKFLVLLFCWAVLCEHN",
    "MVLSEGEWQLVLHVWAKV",
]

results = build_library(
    sequences,
    "yeast",
    target_gc_range=(0.4, 0.6),
    total_length=120,
    adapter_5_prime="GGTCTC",
    adapter_3_prime="GAGACC",
)
```

### Named sequences

```python
sequences = {
    "protein1": "MKKFLVLLFCWAVLCEHN",
    "protein2": "MVLSEGEWQLVLHVWAKV",
}

results = build_library(sequences, "yeast", target_gc_range=(0.45, 0.55))
```

### Saving results

```python
from kea import build_library, save_library

results = build_library(sequences, "yeast")
save_library(results, "optimized_sequences.csv")
```

## Output

`build_library` returns a list of `Sequence` objects, each exposing:

| Attribute                    | Description                                          |
| ---------------------------- | ---------------------------------------------------- |
| `protein_sequence`           | Final protein sequence (including any forced start/stop) |
| `full_dna_sequence`          | Complete DNA sequence with adapters and padding      |
| `coding_sequence`            | DNA sequence encoding the protein                    |
| `gc_content_full_sequence`   | GC content of the complete sequence                  |
| `gc_content_coding_sequence` | GC content of the coding sequence                    |
| `correct_full_translation`   | Whether the full sequence translates correctly       |
| `correct_coding_translation` | Whether the coding sequence translates correctly     |

## Key Parameters

| Parameter             | Default     | Description                                                                 |
| --------------------- | ----------- | --------------------------------------------------------------------------- |
| `protein_sequences`   | *required*  | Input protein sequences (string, list, or dict).                            |
| `codon_frequency_table` | *required* | `"yeast"`, `"s288c"`, `"s288c_unweighted"`, or a custom `{aa: {codon: freq}}` dict. |
| `target_gc_range`     | `(0, 1)`    | Desired GC-content range as `(min, max)`.                                    |
| `total_length`        | `None`      | Target length for all final sequences (padding is added to reach it).       |
| `adapter_5_prime`     | `None`      | 5′ adapter sequence to prepend to each DNA sequence.                         |
| `adapter_3_prime`     | `None`      | 3′ adapter sequence to append to each DNA sequence.                         |
| `pad_location`        | `None`      | Where to add padding: `5`, `3`, or `None` (split across both ends).          |
| `gc_tolerance`        | `0.025`     | Maximum acceptable deviation from the target GC range.                       |

<details>
<summary><b>All parameters</b></summary>

- `force_start_codon`: Ensure sequences begin with a start codon (default: `True`)
- `force_stop_codon`: Ensure sequences end with a stop codon (default: `True`)
- `avoid_adding_start_codons`: Avoid start codons in padding (default: `True`)
- `avoid_adding_stop_codons`: Avoid stop codons in padding (default: `True`)
- `optimization_attempts`: Number of iterations for codon optimization (default: `5000`)
- `gc_finetuning_iterations`: Number of iterations for GC fine-tuning (default: `2000`)
- `padding_attempts`: Maximum attempts to generate padding sequences (default: `10000`)
- `gc_weight`: Weight factor for GC-content optimization (default: `1`)
- `codon_weight`: Weight factor for codon-usage optimization (default: `1`)
- `verify_coding_sequence`: Verify translation of coding sequences (default: `True`)
- `minimum_codon_probability`: Minimum probability threshold for codon selection (default: `None`)
- `show_progress`: Display overall progress bar (default: `True`)
- `verify_unique_protein_sequences`: Check that input protein sequences are unique (default: `True`)
- `verify_start_stop_codons_in_adapters`: Verify adapter sequence requirements (default: `False`)
- `show_optimization_progress`: Display progress during sequence optimization (default: `True`)
- `add_protein_identifiers`: Generate random identifiers for protein sequences (default: `False`)
- `protein_identifier_length`: Length of random protein identifiers (default: `8`)
- `return_best`: Return best optimization result instead of final result (default: `False`)
- `early_stop_threshold`: Stop optimization early when a sequence reaches this fraction of the maximum score (default: `0.95`)

</details>

## Project Status

Kea began as an internal tool to streamline going from protein sequence to oligo
library. It is actively maintained and focused on the features its authors need;
additional functionality (more organisms, restriction-site avoidance, mRNA
structure considerations, and more) is planned. Contributions and issues are
welcome.

## Contributing

Contributions are welcome! Please see the
[contributing guidelines](.github/CONTRIBUTING.md) and
[code of conduct](CODE_OF_CONDUCT.md). To run the test suite locally:

```bash
pip install -e ".[test]"
pytest
```

## License

Released under the [MIT License](LICENSE).

Copyright © 2025-2026 Ryan Emenecker, Holehouse Lab, Washington University School of Medicine.
