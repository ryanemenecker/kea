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

Requires Python 3.8+. `numpy` and `tqdm` are installed automatically. The plotting
helpers in `kea.backend.kea_plotting_utils` additionally need matplotlib:

```bash
pip install "kea[plotting] @ git+https://git@github.com/ryanemenecker/kea.git"
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

### Reproducible Libraries

Pass a `seed` to get the same library back every time:

```python
results = build_library(sequences, "yeast", target_gc_range=(0.45, 0.55), seed=42)
```

### Multiple encodings of the same protein

Request several synonymous coding sequences for each input protein and enforce
an all-pairs minimum nucleotide Hamming distance:

```python
from kea import build_library, nucleotide_hamming_distance

variants = build_library(
    {"MY_PROTEIN": "MKKFLVLLFCWAVLCEHN"},
    "human",
    sequences_per_protein=4,
    minimum_hamming_distance=15,
    hamming_distance_attempts=100,
    avoid_human_splice_sites=True,
    avoid_premature_polyadenylation=True,
    seed=42,
)

for i, left in enumerate(variants):
    for right in variants[:i]:
        assert nucleotide_hamming_distance(
            left.coding_sequence, right.coding_sequence
        ) >= 15
```

The distance is an absolute nucleotide count over the coding sequence, including
forced start/stop codons but excluding adapters and padding. Outputs are named
`MY_PROTEIN_variant_1`, `MY_PROTEIN_variant_2`, and so on. Every variant passes
the same GC, splice, PAS, repeat, translation, and codon-adaptation checks. If
Kea cannot fill the requested set, accepted variants are retained and unfilled
slots are reported with stage `sequence_diversity`; use `on_error="raise"` to
raise `SequenceDiversityError` instead. A search failure is not a mathematical
proof that no set exists; short proteins, proteins dominated by Met/Trp, strict
codon filtering, and competing sequence constraints can all reduce the available
synonymous distance.

## Key Parameters

- `protein_sequences`: Input protein sequences (string, list, or dict)
- `codon_frequency_table`: Codon usage table ("yeast", "s288c", "s288c_unweighted", "human", "human_unweighted", or custom dict)
- `adapter_5_prime`: 5' adapter sequence to add to each DNA sequence (optional)
- `adapter_3_prime`: 3' adapter sequence to add to each DNA sequence (optional)
- `total_length`: Target length for all final DNA sequences (optional)
- `force_start_codon`: Ensure sequences begin with start codon (default: True)
- `force_stop_codon`: Ensure sequences end with stop codon (default: True)
- `target_gc_range`: Desired GC content range as tuple (min, max), both 0-1 (default: None). Leave unset to skip GC targeting entirely. Passing `(0, 1)` explicitly is **not** the same thing — it keeps GC optimization active with a range everything satisfies. Any sequence returned is guaranteed to fall inside this range unless `return_best=True`.
- `pad_location`: Where to add padding (3, 5, or None for both ends)
- `avoid_adding_start_codons`: Avoid start codons in padding (default: True)
- `avoid_adding_stop_codons`: Avoid stop codons in padding (default: True)
- `optimization_attempts`: Number of iterations for codon optimization (default: 5000)
- `gc_finetuning_iterations`: Number of iterations for GC content fine-tuning (default: 2000)
- `padding_attempts`: Maximum attempts to generate padding sequences (default: 10000)
- `gc_weight`: Weight of the GC term in the objective (default: 1). Must be positive when `target_gc_range` is set.
- `codon_weight`: Weight of the codon usage term in the objective (default: 1). The two weights are combined as a weighted average, so only their ratio matters — `(1, 1)` and `(5, 5)` behave identically, while `codon_weight=5, gc_weight=1` favours codon usage 5:1.
- `verify_coding_sequence`: Verify translation of coding sequences (default: True)
- `minimum_codon_probability`: Minimum probability threshold for codon selection (default: None)
- `show_progress`: Display overall progress bar (default: True)
- `gc_tolerance`: Maximum acceptable deviation from target GC content (default: 0.025)
- `verify_unique_protein_sequences`: Check that input protein sequences are unique (default: True)
- `verify_start_stop_codons_in_adapters`: Verify adapter sequence requirements (default: False)
- `show_optimization_progress`: Display progress during sequence optimization (default: True)
- `add_protein_identifiers`: Generate random identifiers for protein sequences (default: False)
- `protein_identifier_length`: Length of random protein identifiers (default: 8)
- `return_best`: Return the best result found even if it falls outside `target_gc_range` (default: False). With the default, a GC range that cannot be reached raises a `ValueError` rather than quietly handing back an off-target sequence.
- `on_error`: `"skip"` (default) records a per-sequence failure and continues; `"raise"` aborts the whole build. See [Handling failures in large libraries](#handling-failures-in-large-libraries).
- `constraint_repair_attempts`: Maximum rounds of synonymous repair per sequence (default: 100)
- `constraint_retry_attempts`: Extra rounds of re-optimization to try when repair fails (default: 2). Set to 0 for single-shot behaviour.
- `escape_gc_boundary`: How hard to push codon adaptation when GC is pinned against a range edge (default: None). `None` escalates per sequence — high-adaptation search first, permissive fallback only where constraints cannot be met with it. `True` forces the high-adaptation search, `False` forces the permissive one. See [Escaping the GC boundary](#escaping-the-gc-boundary).
- `minimum_codon_adaptation`: Reject a finished sequence whose codon adaptation score falls below this, 0-1 (default: None — measured and reported on every sequence, but not enforced).
- `transcribed_5_prime_context` / `transcribed_3_prime_context`: Transcribed vector/UTR bases to include when scanning for splice sites and PAS motifs, so junction-created sites are caught (default: None)
- `seed`: Seed for reproducible libraries (default: None). Codon sampling uses numpy's global RNG and padding uses the stdlib `random` module, so seeding just one of them yourself is not enough — pass this instead.
- `early_stop_threshold`: Stop optimizing a sequence once an in-range candidate reaches this blended codon-usage/GC score, 0-1 (default: None, i.e. disabled). Enabling it makes `optimization_attempts` an upper bound rather than a target; a threshold near 0.95 is reached almost immediately for most proteins and will cut optimization short.
- `sequences_per_protein`: Number of synonymous coding sequences requested for each input protein (default: 1).
- `minimum_hamming_distance`: Required all-pairs coding-sequence distance in differing nucleotide positions (default: 0).
- `hamming_distance_attempts`: Maximum randomized diversification attempts for each additional encoding (default: 100).

## Human transcript constraints

Sequence constraints are opt-in hard requirements. Kea synonymously repairs a
coding sequence after codon/GC optimization and raises an informative error when
the requested constraints cannot all be satisfied.

```python
results = build_library(
    sequences,
    "human",
    avoid_human_splice_sites=True,
    avoid_premature_polyadenylation=True,
    max_homopolymer_length=6,
    max_tandem_repeat_copies=3,
    max_direct_repeat_length=16,
    # Include transcribed vector/UTR bases so junction-created sites are caught.
    transcribed_5_prime_context="GCCACC",
    transcribed_3_prime_context="GCTAGC",
)
```

Human splice avoidance uses the published MaxEntScan models with configurable
thresholds:

- `maxent_donor_threshold=6.0`
- `maxent_acceptor_threshold=7.0`

Premature-polyadenylation avoidance removes `AATAAA`, `ATTAAA`, and 11 common
human PAS variants. Downstream U-rich and GU-rich context is annotated for every
hit. Set `polyadenylation_requires_downstream_element=True` for a less strict
mode that only removes context-supported PAS motifs.

Constraints are enforced on the generated padding as well as the coding sequence,
including sites created across the padding/CDS junction, so a `total_length`
construct does not smuggle back the motifs you excluded.

When repair fails, Kea re-optimizes and tries again (`constraint_retry_attempts`,
default 2) — codon optimization is stochastic and repair is a greedy local
search, so a different starting point often succeeds. Under tight GC plus several
constraints this rescues most otherwise-lost sequences. Violations that no
synonymous encoding can clear are detected by exhaustive enumeration of the
affected codons and fail immediately with a message saying so, rather than
burning retries on something impossible:

```
'ACTB_v3' cannot satisfy the requested constraints under any synonymous encoding,
so retrying will not help. Relax the relevant threshold or change the protein.
Unsatisfiable: human_splice_donor at 102:111 (GTGGTATGT) score=7.64 vs threshold 6.00
```

The complexity controls are independent and disabled when left as `None`:

- `max_homopolymer_length`
- `max_tandem_repeat_copies` and `max_tandem_repeat_unit_length`
- `max_direct_repeat_length`

These are predictive design filters, not guarantees of transcript processing or
expression. For important constructs, include the actual transcribed flanks and
validate splicing experimentally. See [roadmap.md](roadmap.md) for calibration
and planned expression-system-specific features.

## Escaping the GC boundary

When the codon-usage optimum lies outside `target_gc_range` — human usage sits
around 0.59 GC, so a requested (0.30, 0.42) is well below it — a naive optimizer
converges onto the range edge and stops. Every remaining codon swap that would
improve usage pushes GC out of range, so a one-codon-at-a-time search has nowhere
left to go.

Measured on converged sequences: **zero** single-codon improvements remained,
while roughly **3,000 GC-neutral pair moves per sequence** were still available —
improve usage at one position, pay the GC back at another. Kea makes those paired
moves, which drives the count of missed improvements to zero:

```
                      single-codon moves left    pair moves left
before                          0                     19,128
after                           0                          0
```

### How this interacts with constraints

Maximal codon optimization is more repetitive, which makes a sequence harder to
make constraint-compliant — raw optimizer output goes from 2.3 to 7.7 constraint
violations per sequence. Running the high-adaptation search for a whole library
therefore costs yield on tight GC ranges.

So Kea escalates **per sequence** instead of picking one setting for the library:
it runs the high-adaptation search first, and drops to the more permissive one
only for the sequences that cannot satisfy their constraints with it. Nothing is
lost by trying — across 120 proteins, zero sequences built with the boundary
escape that could not also be built without it.

Measured on 40–60 proteins of 180–200 aa with five constraints enabled, averaged
over independent seeds:

| `target_gc_range` | sequences built | codon adaptation | time |
|---|---|---|---|
| (0.30, 0.42) | unchanged | 0.857 → **0.893** | 3.8× |
| (0.28, 0.36) | unchanged | 0.812 → **0.830** | 4.2× |
| (0.20, 0.35) | unchanged | 0.807 → **0.820** | 4.3× |
| (0.45, 0.60) | unchanged | 0.995 → 0.997 | 1.6× |
| (0.35, 0.65) | unchanged | 0.999 → 0.999 | 1.5× |

The gain is largest exactly where it should be — ranges that sit below the
organism's preferred GC (about 0.59 for `"human"`, 0.35 for `"yeast"`) — and
vanishes on ranges that already contain it, where there was no boundary to escape.

`escape_gc_boundary` controls this: `None` (the default) escalates as described,
`True` forces the high-adaptation search only, `False` forces the permissive one
only. Forcing `True` with tight constraints buys a little more adaptation but
loses real yield (46/60 built versus 58/60 at (0.30, 0.42)), which is why it is
not the default.

Simply raising `constraint_retry_attempts` is not a substitute: at (0.28, 0.36),
going from 3 to 20 retries with the high-adaptation search only moved yield from
18/60 to 28/60 while runtime went from 19s to 72s. The sequences that fail are
structurally harder, not unlucky.

## Final quality check

Every sequence passes a single acceptance gate before it is returned, measuring
all targets together on the finished construct — not on an intermediate state.
Each stage optimizes its own objective and is allowed to trade the others off
temporarily; this is where the result is checked as a whole.

```python
print(library[0].quality_report)
# [ok  ] translation: n/a (target: coding sequence translates to the target protein)
# [ok  ] total_length: 400.0000 (target: exactly 400 nt)
# [ok  ] gc_coding: 0.5492 (target: 0.450-0.550)
# [ok  ] gc_full: 0.5020 (target: 0.450-0.550)
# [ok  ] sequence_constraints: 0.0000 (target: 0 violations across padding and coding sequence)
# [ok  ] codon_adaptation: 0.9710 (target: reported, no floor requested)
```

A failed required check makes the sequence a per-sequence failure, handled by
`on_error`. A *missed but permitted* target — asking for best-effort GC with
`return_best=True` — is recorded as a warning instead, so it is visible rather
than silent:

```python
library.with_warnings          # accepted sequences that missed a requested target
library[0].quality_report.warnings
# [warn] gc_coding: 0.6011 (target: 0.800-0.850)
```

Codon adaptation is measured on every sequence and reported. Set
`minimum_codon_adaptation` to make it a hard floor — useful when a tight
`target_gc_range` fights the organism's codon bias and you want a limit on how
much codon optimality gets traded away. With human codon usage and a
`target_gc_range` of (0.30, 0.42), for example, adaptation lands around 0.83–0.89
rather than the ~0.99 achievable at native GC.

## Handling failures in large libraries

Some proteins have no synonymous encoding that satisfies every requested
constraint. By default `build_library` records those and carries on, so one
impossible sequence at position 9,999 does not discard the 9,998 already built.

```python
library = build_library(proteins, "human", avoid_human_splice_sites=True)

print(library.summary())
# Built 9998/10000 sequences.
# 2 failed (constraint_repair: 2).
#   ACTB_v3 [constraint_repair] ConstraintRepairError: ...

for failure in library.failures:
    print(failure.name, failure.stage, failure.reason)
    for violation in failure.remaining_violations:
        print("   ", violation.kind, violation.start, violation.score)
```

The return value is a `LibraryResult`, which *is* a list of `Sequence` objects —
existing code that iterates or indexes it keeps working unchanged. It just also
carries `.failures`, `.n_succeeded`, `.n_failed`, `.failures_by_stage()` and
`.summary()`.

Pass `on_error="raise"` to abort the whole build on the first failure instead.
Configuration mistakes (a bad `pad_location`, `optimization_attempts=0`, an
invalid codon table) always raise regardless of `on_error`, because they apply to
every sequence rather than to one input.

## Saving a library

```python
from kea import build_library, save_library, save_failures

library = build_library(sequences, "human", target_gc_range=(0.45, 0.60))

save_library(library, "library.csv")            # full diagnostics (default)
save_library(library, "short.csv", include_diagnostics=False)
save_failures(library, "failures.csv")          # what didn't build, and why
```

`save_library` writes one row per sequence with the construct itself plus the
measurements you would screen on before ordering. The splice and polyadenylation
screens run on the **complete construct** (adapters + padding + coding sequence)
and are reported whether or not you enabled the matching constraint — so you can
screen a library that was built without those constraints.

## Reading the diagnostics

### Splice-site scores are not probabilities

`donor_max_score` and `acceptor_max_score` are **MaxEntScan scores**: log2 odds
ratios from the Yeo–Burge maximum-entropy model. A score of *s* means the window
is 2^*s* times more likely under the model of real human splice sites than under
a background model of ordinary sequence. They are not percentages and they are not
bounded — a score of 8 is not "8% likely to splice".

Read them by comparison. Measured over 107,400 donor and 105,720 acceptor windows
of ordinary human-codon-optimized coding sequence:

| percentile of ordinary CDS | donor | acceptor |
|---|---|---|
| median | −19.2 | −16.8 |
| 90th | −6.6 | −5.9 |
| 99th | 1.9 | 1.7 |
| 99.9th | 7.9 | 6.5 |
| 99.99th | 10.2 | 9.4 |

For scale, the perfect consensus donor `CAGGTAAGT` scores **10.86**, a strong
pyrimidine-tract acceptor scores **12.75**, and a G-rich window that could never
be an acceptor scores **−20.7**.

So in practice:

- **Below ~0** — unremarkable. Most of any coding sequence sits here.
- **0 to 6** — common. Roughly 1 in 50 windows of ordinary CDS scores above 0.
- **6 to 8** — uncommon; about 2.3 donor sites per kb of ordinary CDS. This is
  where Kea's default `maxent_donor_threshold` sits.
- **Above 8** — rare, about 0.9 donor sites per kb. Worth looking at.
- **Above 10** — top 0.01% of ordinary sequence, comparable to a real consensus
  site. Worth removing or testing.

`donor_sites_over_threshold` counts windows at or above the threshold used for
that build (`splice_donor_threshold` records which). `donor_max_position` is the
predicted junction position in the full construct — the exon/intron boundary for a
donor, the intron/exon boundary for an acceptor — and `donor_max_window` is the
scored window itself, so you can eyeball the `GT`/`AG` dinucleotide.

**Important caveat.** These thresholds are provisional design cutoffs, not
experimentally calibrated ones. A high score means "looks like a splice site to
this model", not "will be spliced"; a low score does not guarantee the transcript
is safe. Only the forward strand is scanned. Treat the numbers as a screening aid
and validate important constructs experimentally.

### Polyadenylation columns are motif counts, not scores

There is no score here — polyadenylation signals are **exact hexamer matches**
against the canonical `AATAAA`, the major variant `ATTAAA`, and 11 further common
human variants.

- `polya_signal_count` — how many of those hexamers appear anywhere in the construct.
- `polya_canonical_aataaa` — how many of them are the canonical `AATAAA`, which is
  the most frequently used signal in real transcripts and the one to prioritise.
- `polya_with_downstream_element` — how many have a U-rich or GU-rich element
  downstream. A hexamer on its own is weak evidence; hexamers are common by chance.
  Real cleavage usually needs that downstream element too, so this column is the
  better risk indicator.
- `polya_signals` — the first ten hits as `HEXAMER@position`.

A construct with `polya_signal_count` of 3 and `polya_with_downstream_element` of 0
is far less concerning than one with a single hit that does have the downstream
context.

### Sequence-complexity columns

These are synthesis and assembly risk indicators rather than expression ones, and
vendors differ in what they accept.

- `max_homopolymer_run` / `max_homopolymer_base` — longest single-base run.
- `longest_direct_repeat_nt` / `longest_direct_repeat` — longest substring that
  occurs more than once anywhere in the construct (occurrences may overlap).
- `max_tandem_repeat_copies` / `max_tandem_repeat_unit` — most consecutive copies
  of a repeating unit of 2–6 nt. Homopolymers are reported separately, above.

### The rest

- `codon_adaptation_index` — 0 to 1, where 1 means every codon is the most frequent
  choice for its amino acid in the chosen table. Scored against the organism's full
  codon usage, not a `minimum_codon_probability`-filtered subset.
- `correct_coding_translation` — the coding sequence translates to the target protein.
  This should always be `True`.
- `correct_full_translation` — the whole construct, read as a scanning ribosome
  would (first ATG to first in-frame stop), gives the target protein. `False` is a
  real warning that an adapter or padding will change the expressed product — for
  example an ATG in the 5' adapter.
- `quality_check_passed` / `quality_warnings` / `quality_failures` — the final
  acceptance gate. A warning names a requested target that was missed but allowed,
  such as a GC range under `return_best=True`.
- `constraint_violations_remaining` — should be 0 for any returned sequence.

## Output

`build_library` returns a `LibraryResult` (a list of Sequence objects). Each
Sequence contains:

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
