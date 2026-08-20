# Kea roadmap

Kea's immediate goal is reliable protein-to-DNA library design with explicit,
testable sequence constraints. It should report infeasible designs rather than
silently returning a construct that violates a requested constraint.

## Implemented

### Human codon tables

- Human highly expressed-gene codon frequencies are available as `"human"`.
- Unweighted human RefSeq codon frequencies are available as
  `"human_unweighted"`.
- These tables describe codon frequencies, not universal decoding rates. Future
  cell-type-specific models should remain distinct from these organism-level
  tables.

### Human splice-site avoidance

- Exact Yeo-Burge MaxEntScan scoring for 9-nt donor and 23-nt acceptor windows.
- Configurable donor and acceptor cutoffs.
- Synonymous single- and double-codon repair with translation and GC-range
  preservation.
- Optional fixed 5' and 3' transcribed context so UTR/CDS and vector/CDS
  junctions are scanned.
- Failure-loud behavior when a site is confined to fixed context or no
  synonymous repair exists.

The current cutoffs are provisional design thresholds, not experimental
guarantees. Kea reports that no *predicted* site exceeds the configured cutoff;
it cannot prove that a transcript will never be spliced.

### Premature polyadenylation avoidance

- Detection and synonymous repair of the canonical `AATAAA`, major `ATTAAA`,
  and 11 common human PAS variants.
- Annotation of downstream U-rich elements and the established GU-rich
  pentamers `GTTGT`, `TGTGT`, and `GTGTT`.
- Conservative mode removes every supported PAS hexamer. An optional
  context-required mode only removes PAS motifs with a downstream U/GU-rich
  element.

### Sequence-complexity controls

- Configurable maximum homopolymer run length.
- Configurable maximum copies of short tandem-repeat units.
- Configurable maximum exact direct-repeat length, including separated and
  overlapping repeats.
- All three use synonymous repair and are disabled unless requested.

### Constraint repair

- Repair runs after codon/GC optimization and is retried from a fresh
  optimization when it fails (`constraint_retry_attempts`, default 2). Codon
  optimization is stochastic and repair is a greedy local search, so a different
  starting point frequently succeeds where the first did not.
- A fix that would push GC outside `target_gc_range` is paired with a
  compensating synonymous edit elsewhere rather than being vetoed. Measured on
  120 x 250 aa proteins at GC (0.30, 0.42) with five constraints enabled:
  85/120 -> 114/120 built, mean CAI 0.8590 -> 0.8578, every result still inside
  the requested range, and 8.6 s -> 6.7 s because fewer sequences exhaust their
  retries.
- Violations that no synonymous encoding can clear are proven by exhaustive
  enumeration of the codons covering the offending window (bounded; wider
  windows report "not proven" rather than guessing). Those fail immediately with
  a message naming the site and score instead of consuming every retry.
- Generated padding is constraint-checked as well as the coding sequence,
  including sites created across the padding/CDS junction. Cloning adapters are
  still excluded by design; pass them as transcribed context to have them
  scanned.

### Synonymous variant sets

- `sequences_per_protein` requests multiple coding sequences for every input
  protein, with stable `_variant_N` names.
- `minimum_hamming_distance` is a hard all-pairs nucleotide-distance requirement
  over coding sequences. `hamming_distance_attempts` bounds the randomized search
  for each additional encoding.
- Diversification starts from the codon-optimized sequence, prioritizes the hard
  distance target while minimizing codon-frequency cost, preserves an already
  satisfied GC range, and then reuses the normal constraint-repair and final-QC
  path.
- Infeasible sets retain accepted variants and report every unfilled slot as a
  `sequence_diversity` failure, or raise `SequenceDiversityError` under
  `on_error="raise"`.

### Final quality check

- Every sequence passes one acceptance gate before it is returned, measuring
  translation, total length, coding and full-sequence GC, remaining constraint
  violations, and codon adaptation together on the finished construct. Each
  stage is free to trade the other targets off temporarily; this is where the
  result is checked as a whole.
- Required checks that fail make the sequence a per-sequence failure. A target
  that was missed but explicitly permitted (best-effort GC under
  `return_best=True`) is recorded as a warning on `sequence.quality_report` and
  surfaced via `LibraryResult.with_warnings`, rather than passing silently as it
  did before.
- Codon adaptation is reported on every sequence and can be enforced with
  `minimum_codon_adaptation`. Measured cost of constraint repair itself is small
  (mean CAI drop 0.0015-0.0018, worst case 0.0095); the large swings come from
  the codon/GC trade-off in the optimizer, not from repair.

### GC centering (evaluated, off by default)

- `gc_centering` tilts the GC term toward the middle of `target_gc_range` instead
  of leaving it flat, so codon usage cannot pin GC against an edge. Implemented
  and tested, but **default 0.0** because measurement did not support enabling it.
- Paired per-protein measurement, 250 aa human sequences, five constraints:

    target_gc_range   gc_centering=0.05
    (0.30, 0.42)      CAI +0.038, ~2 fewer built per 60
    (0.20, 0.35)      CAI +0.028, ~10 fewer built per 35
    (0.62, 0.75)      CAI -0.001 (no effect)
    (0.55, 0.62)      CAI -0.007
    (0.49, 0.51)      CAI -0.014
    (0.45, 0.60)      CAI -0.020

- On the WIDE ranges first tested it buys no additional built sequences -- GC
  compensation in repair already covers that -- and only affects codon adaptation.
  On NARROW ranges below the preferred GC it does improve yield, but only at
  >= 0.10: range (0.28, 0.36) builds 75.0/80 at 0.0, 77.7/80 at 0.10.
- The response is NON-MONOTONIC. Values around 0.03-0.05 are a hazard band where
  the tilt pulls GC off the edge but cannot settle it at the centre: range
  (0.28, 0.36) drops to 55.0/80 at 0.03 and 62.7/80 at 0.05 before recovering to
  77.7/80 at 0.10. Use 0 or >= 0.10, never in between. The originally proposed
  default of 0.05 sat in the worst part of that band.
- The in-range gradient is `gc_centering / half_width`, so the same value behaves
  very differently across range widths: gentle on wide ranges, saturating on
  narrow ones. Worth normalising to absolute GC units if this is revisited.

### Optimizer: boundary local optimum (fixed)

- Diagnosed by exhaustive move census on converged sequences: **0** improving
  single-codon moves remained, but **19,128** GC-neutral pair moves across 6
  sequences were still available. The greedy one-codon-at-a-time search cannot see
  them, because at the GC-range edge every usage-improving single swap leaves the
  range.
- `_improve_sequence` now makes paired swaps: when a swap would leave the range it
  looks up a partner swap that pays the GC back, using an index of swaps grouped
  by GC-count delta and ordered by usage gain. Census after the fix: 0 single and
  **0** pair moves remaining -- a true 2-opt local optimum.
- With no constraints this is a pure win: codon adaptation 0.859 -> 0.910 on
  (0.30, 0.42), 50/50 built, 100% in range, ~1.7x runtime.
- With constraints it costs yield, because maximal codon optimization is more
  repetitive: raw optimizer output goes from 2.3 to 7.7 violations per sequence,
  and (0.28, 0.36) drops from 50.7 to 19.0 built of 60. So `escape_gc_boundary`
  defaults to on only when no constraints are requested; `True` forces it.
- This supersedes `gc_centering` as the answer to the boundary problem: it gains
  more codon adaptation (+0.051 vs +0.041), keeps GC where the caller asked for it
  rather than dragging it to the range centre, and has no non-monotonic hazard band.

### Remaining: constraint-aware optimization

- The yield cost above is the real remaining gap. The optimizer is constraint-blind
  by design (adding an O(n) constraint scan to the O(1) inner loop would undo the
  75x speedup), so it happily produces repetitive sequences that repair then cannot
  fix. A cheap local repetition penalty in the objective -- homopolymer run length
  and short-repeat content are both computable incrementally -- would likely let
  the boundary escape stay on with constraints enabled.

### Library reporting

- `save_library` writes 40 columns per sequence: the construct, GC, codon
  adaptation, translation checks, the final quality-gate outcome, the strongest
  predicted splice donor and acceptor with MaxEntScan score/position/window,
  polyadenylation hexamer counts (total, canonical AATAAA, and how many have a
  downstream U/GU-rich element), and complexity measures. `include_diagnostics=False`
  gives the short form. `save_failures` writes what did not build and why.
- Screens run on the complete construct and are reported whether or not the
  matching constraint was enabled, so a library can be screened after the fact.
- The README documents how to read the numbers, including that MaxEntScan scores
  are log2 odds ratios rather than probabilities, with a measured null distribution
  over 107,400 donor and 105,720 acceptor windows of ordinary human-codon CDS
  (donor median -19.2, p99 1.9, p99.9 7.9; consensus CAGGTAAGT scores 10.86).

### Failure handling

- `build_library` records per-sequence failures and continues (`on_error`,
  default `"skip"`), returning a `LibraryResult` that is still a list of
  Sequence objects but also carries `.failures`, `.n_failed` and `.summary()`.
  Configuration errors always raise, since they apply to every input.

## Near-term priorities

### 1. Calibrate human splice thresholds

- Score annotated human splice junctions, ordinary human CDS windows, and known
  cryptic reporter sites.
- Publish sensitivity/false-positive tradeoffs for recommended donor and
  acceptor cutoffs.
- Benchmark repaired versus unrepaired constructs with junction-spanning
  RT-PCR or targeted RNA sequencing in HEK293 cells and the intended expression
  cell type.
- Add intended-intron annotations so authentic splice junctions can be
  preserved while de novo sites are removed.
- Evaluate optional SpliceAI or Pangolin validation backends. Keep deep models
  optional because of model size, licensing, and synthetic-transcript domain
  shift.

### 2. Improve premature-polyadenylation prediction

- Calibrate the current PAS/DSE heuristic on experimentally mapped human
  cleavage sites and matched decoys.
- Model the expected cleavage interval between the PAS and downstream element,
  rather than treating downstream context as a binary annotation.
- Expose low-, medium-, and high-risk PAS categories without weakening the
  hard-constraint default.

### 3. Add synthesis-vendor complexity profiles

- Named profiles for common limits on homopolymers, tandem repeats, direct
  repeats, local GC extremes, and palindromes.
- Return coordinate-level explanations suitable for vendor troubleshooting.
- Measure the effect of each profile on feasibility and codon adaptation.

### 4. Translation-initiation context

- Kozak-sequence controls.
- Detection of upstream AUGs and unintended upstream open reading frames when a
  5' UTR is provided.
- Local RNA-structure scoring around the 5' UTR/start codon. Avoid treating one
  global minimum-free-energy number as a universal expression objective.

### 5. Assembly and motif constraints

- Integrate required and forbidden restriction sites with the same hard
  constraint interface.
- Detect sites created across adapters, padding, UTRs, and coding-sequence
  junctions.
- General user-supplied forbidden motifs and protected sequence intervals.


### 6. Capacity to check extant libraries for cryptic splice sites and premature polyadenylation
- Enables users to be able to go and check an already-made library for any problems
- Let's us compare to other codon optimization approaches

### 7. GC bias titration across sequence
- Start sequence with lower GC content then go higher after ~10 codons. Still target average.

## Context-dependent features

These should be explicit application profiles rather than universal defaults:

- CpG and UpA modulation for DNA vaccines, viral vectors, or synthetic mRNA.
- Cell-type- and state-specific codon optimality based on tRNA availability and
  measured mRNA stability.
- Codon-pair and codon-context objectives, added only after controlled A/B
  evidence that they improve the intended expression system.
- miRNA, RBP, AU-rich, internal ribosome-entry, and transcription-factor motif
  screens.
- Stop-codon and downstream termination context.
- Viral packaging signals, recombination-prone elements, and vector-specific
  constraints.

## Design principles

1. Translation fidelity and explicitly requested motif exclusions are hard
   constraints; codon preference and other expression proxies are soft
   objectives.
2. Constraints are evaluated on the transcribed sequence, including supplied
   vector/UTR context, not automatically on cloning adapters that may never be
   transcribed.
3. Every repair is synonymous, globally rechecked, and required to preserve an
   already satisfied GC range.
4. Defaults remain off until calibrated. Enabling a constraint must either
   satisfy it or return an informative error with remaining sites.
5. Computational predictions require measured smoke tests before production
   library synthesis.
