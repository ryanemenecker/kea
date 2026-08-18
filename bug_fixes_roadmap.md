# Kea bug fixes roadmap

Tracking doc for the full-codebase audit (2026-08-17). Every item below was reproduced by running
code, not by inspection alone. Line numbers are from the audit snapshot (commit `5a268c3`) and will
drift as fixes land.

Status key: `[ ]` todo · `[~]` in progress · `[x]` done

---

## Phase 1 — Silently wrong output

These produce a library that looks fine and is wrong. Highest priority: you'd order and synthesize
before noticing.

- [x] **1.1 Frame-aware translation.** `kea/backend/translation_utils.py:18,37`
  `translate_sequence` starts at `sequence.find('ATG')` instead of the sequence's own frame.
  Fixes four downstream bugs at once:
  - `force_start_codon=False` always raises `ValueError: Coding sequence does not translate...`
  - proteins with an internal `*` raise the same misleading error
  - `correct_full_translation` is a false negative whenever an adapter/padding contains ATG
  - `force_stop_codon=False` makes `correct_full_translation` always `False` with any 3' flank
  **Fix:** add an explicit `frame_start`/`scan_for_start` argument. Coding-sequence verification
  translates from index 0; full-sequence verification translates from the known offset
  (`len(adapter_5) + len(padding_5)`), not from a text search.
  **Repro:** `build_library('AKKLVAT','yeast',force_start_codon=False)` → ValueError.

- [x] **1.2 `minimum_codon_probability` is not enforced in the output.**
  `kea/backend/optimize_codon_usage.py:245,367`
  `_adjust_gc_content` and `_improve_sequence` draw replacements from the module-global unfiltered
  `aa_to_codons` instead of `codon_table_obj.aa_to_codons`.
  **Fix:** thread the filtered table through both functions.
  **Repro:** `minimum_codon_probability=0.25` + `target_gc_range` → `CGC` (freq 0.032) appears in
  the returned CDS, 3/3 runs.

- [x] **1.3 Undersized `total_length` is silently ignored.** `kea/backend/sequence.py:127`
  `_get_padding_length` returns negative padding; both `> 0` guards fail; no padding is added and
  the over-length sequence is returned.
  **Fix:** raise a clear `ValueError` naming the required minimum length. Validate in
  `build_library` up front against the *longest* input so the whole library fails fast, not midway.
  **Repro:** `build_library('AKKLVATGHW','yeast',total_length=20)` → `len == 36`, padding `-8/-8`.

- [x] **1.4 GC fine-tuning targets the range midpoint, not the range.**
  `kea/backend/optimize_codon_usage.py:655`
  `_adjust_gc_content(final_sequence, codon_table_obj.target_gc, ...)` pulls an already-in-range
  sequence to `sum(gc_range)/2`, which is why every result has GC of exactly the midpoint.
  **Fix:** skip fine-tuning when GC is already inside the range; when outside, move to the nearest
  range *edge* (plus a small margin), not the midpoint.
  **Note:** the audit's claim that this costs ~25% CAI did **not** reproduce — measured −0.028 to
  +0.029 across 12 range×seed combos, improving in 11 of 12. The design flaw is real; the harm is
  mostly loss of diversity, not loss of quality. Fix it, but don't expect a CAI jump.

- [x] **1.5 `optimization_attempts` has no measurable effect.**
  Measured (5 seeds each): 100 → CAI 0.9290 · 500 → 0.9321 · 5000 (default) → 0.9244 ·
  20000 → 0.9281. GC identical at 0.5000 throughout. A 200x increase buys nothing.
  Two causes: the undocumented `early_stop_threshold=0.95` fires almost immediately, and 1.4 then
  overwrites the result. Related: with early stop unreachable, **98% of 3500 `_improve_sequence`
  calls return the sequence unchanged** — the refinement loop never breaks on no-progress.
  **Fix:** depends on 1.4. Add a no-progress break to the phase-2 loop, document
  `early_stop_threshold`, and re-measure the table above to confirm attempts now matter.

- [x] **1.6 `gc_weight >= 2` deletes codon usage from the objective.**
  `kea/backend/optimize_codon_usage.py:212`
  Score is `(usage*(2-gc_priority) + gc*gc_priority)/2` with `gc_priority` capped at 2.0, so at
  `gc_weight >= 2` the usage term is multiplied by zero and ranking inverts:
  | gc_priority | best-codon seq | worse-codon seq | delta |
  |---|---|---|---|
  | 0 | 1.0000 | 0.5783 | +0.4217 |
  | 1 | 0.8889 | 0.7588 | +0.1301 |
  | 2 | 0.7778 | 0.9394 | **−0.1616** (inverted) |
  | 10 | identical to 2 (silently capped) | | |
  Separately `codon_weight` never reaches `_score_sequence` at all — it only affects initial
  sampling weights.
  **Fix:** score with both weights explicitly (`usage*codon_weight + gc*gc_weight`, normalized),
  drop the magic cap, and validate/document the range.

- [x] **1.7 `_score_sequence` GC term is discontinuous at the tolerance edge.**
  `kea/backend/optimize_codon_usage.py:181` — can invert candidate ranking near the boundary.
  **Fix:** make the penalty continuous and monotonic in distance from the range.

- [x] **1.8 `_improve_sequence` alternative-codon cache is keyed wrong.**
  `kea/backend/optimize_codon_usage.py:367` — the list is built excluding the codon at the *first*
  position seen for that amino acid, then reused everywhere. Proven: L at `TTA` then `CTG` → position
  2 is offered `['TTG','CTT','CTC','CTA','CTG']`, i.e. its own current codon as an "alternative",
  and the genuine alternative `TTA` is never considered.
  **Fix:** cache the full codon list per AA and filter out the current codon at use time.

---

## Phase 2 — Crashes on ordinary input

- [x] **2.1 `optimization_attempts <= 3` → opaque `TypeError`.**
  `kea/backend/optimize_codon_usage.py:527` — `int(n_iter*0.3) == 0`, no candidate generated,
  `None` reaches `_adjust_gc_content`. `1/2/3` crash, `4` works.
  **Fix:** floor the exploration phase at >= 1 iteration and validate `n_iter >= 1`.

- [x] **2.2 Padding generation fails for >= ~800 nt.**
  `kea/backend/library_generation_utils.py:104` — the forward-only repair scan re-creates forbidden
  codons behind itself, so the final recheck rejects every attempt. Caps usable `total_length`.
  Measured: 500 OK, 800 ValueError, 2000 ValueError.
  **Fix:** repair iteratively until clean (or re-scan from `i-2` after each edit) instead of
  rejecting the whole attempt.

- [x] **2.3 `save_library` rejects a bare filename.** `kea/kea.py:338`
  `os.path.dirname("out.csv")` is `''` → the README's own example raises
  `ValueError: Path  is not a directory.`
  **Fix:** treat an empty dirname as the cwd.

- [x] **2.4 CSV is written unquoted.** `kea/kea.py:355`
  A comma in a protein name shifts every column (9 header fields vs 11 data fields).
  **Fix:** use the `csv` module.

- [x] **2.5 `pad_location` is unvalidated.** `kea/backend/sequence.py:128`
  Anything other than `None`/`3`/`5` falls through the if/elif chain, returns `(None, None)`, and
  dies with a cryptic `TypeError`.
  **Fix:** validate in `build_library` with a clear message.

- [x] **2.6 Lowercase adapters bypass validation, then corrupt the construct.** `kea/kea.py:251`
  The check is case-sensitive substring matching; adapters are upper-cased *afterward*.
  **Repro:** `adapter_5_prime="atgcc", verify_start_stop_codons_in_adapters=True` → passes,
  construct starts `ATGCC...`, `correct_full_translation=False`.
  **Fix:** upper-case adapters before validating. Also validate they contain only ACGT.

- [x] **2.7 Custom codon tables are never validated.** `kea/backend/codon_table.py:21`
  A dict missing an amino acid dies with a bare `KeyError` deep in the optimizer, and
  `optimize_codon_usage` validates the protein against the *global* `aa_to_codons` rather than the
  supplied table. Also `minimum_codon_probability` is compared against **pre-normalization** values
  (`codon_table.py:76`), so a counts- or per-1000-scale table either no-ops or rejects everything.
  **Fix:** validate the dict up front (all 20 AAs + `*`, known codons, positive values); normalize
  *before* applying the probability filter.

---

## Phase 3 — Parameters that don't do what they say

- [x] **3.1 `verify_coding_sequence` is a complete no-op.** `kea/backend/sequence.py:236`
  Only read at `sequence.py:116`, inside `if coding_sequence:` — a branch `build_library` never
  takes. `add_coding_sequence` re-verifies unconditionally; `False` only changes which line raises.
  **Fix:** honour the flag in `add_coding_sequence`.

- [x] **3.2 `gc_tolerance` never reaches `CodonTable`.** `kea/kea.py:241`
  Only five positional args are passed, so the refinement phase always scores at the hardcoded
  0.025 while the outer loop uses the user's value. Padding separately hardcodes 0.02.
  **Fix:** forward it to both.

- [x] **3.3 `add_padding`'s `target_length` argument is dead.**
  `kea/backend/library_generation_utils.py:137` — never referenced in the body.
  **Note:** the audit claimed this causes spurious "GC not achievable" crashes; 6 trials of the
  stated scenario produced 0 failures. Dead parameter is real, crash consequence unproven.
  **Fix:** use it in the GC algebra or remove it from the signature.

- [x] **3.4 Single-string input ignores `add_protein_identifiers`.** `kea/kea.py:178`
  `str` → always a random ID; `list` → `Protein_0`; `dict` → keys. Inconsistent.
  **Fix:** honour the flag on the `str` path.

- [x] **3.5 `codon_from_fasta(require_start_codon=False)` is ignored.** `kea/backend/kea_utils.py:100`
  An unconditional ATG pre-filter runs before the flag is consulted; verified both settings return
  an identical table. Also double-filters.
  **Fix:** gate the pre-filter on the flag, drop the duplicate pass.

- [x] **3.6 `target_gc_range=(0,1)` is not the documented no-op.** `kea/kea.py:215`
  The docstring says `(0,1)` means no GC targeting; passing it explicitly targets 50% GC. Only
  `None` disables it — and the `None` branch silently zeroes any user-supplied `gc_weight` and
  `gc_finetuning_iterations`.
  **Fix:** decide one semantics and make docs + code agree; warn instead of silently zeroing.

- [x] **3.7 `return_best=False` default contradicts the docs.** `kea/kea.py:141`
  Documented as "the final result after all iterations"; actually turns a missed `target_gc_range`
  into a hard `ValueError`.
  **Fix:** docs, or the default. Prefer fixing the docs and adding a clearer error message.

- [x] **3.8 `early_stop_threshold` is undocumented.** `kea/kea.py:45`
  In the signature but in neither the docstring nor the README — and it is the reason 1.5 happens.

---

## Phase 4 — Packaging and dead modules

- [x] **4.1 `pyproject.toml` declares zero runtime dependencies.** `pyproject.toml:22`
  The block is commented out. `import kea` pulls in numpy, tqdm, matplotlib **and seaborn**, so the
  README's `pip install git+...` yields an install that fails on import in a clean environment.
  **Fix:** declare numpy and tqdm; make matplotlib optional (see 4.3).

- [x] **4.2 `kea_plotting_utils.py` is dead on arrival.** `kea/backend/kea_plotting_utils.py:7`
  `ImportError: cannot import name 'calc_gc_content' from 'kea.backend.kea_utils'` — that function
  does not exist. Behind it, `graph_gc_content` also sets `xlim(0,100)` for 0–1 fractions and labels
  them "%".
  **Fix:** import the real GC helper; fix the axis/label units.

- [x] **4.3 Unused imports make seaborn a hard requirement.** `kea/backend/kea_utils.py:5-8`
  seaborn, matplotlib, numpy and `defaultdict` are imported and never used.
  **Fix:** delete them. Move matplotlib to a lazy import inside the plotting module.

- [x] **4.4 Data files are not packaged.** `pyproject.toml:63`
  `kea/data/*.txt` and `*.fasta` are missing from `package-data`, so they're absent from a wheel.

- [x] **4.5 `generate_barcode_sequences` is a `pass` stub.**
  `kea/backend/library_generation_utils.py:294` — returns `None`.
  **Fix:** implement or remove.

- [~] **4.6 `restriction_utils.py` is broken throughout.** (Unused by `build_library`; header says
  "NEED TO TEST THIS OUT!!")
  - `:125` offset cut positions measured from the start rather than the 3' end of the recognition
    site — wrong for all 72 offset-notation enzymes, silently drops the 18 negative-offset ones
  - `:50` the leading-paren branch is unreachable → 4 double-cut enzymes get a corrupt site
  - `:111` forward strand only (73/251 enzymes are non-palindromic)
  - `:205` `check_restriction_sites` re-runs the full 251-enzyme scan once per requested enzyme —
    measured **53x slowdown** (150 ms vs 2.8 ms)
  **Fix:** low priority since nothing calls it. Either fix properly or mark clearly experimental.

---

## Phase 5 — Robustness, performance, tests

- [x] **5.1 tqdm bar leaks when the optimizer raises.**
  `kea/backend/optimize_codon_usage.py:668` — no `try/finally`; verified 1 instance left open.

- [x] **5.2 No `seed` argument.** Randomness is split across `numpy` and stdlib `random`, so
  reproducing a library requires seeding both. Undocumented. (Seeding both does work — verified.)
  **Fix:** add a `seed` parameter that seeds both.

- [x] **5.3 `generate_random_protein_ids` returns `list(set(...))`.** `kea/backend/kea_utils.py:199`
  ID→protein assignment permutes per process even when seeded.

- [x] **5.4 `check_nucleotide_percent_similarity` edge cases.**
  `kea/backend/library_generation_utils.py:257` — the `N`→4 encoding is dead code (the validator
  rejects `N` first); zero-length input divides by zero; a single sequence reports itself as its own
  nearest neighbour at similarity 0.

- [x] **5.5 `make_codon_table` `'ACGT'` filter is case-sensitive.** `kea/backend/kea_utils.py:51`
  Soft-masked (lowercase) FASTA silently yields an all-zero table.

- [x] **5.6 Performance (both output-identical).**
  - `_adjust_gc_content` rejoins and rescans the whole sequence per candidate swap; an O(1) GC delta
    update gives byte-identical output (27–92x on that function, ~50% of runtime when a GC range is set)
  - `_score_sequence` recomputes `max(codon_frequency_table[aa].values())` per codon; hoisting it to
    a per-AA dict is ~20% of optimizer runtime

- [x] **5.7 Tests.** The only existing test asserts that `import kea` worked.
  Add regression tests for every item above, especially the Phase 1 invariants:
  translation round-trip, exact `total_length`, GC within range, `minimum_codon_probability`
  respected in output, and determinism under a seed.

---

## Verification harness

Measured before and after. Reproduce with `python -m pytest kea/tests/test_kea.py`.

| check | before | after |
|---|---|---|
| `optimization_attempts` 100 → 20000, mean CAI | 0.9290 → 0.9281 (no effect) | 0.9852 → 0.9940 (and 0.92 → 0.99 overall) |
| GC of result with `target_gc_range=(0.40,0.60)` | always exactly 0.5000 | spread across the range, 8/8 in range |
| sub-threshold codons at `minimum_codon_probability=0.25` | 2 per sequence, 3/3 runs | 0 across 15 runs |
| `build_library(..., total_length=20)` on a 36 nt CDS | returns 36 nt silently | `ValueError: ... needs at least 36 nt` |
| max working padding length | ~500 nt | 5000 nt in 2 ms |
| `force_start_codon=False` | always raises | works |
| `save_library(results, "out.csv")` (README example) | `ValueError` | works |
| GC score at 0.4233 vs 0.4267 (closer to range) | 0.9922 vs 0.0667 (inverted) | monotonic, 0 inversions over the full sweep |
| `_adjust_gc_content` on an 864 nt sequence | 218 ms | 8.1 ms (27x, byte-identical output) |
| `check_restriction_sites` with all 251 enzymes | 53x slower than one full scan | one scan |
| test suite | 1 test (asserts `import` worked) | 102 tests, ~1 s |

## Still open

- **4.6 restriction_utils cut coordinates.** Detection, parsing and the O(n^2) scan are
  fixed, and the four leading-paren enzymes (BaeI, BcgI, BsaXI, CspCI) now parse. Two known
  limitations remain, documented at the top of the module: forward-strand-only search
  (73/251 enzymes are non-palindromic), and offset cut positions measured from the start of
  the match rather than the 3' end of the recognition site. Fixing those properly needs
  NEB-derived test vectors rather than guesswork. Nothing in `build_library` uses this module.
- **4.5 `generate_barcode_sequences`** now raises `NotImplementedError` instead of silently
  returning `None`. Implementing it is a feature, not a bug fix.

## Behaviour changes worth knowing about

- `translate_sequence()` now translates from index 0 by default instead of hunting for the
  first ATG, and no longer truncates at the first stop. Use `translate_open_reading_frame()`
  for the old scanning-ribosome behaviour (that is what the whole-construct check uses).
- `build_library(early_stop_threshold=...)` now defaults to `None` (disabled). The old 0.95
  default fired almost immediately and was the reason `optimization_attempts` did nothing.
- `add_padding()` lost its unused `target_length` parameter.
- `_score_sequence()` takes `gc_weight`/`codon_weight` instead of a capped `gc_priority`.
- Runs are slower than before at the same `optimization_attempts`, because the refinement
  phase now actually runs instead of stopping after one iteration. CAI went from ~0.92 to
  ~0.99 as a result.
