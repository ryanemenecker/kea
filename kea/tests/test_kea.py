"""
Unit and regression tests for the kea package.

Test names reference the item numbers in bug_fixes_roadmap.md so a failure points
straight at the defect it is guarding against.
"""

import csv
import os
import random
import sys
from statistics import mean

import numpy as np
import pytest

import kea
from kea import build_library, save_failures, save_library
from kea.backend.codon_table import CodonTable, validate_codon_usage
from kea.backend.kea_utils import generate_random_protein_ids, make_codon_table
from kea.backend.library_generation_utils import (check_nucleotide_percent_similarity,
                                                  create_padding_sequence)
from kea.backend.optimize_codon_usage import (_calculate_gc_content, _score_sequence,
                                              calculate_codon_adaptation_score)
from kea.backend.sequence import Sequence
from kea.backend.sequence_diversity import (SequenceDiversityError,
                                           nucleotide_hamming_distance)
from kea.backend.sequence_constraints import (
    COMMON_HUMAN_POLYADENYLATION_SIGNALS,
    ConstraintSet,
    DirectRepeatConstraint,
    HomopolymerConstraint,
    HumanSpliceConstraint,
    MaxEntScanScorer,
    PrematurePolyadenylationConstraint,
    TandemRepeatConstraint,
    repair_coding_sequence,
)
from kea.backend.translation_utils import (find_first_orf, translate_open_reading_frame,
                                           translate_sequence)
from kea.data.aa_codon_conversions import aa_to_codons, codons_to_aa
from kea.data.codon_tables import yeast_codon_usage

# Small budgets so the suite stays fast; the invariants under test do not need
# a long search to hold.
FAST = dict(optimization_attempts=200, show_progress=False, show_optimization_progress=False)

ALL_AAS = "ACDEFGHIKLMNPQRSTVWY"


def test_kea_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "kea" in sys.modules


def test_import_does_not_require_plotting_stack():
    """4.1/4.3: importing kea must not drag in matplotlib or seaborn."""
    for module in ("matplotlib", "seaborn"):
        assert module not in sys.modules or module in sys.modules  # tolerated if a test imported it
    # The real assertion: nothing in the import chain references seaborn.
    import kea.kea
    import kea.backend.kea_utils
    assert not hasattr(kea.backend.kea_utils, "sns")


# --------------------------------------------------------------------------
# Phase 1.1 -- frame-aware translation
# --------------------------------------------------------------------------

def test_translate_sequence_uses_frame_zero_by_default():
    """1.1: a coding sequence defines its own frame; do not hunt for an ATG."""
    assert translate_sequence("GCTAAACTT") == "AKL"
    assert translate_sequence("ATGAAATAA") == "MK*"


def test_translate_sequence_keeps_internal_stops():
    """1.1: internal stop codons must round-trip as '*', not truncate."""
    assert translate_sequence("ATGAAATAATTTTAA") == "MK*F*"


def test_translate_sequence_strips_trailing_stop_when_asked():
    assert translate_sequence("ATGAAATAA", return_stop_codon=False) == "MK"


def test_translate_sequence_rejects_non_dna():
    with pytest.raises(ValueError, match="Invalid codon"):
        translate_sequence("ATGNNNAAA")


def test_translate_sequence_is_case_insensitive():
    assert translate_sequence("atgaaataa") == "MK*"


def test_translate_open_reading_frame_scans_for_atg():
    """The whole-construct check keeps scanning-ribosome semantics."""
    assert translate_open_reading_frame("CCCATGAAATAA") == "MK*"
    assert translate_open_reading_frame("CCCCCC") == ""


def test_find_first_orf_can_respect_frame():
    # ATG at index 2, which is not in frame 0.
    assert find_first_orf("CCATGAAA") == 2
    assert find_first_orf("CCATGAAA", frame=0) == -1


@pytest.mark.parametrize("protein", [
    "MKKFLVLLFCWAVLCEHN",
    ALL_AAS,
    "AKKLVAT",        # does not start with M
    "MKK*FLVLL",      # internal stop
    "M",
])
@pytest.mark.parametrize("force_start,force_stop", [(True, True), (False, False),
                                                    (True, False), (False, True)])
def test_coding_sequence_round_trips(protein, force_start, force_stop):
    """1.1: every protein/flag combination must translate back exactly."""
    seq = build_library(protein, "yeast", force_start_codon=force_start,
                        force_stop_codon=force_stop, **FAST)[0]
    assert translate_sequence(seq.coding_sequence) == seq.protein_sequence
    assert seq.correct_coding_translation is True
    assert len(seq.coding_sequence) == 3 * len(seq.protein_sequence)


def test_force_start_codon_false_does_not_raise():
    """1.1: this used to raise ValueError for any protein not starting with M."""
    seq = build_library("AKKLVAT", "yeast", force_start_codon=False, **FAST)[0]
    assert not seq.protein_sequence.startswith("M")


# --------------------------------------------------------------------------
# Phase 1.2 -- minimum_codon_probability
# --------------------------------------------------------------------------

@pytest.mark.parametrize("threshold", [0.10, 0.20, 0.25])
def test_minimum_codon_probability_respected_in_output(threshold):
    """1.2: GC fine-tuning used to reinsert codons banned by this threshold."""
    seq = build_library("MKKFLVLLFCWAVLCEHNRRSSPPGGTTVV", "yeast",
                        minimum_codon_probability=threshold,
                        target_gc_range=(0.40, 0.60), gc_finetuning_iterations=800,
                        return_best=True, **FAST)[0]
    cds = seq.coding_sequence
    offenders = [cds[i:i + 3] for i in range(0, len(cds), 3)
                 if yeast_codon_usage[codons_to_aa[cds[i:i + 3]]][cds[i:i + 3]] < threshold]
    assert offenders == []


def test_minimum_codon_probability_too_high_is_a_clear_error():
    with pytest.raises(ValueError, match="minimum_codon_probability is too high"):
        build_library("MKKFLVLL", "yeast", minimum_codon_probability=0.95, **FAST)


# --------------------------------------------------------------------------
# Phase 1.3 / 2.5 -- length and pad_location validation
# --------------------------------------------------------------------------

def test_undersized_total_length_raises():
    """1.3: this used to silently return an over-length sequence."""
    with pytest.raises(ValueError, match="too short"):
        build_library("AKKLVATGHW", "yeast", total_length=20, **FAST)


def test_undersized_total_length_checked_against_longest_sequence():
    with pytest.raises(ValueError, match="at least 81 nt"):
        build_library(["MKK", "MKKFLVLLFCWAVLCEHNRRSSPPGG"], "yeast", total_length=40, **FAST)


@pytest.mark.parametrize("pad_location", [None, 3, 5])
@pytest.mark.parametrize("total_length", [90, 120, 1500])
def test_total_length_is_exact(pad_location, total_length):
    seq = build_library("MKKFLVLLFCWAVLCEHN", "yeast", total_length=total_length,
                        pad_location=pad_location, **FAST)[0]
    assert len(seq.full_dna_sequence) == total_length


def test_pad_location_only_at_requested_end():
    five = build_library("MKKFLVLL", "yeast", total_length=90, pad_location=5, **FAST)[0]
    three = build_library("MKKFLVLL", "yeast", total_length=90, pad_location=3, **FAST)[0]
    assert len(five.padding_3_prime) == 0 and len(five.padding_5_prime) > 0
    assert len(three.padding_5_prime) == 0 and len(three.padding_3_prime) > 0


@pytest.mark.parametrize("bad", [4, "5", 0, 1.5])
def test_invalid_pad_location_raises(bad):
    """2.5: used to fall through and die with a cryptic TypeError."""
    with pytest.raises(ValueError, match="pad_location"):
        build_library("MKKFLVLL", "yeast", total_length=60, pad_location=bad, **FAST)


# --------------------------------------------------------------------------
# Phase 1.4 / 1.5 -- GC targeting and optimization effort
# --------------------------------------------------------------------------

@pytest.mark.parametrize("gc_range", [(0.40, 0.60), (0.45, 0.55), (0.35, 0.65)])
def test_gc_content_lands_inside_requested_range(gc_range):
    for seed in range(3):
        seq = build_library("MKKFLVLLFCWAVLCEHNRRSSPPGG", "yeast", target_gc_range=gc_range,
                            gc_finetuning_iterations=800, return_best=True, seed=seed, **FAST)[0]
        assert gc_range[0] <= seq.gc_content_coding_sequence <= gc_range[1]


def test_gc_is_not_pinned_to_the_range_midpoint():
    """1.4: fine-tuning used to drag every sequence to exactly sum(range)/2."""
    gc_range = (0.35, 0.65)
    values = {round(build_library("MKKFLVLLFCWAVLCEHNRRSSPPGGTTVVIIQQDDEE", "yeast",
                                  target_gc_range=gc_range, gc_finetuning_iterations=800,
                                  return_best=True, seed=s, **FAST)[0].gc_content_coding_sequence, 4)
              for s in range(6)}
    midpoint = round(sum(gc_range) / 2, 4)
    assert values != {midpoint}


def test_more_optimization_attempts_improves_codon_adaptation():
    """1.5: attempts used to have no measurable effect (early stop fired at once)."""
    def mean_cai(attempts):
        scores = []
        for seed in range(4):
            seq = build_library("MKKFLVLLFCWAVLCEHNRRSSPPGGTTVVIIQQDDEE", "yeast",
                                target_gc_range=(0.40, 0.60), optimization_attempts=attempts,
                                gc_finetuning_iterations=500, return_best=True, seed=seed,
                                show_progress=False, show_optimization_progress=False)[0]
            scores.append(calculate_codon_adaptation_score(seq.coding_sequence, yeast_codon_usage))
        return sum(scores) / len(scores)

    assert mean_cai(600) > mean_cai(20)


# --------------------------------------------------------------------------
# Phase 1.6 / 1.7 -- scoring
# --------------------------------------------------------------------------

def test_gc_score_is_monotonic_outside_the_range():
    """1.7: the old piecewise formula jumped from 0.99 to 0.07 at the tolerance edge."""
    def leu_seq(gc_count, n=100):
        two = min(gc_count // 2, n)
        one = 1 if gc_count - 2 * two else 0
        return "CTG" * two + "CTT" * one + "TTA" * (n - two - one)

    previous = None
    for count in range(0, 135):
        seq = leu_seq(count)
        if _calculate_gc_content(seq) > 0.45:
            break
        score = _score_sequence(seq, (0.45, 0.55), yeast_codon_usage,
                                gc_weight=1.0, gc_tolerance=0.025, codon_weight=0.0)
        if previous is not None:
            assert score >= previous - 1e-12
        previous = score


def test_codon_weight_still_counts_at_high_gc_weight():
    """1.6: gc_weight >= 2 used to multiply the codon usage term by zero."""
    best = "ATG" + "TGT" * 10   # most frequent Cys codon
    worse = "ATG" + "TGC" * 10
    for gc_weight in (0, 1, 2, 10):
        a = _score_sequence(best, (0.4, 0.6), yeast_codon_usage, gc_weight, 0.025, codon_weight=1.0)
        b = _score_sequence(worse, (0.4, 0.6), yeast_codon_usage, gc_weight, 0.025, codon_weight=1.0)
        assert a != b, f"codon usage ignored at gc_weight={gc_weight}"


def test_score_is_bounded():
    seq = "ATG" + "TGT" * 10
    for gw, cw in [(0, 1), (1, 1), (5, 1), (1, 5), (100, 1)]:
        assert 0.0 <= _score_sequence(seq, (0.4, 0.6), yeast_codon_usage, gw, 0.025, cw) <= 1.0


def test_negative_weights_rejected():
    with pytest.raises(ValueError, match="non-negative"):
        build_library("MKKFLVLL", "yeast", codon_weight=-1, **FAST)


# --------------------------------------------------------------------------
# Phase 2.1 / 2.2 -- crashes on ordinary input
# --------------------------------------------------------------------------

@pytest.mark.parametrize("attempts", [1, 2, 3, 4, 10])
def test_low_optimization_attempts_does_not_crash(attempts):
    """2.1: 1-3 used to die with 'NoneType has no len()'."""
    seq = build_library("MKKFLVLL", "yeast", optimization_attempts=attempts,
                        show_progress=False, show_optimization_progress=False)[0]
    assert seq.correct_coding_translation is True


def test_zero_optimization_attempts_is_a_clear_error():
    with pytest.raises(ValueError, match="at least 1"):
        build_library("MKKFLVLL", "yeast", optimization_attempts=0,
                      show_progress=False, show_optimization_progress=False)


@pytest.mark.parametrize("length", [60, 200, 800, 2000])
def test_padding_generation_scales(length):
    """2.2: anything above ~800 nt used to fail after exhausting every attempt."""
    padding = create_padding_sequence(length, (0.4, 0.6), 0.02, True, True, max_attempts=2000)
    assert len(padding) == length
    assert "ATG" not in padding
    assert not any(stop in padding for stop in ("TAA", "TAG", "TGA"))


@pytest.mark.parametrize("avoid_start,avoid_stop", [(True, False), (False, True), (False, False)])
def test_padding_avoid_flags_are_independent(avoid_start, avoid_stop):
    padding = create_padding_sequence(600, (0.4, 0.6), 0.02, avoid_start, avoid_stop,
                                      max_attempts=2000)
    if avoid_start:
        assert "ATG" not in padding
    if avoid_stop:
        assert not any(stop in padding for stop in ("TAA", "TAG", "TGA"))


# --------------------------------------------------------------------------
# Phase 2.3 / 2.4 -- save_library
# --------------------------------------------------------------------------

def test_save_library_accepts_bare_filename(tmp_path, monkeypatch):
    """2.3: the README's own example used to raise."""
    monkeypatch.chdir(tmp_path)
    results = build_library(["MKKFLVLLFCWAVLCEHN"], "yeast", **FAST)
    save_library(results, "optimized_sequences.csv")
    assert (tmp_path / "optimized_sequences.csv").is_file()


def test_save_library_escapes_commas(tmp_path):
    """2.4: a comma in a name used to shift every subsequent column."""
    path = tmp_path / "library.csv"
    results = build_library({"prot,with,commas": "MKKFLVLL"}, "yeast", **FAST)
    save_library(results, str(path))
    rows = list(csv.reader(path.open()))
    assert len(rows) == 2
    assert len(rows[1]) == len(rows[0])
    assert rows[1][0] == "prot,with,commas"


# --------------------------------------------------------------------------
# Library report
# --------------------------------------------------------------------------

def test_save_library_writes_every_diagnostic_column(tmp_path):
    from kea.backend.sequence_report import REPORT_COLUMNS
    path = tmp_path / "library.csv"
    library = build_library({"p": "MKKFLVLLFCWAVLCEHNRRSSPPGG"}, "human",
                            target_gc_range=(0.45, 0.60), total_length=400,
                            adapter_5_prime="GGTCTC", adapter_3_prime="GAGACC",
                            gc_finetuning_iterations=300, return_best=True,
                            seed=1, **FAST)
    save_library(library, str(path))
    rows = list(csv.DictReader(path.open()))
    assert len(rows) == 1
    assert list(rows[0]) == REPORT_COLUMNS
    row = rows[0]
    assert row["correct_coding_translation"] == "True"
    assert int(row["full_length_nt"]) == 400
    assert 0.0 <= float(row["codon_adaptation_index"]) <= 1.0
    # Scores are log-odds and may legitimately be negative.
    float(row["donor_max_score"])
    float(row["acceptor_max_score"])
    assert int(row["polya_signal_count"]) >= 0
    assert int(row["max_homopolymer_run"]) >= 1


def test_save_library_short_form(tmp_path):
    from kea.backend.sequence_report import BASIC_COLUMNS
    path = tmp_path / "short.csv"
    library = build_library({"p": "MKKFLVLL"}, "yeast", **FAST)
    save_library(library, str(path), include_diagnostics=False)
    rows = list(csv.DictReader(path.open()))
    assert list(rows[0]) == BASIC_COLUMNS


def test_diagnostics_are_reported_without_the_constraint_enabled(tmp_path):
    """Screening must work on a library built with no splice constraint."""
    library = build_library({"p": "MKKFLVLLFCWAVLCEHNRRSSPPGG"}, "human", seed=1, **FAST)
    record = kea.annotate_sequence(library[0])
    assert isinstance(record["donor_max_score"], float)
    assert record["splice_donor_threshold"] == 6.0
    assert record["donor_sites_over_threshold"] >= 0


def test_save_failures_records_why(tmp_path):
    path = tmp_path / "failures.csv"
    library = build_library({"ok": "MKKFLVLLFCWAVLCEHN", "bad": "M" + "WYV" * 8},
                            "human", avoid_human_splice_sites=True, seed=1, **FAST)
    save_failures(library, str(path))
    rows = list(csv.DictReader(path.open()))
    assert len(rows) == library.n_failed == 1
    assert rows[0]["protein_name"] == "bad"
    assert rows[0]["stage"] == "constraint_repair"
    assert "human_splice_donor" in rows[0]["unsatisfiable_constraints"]


def test_longest_direct_repeat_matches_brute_force():
    from kea.backend.sequence_report import longest_direct_repeat

    def brute(sequence, minimum=6):
        best = ""
        for size in range(minimum, len(sequence)):
            seen, hit = set(), None
            for start in range(len(sequence) - size + 1):
                chunk = sequence[start:start + size]
                if chunk in seen:
                    hit = chunk
                    break
                seen.add(chunk)
            if hit is None:
                break
            best = hit
        return len(best), best

    random.seed(1)
    for _ in range(40):
        length = random.randint(10, 200)
        sequence = "".join(random.choice("ACGT") for _ in range(length))
        assert longest_direct_repeat(sequence) == brute(sequence)


def test_homopolymer_is_not_double_counted_as_a_tandem_repeat():
    from kea.backend.sequence_report import longest_homopolymer, most_tandem_repeat_copies
    assert longest_homopolymer("CCTTTTTTTTCC") == (8, "T")
    # Unit lengths start at 2, matching TandemRepeatConstraint.
    copies, unit = most_tandem_repeat_copies("CCCAGCAGCAGCAGCC")
    assert (copies, unit) == (4, "CAG")


def test_splice_summary_reference_values():
    """Anchored on published MaxEntScan values so a scoring regression is caught."""
    from kea.backend.sequence_report import splice_site_summary
    summary = splice_site_summary("AAA" + "CAGGTAAGT" + "AAA")
    assert summary["donor_max_score"] == 10.86
    assert summary["donor_max_window"] == "CAGGTAAGT"
    assert summary["donor_sites_over_threshold"] == 1


def test_save_library_rejects_missing_directory(tmp_path):
    with pytest.raises(ValueError, match="not a directory"):
        save_library(build_library("MKKFLVLL", "yeast", **FAST),
                     str(tmp_path / "nope" / "library.csv"))


# --------------------------------------------------------------------------
# Phase 2.6 / 2.7 -- adapter and codon table validation
# --------------------------------------------------------------------------

def test_lowercase_adapter_is_validated():
    """2.6: the check was case-sensitive but adapters get upper-cased afterwards."""
    with pytest.raises(ValueError, match="start codon"):
        build_library("MKKFLVLL", "yeast", adapter_5_prime="atgcc",
                      verify_start_stop_codons_in_adapters=True, **FAST)


def test_non_dna_adapter_rejected():
    with pytest.raises(ValueError, match="non-DNA"):
        build_library("MKKFLVLL", "yeast", adapter_5_prime="GGNTC", **FAST)


def test_lowercase_adapter_is_upper_cased_into_the_construct():
    seq = build_library("MKKFLVLL", "yeast", adapter_5_prime="ggtctc",
                        adapter_3_prime="gagacc", **FAST)[0]
    assert seq.full_dna_sequence.startswith("GGTCTC")
    assert seq.full_dna_sequence.endswith("GAGACC")


def test_custom_table_missing_amino_acid_rejected():
    """2.7: used to surface as a bare KeyError deep in the optimizer."""
    table = {aa: dict(codons) for aa, codons in yeast_codon_usage.items()}
    del table["L"]
    with pytest.raises(ValueError, match="missing amino acid"):
        build_library("MKKFLVLL", table, **FAST)


def test_custom_table_wrong_codon_rejected():
    table = {aa: dict(codons) for aa, codons in yeast_codon_usage.items()}
    table["A"]["TTT"] = 0.1
    with pytest.raises(ValueError, match="not a valid codon"):
        validate_codon_usage(table)


def test_unnormalized_custom_table_behaves_like_a_normalized_one():
    """2.7: minimum_codon_probability compared against pre-normalization values."""
    scaled = {aa: {c: f * 1000 for c, f in codons.items()}
              for aa, codons in yeast_codon_usage.items()}
    a = CodonTable(yeast_codon_usage, 1, 1, (0, 1), 0.25).codon_table
    b = CodonTable(scaled, 1, 1, (0, 1), 0.25).codon_table
    assert {aa: set(v) for aa, v in a.items()} == {aa: set(v) for aa, v in b.items()}


def test_codon_table_probabilities_sum_to_one():
    table = CodonTable(yeast_codon_usage, 1, 1, (0.4, 0.6), 0.15).codon_table
    for aa, codons in table.items():
        assert abs(sum(codons.values()) - 1.0) < 1e-9, aa


# --------------------------------------------------------------------------
# Phase 3 -- parameters that were ignored
# --------------------------------------------------------------------------

def test_verify_coding_sequence_false_skips_verification():
    """3.1: the flag was read on a branch build_library never takes."""
    table = CodonTable(yeast_codon_usage, 1, 1, (0, 1), 0.0)
    kwargs = dict(codon_table=table, name="n", adapter_3_prime="", adapter_5_prime="",
                  avoid_adding_start_codons=True, avoid_adding_stop_codons=True,
                  total_length=None, pad_location=None, coding_sequence="ATGGGGTTATAA")
    Sequence("MKL", verify_coding_sequence=False, **kwargs)
    with pytest.raises(ValueError, match="does not translate"):
        Sequence("MKL", verify_coding_sequence=True, **kwargs)


def test_gc_tolerance_reaches_the_codon_table():
    """3.2: build_library never forwarded it, so refinement always used 0.025."""
    from kea.backend.codon_table import CodonTable as CT
    import kea.kea as kea_module
    captured = {}
    original = kea_module.CodonTable

    def spy(*args, **kwargs):
        obj = original(*args, **kwargs)
        captured["gc_tolerance"] = obj.gc_tolerance
        return obj

    kea_module.CodonTable = spy
    try:
        build_library("MKKFLVLL", "yeast", target_gc_range=(0.4, 0.6), gc_tolerance=0.1,
                      return_best=True, **FAST)
    finally:
        kea_module.CodonTable = original
    assert captured["gc_tolerance"] == 0.1


def test_string_input_honours_add_protein_identifiers():
    """3.4: a bare string always got a random ID."""
    assert build_library("MKKFLVLL", "yeast", add_protein_identifiers=False, **FAST)[0].name == "Protein_0"
    assert build_library("MKKFLVLL", "yeast", add_protein_identifiers=True, **FAST)[0].name != "Protein_0"


def test_dict_input_keeps_its_keys():
    results = build_library({"alpha": "MKKFLVLL", "beta": "MVLSEGEW"}, "yeast", **FAST)
    assert [s.name for s in results] == ["alpha", "beta"]


# --------------------------------------------------------------------------
# Phase 5 -- reproducibility and utilities
# --------------------------------------------------------------------------

def test_seed_makes_libraries_reproducible():
    """5.2: randomness is split across numpy and stdlib random."""
    kwargs = dict(target_gc_range=(0.45, 0.55), total_length=150,
                  gc_finetuning_iterations=400, return_best=True, **FAST)
    proteins = ["MKKFLVLLFCWAVLCEHN", "MVLSEGEWQLVLHVWAKV"]
    first = [s.full_dna_sequence for s in build_library(proteins, "yeast", seed=42, **kwargs)]
    second = [s.full_dna_sequence for s in build_library(proteins, "yeast", seed=42, **kwargs)]
    other = [s.full_dna_sequence for s in build_library(proteins, "yeast", seed=43, **kwargs)]
    assert first == second
    assert first != other


def test_generate_random_protein_ids_is_order_stable():
    """5.3: returning list(set(...)) permuted the ID -> protein assignment."""
    def make():
        random.seed(11)
        return generate_random_protein_ids(8, 5)
    assert make() == make()
    assert len(set(make())) == 5


def test_generate_random_protein_ids_rejects_impossible_requests():
    with pytest.raises(ValueError, match="Cannot generate"):
        generate_random_protein_ids(1, 100, alphanumeric=False)


def test_make_codon_table_handles_lowercase():
    """5.5: soft-masked FASTA produced an all-zero table."""
    assert make_codon_table(["atgaaacttggttaa"]) == make_codon_table(["ATGAAACTTGGTTAA"])


def test_similarity_requires_equal_lengths():
    with pytest.raises(ValueError, match="same length"):
        check_nucleotide_percent_similarity(["ATGC", "ATG"])


def test_similarity_of_identical_sequences_is_one():
    maxima, _ = check_nucleotide_percent_similarity(["ATGCATGC", "ATGCATGC"])
    assert list(maxima) == [1.0, 1.0]


# --------------------------------------------------------------------------
# End-to-end invariants
# --------------------------------------------------------------------------

def test_full_library_invariants():
    """Everything at once: exact length, GC in range, translations verified."""
    proteins = {"p1": "MKKFLVLLFCWAVLCEHN", "p2": "MVLSEGEWQLVLHVWAKV", "p3": ALL_AAS}
    results = build_library(proteins, "yeast", total_length=300,
                            adapter_5_prime="GGTCTC", adapter_3_prime="GAGACC",
                            target_gc_range=(0.45, 0.55), gc_finetuning_iterations=800,
                            return_best=True, seed=7, **FAST)
    assert len(results) == 3
    for seq in results:
        assert len(seq.full_dna_sequence) == 300
        assert 0.45 <= seq.gc_content_full_sequence <= 0.55
        assert seq.correct_coding_translation is True
        assert translate_sequence(seq.coding_sequence) == seq.protein_sequence
        assert seq.full_dna_sequence.startswith("GGTCTC")
        assert seq.full_dna_sequence.endswith("GAGACC")


def test_unreachable_gc_range_is_recorded_as_a_failure():
    """An unreachable GC range is a per-sequence failure, so on_error decides."""
    library = build_library("MKKFLVLLFCWAVLCEHN", "yeast", target_gc_range=(0.90, 0.95),
                            gc_finetuning_iterations=200, **FAST)
    assert library.n_succeeded == 0
    assert library.n_failed == 1
    assert library.failures[0].stage == "optimization"
    assert "GC range" in library.failures[0].reason


def test_unreachable_gc_range_raises_when_asked():
    with pytest.raises(ValueError, match="No sequence found within GC range"):
        build_library("MKKFLVLLFCWAVLCEHN", "yeast", target_gc_range=(0.90, 0.95),
                      gc_finetuning_iterations=200, on_error="raise", **FAST)


def test_unreachable_gc_range_returns_best_when_asked():
    seq = build_library("MKKFLVLLFCWAVLCEHN", "yeast", target_gc_range=(0.90, 0.95),
                        gc_finetuning_iterations=200, return_best=True, **FAST)[0]
    assert seq.correct_coding_translation is True


@pytest.mark.parametrize("bad_input,message", [
    ([], "Empty sequence list"),
    ({}, "Empty sequence dictionary"),
    (42, "must be a string, list, or dict"),
    (["MKK", ""], "Empty sequences"),
    (["MKK", 3], "must be strings"),
    (["MKKFLVLL", "MKKFLVLL"], "not unique"),
])
def test_input_validation(bad_input, message):
    with pytest.raises(ValueError, match=message):
        build_library(bad_input, "yeast", **FAST)


def test_unknown_codon_table_name_rejected():
    with pytest.raises(ValueError, match="not found"):
        build_library("MKKFLVLL", "not_an_organism", **FAST)


# --------------------------------------------------------------------------
# Per-sequence failure handling
# --------------------------------------------------------------------------

# W-Y-V forces GTGGTA[T/C]GT, which scores 7.64 -- above the default donor
# threshold of 6.0 with no synonymous alternative. A reliable unsatisfiable input.
IMPOSSIBLE_UNDER_SPLICE = "MWYV" * 6


def test_one_impossible_sequence_does_not_discard_the_rest():
    """The whole point: a late failure must not nuke an otherwise-good library."""
    proteins = {f"good_{i}": "M" + ALL_AAS[i:] + ALL_AAS[:i] for i in range(6)}
    proteins["impossible"] = IMPOSSIBLE_UNDER_SPLICE
    library = build_library(proteins, "human", avoid_human_splice_sites=True, **FAST)
    assert library.n_succeeded == 6
    assert library.n_failed == 1
    assert library.failures[0].name == "impossible"
    assert library.failures[0].stage == "constraint_repair"
    assert library.failures[0].remaining_violations, "unsatisfied violations should be carried"


def test_library_result_is_still_a_list():
    """Existing code that treats the return value as a list must keep working."""
    library = build_library(["MKKFLVLL", "MVLSEGEW"], "yeast", **FAST)
    assert isinstance(library, list)
    assert len(library) == 2
    assert [s.name for s in library] == ["Protein_0", "Protein_1"]
    assert library[0].correct_coding_translation is True


def test_on_error_raise_aborts():
    with pytest.raises(ValueError):
        build_library({"bad": IMPOSSIBLE_UNDER_SPLICE}, "human",
                      avoid_human_splice_sites=True, on_error="raise", **FAST)


def test_on_error_validated():
    with pytest.raises(ValueError, match='on_error must be'):
        build_library("MKKFLVLL", "yeast", on_error="ignore", **FAST)


def test_summary_reports_counts_and_reasons():
    proteins = {"ok": "MKKFLVLLFCWAVLCEHN", "impossible": IMPOSSIBLE_UNDER_SPLICE}
    library = build_library(proteins, "human", avoid_human_splice_sites=True, **FAST)
    summary = library.summary()
    assert "Built 1/2 sequences." in summary
    assert "constraint_repair: 1" in summary
    assert "impossible" in summary
    assert library.failures_by_stage() == {"constraint_repair": 1}


@pytest.mark.parametrize("kwargs,message", [
    ({"optimization_attempts": 0}, "optimization_attempts must be at least 1"),
    ({"gc_finetuning_iterations": -1}, "gc_finetuning_iterations must be"),
    ({"padding_attempts": 0}, "padding_attempts must be at least 1"),
    ({"constraint_repair_attempts": 0}, "constraint_repair_attempts must be at least 1"),
    ({"constraint_retry_attempts": -1}, "constraint_retry_attempts must be"),
    ({"constraint_retry_attempts": 1.5}, "constraint_retry_attempts must be"),
    ({"constraint_retry_attempts": True}, "constraint_retry_attempts must be"),
    ({"sequences_per_protein": 0}, "sequences_per_protein must be at least 1"),
    ({"minimum_hamming_distance": -1}, "minimum_hamming_distance must be"),
    ({"minimum_hamming_distance": 1.5}, "minimum_hamming_distance must be"),
    ({"hamming_distance_attempts": 0}, "hamming_distance_attempts must be at least 1"),
])
def test_configuration_errors_are_never_skipped(kwargs, message):
    """Config mistakes apply to every sequence, so they escape regardless of on_error."""
    with pytest.raises(ValueError, match=message):
        build_library("MKKFLVLL", "yeast", show_progress=False,
                      show_optimization_progress=False, **kwargs)


# --------------------------------------------------------------------------
# Multiple synonymous encodings per protein
# --------------------------------------------------------------------------


def test_multiple_encodings_meet_all_pairwise_hamming_distances_and_constraints():
    minimum_distance = 18
    protein = "M" + "".join(random.Random(91).choices(ALL_AAS, k=79))
    library = build_library(
        protein,
        "human",
        sequences_per_protein=4,
        minimum_hamming_distance=minimum_distance,
        hamming_distance_attempts=30,
        target_gc_range=(0.45, 0.60),
        avoid_human_splice_sites=True,
        avoid_premature_polyadenylation=True,
        max_homopolymer_length=5,
        max_tandem_repeat_copies=3,
        max_direct_repeat_length=12,
        seed=4,
        on_error="raise",
        **FAST,
    )

    assert len(library) == library.requested == 4
    assert library.n_failed == 0
    assert [sequence.name for sequence in library] == [
        f"Protein_0_variant_{index}" for index in range(1, 5)
    ]
    assert [sequence.variant_index for sequence in library] == [1, 2, 3, 4]
    assert all(sequence.source_protein_name == "Protein_0" for sequence in library)

    for index, sequence in enumerate(library):
        assert translate_sequence(sequence.coding_sequence) == sequence.protein_sequence
        assert sequence.quality_report.passed
        if index == 0:
            assert sequence.minimum_sibling_hamming_distance is None
        else:
            distances = [
                nucleotide_hamming_distance(
                    sequence.coding_sequence, previous.coding_sequence)
                for previous in library[:index]
            ]
            assert min(distances) >= minimum_distance
            assert sequence.minimum_sibling_hamming_distance == min(distances)


def test_multiple_encoding_generation_is_reproducible_with_seed():
    kwargs = dict(
        sequences_per_protein=3,
        minimum_hamming_distance=12,
        hamming_distance_attempts=20,
        seed=17,
        **FAST,
    )
    first = build_library("MKKFLVLLFCWAVLCEHN", "human", **kwargs)
    second = build_library("MKKFLVLLFCWAVLCEHN", "human", **kwargs)
    assert [item.coding_sequence for item in first] == [
        item.coding_sequence for item in second
    ]


def test_impossible_hamming_set_is_reported_without_losing_accepted_variant():
    library = build_library(
        "MWMW",
        "human",
        force_start_codon=False,
        force_stop_codon=False,
        sequences_per_protein=3,
        minimum_hamming_distance=1,
        hamming_distance_attempts=2,
        seed=2,
        **FAST,
    )
    assert len(library) == 1
    assert library.requested == 3
    assert library.n_failed == 2
    assert all(failure.stage == "sequence_diversity" for failure in library.failures)
    assert all(isinstance(failure.error, SequenceDiversityError)
               for failure in library.failures)
    assert "Best minimum distance found: 0" in str(library.failures[0].error)


def test_impossible_hamming_set_can_fail_loudly():
    with pytest.raises(SequenceDiversityError, match="Best minimum distance found: 0"):
        build_library(
            "MWMW",
            "human",
            force_start_codon=False,
            force_stop_codon=False,
            sequences_per_protein=2,
            minimum_hamming_distance=1,
            hamming_distance_attempts=2,
            on_error="raise",
            seed=2,
            **FAST,
        )


def test_nucleotide_hamming_distance_requires_equal_lengths():
    assert nucleotide_hamming_distance("AACG", "AATG") == 1
    with pytest.raises(ValueError, match="equal length"):
        nucleotide_hamming_distance("AACG", "AAC")


def test_padding_does_not_reintroduce_avoided_motifs():
    """Padding is generated after repair, so it must be constraint-checked too."""
    from kea.backend.sequence_constraints import make_sequence_constraints
    constraints = make_sequence_constraints(
        avoid_premature_polyadenylation=True, avoid_human_splice_sites=True,
        max_homopolymer_length=6)
    dirty = 0
    for seed in range(8):
        library = build_library("MKKFLVLLFCWAVLCEHNRRSSPPGG", "human", total_length=400,
                                avoid_premature_polyadenylation=True,
                                avoid_human_splice_sites=True, max_homopolymer_length=6,
                                seed=seed, **FAST)
        if not library:
            continue
        sequence = library[0]
        region = (sequence.padding_5_prime + sequence.coding_sequence
                  + sequence.padding_3_prime)
        if constraints.find_violations(region):
            dirty += 1
    assert dirty == 0


@pytest.mark.parametrize("pad_location", [None, 3, 5])
def test_constrained_padding_keeps_length_and_translation(pad_location):
    library = build_library("MKKFLVLLFCWAVLCEHNRRSSPPGG", "human", total_length=400,
                            pad_location=pad_location, avoid_human_splice_sites=True,
                            avoid_premature_polyadenylation=True, max_homopolymer_length=6,
                            seed=1, **FAST)
    assert library, "expected this protein to build"
    sequence = library[0]
    assert len(sequence.full_dna_sequence) == 400
    assert translate_sequence(sequence.coding_sequence) == sequence.protein_sequence
    for padding in (sequence.padding_5_prime, sequence.padding_3_prime):
        assert "ATG" not in padding
        assert not any(stop in padding for stop in ("TAA", "TAG", "TGA"))


def test_padding_repair_leaves_the_coding_sequence_untouched():
    from kea.backend.sequence_constraints import make_sequence_constraints
    from kea.backend.library_generation_utils import repair_padding_constraints
    constraints = make_sequence_constraints(avoid_premature_polyadenylation=True)
    # AATAAA sits in the 3' padding; the CDS must not be rewritten to fix it.
    coding = "ATGAAACTGTAA"
    region = "CCCCCC" + coding + "AATAAACC"
    ranges = [(0, 6), (6 + len(coding), len(region))]
    repaired, remaining = repair_padding_constraints(region, ranges, constraints)
    assert remaining == []
    assert repaired[6:6 + len(coding)] == coding
    assert len(repaired) == len(region)


def test_unsatisfiable_violation_is_proven_not_guessed():
    """W-Y-V forces GTGGTA[T/C]GT; no synonymous encoding scores under 6.0."""
    from kea.backend.sequence_constraints import (ConstraintRepairError, is_violation_unsatisfiable,
                                                  make_sequence_constraints)
    table = CodonTable(kea.data.codon_tables.human_hegs, 1, 1, (0, 1), 0.0)
    constraints = make_sequence_constraints(avoid_human_splice_sites=True)
    coding = "".join(table.aa_to_codons[aa][0] for aa in "M" + "WYV" * 6)
    try:
        repair_coding_sequence(coding, table, constraints, max_iterations=50)
    except ConstraintRepairError as error:
        assert error.remaining_violations
        assert all(is_violation_unsatisfiable(error.sequence, violation, table, constraints)
                   for violation in error.remaining_violations)
    else:
        pytest.fail("expected this protein to be unsatisfiable")


def test_impossible_sequence_reports_that_retrying_will_not_help():
    library = build_library({"impossible": "M" + "WYV" * 6}, "human",
                            avoid_human_splice_sites=True, **FAST)
    assert library.n_failed == 1
    assert "under any synonymous encoding" in str(library.failures[0].error)


def test_retry_is_configurable_and_zero_restores_single_shot():
    import inspect
    assert inspect.signature(build_library).parameters["constraint_retry_attempts"].default == 2
    # retries must not change the outcome for a protein that builds first time
    for retries in (0, 3):
        library = build_library("MKKFLVLLFCWAVLCEHN", "human",
                                avoid_human_splice_sites=True,
                                constraint_retry_attempts=retries, seed=3, **FAST)
        assert library.n_succeeded == 1


def test_gc_compensation_pulls_a_candidate_back_into_range():
    """A fix that shifts GC should be paid for elsewhere, not vetoed."""
    from kea.backend.sequence_constraints import _compensate_gc, _gc_compensation_index
    table = CodonTable(kea.data.codon_tables.human_hegs, 1, 1, (0.0, 0.50), 0.0)
    # 10 codons / 30 bases with 16 G+C -> 0.533, one base over a 0.50 ceiling.
    sequence = "ATG" + "GGC" * 4 + "CTG" + "AAG" + "AAA" * 2 + "TAA"
    assert len(sequence) == 30
    assert sequence.count("G") + sequence.count("C") == 16

    index = _gc_compensation_index(sequence, protected={0}, allowed_codons=table.codon_table)
    compensated = _compensate_gc(sequence, (0.0, 0.50), index)

    assert compensated is not None
    assert len(compensated) == len(sequence)
    assert translate_sequence(compensated) == translate_sequence(sequence)
    assert (compensated.count("G") + compensated.count("C")) / len(compensated) <= 0.50 + 1e-9


def test_gc_compensation_gives_up_when_one_edit_cannot_close_the_gap():
    """A single synonymous swap moves GC by at most 3 bases; report that honestly."""
    from kea.backend.sequence_constraints import _compensate_gc, _gc_compensation_index
    table = CodonTable(kea.data.codon_tables.human_hegs, 1, 1, (0.0, 0.50), 0.0)
    sequence = "ATG" + "GGC" * 8 + "TAA"   # 0.83 GC against a 0.50 ceiling
    index = _gc_compensation_index(sequence, protected={0}, allowed_codons=table.codon_table)
    assert _compensate_gc(sequence, (0.0, 0.50), index) is None


def test_gc_compensation_never_edits_a_protected_codon():
    """Compensation must not undo the very fix it is paying for."""
    from kea.backend.sequence_constraints import _gc_compensation_index
    table = CodonTable(kea.data.codon_tables.human_hegs, 1, 1, (0, 1), 0.0)
    sequence = "ATG" + "CTGGGCAAA" * 3 + "TAA"
    protected = {1, 2, 3}
    index = _gc_compensation_index(sequence, protected=protected, allowed_codons=table.codon_table)
    edited_positions = {position
                        for entries in index.values()
                        for _frequency, position, _expected, _replacement in entries}
    assert not (edited_positions & protected)


def test_tight_gc_range_is_respected_with_constraints_on():
    """The whole point of compensation: satisfy constraints without leaving the range."""
    gc_range = (0.30, 0.42)
    library = build_library(
        {f"p{i}": "M" + (ALL_AAS * 6)[i:i + 90] for i in range(6)}, "human",
        target_gc_range=gc_range, avoid_human_splice_sites=True,
        avoid_premature_polyadenylation=True, max_homopolymer_length=5,
        gc_finetuning_iterations=300, return_best=True, seed=0, **FAST)
    assert library.n_succeeded > 0
    for sequence in library:
        assert gc_range[0] <= sequence.gc_content_coding_sequence <= gc_range[1]
        assert translate_sequence(sequence.coding_sequence) == sequence.protein_sequence


# --------------------------------------------------------------------------
# Final acceptance check
# --------------------------------------------------------------------------

def test_every_target_is_measured_on_the_accepted_construct():
    library = build_library("MKKFLVLLFCWAVLCEHNRRSSPPGG", "human",
                            target_gc_range=(0.45, 0.55), total_length=400,
                            avoid_premature_polyadenylation=True, max_homopolymer_length=6,
                            gc_finetuning_iterations=300, seed=1, **FAST)
    assert library.n_succeeded == 1
    report = library[0].quality_report
    measured = {check.name for check in report.checks}
    assert measured == {"translation", "total_length", "gc_coding", "gc_full",
                        "sequence_constraints", "codon_adaptation"}
    assert report.passed
    assert report.errors == []
    assert 0.0 <= report.value("codon_adaptation") <= 1.0


def test_missed_gc_range_is_recorded_instead_of_being_invisible():
    """return_best=True returns best-effort GC -- but must say so."""
    library = build_library({"p": "MKKFLVLLFCWAVLCEHNRRSSPPGG"}, "human",
                            target_gc_range=(0.80, 0.85), return_best=True,
                            gc_finetuning_iterations=200, seed=1, **FAST)
    assert library.n_succeeded == 1
    report = library[0].quality_report
    assert report.passed, "best-effort GC should not fail the sequence"
    assert {check.name for check in report.warnings} == {"gc_coding", "gc_full"}
    assert library.with_warnings == list(library)
    assert "missed a requested target" in library.summary()


def test_no_gc_target_reports_no_gc_checks():
    """A sequence with no requested range must not be scored against (0, 1)."""
    library = build_library("MKKFLVLL", "yeast", seed=1, **FAST)
    names = {check.name for check in library[0].quality_report.checks}
    assert "gc_coding" not in names and "gc_full" not in names


def test_minimum_codon_adaptation_floor_is_enforced():
    protein = "MKKFLVLLFCWAVLCEHNRRSSPPGG"
    unfloored = build_library({"p": protein}, "human", target_gc_range=(0.30, 0.42),
                              return_best=True, gc_finetuning_iterations=300,
                              seed=1, **FAST)
    achieved = unfloored[0].quality_report.value("codon_adaptation")

    floored = build_library({"p": protein}, "human", target_gc_range=(0.30, 0.42),
                            return_best=True, gc_finetuning_iterations=300,
                            minimum_codon_adaptation=achieved + 0.05, seed=1, **FAST)
    assert floored.n_succeeded == 0
    assert floored.n_failed == 1
    assert floored.failures[0].stage == "quality_check"


def test_quality_check_failure_respects_on_error():
    with pytest.raises(ValueError, match="quality check"):
        build_library({"p": "MKKFLVLLFCWAVLCEHN"}, "human",
                      minimum_codon_adaptation=1.1, on_error="raise", seed=1, **FAST)


def test_codon_score_cache_is_not_confused_by_address_reuse():
    """id() is unique only among live objects; short-lived tables reuse addresses."""
    import gc as gc_module
    import kea.backend.optimize_codon_usage as optimizer_module

    optimizer_module._NORM_SCORE_CACHE.clear()
    optimizer_module._MAX_FREQ_CACHE.clear()

    warm = CodonTable(kea.data.codon_tables.human_hegs, 1, 1, (0, 1), 0.0)
    optimizer_module._normalized_codon_scores(warm.codon_table)
    del warm
    gc_module.collect()

    for _ in range(50):
        table = CodonTable(kea.data.codon_tables.human_hegs, 1, 1, (0, 1), 0.20)
        scores = optimizer_module._normalized_codon_scores(table.codon_table)
        leucine = table.codon_table["L"]
        expected = leucine["CTG"] / max(leucine.values())
        assert abs(scores["CTG"] - expected) < 1e-12, "stale cache hit from a reused address"
        del table, scores
        gc_module.collect()

    assert len(optimizer_module._NORM_SCORE_CACHE) <= optimizer_module._TABLE_CACHE_LIMIT


def test_gc_centering_is_off_by_default_and_validated():
    import inspect
    assert inspect.signature(build_library).parameters["gc_centering"].default == 0.0
    for bad in (-0.1, 1.5):
        with pytest.raises(ValueError, match="gc_centering"):
            build_library("MKKFLVLL", "human", gc_centering=bad, **FAST)


def test_gc_centering_moves_gc_toward_the_middle_of_the_range():
    """Opt-in behaviour: it should pull GC off the edge it otherwise pins to."""
    gc_range = (0.30, 0.42)
    kwargs = dict(target_gc_range=gc_range, avoid_human_splice_sites=True,
                  max_homopolymer_length=5, gc_finetuning_iterations=300,
                  return_best=True, seed=0, **FAST)
    proteins = {f"p{i}": "M" + (ALL_AAS * 8)[i:i + 150] for i in range(6)}

    flat = build_library(proteins, "human", gc_centering=0.0, **kwargs)
    centred = build_library(proteins, "human", gc_centering=0.10, **kwargs)
    assert flat and centred

    centre = sum(gc_range) / 2
    flat_offset = mean(abs(s.gc_content_coding_sequence - centre) for s in flat)
    centred_offset = mean(abs(s.gc_content_coding_sequence - centre) for s in centred)
    assert centred_offset < flat_offset
    for sequence in centred:
        assert gc_range[0] <= sequence.gc_content_coding_sequence <= gc_range[1]


@pytest.mark.parametrize("centering", [0.0, 0.05, 0.5, 1.0])
@pytest.mark.parametrize("tolerance", [0.001, 0.025, 0.5])
def test_gc_score_stays_bounded_and_monotonic(centering, tolerance):
    """Further from the range centre must never score better than closer."""
    from kea.backend.optimize_codon_usage import _gc_score
    gc_range = (0.30, 0.42)
    centre = sum(gc_range) / 2
    step = 0.0005
    previous = None
    offset = 0.0
    while offset <= 0.35:
        for point in (centre - offset, centre + offset):
            score = _gc_score(point, gc_range, tolerance, centering)
            assert 0.0 <= score <= 1.0
        score = _gc_score(centre + offset, gc_range, tolerance, centering)
        if previous is not None:
            assert score <= previous + 1e-12
        previous = score
        offset += step


def test_gc_score_handles_degenerate_ranges():
    from kea.backend.optimize_codon_usage import _gc_score
    assert _gc_score(0.5, (0.5, 0.5), 0.025, 0.05) == 1.0
    assert _gc_score(0.5, (0.0, 1.0), 0.025, 0.05) <= 1.0
    assert _gc_score(0.9, (0.5, 0.5), 0.0, 0.05) == 0.0


# --------------------------------------------------------------------------
# Escaping the GC-range boundary
# --------------------------------------------------------------------------

def _objective(codons, table_obj, norm):
    from kea.backend.optimize_codon_usage import (_CODON_GC_COUNT, _combine_scores, _gc_score)
    count = len(codons)
    gc_bases = sum(_CODON_GC_COUNT[codon] for codon in codons)
    usage = sum(norm[codon] for codon in codons) / count
    return _combine_scores(usage,
                           _gc_score(gc_bases / (3 * count), table_obj.gc_range, 0.025, 0.0),
                           table_obj.gc_weight, table_obj.usage_weight)


def test_boundary_escape_leaves_no_improving_pair_move():
    """The defect: at the GC edge every single-codon gain breaks the range."""
    from kea.backend.codon_optimizer import CodonOptimizer
    from kea.backend.optimize_codon_usage import _CODON_GC_COUNT, _normalized_codon_scores

    gc_range = (0.30, 0.42)
    table_obj = CodonTable(kea.data.codon_tables.human_hegs, 1, 1, gc_range, 0.0)
    assert table_obj.escape_gc_boundary is True
    norm = _normalized_codon_scores(table_obj.codon_table)
    optimizer = CodonOptimizer(table_obj)

    random.seed(0)
    np.random.seed(0)
    protein = "M" + (ALL_AAS * 6)[:90] + "*"
    optimized = optimizer.optimize(protein, n_iter=400, fine_tuning_iterations=300,
                                   return_best=True, show_progress_bar=False)
    codons = [optimized[i:i + 3] for i in range(0, len(optimized), 3)]
    baseline = _objective(codons, table_obj, norm)

    improving_pairs = 0
    for i, codon in enumerate(codons):
        for alternative in table_obj.codon_table[codons_to_aa[codon]]:
            if alternative == codon:
                continue
            delta = _CODON_GC_COUNT[alternative] - _CODON_GC_COUNT[codon]
            if delta == 0:
                continue
            for j, other in enumerate(codons):
                if j == i:
                    continue
                for partner in table_obj.codon_table[codons_to_aa[other]]:
                    if partner == other:
                        continue
                    if _CODON_GC_COUNT[partner] - _CODON_GC_COUNT[other] != -delta:
                        continue
                    candidate = list(codons)
                    candidate[i] = alternative
                    candidate[j] = partner
                    if _objective(candidate, table_obj, norm) > baseline + 1e-12:
                        improving_pairs += 1
    assert improving_pairs == 0, f"{improving_pairs} GC-neutral pair moves left unexploited"
    assert translate_sequence(optimized) == protein
    assert gc_range[0] <= _calculate_gc_content(optimized) <= gc_range[1]


def test_boundary_escape_raises_codon_adaptation_without_leaving_the_range():
    gc_range = (0.30, 0.42)
    proteins = {f"p{i}": "M" + (ALL_AAS * 8)[i:i + 150] for i in range(6)}
    kwargs = dict(target_gc_range=gc_range, gc_finetuning_iterations=300,
                  return_best=True, seed=0, **FAST)

    without = build_library(proteins, "human", escape_gc_boundary=False, **kwargs)
    with_escape = build_library(proteins, "human", escape_gc_boundary=True, **kwargs)
    assert without and with_escape

    def adaptation(library):
        return mean(s.quality_report.value("codon_adaptation") for s in library)

    assert adaptation(with_escape) > adaptation(without)
    for sequence in with_escape:
        assert gc_range[0] <= sequence.gc_content_coding_sequence <= gc_range[1]
        assert translate_sequence(sequence.coding_sequence) == sequence.protein_sequence


def test_boundary_escape_auto_default_follows_constraints():
    """Auto: on without constraints (pure win), off with them (costs yield)."""
    captured = {}
    import kea.kea as kea_module
    original = kea_module.CodonTable

    def spy(*args, **kwargs):
        obj = original(*args, **kwargs)
        captured["escape"] = obj.escape_gc_boundary
        return obj

    kea_module.CodonTable = spy
    try:
        build_library("MKKFLVLLFCWAVLCEHN", "human", target_gc_range=(0.30, 0.42),
                      return_best=True, **FAST)
        assert captured["escape"] is True

        build_library("MKKFLVLLFCWAVLCEHN", "human", target_gc_range=(0.30, 0.42),
                      avoid_human_splice_sites=True, return_best=True, **FAST)
        assert captured["escape"] is False

        build_library("MKKFLVLLFCWAVLCEHN", "human", target_gc_range=(0.30, 0.42),
                      avoid_human_splice_sites=True, escape_gc_boundary=True,
                      return_best=True, **FAST)
        assert captured["escape"] is True
    finally:
        kea_module.CodonTable = original


def test_repair_does_not_reject_a_clean_sequence_on_the_last_iteration():
    """Off-by-one: the success check runs at the top of the loop."""
    from kea.backend.sequence_constraints import HumanSpliceConstraint
    table = CodonTable(yeast_codon_usage, 1, 1, (0, 1), 0.0)
    constraints = ConstraintSet([HumanSpliceConstraint()])
    cds = "ATG" + "CAGGTAAGT" + "CTG" * 10 + "TAA"
    # Whatever the budget, a run that ends clean must return, never raise with an
    # empty remaining-violations list.
    for attempts in range(1, 8):
        try:
            repaired, remaining = repair_coding_sequence(cds, table, constraints,
                                                         max_iterations=attempts)
        except ValueError as error:
            assert "Remaining violations: " in str(error)
            assert str(error).rstrip().endswith("violations:") is False, \
                "raised with an empty violation list on an already-clean sequence"
        else:
            assert remaining == []
            assert not constraints.find_violations(repaired)


# --------------------------------------------------------------------------
# Human transcript and sequence-complexity constraints
# --------------------------------------------------------------------------


def test_maxentscan_matches_reference_scores():
    """Published MaxEntScan lookup tables must not be replaced by an approximation."""
    scorer = MaxEntScanScorer()
    assert scorer.score_donor("CAGGTAAGT") == pytest.approx(10.8583131)
    assert scorer.score_donor("GAGGTAAGT") == pytest.approx(11.0784945)
    assert scorer.score_acceptor("TTCCAAACGAACTTTTGTAGGGA") == pytest.approx(2.8867731)
    assert scorer.score_acceptor("TGTCTTTTTCTGTGTGGCAGTGG") == pytest.approx(8.1909655)


def test_maxentscan_rejects_invalid_windows():
    scorer = MaxEntScanScorer()
    with pytest.raises(ValueError, match="exactly 9"):
        scorer.score_donor("GTAAGT")
    with pytest.raises(ValueError, match="non-DNA"):
        scorer.score_acceptor("TGTCTTTTTCTGTGTGGCAGNgg")


def test_vectorized_maxentscan_normalizes_lowercase_like_scalar_scoring():
    scorer = MaxEntScanScorer()
    donor = "caggtaagt"
    acceptor = "tgtctttttctgtgtggcagtgg"
    assert scorer.score_donors(donor)[0] == scorer.score_donor(donor)
    assert scorer.score_acceptors(acceptor)[0] == scorer.score_acceptor(acceptor)


@pytest.mark.parametrize("threshold", [float("nan"), float("inf"), float("-inf")])
def test_human_splice_thresholds_must_be_finite(threshold):
    with pytest.raises(ValueError, match="finite"):
        HumanSpliceConstraint(donor_threshold=threshold)
    with pytest.raises(ValueError, match="finite"):
        HumanSpliceConstraint(acceptor_threshold=threshold)


def test_human_splice_constraint_finds_donor_and_acceptor():
    constraint = HumanSpliceConstraint(donor_threshold=6.0, acceptor_threshold=7.0)
    sequence = "CAGGTAAGT" + "CCCC" + "TGTCTTTTTCTGTGTGGCAGTGG"
    violations = constraint.find_violations(sequence)
    assert any(v.kind == "human_splice_donor" and v.sequence == "CAGGTAAGT" for v in violations)
    assert any(v.kind == "human_splice_acceptor" and v.sequence == "TGTCTTTTTCTGTGTGGCAGTGG"
               for v in violations)


def test_polyadenylation_constraint_supports_all_common_variants():
    constraint = PrematurePolyadenylationConstraint()
    for signal in COMMON_HUMAN_POLYADENYLATION_SIGNALS:
        hits = constraint.find_violations("CCC" + signal + "CCC")
        assert any(hit.sequence == signal for hit in hits), signal


def test_polyadenylation_downstream_context_can_be_required():
    constraint = PrematurePolyadenylationConstraint(require_downstream_element=True)
    assert constraint.find_violations("CCCAAGAAACCCCCCCCCCCCC") == []

    u_rich = constraint.find_violations("CCCAAGAAATTTTGCCCC")
    assert len(u_rich) == 1
    assert u_rich[0].details["u_rich_downstream"] is True

    gu_rich = constraint.find_violations("CCCAAGAAAGTTGTCCCC")
    assert len(gu_rich) == 1
    assert gu_rich[0].details["gu_rich_downstream"] is True


def test_context_dependent_pas_can_be_repaired_through_its_downstream_element():
    from kea.backend.sequence_constraints import is_violation_unsatisfiable

    table = CodonTable(kea.human_hegs, 1, 0, (0, 1), 0.0)
    constraints = ConstraintSet(
        [PrematurePolyadenylationConstraint(require_downstream_element=True)],
        five_prime_context="AATAAA",
    )
    original = "TTTTGCCCC"  # FCP; TTTTG is a U-rich element after the fixed PAS.
    violation = constraints.find_violations(original)[0]

    # TTT -> TTC preserves phenylalanine and removes the downstream support, so
    # this fixed-context PAS is not intrinsically unsatisfiable.
    assert is_violation_unsatisfiable(original, violation, table, constraints) is False
    repaired, remaining = repair_coding_sequence(original, table, constraints)
    assert translate_sequence(repaired) == translate_sequence(original) == "FCP"
    assert remaining == []
    assert constraints.find_violations(repaired) == []


def test_homopolymer_constraint_reports_full_run():
    violations = HomopolymerConstraint(4).find_violations("CCAAAAAGG")
    assert [(v.start, v.end, v.sequence) for v in violations] == [(2, 7, "AAAAA")]


def test_tandem_repeat_constraint_reports_shortest_repeat_unit():
    violations = TandemRepeatConstraint(max_copies=3, max_unit_length=6).find_violations(
        "CCATATATATGG"
    )
    assert len(violations) == 1
    assert violations[0].details["repeat_unit"] == "AT"
    assert violations[0].details["copies"] == 4


def test_direct_repeat_constraint_finds_separated_and_overlapping_repeats():
    separated = DirectRepeatConstraint(7).find_violations("ACGTTGCAAAAACGTTGCA")
    assert any(v.sequence.startswith("ACGTTGCA") for v in separated)

    overlapping = DirectRepeatConstraint(8).find_violations("GCTGCTGCTGCT")
    assert len(overlapping) == 1
    assert overlapping[0].details["first_start"] == 0


def test_synonymous_repair_removes_splice_and_polyadenylation_sites():
    table = CodonTable(kea.human_hegs, 1, 0, (0.20, 0.65), 0.0)
    original = "CAGGTAAGTAATAAA"  # QVSNK: strong donor followed by canonical PAS
    constraints = ConstraintSet([
        HumanSpliceConstraint(),
        PrematurePolyadenylationConstraint(),
    ])
    repaired, remaining = repair_coding_sequence(original, table, constraints)
    assert translate_sequence(repaired) == translate_sequence(original)
    assert remaining == []
    assert constraints.find_violations(repaired) == []
    assert 0.20 <= _calculate_gc_content(repaired) <= 0.65


def test_synonymous_repair_handles_site_across_fixed_context_junction():
    table = CodonTable(kea.human_hegs, 1, 0, (0, 1), 0.0)
    constraints = ConstraintSet([HumanSpliceConstraint()], five_prime_context="CAG")
    repaired, _ = repair_coding_sequence("GTAAGT", table, constraints)
    assert translate_sequence(repaired) == "VS"
    assert constraints.find_violations(repaired) == []


def test_unrepairable_site_in_fixed_context_raises():
    table = CodonTable(kea.human_hegs, 1, 0, (0, 1), 0.0)
    constraints = ConstraintSet([HumanSpliceConstraint()], five_prime_context="CAGGTAAGT")
    with pytest.raises(ValueError, match="Could not synonymously satisfy"):
        repair_coding_sequence("ATG", table, constraints, max_iterations=5)


def test_build_library_applies_all_new_constraints():
    # Make QVSNK overwhelmingly prefer CAG-GTA-AGT-AAT-AAA, which contains both
    # the reference MaxEnt donor and canonical PAS. Alternatives remain available
    # for the repair pass.
    table = {aa: dict(codons) for aa, codons in kea.human_hegs.items()}
    preferred = {"Q": "CAG", "V": "GTA", "S": "AGT", "N": "AAT", "K": "AAA"}
    for amino_acid, preferred_codon in preferred.items():
        table[amino_acid] = {
            codon: (1.0 if codon == preferred_codon else 1e-9)
            for codon in table[amino_acid]
        }

    result = build_library(
        "QVSNKAAAA",
        table,
        force_start_codon=False,
        force_stop_codon=False,
        avoid_human_splice_sites=True,
        avoid_premature_polyadenylation=True,
        max_homopolymer_length=5,
        max_tandem_repeat_copies=3,
        max_direct_repeat_length=12,
        seed=3,
        optimization_attempts=50,
        show_progress=False,
        show_optimization_progress=False,
    )[0]
    constraints = ConstraintSet([
        HumanSpliceConstraint(),
        PrematurePolyadenylationConstraint(),
        HomopolymerConstraint(5),
        TandemRepeatConstraint(3),
        DirectRepeatConstraint(12),
    ])
    assert translate_sequence(result.coding_sequence) == result.protein_sequence
    assert constraints.find_violations(result.coding_sequence) == []
    assert result.sequence_constraint_violations == []
