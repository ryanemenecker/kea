"""
Unit and regression tests for the kea package.
"""

import sys

import matplotlib

matplotlib.use("Agg")  # headless backend for plotting tests

import pytest

from kea import build_library, save_library
from kea.backend.translation_utils import find_first_orf, translate_sequence
from kea.backend.optimize_codon_usage import (
    _calculate_gc_content,
    calculate_codon_adaptation_score,
)
from kea.backend.kea_utils import (
    generate_random_protein_ids,
    make_codon_table,
    read_fasta,
    codon_from_fasta,
)
from kea.backend.library_generation_utils import check_nucleotide_percent_similarity
from kea.backend.codon_table import CodonTable
from kea.backend.sequence import Sequence
from kea.backend.restriction_utils import (
    parse_restriction_site,
    site_to_regex,
    identify_restriction_sites,
    check_restriction_sites,
)
from kea.data.codon_tables import all_codon_tables, yeast_codon_usage


# Keep optimization cheap so the suite runs quickly.
FAST = dict(
    optimization_attempts=200,
    gc_finetuning_iterations=50,
    show_progress=False,
    show_optimization_progress=False,
)


# ---------------------------------------------------------------------------
# Package import
# ---------------------------------------------------------------------------
def test_kea_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "kea" in sys.modules


def test_public_api_available():
    assert callable(build_library)
    assert callable(save_library)


# ---------------------------------------------------------------------------
# build_library: input handling
# ---------------------------------------------------------------------------
def test_build_library_single_string():
    results = build_library("MKKFLVLL", "yeast", **FAST)
    assert len(results) == 1
    seq_obj = results[0]
    # default forces start (M already present) and stop codon (*)
    assert seq_obj.protein_sequence == "MKKFLVLL*"
    assert seq_obj.correct_coding_translation is True
    assert len(seq_obj.coding_sequence) % 3 == 0


def test_build_library_list():
    seqs = ["MKKFLVLL", "MVLSEGEW"]
    results = build_library(seqs, "yeast", **FAST)
    assert len(results) == 2
    assert [r.name for r in results] == ["Protein_0", "Protein_1"]


def test_build_library_dict_preserves_names():
    seqs = {"alpha": "MKKFLVLL", "beta": "MVLSEGEW"}
    results = build_library(seqs, "yeast", **FAST)
    assert [r.name for r in results] == ["alpha", "beta"]


def test_build_library_translates_correctly():
    results = build_library("MWADEFGH", "yeast", **FAST)
    coding = results[0].coding_sequence
    assert translate_sequence(coding, return_stop_codon=True) == "MWADEFGH*"


# ---------------------------------------------------------------------------
# build_library: validation / error handling
# ---------------------------------------------------------------------------
def test_build_library_empty_list_raises():
    with pytest.raises(ValueError):
        build_library([], "yeast", **FAST)


def test_build_library_empty_dict_raises():
    with pytest.raises(ValueError):
        build_library({}, "yeast", **FAST)


def test_build_library_empty_sequence_raises():
    with pytest.raises(ValueError):
        build_library([""], "yeast", **FAST)


def test_build_library_non_string_sequence_raises():
    with pytest.raises(ValueError):
        build_library([123], "yeast", **FAST)


def test_build_library_bad_input_type_raises():
    with pytest.raises(ValueError):
        build_library(123, "yeast", **FAST)


def test_build_library_unknown_codon_table_raises():
    with pytest.raises(ValueError):
        build_library("MKKF", "not_a_real_table", **FAST)


def test_build_library_duplicate_sequences_raises():
    with pytest.raises(ValueError):
        build_library(["MKKF", "MKKF"], "yeast", **FAST)


def test_build_library_duplicate_allowed_when_disabled():
    results = build_library(
        ["MKKF", "MKKF"],
        "yeast",
        verify_unique_protein_sequences=False,
        **FAST,
    )
    assert len(results) == 2


def test_build_library_invalid_gc_range_raises():
    with pytest.raises(ValueError):
        build_library("MKKF", "yeast", target_gc_range=(0.7, 0.3), **FAST)


def test_build_library_bad_min_codon_probability_raises():
    with pytest.raises(ValueError):
        build_library("MKKF", "yeast", minimum_codon_probability=2.0, **FAST)


# ---------------------------------------------------------------------------
# build_library: GC range, length and adapters
# ---------------------------------------------------------------------------
def test_build_library_respects_gc_range():
    results = build_library(
        "MKKFLVLLFCWAVLCEHN",
        "yeast",
        target_gc_range=(0.4, 0.6),
        gc_tolerance=0.05,
        **FAST,
    )
    gc = results[0].gc_content_coding_sequence
    assert 0.4 - 0.05 <= gc <= 0.6 + 0.05


def test_build_library_total_length_and_adapters():
    results = build_library(
        ["MKKFLVLLFCWAVLCEHN", "MVLSEGEWQLVLHVWAKV"],
        "yeast",
        target_gc_range=(0.4, 0.6),
        total_length=120,
        adapter_5_prime="GGTCTC",
        adapter_3_prime="GAGACC",
        gc_tolerance=0.05,
        **FAST,
    )
    for r in results:
        assert len(r.full_dna_sequence) == 120
        assert r.full_dna_sequence.startswith("GGTCTC")
        assert r.full_dna_sequence.endswith("GAGACC")
        assert r.correct_full_translation is True


def test_build_library_adapter_codon_check():
    # 3' adapter contains a stop codon (TGA) -> should raise when checked
    with pytest.raises(ValueError):
        build_library(
            "MKKF",
            "yeast",
            adapter_3_prime="TGA",
            verify_start_stop_codons_in_adapters=True,
            **FAST,
        )


def test_build_library_no_force_start_stop():
    results = build_library(
        "MKKF",
        "yeast",
        force_start_codon=False,
        force_stop_codon=False,
        **FAST,
    )
    assert results[0].protein_sequence == "MKKF"


# ---------------------------------------------------------------------------
# save_library
# ---------------------------------------------------------------------------
def test_save_library_writes_file(tmp_path):
    results = build_library("MKKFLVLL", "yeast", **FAST)
    out = tmp_path / "library.csv"
    save_library(results, str(out))
    assert out.exists()
    lines = out.read_text().strip().split("\n")
    assert lines[0].startswith("Protein Name,")
    assert len(lines) == 2  # header + one sequence


def test_save_library_bare_filename(tmp_path, monkeypatch):
    # Regression: a bare filename (no directory component) must be accepted.
    monkeypatch.chdir(tmp_path)
    results = build_library("MKKFLVLL", "yeast", **FAST)
    save_library(results, "library.csv")
    assert (tmp_path / "library.csv").exists()


def test_save_library_bad_directory_raises():
    results = build_library("MKKFLVLL", "yeast", **FAST)
    with pytest.raises(ValueError):
        save_library(results, "/this/dir/does/not/exist/out.csv")


# ---------------------------------------------------------------------------
# translation utilities
# ---------------------------------------------------------------------------
def test_find_first_orf():
    assert find_first_orf("GGGATGCCC") == 3
    assert find_first_orf("CCCGGG") == -1


def test_translate_sequence_basic():
    # ATG AAA TAA -> M K *
    assert translate_sequence("ATGAAATAA", return_stop_codon=True) == "MK*"
    assert translate_sequence("ATGAAATAA", return_stop_codon=False) == "MK"


def test_translate_sequence_no_start():
    assert translate_sequence("CCCGGG") == ""


def test_translate_returns_nucleotides():
    nt = translate_sequence(
        "ATGAAATAA", return_stop_codon=False, return_nucleotide_sequence=True
    )
    assert nt == "ATGAAA"


# ---------------------------------------------------------------------------
# GC content + adaptation score
# ---------------------------------------------------------------------------
def test_calculate_gc_content():
    assert _calculate_gc_content("GGCC") == 1.0
    assert _calculate_gc_content("ATAT") == 0.0
    assert _calculate_gc_content("ATGC") == 0.5
    assert _calculate_gc_content("") == 0


def test_codon_adaptation_score_range():
    # all-optimal codons should score near 1 for yeast
    score = calculate_codon_adaptation_score("ATGTGG", yeast_codon_usage)
    assert 0.0 <= score <= 1.0


def test_codon_adaptation_score_bad_length_raises():
    with pytest.raises(ValueError):
        calculate_codon_adaptation_score("ATGA", yeast_codon_usage)


# ---------------------------------------------------------------------------
# random protein ids
# ---------------------------------------------------------------------------
def test_generate_random_protein_ids_unique():
    ids = generate_random_protein_ids(8, 10)
    assert len(ids) == 10
    assert len(set(ids)) == 10
    assert all(len(i) == 8 for i in ids)


def test_generate_random_protein_ids_impossible_raises():
    with pytest.raises(ValueError):
        generate_random_protein_ids(1, 1000, alphanumeric=False)  # only 10 digits


# ---------------------------------------------------------------------------
# make_codon_table
# ---------------------------------------------------------------------------
def test_make_codon_table_frequencies_sum_to_one():
    # M (ATG) appears twice, K split AAA/AAG
    table = make_codon_table("ATGATGAAAAAG")
    assert table["M"]["ATG"] == 1.0
    assert pytest.approx(sum(table["K"].values())) == 1.0
    assert table["K"]["AAA"] == 0.5
    assert table["K"]["AAG"] == 0.5


# ---------------------------------------------------------------------------
# nucleotide similarity
# ---------------------------------------------------------------------------
def test_nucleotide_similarity_identical():
    max_sim, idx = check_nucleotide_percent_similarity(["ATGC", "ATGC"])
    assert max_sim[0] == pytest.approx(1.0)


def test_nucleotide_similarity_length_mismatch_raises():
    with pytest.raises(ValueError):
        check_nucleotide_percent_similarity(["ATGC", "ATG"])


def test_nucleotide_similarity_invalid_nt_raises():
    with pytest.raises(ValueError):
        check_nucleotide_percent_similarity(["ATGX", "ATGC"])


def test_nucleotide_similarity_matrix_shape():
    matrix = check_nucleotide_percent_similarity(
        ["ATGC", "ATGG", "TTTT"], return_max=False
    )
    assert matrix.shape == (3, 3)


# ---------------------------------------------------------------------------
# CodonTable
# ---------------------------------------------------------------------------
def test_codon_table_construction():
    ct = CodonTable(yeast_codon_usage, 1, 1, (0.4, 0.6))
    assert ct.target_gc == pytest.approx(0.5)
    # every amino acid maps to at least one codon
    assert all(len(codons) > 0 for codons in ct.aa_to_codons.values())
    # weights normalised per amino acid
    for _aa, (codons, weights) in ct.aa_weights.items():
        assert pytest.approx(weights.sum()) == 1.0


def test_codon_table_high_min_probability_raises():
    with pytest.raises(ValueError):
        CodonTable(yeast_codon_usage, 1, 1, (0.4, 0.6), minimum_codon_probability=1.1)


def test_all_codon_tables_keys():
    assert set(all_codon_tables) == {"yeast", "s288c", "s288c_unweighted"}


# ---------------------------------------------------------------------------
# Sequence object
# ---------------------------------------------------------------------------
def _yeast_codon_table():
    return CodonTable(yeast_codon_usage, 1, 0, (0, 1))


def test_sequence_add_coding_sequence_and_getters():
    ct = _yeast_codon_table()
    s = Sequence("MK", ct, "p1", "", "", True, True, None, None)
    assert s.protein_sequence == "MK*"
    s.add_coding_sequence("ATGAAATAA")  # M K *
    assert s.get_coding_sequence() == "ATGAAATAA"
    assert s.get_dna_sequence() == "ATGAAATAA"
    assert s.correct_coding_translation is True
    assert s.correct_full_translation is True
    assert s.get_gc_content() == pytest.approx(s.get_gc_content(full_sequence=True))
    assert "Protein sequence: MK*" in str(s)


def test_sequence_rejects_non_divisible_coding():
    ct = _yeast_codon_table()
    s = Sequence("MK", ct, "p1", "", "", True, True, None, None)
    with pytest.raises(ValueError):
        s.add_coding_sequence("ATGAA")  # length not divisible by 3


def test_sequence_rejects_mismatched_provided_coding():
    ct = _yeast_codon_table()
    with pytest.raises(ValueError):
        Sequence(
            "MK", ct, "p1", "", "", True, True, None, None,
            coding_sequence="ATGTGGTAA",  # translates to MW*, not MK*
        )


# ---------------------------------------------------------------------------
# fasta + codon usage from file
# ---------------------------------------------------------------------------
def _write_fasta(tmp_path):
    fa = tmp_path / "seqs.fasta"
    fa.write_text(
        ">seq1\nATGAAATAA\n>seq2\nATG\nAAA\nAAG\n"  # wrapped lines for seq2
    )
    return fa


def test_read_fasta(tmp_path):
    fa = _write_fasta(tmp_path)
    seqs = read_fasta(str(fa))
    assert seqs == ["ATGAAATAA", "ATGAAAAAG"]


def test_codon_from_fasta(tmp_path):
    fa = _write_fasta(tmp_path)
    table = codon_from_fasta(str(fa))
    assert table["M"]["ATG"] == 1.0
    assert pytest.approx(sum(table["K"].values())) == 1.0


def test_codon_from_fasta_missing_file_returns_none():
    assert codon_from_fasta("/no/such/file.fasta") is None


# ---------------------------------------------------------------------------
# restriction enzyme utilities
# ---------------------------------------------------------------------------
def test_parse_restriction_site_simple():
    info = parse_restriction_site("G/AATTC")
    assert info["recognition_site"] == "GAATTC"
    assert info["top_cut"] == 1
    assert info["pattern"] == "GAATTC"


def test_parse_restriction_site_with_positions():
    info = parse_restriction_site("GAAGAC(2/6)")
    assert info["recognition_site"] == "GAAGAC"
    assert info["top_cut"] == 2
    assert info["bottom_cut"] == 6


def test_site_to_regex_degenerate():
    assert site_to_regex("GGNCC") == "GG[ACGT]CC"


def test_identify_restriction_sites_finds_ecori():
    results = identify_restriction_sites("AAAAGAATTCAAAA")
    assert "EcoRI" in results
    assert results["EcoRI"][0]["start"] == 4


def test_check_restriction_sites_filters():
    seq = "AAAAGAATTCAAAA"
    results = check_restriction_sites(seq, enzymes=["EcoRI"])
    assert list(results.keys()) == ["EcoRI"]
    # an enzyme that does not cut returns nothing
    assert check_restriction_sites(seq, enzymes=["NotI"]) == {}


# ---------------------------------------------------------------------------
# plotting (smoke tests - just ensure a figure is produced)
# ---------------------------------------------------------------------------
def test_graph_gc_content_returns_figure():
    from matplotlib.figure import Figure
    from kea.backend.kea_plotting_utils import graph_gc_content

    fig = graph_gc_content(["ATGCGC", "ATATAT", "GGGCCC"])
    assert isinstance(fig, Figure)


def test_graph_codon_usage_returns_figure():
    from matplotlib.figure import Figure
    from kea.backend.kea_plotting_utils import graph_codon_usage

    fig = graph_codon_usage(["ATGAAAAAG", "ATGAAATGG"], yeast_codon_usage)
    assert isinstance(fig, Figure)
