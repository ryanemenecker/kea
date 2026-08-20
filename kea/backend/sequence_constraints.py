"""Sequence-level constraints and synonymous repair for coding sequences.

The codon optimizer handles codon usage and GC content.  This module handles
features that depend on neighboring bases, such as splice sites, premature
polyadenylation signals, homopolymers, and repeated sequence.  Constraints are
hard requirements: when enabled, Kea either synonymously repairs every detected
violation or raises a clear error.
"""

from dataclasses import dataclass, field
from functools import lru_cache
import importlib.resources
import itertools
import math
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from kea.data.aa_codon_conversions import codons_to_aa


DNA_ALPHABET = frozenset("ACGT")


class ConstraintRepairError(ValueError):
    """Raised when no synonymous rewrite satisfies every requested constraint.

    Subclasses ValueError so existing `except ValueError` handlers keep working,
    while carrying the unsatisfied violations so callers can report which
    constraint blocked which sequence instead of just a formatted string.
    """

    def __init__(self, message, remaining_violations=None, sequence=None):
        super().__init__(message)
        self.remaining_violations = list(remaining_violations or [])
        self.sequence = sequence


def _normalize_dna(sequence: str, label: str = "sequence") -> str:
    if not isinstance(sequence, str):
        raise ValueError(f"{label} must be a string")
    sequence = sequence.upper()
    invalid = sorted(set(sequence) - DNA_ALPHABET)
    if invalid:
        raise ValueError(
            f"{label} contains non-DNA characters: {', '.join(invalid)}. "
            "Only A, C, G and T are supported."
        )
    return sequence


@dataclass(frozen=True)
class ConstraintViolation:
    """A sequence interval that violates one hard design constraint.

    Coordinates are zero-based half-open offsets in the complete sequence seen
    by the constraint (5' context + coding sequence + 3' context).
    """

    kind: str
    start: int
    end: int
    sequence: str
    score: float = 1.0
    threshold: float = 0.0
    details: Dict[str, object] = field(default_factory=dict)

    @property
    def excess(self) -> float:
        """Amount by which the violation exceeds its configured threshold."""
        return max(0.0, float(self.score) - float(self.threshold))


class SequenceConstraint:
    """Small interface implemented by all sequence constraints."""

    name = "sequence_constraint"

    def find_violations(self, sequence: str) -> List[ConstraintViolation]:
        raise NotImplementedError


class MaxEntScanScorer:
    """Exact Yeo-Burge MaxEntScan donor and acceptor scoring.

    Donors use the published 9-nt window (3 exonic + 6 intronic bases).
    Acceptors use the published 23-nt window (20 intronic + 3 exonic bases).
    The fixed lookup tables are loaded lazily from ``kea.data``.
    """

    _BGD_5 = {"A": 0.27, "C": 0.23, "G": 0.23, "T": 0.27}
    _CONS1_5 = {"A": 0.004, "C": 0.0032, "G": 0.9896, "T": 0.0032}
    _CONS2_5 = {"A": 0.0034, "C": 0.0039, "G": 0.0042, "T": 0.9884}

    _BGD_3 = {"A": 0.27, "C": 0.23, "G": 0.23, "T": 0.27}
    _CONS1_3 = {"A": 0.9903, "C": 0.0032, "G": 0.0034, "T": 0.0030}
    _CONS2_3 = {"A": 0.0027, "C": 0.0037, "G": 0.9905, "T": 0.0030}
    _BASE_TO_INT = {"A": 0, "C": 1, "G": 2, "T": 3}

    @staticmethod
    def _hash_sequence(sequence: str) -> int:
        value = 0
        for base in sequence:
            value = value * 4 + MaxEntScanScorer._BASE_TO_INT[base]
        return value

    @staticmethod
    @lru_cache(maxsize=1)
    def _load_tables() -> Tuple[np.ndarray, Tuple[np.ndarray, ...]]:
        # Use the modern Traversable API where available, with the Python 3.8
        # resource API as a compatibility fallback. Arrays are copied before the
        # handle closes so this also works for non-filesystem importers.
        if hasattr(importlib.resources, "files"):
            handle_context = importlib.resources.files("kea.data").joinpath("maxent_scan.npz").open("rb")
        else:  # pragma: no cover - exercised only on Python 3.8
            handle_context = importlib.resources.open_binary("kea.data", "maxent_scan.npz")
        with handle_context as handle:
            with np.load(handle) as archive:
                donor = archive["score5"].copy()
                acceptor = tuple(archive[f"score3_{index}"].copy() for index in range(9))
        return donor, acceptor

    def score_donor(self, sequence: str) -> float:
        sequence = _normalize_dna(sequence, "MaxEntScan donor window")
        if len(sequence) != 9:
            raise ValueError("MaxEntScan donor windows must be exactly 9 nt")

        donor_table, _ = self._load_tables()
        first, second = sequence[3], sequence[4]
        consensus = (
            self._CONS1_5[first]
            * self._CONS2_5[second]
            / (self._BGD_5[first] * self._BGD_5[second])
        )
        rest = sequence[:3] + sequence[5:]
        return math.log2(consensus * donor_table[self._hash_sequence(rest)])

    # Byte -> base-code lookup. Anything that is not ACGT maps to -1 so callers
    # can detect it; sequences are validated before they reach here.
    _BYTE_TO_CODE = np.full(256, -1, dtype=np.int64)
    for _code, _base in enumerate("ACGT"):
        _BYTE_TO_CODE[ord(_base)] = _code
    del _code, _base

    # Donor: positions 3 and 4 are the consensus GT, the other 7 index score5.
    _DONOR_REST_COLUMNS = np.array([0, 1, 2, 5, 6, 7, 8])
    # Acceptor: positions 18 and 19 are the consensus AG, the other 21 form "rest".
    _ACCEPTOR_REST_COLUMNS = np.array(list(range(18)) + [20, 21, 22])
    # (table index, slice of "rest", divide?) -- mirrors the reference score3 product.
    _ACCEPTOR_TABLE_SLICES = (
        (0, 0, 7, False), (1, 7, 14, False), (2, 14, 21, False),
        (3, 4, 11, False), (4, 11, 18, False),
        (5, 4, 7, True), (6, 7, 11, True), (7, 11, 14, True), (8, 14, 18, True),
    )

    @staticmethod
    def _encode(sequence: str) -> np.ndarray:
        """Encode an ACGT string as base codes 0-3."""
        codes = MaxEntScanScorer._BYTE_TO_CODE[
            np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)]
        if codes.size and codes.min() < 0:
            raise ValueError("Sequence contains non-DNA characters; only A, C, G and T "
                             "are supported.")
        return codes

    @staticmethod
    def _hash_columns(windows: np.ndarray, columns) -> np.ndarray:
        """Positional base-4 hash of the given columns, for every window at once."""
        selected = windows[:, columns]
        weights = 4 ** np.arange(selected.shape[1] - 1, -1, -1)
        return selected @ weights

    @classmethod
    def _factor(cls, table: dict, codes: np.ndarray) -> np.ndarray:
        return np.array([table[base] for base in "ACGT"])[codes]

    def score_donors(self, sequence: str) -> np.ndarray:
        """Score every 9-nt donor window in one vectorized pass.

        Returns an array of length max(0, len(sequence) - 8); entry i is the score
        of sequence[i:i+9]. Mathematically identical to calling score_donor on each
        window, but the per-window Python loop and its 1.2M-call hashing were the
        single largest cost in a constrained build.
        """
        sequence = _normalize_dna(sequence, "MaxEntScan donor sequence")
        if len(sequence) < 9:
            return np.empty(0, dtype=np.float64)
        codes = self._encode(sequence)
        windows = np.lib.stride_tricks.sliding_window_view(codes, 9)

        donor_table, _ = self._load_tables()
        first, second = windows[:, 3], windows[:, 4]
        consensus = (self._factor(self._CONS1_5, first) * self._factor(self._CONS2_5, second)
                     / (self._factor(self._BGD_5, first) * self._factor(self._BGD_5, second)))
        return np.log2(consensus * donor_table[self._hash_columns(windows, self._DONOR_REST_COLUMNS)])

    def score_acceptors(self, sequence: str) -> np.ndarray:
        """Score every 23-nt acceptor window in one vectorized pass.

        Returns an array of length max(0, len(sequence) - 22); entry i is the score
        of sequence[i:i+23]. Identical results to score_acceptor per window.
        """
        sequence = _normalize_dna(sequence, "MaxEntScan acceptor sequence")
        if len(sequence) < 23:
            return np.empty(0, dtype=np.float64)
        codes = self._encode(sequence)
        windows = np.lib.stride_tricks.sliding_window_view(codes, 23)

        _, tables = self._load_tables()
        first, second = windows[:, 18], windows[:, 19]
        consensus = (self._factor(self._CONS1_3, first) * self._factor(self._CONS2_3, second)
                     / (self._factor(self._BGD_3, first) * self._factor(self._BGD_3, second)))

        # Accumulate the table product first and apply the consensus factor last,
        # matching score_acceptor's operation order exactly so both paths agree to
        # the last bit rather than to within float rounding.
        rest = windows[:, self._ACCEPTOR_REST_COLUMNS]
        score = np.ones(windows.shape[0], dtype=np.float64)
        for table_index, start, stop, divide in self._ACCEPTOR_TABLE_SLICES:
            values = tables[table_index][self._hash_columns(rest, np.arange(start, stop))]
            if divide:
                score /= values
            else:
                score *= values
        return np.log2(consensus * score)

    def score_acceptor(self, sequence: str) -> float:
        sequence = _normalize_dna(sequence, "MaxEntScan acceptor window")
        if len(sequence) != 23:
            raise ValueError("MaxEntScan acceptor windows must be exactly 23 nt")

        _, tables = self._load_tables()
        first, second = sequence[18], sequence[19]
        consensus = (
            self._CONS1_3[first]
            * self._CONS2_3[second]
            / (self._BGD_3[first] * self._BGD_3[second])
        )
        rest = sequence[:18] + sequence[20:]
        score = 1.0
        score *= tables[0][self._hash_sequence(rest[:7])]
        score *= tables[1][self._hash_sequence(rest[7:14])]
        score *= tables[2][self._hash_sequence(rest[14:])]
        score *= tables[3][self._hash_sequence(rest[4:11])]
        score *= tables[4][self._hash_sequence(rest[11:18])]
        score /= tables[5][self._hash_sequence(rest[4:7])]
        score /= tables[6][self._hash_sequence(rest[7:11])]
        score /= tables[7][self._hash_sequence(rest[11:14])]
        score /= tables[8][self._hash_sequence(rest[14:18])]
        return math.log2(consensus * score)


class HumanSpliceConstraint(SequenceConstraint):
    """Reject high-scoring human splice donors and acceptors."""

    name = "human_splice_site"

    def __init__(self, donor_threshold: float = 6.0, acceptor_threshold: float = 7.0):
        if isinstance(donor_threshold, bool) or not isinstance(donor_threshold, (int, float)):
            raise ValueError("donor_threshold must be numeric")
        if isinstance(acceptor_threshold, bool) or not isinstance(acceptor_threshold, (int, float)):
            raise ValueError("acceptor_threshold must be numeric")
        self.donor_threshold = float(donor_threshold)
        self.acceptor_threshold = float(acceptor_threshold)
        if not math.isfinite(self.donor_threshold):
            raise ValueError("donor_threshold must be finite")
        if not math.isfinite(self.acceptor_threshold):
            raise ValueError("acceptor_threshold must be finite")
        self.scorer = MaxEntScanScorer()

    def find_violations(self, sequence: str) -> List[ConstraintViolation]:
        sequence = _normalize_dna(sequence)
        violations = []

        # Score every window in two vectorized passes, then build violation records
        # only for the handful that cross a threshold. Scoring windows one at a time
        # in Python was the largest single cost in a constrained build.
        donor_scores = self.scorer.score_donors(sequence)
        for start in np.flatnonzero(donor_scores >= self.donor_threshold):
            start = int(start)
            violations.append(
                ConstraintViolation(
                    kind="human_splice_donor",
                    start=start,
                    end=start + 9,
                    sequence=sequence[start:start + 9],
                    score=float(donor_scores[start]),
                    threshold=self.donor_threshold,
                    details={"junction": start + 3},
                )
            )

        acceptor_scores = self.scorer.score_acceptors(sequence)
        for start in np.flatnonzero(acceptor_scores >= self.acceptor_threshold):
            start = int(start)
            violations.append(
                ConstraintViolation(
                    kind="human_splice_acceptor",
                    start=start,
                    end=start + 23,
                    sequence=sequence[start:start + 23],
                    score=float(acceptor_scores[start]),
                    threshold=self.acceptor_threshold,
                    details={"junction": start + 20},
                )
            )

        return violations


# The 13 well-established human PAS hexamers used in PolyA_DB-style analyses.
COMMON_HUMAN_POLYADENYLATION_SIGNALS = (
    "AATAAA",
    "ATTAAA",
    "TATAAA",
    "AGTAAA",
    "AAGAAA",
    "AATATA",
    "AATACA",
    "CATAAA",
    "GATAAA",
    "AATGAA",
    "TTTAAA",
    "ACTAAA",
    "AATAGA",
)


class PrematurePolyadenylationConstraint(SequenceConstraint):
    """Detect human PAS hexamers and annotate downstream U/GU-rich context.

    The canonical AATAAA and major ATTAAA signal are always considered strong.
    All configured PAS variants are hard constraints by default.  Setting
    ``require_downstream_element=True`` reduces false positives by requiring a
    U-rich or GU-rich element within 60 nt downstream of the PAS.
    """

    name = "premature_polyadenylation"
    _STRONG_SIGNALS = frozenset(("AATAAA", "ATTAAA"))
    _GU_RICH_PENTAMERS = frozenset(("GTTGT", "TGTGT", "GTGTT"))

    def __init__(
        self,
        signals: Sequence[str] = COMMON_HUMAN_POLYADENYLATION_SIGNALS,
        require_downstream_element: bool = False,
        downstream_window: int = 60,
    ):
        normalized = tuple(_normalize_dna(signal, "polyadenylation signal") for signal in signals)
        if not normalized or any(len(signal) != 6 for signal in normalized):
            raise ValueError("polyadenylation signals must be non-empty 6-nt DNA sequences")
        if not isinstance(require_downstream_element, bool):
            raise ValueError("require_downstream_element must be True or False")
        if not isinstance(downstream_window, int) or downstream_window < 5:
            raise ValueError("downstream_window must be an integer of at least 5 nt")
        self.signals = normalized
        self.require_downstream_element = require_downstream_element
        self.downstream_window = downstream_window

    @classmethod
    def _downstream_feature_intervals(
        cls, downstream: str
    ) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
        """Return relative intervals for every supported downstream element."""
        u_rich = []
        gu_rich = []
        for start, window in enumerate(_windows(downstream, 5)):
            interval = (start, start + 5)
            if window.count("T") >= 4:
                u_rich.append(interval)
            if window in cls._GU_RICH_PENTAMERS:
                gu_rich.append(interval)
        return u_rich, gu_rich

    @classmethod
    def _downstream_features(cls, downstream: str) -> Tuple[bool, bool]:
        u_rich, gu_rich = cls._downstream_feature_intervals(downstream)
        return bool(u_rich), bool(gu_rich)

    def find_violations(self, sequence: str) -> List[ConstraintViolation]:
        sequence = _normalize_dna(sequence)
        violations = []
        for signal in self.signals:
            search_from = 0
            while True:
                start = sequence.find(signal, search_from)
                if start < 0:
                    break
                end = start + len(signal)
                downstream = sequence[end:end + self.downstream_window]
                u_rich_intervals, gu_rich_intervals = self._downstream_feature_intervals(
                    downstream
                )
                u_rich = bool(u_rich_intervals)
                gu_rich = bool(gu_rich_intervals)
                if not self.require_downstream_element or u_rich or gu_rich:
                    score = (3.0 if signal in self._STRONG_SIGNALS else 1.0)
                    score += float(u_rich) + float(gu_rich)
                    downstream_intervals = tuple(
                        (end + relative_start, end + relative_end)
                        for relative_start, relative_end in (
                            u_rich_intervals + gu_rich_intervals
                        )
                    )
                    violations.append(
                        ConstraintViolation(
                            kind="premature_polyadenylation_signal",
                            start=start,
                            end=end,
                            sequence=signal,
                            score=score,
                            details={
                                "strong_signal": signal in self._STRONG_SIGNALS,
                                "u_rich_downstream": u_rich,
                                "gu_rich_downstream": gu_rich,
                                "downstream_sequence": downstream,
                                "downstream_element_intervals": downstream_intervals,
                                # In context-required mode, clearing every downstream
                                # element is a valid alternative to changing the PAS
                                # itself. Expose those bases to synonymous repair.
                                "repair_intervals": (
                                    downstream_intervals
                                    if self.require_downstream_element else ()
                                ),
                            },
                        )
                    )
                search_from = start + 1
        return sorted(violations, key=lambda violation: (violation.start, violation.kind))


def _windows(sequence: str, size: int) -> Iterable[str]:
    for start in range(max(0, len(sequence) - size + 1)):
        yield sequence[start:start + size]


class HomopolymerConstraint(SequenceConstraint):
    """Limit runs of one repeated nucleotide."""

    name = "homopolymer"

    def __init__(self, max_run_length: int):
        if not isinstance(max_run_length, int) or isinstance(max_run_length, bool) or max_run_length < 1:
            raise ValueError("max_homopolymer_length must be a positive integer")
        self.max_run_length = max_run_length

    def find_violations(self, sequence: str) -> List[ConstraintViolation]:
        sequence = _normalize_dna(sequence)
        violations = []
        start = 0
        while start < len(sequence):
            end = start + 1
            while end < len(sequence) and sequence[end] == sequence[start]:
                end += 1
            run_length = end - start
            if run_length > self.max_run_length:
                violations.append(
                    ConstraintViolation(
                        kind="homopolymer",
                        start=start,
                        end=end,
                        sequence=sequence[start:end],
                        score=float(run_length),
                        threshold=float(self.max_run_length),
                        details={"base": sequence[start], "run_length": run_length},
                    )
                )
            start = end
        return violations


class TandemRepeatConstraint(SequenceConstraint):
    """Limit consecutive copies of short repeat units."""

    name = "tandem_repeat"

    def __init__(self, max_copies: int, max_unit_length: int = 6):
        if not isinstance(max_copies, int) or isinstance(max_copies, bool) or max_copies < 1:
            raise ValueError("max_tandem_repeat_copies must be a positive integer")
        if not isinstance(max_unit_length, int) or max_unit_length < 2:
            raise ValueError("max_tandem_repeat_unit_length must be an integer of at least 2")
        self.max_copies = max_copies
        self.max_unit_length = max_unit_length

    def find_violations(self, sequence: str) -> List[ConstraintViolation]:
        sequence = _normalize_dna(sequence)
        violations = []
        start = 0
        while start < len(sequence):
            found = None
            for unit_length in range(2, self.max_unit_length + 1):
                if start + unit_length * (self.max_copies + 1) > len(sequence):
                    continue
                unit = sequence[start:start + unit_length]
                copies = 1
                while sequence[
                    start + copies * unit_length:start + (copies + 1) * unit_length
                ] == unit:
                    copies += 1
                if copies > self.max_copies:
                    found = (unit, unit_length, copies)
                    break
            if found is None:
                start += 1
                continue
            unit, unit_length, copies = found
            end = start + unit_length * copies
            violations.append(
                ConstraintViolation(
                    kind="tandem_repeat",
                    start=start,
                    end=end,
                    sequence=sequence[start:end],
                    score=float(copies),
                    threshold=float(self.max_copies),
                    details={"repeat_unit": unit, "unit_length": unit_length, "copies": copies},
                )
            )
            start = end
        return violations


class DirectRepeatConstraint(SequenceConstraint):
    """Limit exact repeated substrings, including separated direct repeats."""

    name = "direct_repeat"

    def __init__(self, max_repeat_length: int):
        if not isinstance(max_repeat_length, int) or isinstance(max_repeat_length, bool) or max_repeat_length < 2:
            raise ValueError("max_direct_repeat_length must be an integer of at least 2")
        self.max_repeat_length = max_repeat_length

    def find_violations(self, sequence: str) -> List[ConstraintViolation]:
        sequence = _normalize_dna(sequence)
        window_size = self.max_repeat_length + 1
        # Overlapping occurrences count too (for example, a long tandem repeat),
        # so two full non-overlapping windows are not required.
        if len(sequence) <= window_size:
            return []

        first_occurrence = {}
        violations = []
        covered_until = -1
        for start in range(len(sequence) - window_size + 1):
            word = sequence[start:start + window_size]
            if word not in first_occurrence:
                first_occurrence[word] = start
                continue
            first = first_occurrence[word]
            if start < covered_until:
                continue
            repeat_length = window_size
            while (
                first + repeat_length < len(sequence)
                and start + repeat_length < len(sequence)
                and sequence[first + repeat_length] == sequence[start + repeat_length]
            ):
                repeat_length += 1
            end = start + repeat_length
            violations.append(
                ConstraintViolation(
                    kind="direct_repeat",
                    start=start,
                    end=end,
                    sequence=sequence[start:end],
                    score=float(repeat_length),
                    threshold=float(self.max_repeat_length),
                    details={"first_start": first, "repeat_length": repeat_length},
                )
            )
            covered_until = end
        return violations


class ConstraintSet:
    """Apply constraints to a CDS together with fixed transcribed context."""

    def __init__(
        self,
        constraints: Sequence[SequenceConstraint],
        five_prime_context: str = "",
        three_prime_context: str = "",
    ):
        self.constraints = tuple(constraints)
        self.five_prime_context = _normalize_dna(five_prime_context, "5' transcribed context")
        self.three_prime_context = _normalize_dna(three_prime_context, "3' transcribed context")

    @property
    def coding_offset(self) -> int:
        return len(self.five_prime_context)

    def complete_sequence(self, coding_sequence: str) -> str:
        return self.five_prime_context + coding_sequence + self.three_prime_context

    def find_violations(self, coding_sequence: str) -> List[ConstraintViolation]:
        complete = self.complete_sequence(_normalize_dna(coding_sequence, "coding sequence"))
        violations = []
        for constraint in self.constraints:
            violations.extend(constraint.find_violations(complete))
        return sorted(violations, key=lambda item: (item.start, item.kind, -item.excess))


def make_sequence_constraints(
    avoid_human_splice_sites: bool = False,
    maxent_donor_threshold: float = 6.0,
    maxent_acceptor_threshold: float = 7.0,
    avoid_premature_polyadenylation: bool = False,
    polyadenylation_requires_downstream_element: bool = False,
    max_homopolymer_length: Optional[int] = None,
    max_tandem_repeat_copies: Optional[int] = None,
    max_tandem_repeat_unit_length: int = 6,
    max_direct_repeat_length: Optional[int] = None,
    transcribed_5_prime_context: str = "",
    transcribed_3_prime_context: str = "",
) -> ConstraintSet:
    """Build the public Kea constraint collection from user-facing options."""
    constraints = []
    if avoid_human_splice_sites:
        constraints.append(HumanSpliceConstraint(maxent_donor_threshold, maxent_acceptor_threshold))
    if avoid_premature_polyadenylation:
        constraints.append(
            PrematurePolyadenylationConstraint(
                require_downstream_element=polyadenylation_requires_downstream_element
            )
        )
    if max_homopolymer_length is not None:
        constraints.append(HomopolymerConstraint(max_homopolymer_length))
    if max_tandem_repeat_copies is not None:
        constraints.append(
            TandemRepeatConstraint(max_tandem_repeat_copies, max_tandem_repeat_unit_length)
        )
    if max_direct_repeat_length is not None:
        constraints.append(DirectRepeatConstraint(max_direct_repeat_length))
    return ConstraintSet(constraints, transcribed_5_prime_context, transcribed_3_prime_context)


def _codon_preference_score(sequence: str, codon_frequency_table: Dict[str, Dict[str, float]]) -> float:
    if not sequence:
        return 0.0
    score = 0.0
    codon_count = len(sequence) // 3
    for start in range(0, len(sequence), 3):
        codon = sequence[start:start + 3]
        amino_acid = codons_to_aa[codon]
        frequencies = codon_frequency_table[amino_acid]
        maximum = max(frequencies.values())
        if maximum > 0:
            score += frequencies.get(codon, 0.0) / maximum
    return score / codon_count


def _gc_content(sequence: str) -> float:
    if not sequence:
        return 0.0
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


def _gc_distance(sequence: str, gc_range: Tuple[float, float]) -> float:
    gc_content = _gc_content(sequence)
    if gc_content < gc_range[0]:
        return gc_range[0] - gc_content
    if gc_content > gc_range[1]:
        return gc_content - gc_range[1]
    return 0.0


def _violation_rank(
    sequence: str,
    violations: Sequence[ConstraintViolation],
    codon_frequency_table: Dict[str, Dict[str, float]],
    gc_range: Tuple[float, float],
) -> Tuple[float, ...]:
    excesses = [violation.excess for violation in violations]
    return (
        float(len(violations)),
        max(excesses, default=0.0),
        sum(excesses),
        _gc_distance(sequence, gc_range),
        -_codon_preference_score(sequence, codon_frequency_table),
    )


def _mutable_codon_positions(
    violation: ConstraintViolation,
    coding_offset: int,
    coding_length: int,
) -> List[int]:
    intervals = [(violation.start, violation.end)]
    # A direct repeat has two disjoint mutable copies. The violation's primary
    # interval is the later copy; include the first copy as well, which matters
    # when one occurrence lies in fixed transcribed context.
    if violation.kind == "direct_repeat" and "first_start" in violation.details:
        first_start = int(violation.details["first_start"])
        repeat_length = int(violation.details["repeat_length"])
        intervals.append((first_start, first_start + repeat_length))

    # Some constraints can be cleared outside their primary reporting interval.
    # For a context-dependent PAS, either the hexamer or all supporting downstream
    # elements may be changed. The intervals are absolute coordinates in the same
    # complete sequence as violation.start/end.
    intervals.extend(violation.details.get("repair_intervals", ()))

    positions = set()
    for interval_start, interval_end in intervals:
        overlap_start = max(interval_start, coding_offset)
        overlap_end = min(interval_end, coding_offset + coding_length)
        if overlap_start >= overlap_end:
            continue
        first = (overlap_start - coding_offset) // 3
        last = (overlap_end - 1 - coding_offset) // 3
        positions.update(range(first, last + 1))
    return sorted(positions)


def is_violation_unsatisfiable(
    coding_sequence: str,
    violation: ConstraintViolation,
    codon_table_obj,
    constraint_set: "ConstraintSet",
    max_combinations: int = 20000,
) -> bool:
    """
    Prove that no synonymous rewrite can clear this violation.

    A fixed-width window (a 9-nt splice donor, say) is covered by only a handful of
    codons, so every synonymous encoding of those codons can be enumerated. If none
    of them removes the violation, the protein simply cannot be encoded under the
    current settings and re-optimizing is wasted work -- the caller should be told
    to relax the threshold rather than to try again.

    Returns False when the search space is too large to enumerate (see
    max_combinations), so a False result means "not proven impossible", not
    "definitely fixable".
    """
    positions = _mutable_codon_positions(
        violation, constraint_set.coding_offset, len(coding_sequence))
    if not positions:
        # Entirely inside fixed context; no synonymous rewrite can reach it.
        return True

    codons = [coding_sequence[index * 3:index * 3 + 3] for index in range(len(coding_sequence) // 3)]

    choices = []
    for index in positions:
        aa = codons_to_aa.get(codons[index])
        if aa is None:
            return False
        choices.append(tuple(codon_table_obj.codon_table[aa]))

    total = 1
    for options in choices:
        total *= len(options)
        if total > max_combinations:
            return False

    for combination in itertools.product(*choices):
        candidate_codons = list(codons)
        for index, replacement in zip(positions, combination):
            candidate_codons[index] = replacement
        candidate = constraint_set.complete_sequence(''.join(candidate_codons))
        # Ask only whether THIS window is still violated. Other violations in the
        # complete sequence belong to other codons and are not this enumeration's
        # business. Reconstructing the complete CDS also keeps disjoint repair
        # intervals (direct repeats and PAS downstream elements) at their real
        # coordinates instead of accidentally closing the bases between them.
        still_violated = any(
            item.kind == violation.kind and item.start == violation.start
            for constraint in constraint_set.constraints
            for item in constraint.find_violations(candidate)
        )
        if not still_violated:
            return False
    return True


# G+C base count per codon, for exact integer GC bookkeeping during compensation.
_CODON_GC = {codon: codon.count('G') + codon.count('C') for codon in codons_to_aa}


def _gc_compensation_index(
    sequence: str,
    protected: Iterable[int],
    allowed_codons: Dict[str, Dict[str, float]],
) -> Dict[int, List[Tuple[float, int, str, str]]]:
    """
    Index of synonymous edits by the GC-count change they produce.

    Built once per repair round from the current sequence. Codons covered by the
    violation being repaired are excluded, so a compensating edit can never undo
    the fix it is paying for. Entries are ordered by the replacement codon's
    frequency, so compensation costs as little codon adaptation as possible.
    """
    protected = set(protected)
    index: Dict[int, List[Tuple[float, int, str, str]]] = {}
    for position in range(len(sequence) // 3):
        if position in protected:
            continue
        start = position * 3
        current = sequence[start:start + 3]
        amino_acid = codons_to_aa.get(current)
        if amino_acid is None:
            continue
        current_gc = _CODON_GC[current]
        for replacement, frequency in allowed_codons[amino_acid].items():
            if replacement == current:
                continue
            delta = _CODON_GC[replacement] - current_gc
            if delta:
                index.setdefault(delta, []).append((frequency, position, current, replacement))
    for entries in index.values():
        entries.sort(reverse=True)
    return index


def _compensate_gc(candidate: str, gc_range: Tuple[float, float], index) -> Optional[str]:
    """
    Pull `candidate` back inside the GC range with one extra synonymous edit.

    The optimizer typically hands repair a sequence pinned against a GC boundary,
    which leaves no room for a fix that shifts GC. Rather than vetoing that fix,
    pay the GC back somewhere else in the sequence. Returns None when no single
    edit can close the gap.
    """
    length = len(candidate)
    gc_count = candidate.count('G') + candidate.count('C')
    lowest = gc_range[0] * length
    highest = gc_range[1] * length
    if lowest <= gc_count <= highest:
        return candidate

    # Deltas that would land the count back inside the range. Codon swaps move it
    # by at most 3 bases, so this is a very small search.
    delta_min = math.ceil(lowest - gc_count - 1e-9)
    delta_max = math.floor(highest - gc_count + 1e-9)
    for delta in range(delta_min, delta_max + 1):
        if delta == 0:
            continue
        for _frequency, position, expected, replacement in index.get(delta, ()):
            start = position * 3
            # The index was built before this candidate's own edit; skip entries
            # whose codon has since changed.
            if candidate[start:start + 3] != expected:
                continue
            return candidate[:start] + replacement + candidate[start + 3:]
    return None


def _single_synonymous_edits(
    sequence: str,
    positions: Iterable[int],
    allowed_codons: Dict[str, Dict[str, float]],
) -> List[Tuple[int, str, str]]:
    edits = []
    for position in sorted(set(positions)):
        start = position * 3
        current = sequence[start:start + 3]
        amino_acid = codons_to_aa[current]
        for replacement in allowed_codons[amino_acid]:
            if replacement != current:
                edits.append((position, current, replacement))
    return edits


def _apply_edit(sequence: str, edit: Tuple[int, str, str]) -> str:
    position, _current, replacement = edit
    start = position * 3
    return sequence[:start] + replacement + sequence[start + 3:]


def repair_coding_sequence(
    coding_sequence: str,
    codon_table_obj,
    constraint_set: ConstraintSet,
    max_iterations: int = 100,
) -> Tuple[str, List[ConstraintViolation]]:
    """Synonymously repair every enabled sequence constraint.

    Single-codon changes are tried first.  If no single edit improves the hard
    constraint rank, all two-codon combinations within the current violation are
    considered.  Codon usage is only used to break ties after hard constraints,
    and an initially satisfied GC range is preserved.
    """
    coding_sequence = _normalize_dna(coding_sequence, "coding sequence")
    if len(coding_sequence) % 3:
        raise ValueError("Coding sequence length must be divisible by 3 for synonymous repair")
    if not isinstance(max_iterations, int) or max_iterations < 1:
        raise ValueError("constraint_repair_attempts must be a positive integer")
    if not constraint_set.constraints:
        return coding_sequence, []

    current = coding_sequence
    gc_range = tuple(codon_table_obj.gc_range)
    codon_table = codon_table_obj.codon_table
    preserve_in_range = _gc_distance(current, gc_range) == 0.0

    for _ in range(max_iterations):
        violations = constraint_set.find_violations(current)
        if not violations:
            return current, []
        current_rank = _violation_rank(current, violations, codon_table, gc_range)

        # Work on the most severe violation first.  Its full scoring window is
        # mutable, not just the GT/AG or PAS core, because flanking bases are part
        # of MaxEntScan and downstream context can be sequence-dependent.
        violation = max(violations, key=lambda item: (item.excess, item.score))
        positions = _mutable_codon_positions(
            violation, constraint_set.coding_offset, len(current)
        )
        if not positions:
            break

        edits = _single_synonymous_edits(current, positions, codon_table)
        best_sequence = None
        best_rank = current_rank

        # Built once per round, excluding the codons under repair.
        compensation_index = (_gc_compensation_index(current, positions, codon_table)
                              if preserve_in_range else None)

        def consider(candidate: str) -> None:
            nonlocal best_sequence, best_rank
            if preserve_in_range and _gc_distance(candidate, gc_range) > 0.0:
                # Do not veto a fix just because it shifts GC -- try to pay the GC
                # back elsewhere first. Any new violation the compensating edit
                # introduces is caught by the rank comparison below.
                candidate = _compensate_gc(candidate, gc_range, compensation_index)
                if candidate is None:
                    return
            candidate_violations = constraint_set.find_violations(candidate)
            rank = _violation_rank(candidate, candidate_violations, codon_table, gc_range)
            if rank < best_rank:
                best_sequence = candidate
                best_rank = rank

        for edit in edits:
            consider(_apply_edit(current, edit))

        if best_sequence is None:
            # A pair of individually neutral or unfavorable changes can remove a
            # long acceptor/repeat window while preserving GC.  Bound the search
            # to codons overlapping this violation so it remains inexpensive.
            for first_index, first_edit in enumerate(edits):
                intermediate = _apply_edit(current, first_edit)
                for second_edit in edits[first_index + 1:]:
                    if second_edit[0] == first_edit[0]:
                        continue
                    consider(_apply_edit(intermediate, second_edit))

        if best_sequence is None:
            break
        current = best_sequence

    # The success check happens at the top of the loop, so a repair that lands its
    # final edit on the last permitted iteration falls through to here with nothing
    # actually wrong. Re-check before failing, otherwise a clean sequence is thrown
    # away with a self-contradicting "Remaining violations:" list that is empty.
    remaining = constraint_set.find_violations(current)
    if not remaining:
        return current, []

    preview = "; ".join(
        f"{item.kind} {item.start}:{item.end} {item.sequence} score={item.score:.2f}"
        for item in remaining[:5]
    )
    raise ConstraintRepairError(
        "Could not synonymously satisfy the requested sequence constraints after "
        f"{max_iterations} repair iterations. Remaining violations: {preview}",
        remaining_violations=remaining,
        sequence=current,
    )
