"""Synonymous sequence diversification with hard Hamming-distance targets."""

import math
from typing import Optional, Sequence

import numpy as np

from kea.data.aa_codon_conversions import codons_to_aa


class SequenceDiversityError(ValueError):
    """Raised when the requested set of sufficiently distinct encodings is not found."""

    def __init__(self, message, minimum_distance, attempts, best_distance=None):
        super().__init__(message)
        self.minimum_distance = minimum_distance
        self.attempts = attempts
        self.best_distance = best_distance


def nucleotide_hamming_distance(first: str, second: str) -> int:
    """Return the number of nucleotide positions at which two sequences differ."""
    if not isinstance(first, str) or not isinstance(second, str):
        raise ValueError("Hamming distance requires two nucleotide strings")
    if len(first) != len(second):
        raise ValueError("Hamming distance requires sequences of equal length")
    return sum(left != right for left, right in zip(first.upper(), second.upper()))


def minimum_hamming_distance(sequence: str, references: Sequence[str]) -> Optional[int]:
    """Minimum Hamming distance from ``sequence`` to a non-empty reference set."""
    if not references:
        return None
    return min(nucleotide_hamming_distance(sequence, reference) for reference in references)


def _distance_vector(codons, reference_codons):
    return [
        sum(
            nucleotide_hamming_distance(codon, reference[position])
            for position, codon in enumerate(codons)
        )
        for reference in reference_codons
    ]


def _capped_distance_rank(distances, target):
    capped = [min(distance, target) for distance in distances]
    return min(capped), sum(capped)


def _distance_after_swap(distances, references, position, current, replacement):
    projected = []
    for distance, reference in zip(distances, references):
        projected.append(
            distance
            - nucleotide_hamming_distance(current, reference[position])
            + nucleotide_hamming_distance(replacement, reference[position])
        )
    return projected


def _gc_gap(gc_count, low_count, high_count):
    if gc_count < low_count:
        return low_count - gc_count
    if gc_count > high_count:
        return gc_count - high_count
    return 0


def _restore_gc_range(codons, references, distances, target, codon_table, gc_range):
    """Restore GC with synonymous single-codon edits that retain every distance."""
    length = len(codons) * 3
    low_count = math.ceil(gc_range[0] * length - 1e-9)
    high_count = math.floor(gc_range[1] * length + 1e-9)
    gc_count = sum(codon.count("G") + codon.count("C") for codon in codons)

    while _gc_gap(gc_count, low_count, high_count):
        current_gap = _gc_gap(gc_count, low_count, high_count)
        best_key = None
        best_moves = []
        for position, current in enumerate(codons):
            amino_acid = codons_to_aa[current]
            current_gc = current.count("G") + current.count("C")
            for replacement, frequency in codon_table[amino_acid].items():
                if replacement == current:
                    continue
                projected = _distance_after_swap(
                    distances, references, position, current, replacement
                )
                if min(projected) < target:
                    continue
                replacement_gc = replacement.count("G") + replacement.count("C")
                projected_gc = gc_count - current_gc + replacement_gc
                gap = _gc_gap(projected_gc, low_count, high_count)
                if gap >= current_gap:
                    continue
                key = (gap, -frequency)
                if best_key is None or key < best_key:
                    best_key = key
                    best_moves = [(position, replacement, projected, projected_gc)]
                elif key == best_key:
                    best_moves.append((position, replacement, projected, projected_gc))

        if not best_moves:
            return None
        position, replacement, distances, gc_count = best_moves[
            int(np.random.randint(len(best_moves)))
        ]
        codons[position] = replacement

    return codons, distances


def diversify_coding_sequence(
    sequence: str,
    codon_table_obj,
    references: Sequence[str],
    minimum_distance: int,
) -> Optional[str]:
    """Create a synonymous encoding separated from every accepted reference.

    The search starts from an optimized coding sequence and changes only as many
    codons as needed to meet the hard distance target. Diversity is optimized
    lexicographically before codon frequency; an already satisfied GC range is
    restored before returning. ``None`` means this randomized search reached a
    local infeasibility, not that no encoding exists.
    """
    if minimum_distance <= 0 or not references:
        return sequence
    if len(sequence) % 3:
        raise ValueError("Coding sequence length must be divisible by 3")
    if any(len(reference) != len(sequence) for reference in references):
        raise ValueError("Diversity references must match the coding-sequence length")

    codons = [sequence[start:start + 3] for start in range(0, len(sequence), 3)]
    reference_codons = [
        [reference[start:start + 3] for start in range(0, len(reference), 3)]
        for reference in references
    ]
    codon_table = codon_table_obj.codon_table
    distances = _distance_vector(codons, reference_codons)
    if min(distances) >= minimum_distance:
        return sequence

    mutable_positions = [
        position for position, codon in enumerate(codons)
        if len(codon_table[codons_to_aa[codon]]) > 1
    ]
    if not mutable_positions:
        return None

    original_gc = sequence.count("G") + sequence.count("C")
    gc_count = original_gc
    max_rounds = max(2, len(mutable_positions) * 2)
    for _ in range(max_rounds):
        changed = False
        for position in np.random.permutation(mutable_positions):
            position = int(position)
            current = codons[position]
            current_rank = _capped_distance_rank(distances, minimum_distance)
            amino_acid = codons_to_aa[current]
            best_rank = current_rank
            best_options = []
            for replacement, frequency in codon_table[amino_acid].items():
                if replacement == current:
                    continue
                projected = _distance_after_swap(
                    distances, reference_codons, position, current, replacement
                )
                rank = _capped_distance_rank(projected, minimum_distance)
                if rank > best_rank:
                    best_rank = rank
                    best_options = [(replacement, projected, frequency)]
                elif rank == best_rank and rank > current_rank:
                    best_options.append((replacement, projected, frequency))

            if not best_options:
                continue

            current_gc = current.count("G") + current.count("C")
            # Among equally useful diversity moves, minimize drift from the optimized
            # sequence's GC count and then sample according to allowed codon usage.
            drifts = [
                abs(gc_count - current_gc
                    + replacement.count("G") + replacement.count("C") - original_gc)
                for replacement, _projected, _frequency in best_options
            ]
            minimum_drift = min(drifts)
            choices = [option for option, drift in zip(best_options, drifts)
                       if drift == minimum_drift]
            weights = np.array([frequency for _replacement, _projected, frequency in choices],
                               dtype=float)
            weights = weights / weights.sum() if weights.sum() > 0 else None
            selected = int(np.random.choice(len(choices), p=weights))
            replacement, distances, _frequency = choices[selected]
            gc_count += (replacement.count("G") + replacement.count("C")
                         - current.count("G") - current.count("C"))
            codons[position] = replacement
            changed = True

            if min(distances) >= minimum_distance:
                break

        if min(distances) >= minimum_distance:
            break
        if not changed:
            return None

    if min(distances) < minimum_distance:
        return None

    gc_range = tuple(codon_table_obj.gc_range)
    original_gc_fraction = original_gc / len(sequence)
    if gc_range[0] <= original_gc_fraction <= gc_range[1]:
        restored = _restore_gc_range(
            codons, reference_codons, distances, minimum_distance,
            codon_table, gc_range,
        )
        if restored is None:
            return None
        codons, _distances = restored

    return "".join(codons)
