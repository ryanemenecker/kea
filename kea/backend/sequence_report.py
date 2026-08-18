'''
Per-sequence diagnostics for a finished library.

Everything here is measured on the finished construct and is reported whether or
not the matching constraint was enabled during the build, so a library built
without splice avoidance can still be screened for splice-like sites afterwards.

None of these numbers are probabilities. MaxEntScan scores are log2 odds ratios
(see score_donor/score_acceptor), and polyadenylation signals are exact hexamer
matches. See the README for how to read them.
'''

from kea.backend.sequence_constraints import (COMMON_HUMAN_POLYADENYLATION_SIGNALS,
                                              MaxEntScanScorer,
                                              PrematurePolyadenylationConstraint)

# Defaults used when a report is generated for a library that did not enable the
# splice constraint, so the counts still mean something specific.
DEFAULT_DONOR_THRESHOLD = 6.0
DEFAULT_ACCEPTOR_THRESHOLD = 7.0

_SCORER = MaxEntScanScorer()


def longest_homopolymer(sequence):
    '''Length of the longest single-base run, and the base.'''
    if not sequence:
        return 0, ''
    best_length, best_base = 1, sequence[0]
    run_length = 1
    for index in range(1, len(sequence)):
        if sequence[index] == sequence[index - 1]:
            run_length += 1
            if run_length > best_length:
                best_length, best_base = run_length, sequence[index]
        else:
            run_length = 1
    return best_length, best_base


def longest_direct_repeat(sequence, minimum=6):
    '''
    Length of the longest substring that occurs more than once (occurrences may
    overlap), and the substring itself.

    If a repeat of length L exists then so does one of length L-1, so the property
    is monotone and can be searched rather than scanned. Real sequences repeat over
    short stretches, so the search gallops up from `minimum` before bisecting --
    starting a plain bisection at n/2 would hash half-length substrings on every
    sequence for nothing.
    '''
    length = len(sequence)
    if length < minimum + 1:
        return 0, ''

    def duplicate_of_length(size):
        if size > length - 1:
            return None
        seen = set()
        for start in range(length - size + 1):
            chunk = sequence[start:start + size]
            if chunk in seen:
                return chunk
            seen.add(chunk)
        return None

    if duplicate_of_length(minimum) is None:
        return 0, ''

    # Gallop: find a length that fails, keeping the successful one as the floor.
    low, best = minimum, duplicate_of_length(minimum)
    high = minimum * 2
    while high <= length - 1:
        found = duplicate_of_length(high)
        if found is None:
            break
        best, low = found, high
        high *= 2
    high = min(high, length - 1)

    # Bisect between the last success and the first failure.
    while low + 1 < high:
        middle = (low + high) // 2
        found = duplicate_of_length(middle)
        if found is None:
            high = middle
        else:
            best, low = found, middle
    return len(best), best


def most_tandem_repeat_copies(sequence, max_unit_length=6):
    '''
    Largest number of consecutive copies of any repeat unit, and the unit.

    Unit lengths start at 2, matching TandemRepeatConstraint. A unit of length 1
    is a homopolymer, which is measured separately by longest_homopolymer --
    counting it here too would report the same feature twice.
    '''
    best_copies, best_unit = 0, ''
    length = len(sequence)
    for unit_length in range(2, max_unit_length + 1):
        start = 0
        while start + unit_length <= length:
            unit = sequence[start:start + unit_length]
            copies = 1
            while (sequence[start + copies * unit_length:
                            start + (copies + 1) * unit_length] == unit):
                copies += 1
            if copies > best_copies:
                best_copies, best_unit = copies, unit
            start += 1
    return best_copies, best_unit


def splice_site_summary(sequence, donor_threshold=DEFAULT_DONOR_THRESHOLD,
                        acceptor_threshold=DEFAULT_ACCEPTOR_THRESHOLD):
    '''
    Strongest predicted splice donor and acceptor, plus counts over threshold.

    Scores are MaxEntScan log2 odds ratios on the forward strand. Positions are
    zero-based offsets of the scoring window within `sequence`.
    '''
    summary = {}
    for label, scores, width, threshold, junction_offset in (
            ("donor", _SCORER.score_donors(sequence), 9, donor_threshold, 3),
            ("acceptor", _SCORER.score_acceptors(sequence), 23, acceptor_threshold, 20)):
        if len(scores) == 0:
            summary.update({
                f"{label}_max_score": '', f"{label}_max_position": '',
                f"{label}_max_window": '', f"{label}_sites_over_threshold": 0,
            })
            continue
        best = int(scores.argmax())
        summary.update({
            f"{label}_max_score": round(float(scores[best]), 2),
            f"{label}_max_position": best + junction_offset,
            f"{label}_max_window": sequence[best:best + width],
            f"{label}_sites_over_threshold": int((scores >= threshold).sum()),
        })
    return summary


def polyadenylation_summary(sequence):
    '''
    Every supported polyadenylation hexamer in the sequence.

    Reports both the plain motif count and how many have a downstream U-rich or
    GU-rich element, which is the context that makes a hexamer more likely to be
    used as a real cleavage signal.
    '''
    plain = PrematurePolyadenylationConstraint()
    with_context = PrematurePolyadenylationConstraint(require_downstream_element=True)
    hits = plain.find_violations(sequence)
    supported = with_context.find_violations(sequence)
    return {
        "polya_signal_count": len(hits),
        "polya_with_downstream_element": len(supported),
        "polya_signals": "; ".join(f"{item.sequence}@{item.start}" for item in hits[:10]),
        "polya_canonical_aataaa": sum(1 for item in hits if item.sequence == "AATAAA"),
    }


def annotate_sequence(sequence_object, donor_threshold=None, acceptor_threshold=None):
    '''
    Full diagnostic record for one finished Sequence.

    Splice and polyadenylation screens run on the complete construct
    (adapters + padding + coding sequence), because that is what gets
    synthesized. Complexity measures are reported for the construct too.

    Returns
    -------
    dict
        Column name -> value, ready to write as one CSV row.
    '''
    coding = sequence_object.coding_sequence or ''
    full = sequence_object.full_dna_sequence or coding

    constraint_set = getattr(sequence_object, 'constraint_set', None)
    if donor_threshold is None or acceptor_threshold is None:
        donor, acceptor = DEFAULT_DONOR_THRESHOLD, DEFAULT_ACCEPTOR_THRESHOLD
        for constraint in getattr(constraint_set, 'constraints', ()):  # pragma: no branch
            if hasattr(constraint, 'donor_threshold'):
                donor = constraint.donor_threshold
                acceptor = constraint.acceptor_threshold
        donor_threshold = donor if donor_threshold is None else donor_threshold
        acceptor_threshold = acceptor if acceptor_threshold is None else acceptor_threshold

    report = getattr(sequence_object, 'quality_report', None)
    homopolymer_length, homopolymer_base = longest_homopolymer(full)
    repeat_length, repeat_unit = longest_direct_repeat(full)
    tandem_copies, tandem_unit = most_tandem_repeat_copies(full)

    record = {
        "protein_name": sequence_object.name,
        "protein_sequence": sequence_object.protein_sequence,
        "protein_length": len(sequence_object.protein_sequence or ''),
        "full_dna_sequence": full,
        "coding_sequence": coding,
        "full_length_nt": len(full),
        "coding_length_nt": len(coding),
        "adapter_5_prime": sequence_object.adapter_5_prime,
        "adapter_3_prime": sequence_object.adapter_3_prime,
        "padding_5_prime_nt": len(sequence_object.padding_5_prime or ''),
        "padding_3_prime_nt": len(sequence_object.padding_3_prime or ''),
        "gc_full": round(sequence_object.gc_content_full_sequence, 4)
                   if sequence_object.gc_content_full_sequence is not None else '',
        "gc_coding": round(sequence_object.gc_content_coding_sequence, 4)
                     if sequence_object.gc_content_coding_sequence is not None else '',
        "correct_coding_translation": sequence_object.correct_coding_translation,
        "correct_full_translation": sequence_object.correct_full_translation,
        "splice_donor_threshold": donor_threshold,
        "splice_acceptor_threshold": acceptor_threshold,
        "max_homopolymer_run": homopolymer_length,
        "max_homopolymer_base": homopolymer_base,
        "longest_direct_repeat_nt": repeat_length,
        "longest_direct_repeat": repeat_unit,
        "max_tandem_repeat_copies": tandem_copies,
        "max_tandem_repeat_unit": tandem_unit,
        "constraint_violations_remaining":
            len(getattr(sequence_object, 'sequence_constraint_violations', []) or []),
    }
    record.update(splice_site_summary(full, donor_threshold, acceptor_threshold))
    record.update(polyadenylation_summary(full))

    if report is not None:
        record["quality_check_passed"] = report.passed
        record["quality_warnings"] = "; ".join(check.name for check in report.warnings)
        record["quality_failures"] = "; ".join(check.name for check in report.errors)
        record["codon_adaptation_index"] = (round(report.value("codon_adaptation"), 4)
                                            if report.value("codon_adaptation") is not None else '')
    else:
        record["quality_check_passed"] = ''
        record["quality_warnings"] = ''
        record["quality_failures"] = ''
        record["codon_adaptation_index"] = ''
    return record


# Column order for the CSV. Identity first, then sequences, then the measurements
# a user actually screens on.
REPORT_COLUMNS = [
    "protein_name", "protein_sequence", "protein_length",
    "full_dna_sequence", "coding_sequence",
    "full_length_nt", "coding_length_nt",
    "adapter_5_prime", "adapter_3_prime",
    "padding_5_prime_nt", "padding_3_prime_nt",
    "gc_full", "gc_coding", "codon_adaptation_index",
    "correct_coding_translation", "correct_full_translation",
    "quality_check_passed", "quality_warnings", "quality_failures",
    "donor_max_score", "donor_max_position", "donor_max_window",
    "donor_sites_over_threshold", "splice_donor_threshold",
    "acceptor_max_score", "acceptor_max_position", "acceptor_max_window",
    "acceptor_sites_over_threshold", "splice_acceptor_threshold",
    "polya_signal_count", "polya_canonical_aataaa",
    "polya_with_downstream_element", "polya_signals",
    "max_homopolymer_run", "max_homopolymer_base",
    "longest_direct_repeat_nt", "longest_direct_repeat",
    "max_tandem_repeat_copies", "max_tandem_repeat_unit",
    "constraint_violations_remaining",
]

# The lean set, matching what save_library wrote before diagnostics existed.
BASIC_COLUMNS = [
    "protein_name", "protein_sequence", "full_dna_sequence", "coding_sequence",
    "gc_full", "gc_coding", "correct_full_translation", "correct_coding_translation",
    "codon_adaptation_index",
]
