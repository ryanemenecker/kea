'''
Result container for build_library.

A library build over thousands of sequences should not lose 9,999 good
constructs because the 10,000th protein could not satisfy a constraint. These
types let build_library keep going and report exactly what failed and why, while
still behaving like the plain list it used to return.
'''

from dataclasses import dataclass, field
from typing import Any, List, Optional


@dataclass
class QualityCheck:
    '''
    One target measured on a finished sequence.

    severity is 'error' for anything that must hold before a sequence is accepted,
    and 'warning' for a target that was requested but missed in a way the caller
    explicitly allowed (asking for the best available GC with return_best=True,
    for instance).
    '''
    name: str
    passed: bool
    value: Optional[float]
    target: str
    severity: str = "error"

    def __str__(self) -> str:
        mark = "ok  " if self.passed else ("FAIL" if self.severity == "error" else "warn")
        shown = "n/a" if self.value is None else f"{self.value:.4f}"
        return f"[{mark}] {self.name}: {shown} (target: {self.target})"


@dataclass
class SequenceQualityReport:
    '''
    Every target checked together on the finished sequence.

    Individual stages each enforce their own piece, so it was previously possible
    for a sequence to be handed back having quietly missed a requested GC range
    with nothing recording that. This is the single place where all the targets
    are measured on the final construct.
    '''
    checks: List[QualityCheck] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        '''True when no error-severity check failed.'''
        return not self.errors

    @property
    def errors(self) -> List[QualityCheck]:
        return [check for check in self.checks
                if not check.passed and check.severity == "error"]

    @property
    def warnings(self) -> List[QualityCheck]:
        return [check for check in self.checks
                if not check.passed and check.severity == "warning"]

    def value(self, name: str) -> Optional[float]:
        for check in self.checks:
            if check.name == name:
                return check.value
        return None

    def __str__(self) -> str:
        return "\n".join(str(check) for check in self.checks)


class QualityCheckError(ValueError):
    '''Raised when a finished sequence fails a required target.'''

    def __init__(self, message, report=None):
        super().__init__(message)
        self.report = report


@dataclass
class FailedSequence:
    '''
    One input sequence that could not be built.

    Attributes
    ----------
    name : str
        Identifier the sequence was given.
    protein_sequence : str
        The input protein sequence, as supplied.
    error : Exception
        The exception that stopped this sequence.
    stage : str
        Which step failed: 'sequence_setup', 'optimization',
        'constraint_repair', 'sequence_diversity', 'padding', or 'quality_check'.
    remaining_violations : list
        Unsatisfied ConstraintViolation records, when the failure was a
        constraint repair. Empty otherwise.
    '''
    name: str
    protein_sequence: str
    error: Exception
    stage: str
    remaining_violations: List[Any] = field(default_factory=list)

    @property
    def reason(self) -> str:
        '''Short human-readable reason.'''
        return f"{type(self.error).__name__}: {self.error}"

    def __str__(self) -> str:
        return f"{self.name} [{self.stage}] {self.reason}"


def evaluate_sequence(sequence_object,
                      target_gc_range=None,
                      total_length=None,
                      constraint_set=None,
                      codon_usage=None,
                      minimum_codon_adaptation=None,
                      gc_is_required=True):
    '''
    Measure every requested target on a finished sequence.

    Called once, after optimization, constraint repair and padding, so the numbers
    describe the construct that would actually be ordered rather than an
    intermediate state.

    Parameters
    ----------
    sequence_object : Sequence
        The finished sequence.
    target_gc_range : tuple(float, float), optional
        Requested GC range. Checked against both the coding sequence and the full
        construct.
    total_length : int, optional
        Requested total length.
    constraint_set : ConstraintSet, optional
        Constraints to re-verify on padding + coding sequence.
    codon_usage : dict, optional
        Codon frequency table to score codon adaptation against.
    minimum_codon_adaptation : float, optional
        If set, codon adaptation below this fails the sequence.
    gc_is_required : bool, default=True
        If False, a missed GC range is a warning rather than an error. Used when
        the caller passed return_best=True and explicitly accepted best-effort GC.

    Returns
    -------
    SequenceQualityReport
    '''
    # Imported here to keep this module free of import cycles.
    from kea.backend.optimize_codon_usage import calculate_codon_adaptation_score
    from kea.backend.translation_utils import translate_sequence

    checks = []

    translated = translate_sequence(sequence_object.coding_sequence)
    checks.append(QualityCheck(
        name="translation",
        passed=translated == sequence_object.protein_sequence,
        value=None,
        target="coding sequence translates to the target protein",
    ))

    if total_length is not None:
        actual = len(sequence_object.full_dna_sequence or "")
        checks.append(QualityCheck(
            name="total_length",
            passed=actual == total_length,
            value=float(actual),
            target=f"exactly {total_length} nt",
        ))

    if target_gc_range is not None:
        low, high = target_gc_range
        severity = "error" if gc_is_required else "warning"
        for name, value in (("gc_coding", sequence_object.gc_content_coding_sequence),
                            ("gc_full", sequence_object.gc_content_full_sequence)):
            if value is None:
                continue
            checks.append(QualityCheck(
                name=name,
                passed=low <= value <= high,
                value=float(value),
                target=f"{low:.3f}-{high:.3f}",
                severity=severity,
            ))

    if constraint_set is not None and constraint_set.constraints:
        region = ((sequence_object.padding_5_prime or "")
                  + sequence_object.coding_sequence
                  + (sequence_object.padding_3_prime or ""))
        remaining = constraint_set.find_violations(region)
        checks.append(QualityCheck(
            name="sequence_constraints",
            passed=not remaining,
            value=float(len(remaining)),
            target="0 violations across padding and coding sequence",
        ))

    if codon_usage is not None:
        score = calculate_codon_adaptation_score(sequence_object.coding_sequence, codon_usage)
        if minimum_codon_adaptation is None:
            passed, target = True, "reported, no floor requested"
        else:
            passed, target = score >= minimum_codon_adaptation, f">= {minimum_codon_adaptation:.4f}"
        checks.append(QualityCheck(
            name="codon_adaptation",
            passed=passed,
            value=float(score),
            target=target,
        ))

    return SequenceQualityReport(checks=checks)


class LibraryResult(list):
    '''
    The successfully built Sequence objects, with the failures attached.

    Behaves exactly like the list of Sequence objects build_library has always
    returned, so existing code keeps working:

        for seq in build_library(...):
            ...

    Failures are available alongside:

        library = build_library(proteins, "human", on_error="skip")
        print(library.summary())
        for failure in library.failures:
            print(failure.name, failure.reason)
    '''

    def __init__(self, sequences=(), failures=None, requested=None):
        super().__init__(sequences)
        self.failures: List[FailedSequence] = list(failures or [])
        # How many inputs were requested; may exceed len(self) + len(failures)
        # if a build was interrupted.
        self.requested: Optional[int] = requested

    @property
    def n_succeeded(self) -> int:
        return len(self)

    @property
    def n_failed(self) -> int:
        return len(self.failures)

    def failures_by_stage(self) -> dict:
        '''Count of failures grouped by the stage that raised.'''
        counts: dict = {}
        for failure in self.failures:
            counts[failure.stage] = counts.get(failure.stage, 0) + 1
        return counts

    @property
    def with_warnings(self) -> list:
        '''
        Accepted sequences that missed a requested target without failing.

        The usual case is a GC range that could not be reached under
        return_best=True: the sequence is still returned, but it did not hit what
        was asked for and you probably want to know.
        '''
        return [sequence for sequence in self
                if getattr(sequence, "quality_report", None) is not None
                and sequence.quality_report.warnings]

    def summary(self) -> str:
        '''One-paragraph description of what succeeded and what did not.'''
        requested = self.requested if self.requested is not None else self.n_succeeded + self.n_failed
        lines = [f"Built {self.n_succeeded}/{requested} sequences."]
        warned = self.with_warnings
        if warned:
            reasons = sorted({check.name
                              for sequence in warned
                              for check in sequence.quality_report.warnings})
            lines.append(f"{len(warned)} built but missed a requested target "
                         f"({', '.join(reasons)}).")
        if self.failures:
            by_stage = ", ".join(f"{stage}: {count}"
                                 for stage, count in sorted(self.failures_by_stage().items()))
            lines.append(f"{self.n_failed} failed ({by_stage}).")
            for failure in self.failures[:10]:
                lines.append(f"  {failure}")
            if self.n_failed > 10:
                lines.append(f"  ... and {self.n_failed - 10} more")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return (f"<LibraryResult {self.n_succeeded} built, {self.n_failed} failed>")
