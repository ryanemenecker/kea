'''
Class that takes care of the codon table.
'''

import numpy as np

from kea.data.aa_codon_conversions import aa_to_codons as CANONICAL_AA_TO_CODONS


def validate_codon_usage(codon_usage):
    '''
    Check a codon usage table before anything downstream indexes into it.

    build_library accepts a user-supplied dict. Without this check a table that is
    missing an amino acid, or that lists a codon under the wrong one, surfaces much
    later as a bare KeyError from inside the optimizer.

    Raises
    ------
    ValueError
        With a message naming the specific problem.
    '''
    if not isinstance(codon_usage, dict) or not codon_usage:
        raise ValueError("Codon frequency table must be a non-empty dictionary of "
                         "{amino_acid: {codon: frequency}}")

    missing = sorted(set(CANONICAL_AA_TO_CODONS) - set(codon_usage))
    if missing:
        raise ValueError(f"Codon frequency table is missing amino acid(s): {', '.join(missing)}. "
                         "It must cover all 20 amino acids plus '*' (stop).")

    for aa, codons in codon_usage.items():
        if aa not in CANONICAL_AA_TO_CODONS:
            raise ValueError(f"Codon frequency table has an unknown amino acid key: {aa!r}")
        if not isinstance(codons, dict) or not codons:
            raise ValueError(f"Codon frequency table entry for {aa!r} must be a non-empty "
                             "dictionary of {codon: frequency}")
        valid = set(CANONICAL_AA_TO_CODONS[aa])
        for codon, freq in codons.items():
            if codon not in valid:
                raise ValueError(f"Codon {codon!r} is not a valid codon for amino acid {aa!r}. "
                                 f"Expected one of: {', '.join(sorted(valid))}")
            if isinstance(freq, bool) or not isinstance(freq, (int, float)):
                raise ValueError(f"Frequency for codon {codon!r} ({aa!r}) must be a number, "
                                 f"got {type(freq).__name__}")
            if freq < 0:
                raise ValueError(f"Frequency for codon {codon!r} ({aa!r}) must be non-negative")
        if sum(codons.values()) <= 0:
            raise ValueError(f"All codon frequencies for amino acid {aa!r} are zero")


class CodonTable:
    """Codon table with GC content preferences."""
    
    def __init__(self, codon_usage, 
                 usage_weight, gc_weight, gc_range,
                 minimum_codon_probability=0.0, gc_tolerance=0.025,
                 gc_centering=0.0, escape_gc_boundary=True):
        """
        Initialize codon table with weights and GC preferences.
        
        Parameters:
        ...existing parameters...
        gc_tolerance : float
            Allowable deviation from GC range (default 0.025 or 2.5%)
        gc_centering : float
            How strongly to prefer the middle of gc_range over its edges, 0-1.
            0 (the default) makes the GC term flat across the whole range, which
            lets codon usage pin GC against an edge and leaves constraint repair
            no room to move it.
        escape_gc_boundary : bool
            Allow the refinement search to make paired swaps that improve codon
            usage at one position and pay the GC back at another. Without it the
            search stalls at the GC-range edge whenever the codon-usage optimum
            lies outside the range.
        """
        validate_codon_usage(codon_usage)
        # Outside 0-1 the GC objective inverts: a negative value rewards the range
        # edges over its centre, and above 1 makes the out-of-range branch negative.
        if not 0.0 <= gc_centering <= 1.0:
            raise ValueError("gc_centering must be between 0 and 1")
        self.codon_usage = codon_usage
        self.usage_weight = usage_weight
        self.gc_weight = gc_weight
        self.gc_range = gc_range
        self.minimum_codon_probability = minimum_codon_probability
        self.gc_tolerance = gc_tolerance
        self.gc_centering = gc_centering
        self.escape_gc_boundary = escape_gc_boundary
        self.target_gc = sum(gc_range) / 2
        
        # Calculate everything once during initialization
        self.codon_table = self.make_codon_table()
        self.aa_to_codons = self.get_aa_to_codons()
        self.codons_to_aa = self.get_codons_to_aa()
        self.codon_gc_content = self.get_codon_gc_content()
        self.aa_weights = self.get_aa_weights()
        self.target_gc = self.get_target_gc()

    def get_aa_to_codons(self):
        '''
        Function to get the amino acid to codon mapping.

        Derived from self.codon_table so that it always agrees with it. Building it
        separately from the raw codon_usage meant the two could disagree about which
        codons minimum_codon_probability had removed.
        '''
        return {aa: list(codons) for aa, codons in self.codon_table.items()}

    def get_codons_to_aa(self):
        '''
        Function to get the codon to amino acid mapping.
        '''
        codons_to_aa = {}
        for aa in self.codon_usage:
            for codon in self.codon_usage[aa]:    
                codons_to_aa[codon] = aa
        return codons_to_aa

    def make_codon_table(self, minimum_codon_probability=None, force_minimum_one_codon=True):
        '''
        Function to make a codon table from a codon usage dictionary.
        '''
        if minimum_codon_probability is None:
            minimum_codon_probability = self.minimum_codon_probability
        
        # dict to hold codon table.
        codon_table = {}

        # make dictionary where the key is the amino acid and the value
        # is a dictionary holding each codon and its corresponding probabilities.
        for aa, raw in self.codon_usage.items():
            # Normalize BEFORE filtering. minimum_codon_probability is documented as
            # a probability, so it has to be compared against relative frequencies.
            # Filtering the raw values first meant a table given as counts or as
            # per-1000 values kept everything (all values above the threshold) while
            # a table on a small scale lost everything.
            total = sum(raw.values())
            if total > 0:
                normalized = {codon: value / total for codon, value in raw.items()}
            else:
                normalized = {codon: 0.0 for codon in raw}

            kept = {codon: freq for codon, freq in normalized.items()
                    if freq >= minimum_codon_probability}

            if force_minimum_one_codon and not kept:
                raise ValueError(f"No codons found for amino acid {aa}! The minimum_codon_probability is too high.")

            # Renormalize across the surviving codons so they still sum to 1.
            kept_total = sum(kept.values())
            if kept_total > 0:
                kept = {codon: freq / kept_total for codon, freq in kept.items()}

            codon_table[aa] = kept

        return codon_table
    
    def get_codon_gc_content(self):
        '''
        Function to get the GC content of each codon.
        '''
        codon_gc_content = {codon: (codon.count('G') + codon.count('C'))/3
                            for aa in self.aa_to_codons 
                            for codon in self.aa_to_codons[aa]} 
        return codon_gc_content
    
    def get_aa_weights(self, gc_range=None,
                       usage_weight=None, gc_weight=None):
        '''
        Function to get the weights of each amino acid.
        '''
        if gc_range is None:
            gc_range = self.gc_range
        if usage_weight is None:
            usage_weight = self.usage_weight
        if gc_weight is None:
            gc_weight = self.gc_weight

        # Pre-calculate target GC
        target_gc = sum(gc_range) / 2
        
        # Pre-calculate weights for each amino acid
        aa_weights = {}
        for aa in self.aa_to_codons:  
            possible_codons = self.aa_to_codons[aa]
            # Read frequencies from the normalized/filtered table, not the raw one,
            # so a table supplied as counts or per-1000 behaves the same as a
            # normalized one.
            usage_weights = [self.codon_table[aa].get(codon, 0.01) for codon in possible_codons]
            
            # Stronger GC weighting function
            gc_weights = []
            for codon in possible_codons:
                codon_gc = self.codon_gc_content[codon]  
                # Sharper penalty for being outside target range
                if codon_gc < gc_range[0]:
                    gc_factor = 1 - 2 * (gc_range[0] - codon_gc)
                elif codon_gc > gc_range[1]:
                    gc_factor = 1 - 2 * (codon_gc - gc_range[1])
                else:
                    # Bonus for being in range with peak at target
                    gc_factor = 1 + (1 - abs(codon_gc - target_gc) * 2)
                
                gc_weights.append(max(0.01, gc_factor))
            
            # Combine weights 
            combined_weights = np.array([uw**usage_weight * gw**(gc_weight) for uw, gw in zip(usage_weights, gc_weights)])
            # Handle case where all weights are zero
            if combined_weights.sum() > 0:
                aa_weights[aa] = (possible_codons, combined_weights/combined_weights.sum())
            else:
                # Provide equal weights if all are zero
                equal_weights = np.ones(len(possible_codons)) / len(possible_codons)
                aa_weights[aa] = (possible_codons, equal_weights)
                
        return aa_weights
    
    def get_target_gc(self):
        '''
        Function to get the target GC content.
        '''
        return sum(self.gc_range) / 2
    
    def __str__(self):
        '''
        String representation of the Sequence object.
        '''
        parts = []
        parts.append(self.codon_usage)
        parts.append(f"Usage weight: {self.usage_weight}")
        parts.append(f"GC weight: {self.gc_weight}")
        parts.append(f"GC range: {self.gc_range}")
        parts.append(f"Minimum codon probability: {self.minimum_codon_probability}")
        return "\n".join(parts)