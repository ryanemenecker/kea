'''
Class that takes care of the codon table.
'''

import numpy as np

class CodonTable:
    def __init__(self, codon_usage, 
                 usage_weight, gc_weight, gc_range,
                 minimum_codon_probability=0.0):
        '''
        Initialize the CodonTable object with a codon usage dictionary.
        '''
        self.codon_usage = codon_usage
        self.usage_weight = usage_weight
        self.gc_weight = gc_weight
        self.gc_range = gc_range
        self.minimum_codon_probability = minimum_codon_probability
        
        # Calculate everything once during initialization
        self.codon_table = self.make_codon_table()
        self.aa_to_codons = self.get_aa_to_codons()
        self.codons_to_aa = self.get_codons_to_aa()
        self.codon_gc_content = self.get_codon_gc_content()
        self.aa_weights = self.get_aa_weights()
        self.target_gc = self.get_target_gc()

    def get_aa_to_codons(self, minimum_codon_probability=None):
        '''
        Function to get the amino acid to codon mapping.
        '''
        if minimum_codon_probability is None:
            minimum_codon_probability = self.minimum_codon_probability
        aa_to_codons = {}
        for aa in self.codon_usage:
            aa_to_codons[aa] = []
            for codon in self.codon_usage[aa]:
                if self.codon_usage[aa][codon] >= minimum_codon_probability:
                    aa_to_codons[aa].append(codon)
        return aa_to_codons

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
        for aa in self.codon_usage:
            codon_table[aa] = {}
            for codon in self.codon_usage[aa]:
                if self.codon_usage[aa][codon] >= minimum_codon_probability:
                    codon_table[aa][codon] = self.codon_usage[aa][codon]  

        # now make sure that the probabilities sum to 1.
        for aa in codon_table:
            if force_minimum_one_codon and len(codon_table[aa]) == 0:
                raise ValueError(f"No codons found for amino acid {aa}! The minimum_codon_probability is too high.")
            total = sum(codon_table[aa].values())
            if total > 0:
                for codon in codon_table[aa]:
                    codon_table[aa][codon] /= total
            else:
                for codon in codon_table[aa]:
                    codon_table[aa][codon] = 0.0

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
            usage_weights = [self.codon_usage[aa].get(codon, 0.01) for codon in possible_codons]
            
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
            
            # Combine weights with stronger GC emphasis
            combined_weights = np.array([uw**usage_weight * gw**(gc_weight*2) for uw, gw in zip(usage_weights, gc_weights)])
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