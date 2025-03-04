"""
Codon optimization class implementation with improved structure and performance.
"""

from typing import Dict, List, Tuple, Optional, Union, Set
import random
import numpy as np
from tqdm import tqdm
from kea.data.aa_codon_conversions import aa_to_codons, codons_to_aa
from kea.backend.optimize_codon_usage import (
    _calculate_gc_content, _score_sequence, _adjust_gc_content,
    _improve_sequence, _build_sequence, calculate_codon_adaptation_score,
    optimize_codon_usage  # Add this import
)

class CodonOptimizer:
    """Class for optimizing codon usage with configurable parameters."""
    
    def __init__(self, codon_table_obj):
        """
        Initialize the optimizer with a codon table object.
        
        Parameters:
        -----------
        codon_table_obj : object
            An object containing codon usage frequencies and GC preferences
        """
        self.codon_table_obj = codon_table_obj
        self.validate_codon_table()
        
    def validate_codon_table(self):
        """Validate that the codon table object has the required attributes."""
        required_attrs = ['codon_table', 'gc_range', 'target_gc', 'aa_weights', 'gc_weight']
        missing = [attr for attr in required_attrs if not hasattr(self.codon_table_obj, attr)]
        if missing:
            raise ValueError(f"Codon table object missing required attributes: {', '.join(missing)}")
            
    def optimize(self, amino_acid_sequence: str, 
                 n_iter: int = 10000,
                 fine_tuning_iterations: int = 500,
                 return_best: bool = True,
                 early_stop_threshold: float = 0.95,
                 show_progress_bar: bool = True) -> str:
        """
        Optimize codon usage for a given amino acid sequence.
        
        Parameters and return value same as optimize_codon_usage function.
        """
        # Call the original function with our parameters
        return optimize_codon_usage(
            amino_acid_sequence,
            self.codon_table_obj,
            n_iter=n_iter,
            fine_tuning_iterations=fine_tuning_iterations,
            return_best=return_best,
            early_stop_threshold=early_stop_threshold,
            show_progress_bar=show_progress_bar
        )
        
    def calculate_adaptation_score(self, sequence: str) -> float:
        """Calculate codon adaptation score for a sequence."""
        return calculate_codon_adaptation_score(sequence, self.codon_table_obj.codon_table)
