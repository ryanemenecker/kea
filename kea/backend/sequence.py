'''
Class that holds sequence information and manages DNA sequence generation from protein sequences.
This includes handling of protein sequences, DNA sequences, adapters, padding, and verification
of sequence properties.
'''
from kea.backend.translation_utils import (find_first_orf, translate_sequence,
                                           translate_open_reading_frame)
from kea.backend.library_generation_utils import add_padding, repair_padding_constraints
from kea.backend.sequence_constraints import ConstraintRepairError

class Sequence:
    def __init__(self, 
                 protein_sequence, 
                 codon_table,
                 name,
                 adapter_3_prime,
                 adapter_5_prime,
                 avoid_adding_start_codons,
                 avoid_adding_stop_codons,
                 total_length,
                 pad_location,
                 force_start_codon=True,
                 force_stop_codon=True,
                 coding_sequence=None,
                 verify_coding_sequence=True,
                 return_best=False,
                 padding_attempts=5000,
                 constraint_set=None):
        '''
        Initialize a Sequence object for DNA sequence generation and management.

        Parameters
        ----------
        protein_sequence : str
            The target protein sequence to be encoded
        codon_table : CodonTable
            CodonTable object containing codon usage information and weights
        name : str
            Identifier for the sequence
        adapter_3_prime : str
            3' adapter sequence to be added to the DNA sequence
        adapter_5_prime : str
            5' adapter sequence to be added to the DNA sequence
        avoid_adding_start_codons : bool
            If True, avoid start codons in padding sequences
        avoid_adding_stop_codons : bool
            If True, avoid stop codons in padding sequences
        total_length : int or None
            Target total length for the final DNA sequence
        pad_location : int or None
            Location for padding (3, 5, or None for both ends)
        force_start_codon : bool, default=True
            If True, ensure protein sequence starts with methionine
        force_stop_codon : bool, default=True
            If True, ensure protein sequence ends with stop codon
        coding_sequence : str, optional
            Pre-defined coding sequence (if None, will be generated)
        verify_coding_sequence : bool, default=True
            If True, verify translation of coding sequence
        return_best : bool, default=False
            If True, return best sequence based on optimization criteria
        padding_attempts : int, default=5000
            Number of attempts for generating padding sequences
        constraint_set : ConstraintSet, optional
            Sequence constraints to enforce on generated padding as well as on the
            coding sequence. If None, padding is only checked for GC content and
            start/stop codons.

        Attributes
        ----------
        protein_sequence : str
            The final protein sequence (including any forced start/stop)
        full_dna_sequence : str
            Complete DNA sequence including adapters and padding
        coding_sequence : str
            DNA sequence encoding the protein
        gc_content_full_sequence : float
            GC content of the complete sequence
        gc_content_coding_sequence : float
            GC content of the coding sequence only
        correct_full_translation : bool
            Whether full sequence translates correctly
        correct_coding_translation : bool
            Whether coding sequence translates correctly
        '''
        # Store the raw protein sequence first
        self._raw_protein_sequence = protein_sequence.upper()
        self.codon_table = codon_table
        self.name = name
        self.adapter_3_prime = adapter_3_prime.upper()
        self.adapter_5_prime = adapter_5_prime.upper()
        self.total_length = total_length
        self.force_start_codon = force_start_codon
        self.force_stop_codon = force_stop_codon
        self.pad_location = pad_location
        self.avoid_adding_start_codons = avoid_adding_start_codons
        self.avoid_adding_stop_codons = avoid_adding_stop_codons
        self.coding_sequence = coding_sequence
        self.verify_coding_sequence=verify_coding_sequence
        self.return_best = return_best
        self.padding_attempts = padding_attempts
        # Constraints are enforced on the padding as well as the coding sequence;
        # None means no constraints were requested.
        self.constraint_set = constraint_set
        self.sequence_constraint_violations = []
        # Populated by build_library's final acceptance check.
        self.quality_report = None
        # Generate the modified protein sequence
        self.protein_sequence = self.get_protein_sequence()
        
        # Initialize other attributes
        self.translated_coding_sequence = None
        self.translated_full_sequence = None
        self.full_dna_sequence = None
        self.gc_content_full_sequence = None
        self.gc_content_coding_sequence = None
        self.correct_full_translation = None
        self.correct_coding_translation = None
        # Populated by build_library after optional synonymous sequence repair.
        # A successfully returned constrained sequence has an empty list.
        self.sequence_constraint_violations = []
        self.padding_3_prime = ""
        self.padding_5_prime = ""

        # figure out padding length.     
        self.padding_5_prime_length=None 
        self.padding_3_prime_length = None 
        
        # Validate the coding sequence if provided. add_coding_sequence honours
        # verify_coding_sequence itself, so there is no separate check here.
        if coding_sequence:
            self.add_coding_sequence(coding_sequence)
    
    def _get_padding_length(self):
        '''
        Get the padding length.
        '''
        if self.total_length==None:
            return 0, 0
        total_padding = self.total_length - len(self.coding_sequence) - len(self.adapter_5_prime) - len(self.adapter_3_prime)
        # A negative value here means total_length cannot accommodate this sequence.
        # Both padding branches below would silently skip, returning an over-length
        # sequence, so fail loudly instead.
        if total_padding < 0:
            raise ValueError(
                f"total_length={self.total_length} is too short for sequence '{self.name}': "
                f"it needs at least {self.total_length - total_padding} nt "
                f"(coding sequence + adapters).")
        if self.pad_location == None:
            self.padding_5_prime_length = total_padding // 2
            self.padding_3_prime_length = total_padding - self.padding_5_prime_length
        elif self.pad_location == 5:
            self.padding_5_prime_length = total_padding
            self.padding_3_prime_length = 0
        elif self.pad_location == 3:
            self.padding_5_prime_length = 0
            self.padding_3_prime_length = total_padding
        else:
            raise ValueError("pad_location must be 3 (3' end), 5 (5' end), or None (both ends)")
        return self.padding_5_prime_length, self.padding_3_prime_length
    
    def generate_padding(self):
        '''
        Generate padding sequences to reach target length.
        
        Generates padding sequences that maintain GC content and avoid
        unwanted features (like start/stop codons if specified).
        Padding is added according to pad_location setting.

        Returns
        -------
        self
            Returns self for method chaining
        '''
        self.padding_5_prime_length, self.padding_3_prime_length = self._get_padding_length()
        if self.padding_5_prime_length > 0:
            cur_padding = add_padding(self.gc_content_full_sequence,
                                               self.padding_5_prime_length,
                                               len(self.full_dna_sequence),
                                               self.codon_table.gc_range,
                                               avoid_adding_start_codons=self.avoid_adding_start_codons,
                                               avoid_adding_stop_codons=self.avoid_adding_stop_codons,
                                               return_best=self.return_best,
                                               num_attempts=self.padding_attempts)
            self._add_5_prime_padding(cur_padding)
        if self.padding_3_prime_length >0:
            cur_padding = add_padding(self.gc_content_full_sequence,
                                               self.padding_3_prime_length,
                                               len(self.full_dna_sequence),
                                               self.codon_table.gc_range,
                                               avoid_adding_start_codons=self.avoid_adding_start_codons,
                                               avoid_adding_stop_codons=self.avoid_adding_stop_codons,
                                               return_best=self.return_best,
                                               num_attempts=self.padding_attempts)
            self._add_3_prime_padding(cur_padding)

        self._enforce_constraints_on_padding()
        return self

    def _enforce_constraints_on_padding(self):
        '''
        Rewrite the padding so it does not re-introduce constraints the caller
        asked to avoid.

        Padding is generated after the coding sequence has been repaired and knows
        only about GC content and start/stop codons, so it routinely creates fresh
        splice sites and polyadenylation signals -- including at the junction with
        the coding sequence. Only padding bases are rewritten; the coding sequence
        is untouched.

        Raises
        ------
        ConstraintRepairError
            If some violation touching the padding cannot be removed.
        '''
        if self.constraint_set is None:
            return self
        if not self.padding_5_prime and not self.padding_3_prime:
            return self

        forbidden = set()
        if self.avoid_adding_start_codons:
            forbidden.update(('ATG',))
        if self.avoid_adding_stop_codons:
            forbidden.update(('TAA', 'TAG', 'TGA'))

        five_length = len(self.padding_5_prime)
        three_length = len(self.padding_3_prime)
        region = self.padding_5_prime + self.coding_sequence + self.padding_3_prime
        coding_end = five_length + len(self.coding_sequence)

        repaired, remaining = repair_padding_constraints(
            region,
            [(0, five_length), (coding_end, len(region))],
            self.constraint_set,
            forbidden=forbidden,
        )

        if remaining:
            preview = "; ".join(
                f"{item.kind} {item.start}:{item.end} {item.sequence} score={item.score:.2f}"
                for item in remaining[:5]
            )
            raise ConstraintRepairError(
                f"Could not generate padding for '{self.name}' that satisfies the requested "
                f"sequence constraints. Remaining violations: {preview}",
                remaining_violations=remaining,
                sequence=repaired,
            )

        # The coding sequence must come back byte-identical.
        assert repaired[five_length:coding_end] == self.coding_sequence

        self.padding_5_prime = repaired[:five_length]
        self.padding_3_prime = repaired[coding_end:] if three_length else ""
        self._update_full_dna_sequence()
        return self

    def get_protein_sequence(self, 
                             force_start_codon=None,
                             force_stop_codon=None):
        '''
        Get the protein sequence with optional start/stop codon modifications.

        Parameters
        ----------
        force_start_codon : bool, optional
            Override instance setting for forcing start codon
        force_stop_codon : bool, optional
            Override instance setting for forcing stop codon

        Returns
        -------
        str
            Modified protein sequence
        '''
        if force_start_codon is None:
            force_start_codon = self.force_start_codon
        if force_stop_codon is None:
            force_stop_codon = self.force_stop_codon

        sequence = self._raw_protein_sequence
        
        if force_start_codon:
            if not sequence.startswith('M'):
                sequence = 'M' + sequence
        if force_stop_codon:
            if not sequence.endswith('*'):
                sequence = sequence + '*'
        return sequence

    def add_coding_sequence(self, coding_sequence):
        '''
        Add and validate a coding DNA sequence.

        Parameters
        ----------
        coding_sequence : str
            DNA sequence encoding the protein sequence

        Returns
        -------
        self
            Returns self for method chaining

        Raises
        ------
        ValueError
            If sequence is invalid or doesn't translate correctly
        '''
        if not isinstance(coding_sequence, str):
            raise ValueError("Coding sequence must be a string")

        if len(coding_sequence) % 3 != 0:
            raise ValueError("Coding sequence length must be divisible by 3")

        coding_sequence = coding_sequence.upper()

        # Translate from the start of the coding sequence. The coding sequence
        # defines its own reading frame, so searching for an ATG here would break
        # any protein that does not begin with methionine.
        translated = translate_sequence(coding_sequence, return_stop_codon=True)
        if self.verify_coding_sequence and translated != self.protein_sequence:
            raise ValueError("Coding sequence does not translate to the target protein sequence")

        self.coding_sequence = coding_sequence
        self.gc_content_coding_sequence = self._calculate_gc_content(coding_sequence)
        self.translated_coding_sequence = translated
        self._update_full_dna_sequence()
        self.correct_coding_translation = self.verify_coding_sequence_translation()
        return self
    
    def _add_5_prime_padding(self, padding):
        '''
        Add padding sequence.
        
        Parameters:
        -----------
        padding : str
            The padding DNA sequence
        '''
        self.padding_5_prime = padding
        self._update_full_dna_sequence()
        return self
    
    def _add_3_prime_padding(self, padding):
        '''
        Add padding sequence.
        
        Parameters:
        -----------
        padding : str
            The padding DNA sequence
        '''
        self.padding_3_prime = padding
        self._update_full_dna_sequence()
        return self
                 

    def _update_full_dna_sequence(self):
        '''
        Update the complete DNA sequence and related properties.

        Internal method that updates the full DNA sequence by combining
        adapters, padding, and coding sequence. Also updates GC content
        and translation verification.

        Returns
        -------
        self
            Returns self for method chaining
        '''
        if self.coding_sequence:
            self.full_dna_sequence = self.adapter_5_prime + self.padding_5_prime + self.coding_sequence + self.padding_3_prime + self.adapter_3_prime 
            self.gc_content_full_sequence = self._calculate_gc_content(self.full_dna_sequence)
            # The whole construct is checked with the scanning-ribosome model: start
            # at the first ATG and read to the first stop. A False here is a real
            # warning that an adapter or padding will change the expressed product.
            self.translated_full_sequence = translate_open_reading_frame(self.full_dna_sequence,
                                                                         return_stop_codon=True)
            self.correct_full_translation = self.verify_translation_product()
        return self
    
    def _calculate_gc_content(self, target_sequence):
        '''
        Calculate GC content of a DNA sequence.

        Parameters
        ----------
        target_sequence : str
            DNA sequence to analyze

        Returns
        -------
        float
            GC content as fraction between 0 and 1
        '''
        if not target_sequence:
            return 0.0
        gc_count = target_sequence.count('G') + target_sequence.count('C')
        return gc_count / len(target_sequence)

    def get_coding_sequence(self):
        '''
        Get the coding sequence from the DNA sequence.
        
        Returns:
        --------
        str
            The coding DNA sequence
        '''
        return self.coding_sequence
    
    def get_dna_sequence(self):
        '''
        Get the complete DNA sequence.
        
        Returns:
        --------
        str
            The complete DNA sequence with adapters and padding
        '''
        return self.full_dna_sequence
    
    def get_gc_content(self, full_sequence=True):
        '''
        Get the GC content of the sequence.

        Parameters
        ----------
        full_sequence : bool, default=True
            If True, return GC content for complete sequence
            If False, return GC content for coding sequence only

        Returns
        -------
        float
            GC content as fraction between 0 and 1
        '''
        if full_sequence:
            return self.gc_content_full_sequence
        else:
            return self.gc_content_coding_sequence
        
    def verify_translation_product(self):
        '''
        Verify that the complete sequence translates to target protein.

        Returns
        -------
        bool
            True if translation matches target protein sequence
        '''
        if self.translated_full_sequence == self.protein_sequence:
            return True
        return False
    
    def verify_coding_sequence_translation(self):
        '''
        Verify that coding sequence translates to target protein.

        Returns
        -------
        bool
            True if coding sequence translates to target protein sequence
        '''
        if not self.translated_coding_sequence or not self.protein_sequence:
            return False
        return self.translated_coding_sequence == self.protein_sequence

    def __str__(self):
        '''
        Get string representation of the Sequence object.

        Returns
        -------
        str
            Multi-line string containing sequence information
        '''
        parts = []
        if self.protein_sequence:
            parts.append(f"Protein sequence: {self.protein_sequence}")
        if self.name:
            parts.append(f"Name: {self.name}")
        if self.correct_coding_translation==True:
            parts.append("Coding sequence translation: Correct")
        elif self.correct_coding_translation==False:
            parts.append("Coding sequence translation: Incorrect")
        else:
            parts.append("Coding sequence translation: Not available")
        if self.correct_full_translation==True:
            parts.append("Full sequence translation: Correct")
        elif self.correct_full_translation==False:
            parts.append("Full sequence translation: Incorrect")
        else:
            parts.append("Full sequence translation: Not available")
        if self.full_dna_sequence:
            parts.append(f"Full DNA sequence: {self.full_dna_sequence}")
            parts.append(f"GC content full sequence: {self.gc_content_full_sequence:.2f}")
        if self.coding_sequence:
            parts.append(f"Coding sequence: {self.coding_sequence}")
            parts.append(f"GC content coding sequence: {self.gc_content_coding_sequence:.2f}")
        return "\n".join(parts)
