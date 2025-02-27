'''
Utilities for restriction enzyme analysis, including site identification and manipulation.

NEED TO TEST THIS OUT!!

'''

import re
import numpy as np
from kea.data.restriction_enzymes import restriction_enzymes, nt_code

def parse_restriction_site(site_notation):
    """
    Parse a restriction enzyme cut site notation.
    
    Parameters
    ----------
    site_notation : str
        The restriction site notation (e.g., "G/GATCC", "GAAGAC(2/6)")
        
    Returns
    -------
    dict
        Dictionary containing recognition sequence, cut positions, and regex pattern
    """
    result = {'original': site_notation}
    
    # Handle cut site with position numbers
    if '(' in site_notation and ')' in site_notation:
        # Extract the recognition sequence and cut positions
        pattern = r'([^(]+)\(([^/]+)/([^)]+)\)'
        match = re.search(pattern, site_notation)
        if match:
            site = match.group(1)
            top_cut = match.group(2)
            bottom_cut = match.group(3)
            
            result['recognition_site'] = site
            # Handle 'none' values for nicking enzymes
            result['top_cut'] = None if top_cut.lower() == 'none' else int(top_cut)
            result['bottom_cut'] = None if bottom_cut.lower() == 'none' else int(bottom_cut)
    # Handle standard notation with slash
    elif '/' in site_notation:
        parts = site_notation.split('/')
        if len(parts) == 2:
            result['recognition_site'] = parts[0] + parts[1]
            result['top_cut'] = len(parts[0])
            result['bottom_cut'] = len(parts[0])
    # Handle complex double-cut enzymes
    elif site_notation.startswith('('):
        pattern = r'\((\d+)/(\d+)\)([^(]+)\((\d+)/(\d+)\)'
        match = re.search(pattern, site_notation)
        if match:
            result['recognition_site'] = match.group(3)
            # For these complex sites, store both cut positions
            result['first_cut'] = (int(match.group(1)), int(match.group(2)))
            result['second_cut'] = (int(match.group(4)), int(match.group(5)))
    else:
        # No cut site specified, just a recognition site
        result['recognition_site'] = site_notation
    
    # Generate regex pattern
    result['pattern'] = site_to_regex(result.get('recognition_site', ''))
    
    return result

def site_to_regex(site):
    """Convert a recognition site with degenerate bases to a regex pattern."""
    pattern = ''
    for char in site:
        if char in 'ACGT':
            pattern += char
        elif char in nt_code:
            # Ensure proper handling of the nt_code dictionary
            bases = nt_code.get(char, '')
            pattern += f'[{bases}]'
        else:
            # Skip non-nucleotide characters
            continue
    return pattern

def identify_restriction_sites(sequence, min_buffer=0):
    """
    Identifies restriction enzyme cut sites in a given DNA sequence.
    
    Parameters
    ----------
    sequence : str
        The DNA sequence to check for restriction sites
    min_buffer : int, optional
        Minimum number of nucleotides required before/after the cut site
        
    Returns
    -------
    dict
        Dictionary with enzyme names as keys and lists of cut positions as values
    """
    sequence = sequence.upper()
    results = {}
    
    for enzyme, site_notation in restriction_enzymes.items():
        # Parse the restriction site
        site_info = parse_restriction_site(site_notation)
        
        # Skip if we couldn't parse the site properly
        if 'pattern' not in site_info or not site_info['pattern']:
            continue
            
        # Find all matches
        matches = []
        for match in re.finditer(site_info['pattern'], sequence, re.IGNORECASE):
            start = match.start()
            end = match.end()
            match_info = {
                'start': start,
                'end': end,
                'match_sequence': sequence[start:end]
            }
            
            # Calculate cut positions
            valid_match = True
            
            # Standard cut sites
            if 'top_cut' in site_info and site_info['top_cut'] is not None:
                top_pos = start + site_info['top_cut']
                match_info['top_cut'] = top_pos
                
                # Check buffer
                if top_pos < min_buffer or top_pos > len(sequence) - min_buffer:
                    valid_match = False
                    
            if 'bottom_cut' in site_info and site_info['bottom_cut'] is not None:
                bottom_pos = start + site_info['bottom_cut'] 
                match_info['bottom_cut'] = bottom_pos
                
                # Check buffer
                if bottom_pos < min_buffer or bottom_pos > len(sequence) - min_buffer:
                    valid_match = False
            
            # Handle complex double-cut enzymes
            if 'first_cut' in site_info and 'second_cut' in site_info:
                first_top, first_bottom = site_info['first_cut']
                second_top, second_bottom = site_info['second_cut']
                
                match_info['first_top_cut'] = start - first_top
                match_info['first_bottom_cut'] = start - first_bottom
                match_info['second_top_cut'] = end + second_top
                match_info['second_bottom_cut'] = end + second_bottom
                
                # Check all positions for buffer
                cut_positions = [start - first_top, start - first_bottom, 
                                end + second_top, end + second_bottom]
                for pos in cut_positions:
                    if pos < min_buffer or pos > len(sequence) - min_buffer:
                        valid_match = False
            
            if valid_match:
                matches.append(match_info)
        
        if matches:
            results[enzyme] = matches
            
    return results

def check_restriction_sites(sequence, enzymes=None, min_buffer=0):
    """
    Check a sequence for specified restriction enzyme sites.
    
    Parameters
    ----------
    sequence : str
        DNA sequence to analyze
    enzymes : list, optional
        List of enzyme names to check (None for all enzymes)
    min_buffer : int, optional
        Minimum nucleotides needed before/after the cut site
        
    Returns
    -------
    dict
        Dictionary with enzyme names as keys and match information
    """
    if enzymes is None:
        return identify_restriction_sites(sequence, min_buffer)
    
    sequence = sequence.upper()
    results = {}
    
    # Filter the enzymes to only those requested
    filtered_enzymes = {name: site for name, site in restriction_enzymes.items() 
                        if name in enzymes}
    
    # Use the same logic as identify_restriction_sites but with filtered enzymes
    for enzyme, site_notation in filtered_enzymes.items():
        site_info = parse_restriction_site(site_notation)
        
        # Skip if we couldn't parse the site properly
        if 'pattern' not in site_info or not site_info['pattern']:
            continue
            
        # Find all matches (same code as in identify_restriction_sites)
        # ... (rest of the matching logic)
        
        # For brevity, I'll use the identify_restriction_sites function
        enzyme_result = identify_restriction_sites(sequence, min_buffer)
        if enzyme in enzyme_result:
            results[enzyme] = enzyme_result[enzyme]
    
    return results
