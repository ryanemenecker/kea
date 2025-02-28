'''
Codon usage frequencies.
'''

# From Genscript. Unsure exactly how they calculated this. 
yeast_codon_usage = {
    '*': {'TAA': 0.28, 'TAG': 0.22, 'TGA': 0.5},
    'A': {'GCT': 0.177, 'GCC': 0.222, 'GCA': 0.248, 'GCG': 0.353},
    'C': {'TGT': 0.651, 'TGC': 0.349},
    'D': {'GAT': 0.638, 'GAC': 0.362},
    'E': {'GAA': 0.685, 'GAG': 0.315},
    'F': {'TTT': 0.59, 'TTC': 0.41},
    'G': {'GGT': 0.351, 'GGC': 0.214, 'GGA': 0.284, 'GGG': 0.151},
    'H': {'CAT': 0.62, 'CAC': 0.38},
    'I': {'ATT': 0.512, 'ATC': 0.388, 'ATA': 0.1},
    'K': {'AAA': 0.584, 'AAG': 0.416},
    'L': {'TTA': 0.283, 'TTG': 0.114, 'CTT': 0.137, 'CTC': 0.069, 'CTA': 0.127, 'CTG': 0.269},
    'M': {'ATG': 1.0},
    'N': {'AAT': 0.476, 'AAC': 0.524},
    'P': {'CCT': 0.263, 'CCC': 0.217, 'CCA': 0.301, 'CCG': 0.218},
    'Q': {'CAA': 0.652, 'CAG': 0.348},
    'R': {'CGT': 0.11, 'CGC': 0.032, 'CGA': 0.1, 'CGG': 0.052, 'AGA': 0.544, 'AGG': 0.162},
    'S': {'TCT': 0.283, 'TCC': 0.151, 'TCA': 0.164, 'TCG': 0.08, 'AGT': 0.163, 'AGC': 0.159},
    'T': {'ACT': 0.237, 'ACC': 0.345, 'ACA': 0.281, 'ACG': 0.137},
    'V': {'GTT': 0.34, 'GTC': 0.192, 'GTA': 0.145, 'GTG': 0.323},
    'W': {'TGG': 1.0},
    'Y': {'TAT': 0.597, 'TAC': 0.403}
}

# from refseq, S288C specific. Is the normalized values for transcriptome-weighted codon usage.
# normalized by taking the transcriptome weighted value for each 
s288c_codon_usage = {
    '*': {'TAA': 0.4732, 'TAG': 0.2293, 'TGA': 0.2976},
    'A': {'GCT': 0.3693, 'GCC': 0.2211, 'GCA': 0.2967, 'GCG': 0.1129},
    'C': {'TGT': 0.6213, 'TGC': 0.3787},
    'D': {'GAT': 0.6512, 'GAC': 0.3488},
    'E': {'GAA': 0.7006, 'GAG': 0.2994},
    'F': {'TTT': 0.5947, 'TTC': 0.4053},
    'G': {'GGT': 0.4544, 'GGC': 0.1975, 'GGA': 0.2257, 'GGG': 0.1224},
    'H': {'CAT': 0.6419, 'CAC': 0.3581},
    'I': {'ATT': 0.4602, 'ATC': 0.2596, 'ATA': 0.2802},
    'K': {'AAA': 0.5843, 'AAG': 0.4157},
    'L': {'TTA': 0.2773, 'TTG': 0.279, 'CTT': 0.13, 'CTC': 0.0585, 'CTA': 0.1426, 'CTG': 0.1125},
    'M': {'ATG': 1.0},
    'N': {'AAT': 0.5961, 'AAC': 0.4039},
    'P': {'CCT': 0.3104, 'CCC': 0.1583, 'CCA': 0.4069, 'CCG': 0.1244},
    'Q': {'CAA': 0.6854, 'CAG': 0.3146},
    'R': {'CGT': 0.1414, 'CGC': 0.0597, 'CGA': 0.0703, 'CGG': 0.0417, 'AGA': 0.4736, 'AGG': 0.2133},
    'S': {'TCT': 0.2601, 'TCC': 0.1566, 'TCA': 0.2123, 'TCG': 0.0972, 'AGT': 0.1626, 'AGC': 0.1112},
    'T': {'ACT': 0.3425, 'ACC': 0.2112, 'ACA': 0.308, 'ACG': 0.1382},
    'V': {'GTT': 0.3866, 'GTC': 0.2022, 'GTA': 0.2176, 'GTG': 0.1936},
    'W': {'TGG': 1.0},
    'Y': {'TAT': 0.5665, 'TAC': 0.4335}
}

# calculated by raw numbers of codons, not transcriptome weighted. 
s288c_unweighted={
    '*': {'TAA': 0.4718, 'TGA': 0.2974, 'TAG': 0.2308},
    'A': {'GCT': 0.3692, 'GCC': 0.2212, 'GCA': 0.2967, 'GCG': 0.1128},
    'C': {'TGT': 0.6211, 'TGC': 0.3789},
    'D': {'GAT': 0.6512, 'GAC': 0.3488},
    'E': {'GAA': 0.7006, 'GAG': 0.2994},
    'F': {'TTT': 0.5947, 'TTC': 0.4053},
    'G': {'GGT': 0.4545, 'GGC': 0.1974, 'GGA': 0.2257, 'GGG': 0.1224},
    'H': {'CAT': 0.642, 'CAC': 0.358},
    'I': {'ATT': 0.4602, 'ATC': 0.2595, 'ATA': 0.2803},
    'K': {'AAA': 0.5843, 'AAG': 0.4157},
    'L': {'TTA': 0.2773, 'TTG': 0.2791, 'CTT': 0.13, 'CTC': 0.0585, 'CTA': 0.1426, 'CTG': 0.1125},
    'M': {'ATG': 1.0},
    'N': {'AAT': 0.596, 'AAC': 0.404},
    'P': {'CCT': 0.3104, 'CCC': 0.1583, 'CCA': 0.407, 'CCG': 0.1244},
    'Q': {'CAA': 0.6854, 'CAG': 0.3146},
    'R': {'CGT': 0.1414, 'CGC': 0.0597, 'CGA': 0.0703, 'CGG': 0.0416, 'AGA': 0.4737, 'AGG': 0.2133},
    'S': {'TCT': 0.2601, 'TCC': 0.1566, 'TCA': 0.2123, 'TCG': 0.0971, 'AGT': 0.1626, 'AGC': 0.1112},
    'T': {'ACT': 0.3425, 'ACC': 0.2112, 'ACA': 0.308, 'ACG': 0.1383},
    'V': {'GTT': 0.3866, 'GTC': 0.2022, 'GTA': 0.2176, 'GTG': 0.1935},
    'W': {'TGG': 1.0},
    'Y': {'TAT': 0.5665, 'TAC': 0.4335}
}

# add codon tables here. 
all_codon_tables={
    'yeast': yeast_codon_usage,
    's288c': s288c_codon_usage,
    's288c_unweighted': s288c_unweighted
}