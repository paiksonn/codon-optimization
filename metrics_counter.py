from collections import Counter
from Bio.Seq import Seq
from math import log
from math import exp
import sequence_preparations

sequence = sequence_preparations.sequence
seq_len = len(sequence)
codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]

codon_freq = sequence_preparations.calc_codon_freq(sequence)
nucl_freq = sequence_preparations.calc_nucl_freq(sequence) 

norm_codon_freq = sequence_preparations.calc_norm_codon_freq(sequence)
norm_nucl_freq = sequence_preparations.calc_norm_nucl_freq(sequence)

positional_nucl_freq = sequence_preparations.calc_positional_nucl_freq(sequence)
norm_positional_nucl_freq = sequence_preparations.calc_norm_positional_nucl_freq(sequence)

codon_pairs = sequence_preparations.count_codon_pair(sequence)

aa_table = {
    'F': ['TTT', 'TTC'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'I': ['ATT', 'ATC', 'ATA'],
    'M': ['ATG'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'Y': ['TAT', 'TAC'],
    'H': ['CAT', 'CAC'],
    'Q': ['CAA', 'CAG'],
    'N': ['AAT', 'AAC'],
    'K': ['AAA', 'AAG'],
    'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'],
    'C': ['TGT', 'TGC'],
    'W': ['TGG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG']
}

codon_to_aa = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }

# Relative Codon Adaptation index (RCA)
def RCA_calc(sequence: str, seq_len: int) -> float:
    RCA = 1
    for l in range(0, seq_len, 3):
        codon = sequence[l:l+3]
        f = codon_freq[codon]
        f1 = positional_nucl_freq[1][codon[0]]
        f2 = positional_nucl_freq[2][codon[1]]
        f3 = positional_nucl_freq[3][codon[2]]
        RCAxyz = f / (f1 * f2 * f3)
        RCA *= RCAxyz ** (1 / (seq_len / 3))
    
    return RCA


# Relative Codon Bias Strength index (RCBS)
def RCBS_calc(sequence: str, seq_len: int) -> float:
    RCBS = 1
    for l in range(0, seq_len, 3):
        codon = sequence[l:l+3]
        f = norm_codon_freq[codon]
        f1 = norm_positional_nucl_freq[0][codon[0]]
        f2 = norm_positional_nucl_freq[1][codon[1]]
        f3 = norm_positional_nucl_freq[2][codon[2]]
        RCBxyz = f / (f1 * f2 * f3)
        RCBS *= RCBxyz ** (1 / (seq_len / 3))
    
    return RCBS - 1


# GC frequency 
def GC_freq_calc() -> float:
    G_freq = nucl_freq['G'] 
    C_freq = nucl_freq['C'] 
    A_freq = nucl_freq['A']
    T_freq = nucl_freq['T']
    GC_freq = (G_freq + C_freq) / (G_freq + C_freq + A_freq + T_freq)
    return GC_freq * 100


# RSCU
def calculate_RSCU(codon_count):
    total_count = sum(codon_count.values())
    codon_RSCU = {}
    for codon, count in codon_count.items():
        expected_count = total_count / len(codon_count)
        if expected_count == 0:
            codon_RSCU[codon] = 0
        else:
            codon_RSCU[codon] = count / expected_count
    return codon_RSCU


# Average ration of Relative Synonymouse Codon Usage (ARSCU)
def ARSCU_calc(codons: list) -> float:
    codon_count = Counter(codons)
    GC_codons = [codon for codon in codon_count if codon.endswith(('G', 'C'))]
    AT_codons = [codon for codon in codon_count if codon.endswith(('A', 'T'))]
    GC_RSCU = calculate_RSCU(Counter(GC_codons))
    AT_RSCU = calculate_RSCU(Counter(AT_codons))
    GC_ending_RSCU_sum = sum(GC_RSCU[codon] for codon in GC_codons)
    AT_ending_RSCU_sum = sum(AT_RSCU[codon] for codon in AT_codons)
    ARSCU = GC_ending_RSCU_sum / (AT_ending_RSCU_sum / len(AT_codons))
    return ARSCU


# translate nucleotides to amino acids
def count_amino_acids(seq: str, pseudocount: float = 0.1) -> dict:
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    aa_seq = ''.join([codon_to_aa[c] for c in codons])
    aa_counts = {aa: aa_seq.count(aa) for aa in aa_table.keys()}
    
    # add pseudocount to each amino acid count
    for aa in aa_counts:
        aa_counts[aa] += pseudocount
        
    return aa_counts


# Codon Pair Score (CPS)
def CPS_calc(sequence: str, codon_pairs: dict, codon_freq: dict, codon_to_aa: dict) -> float:
    # count amino acids
    aa_counts = count_amino_acids(sequence)
    
    # calculate total number of codon pairs and amino acid pairs
    total_codon_pairs = sum(codon_pairs.values())
    total_aa_pairs = sum(aa_counts.values()) - 1
    
    # initialize variables for numerator and denominator
    numerator = 0
    denominator = 0
    
    # iterate over codon pairs and calculate CPS components
    for codon_pair, count in codon_pairs.items():
        codon_a, codon_b = codon_pair[:3], codon_pair[3:]
        count_a = codon_freq[codon_a]
        count_b = codon_freq[codon_b]
        aa_a = codon_to_aa[codon_a]
        aa_b = codon_to_aa[codon_b]
        count_aa_pair = aa_counts[aa_a] * aa_counts[aa_b]
        
        # calculate CPS components and add to numerator and denominator
        component = count / total_codon_pairs
        expected_count = (count_a * count_b / count_aa_pair) * total_codon_pairs
        numerator += component * log(count / expected_count)
        denominator += component
    
    # calculate CPS
    CPS = numerator / denominator
    
    return CPS


# Codon Pair Bias (CPB)
def CPB_calc(sequence: str) -> float:
    # calculate CPS for each codon pair in the sequence
    CPS_list = []
    for i in range(0, len(sequence)-6, 3):
        CPS = CPS_calc(sequence[i:i+6], codon_pairs, codon_freq, codon_to_aa)
        if CPS is not None:
            CPS_list.append(CPS)
    
    # calculate CPB
    if len(CPS_list) == 0:
        return None  # or raise an exception if appropriate
    else:
        CPB = sum(CPS_list) / ((len(sequence) - 6) / 3)
        return CPB


print('RCA metric:', RCA_calc(sequence, seq_len))
print('RCBS metric:', RCBS_calc(sequence, seq_len))
print('CPS metric:', CPS_calc(sequence, codon_pairs, codon_freq, codon_to_aa))
print('CPB metric:', CPB_calc(sequence))
print('GC frequency metric:', GC_freq_calc())
print('ARSCU metric:', ARSCU_calc(codons))
