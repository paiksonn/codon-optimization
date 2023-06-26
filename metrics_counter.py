from collections import Counter
from math import log
from math import exp
from Bio.SeqUtils import gc_fraction
from Bio.Data import CodonTable
import sequence_preparations
import codonpair
import python_codon_tables as pct

sequence = sequence_preparations.sequence  # read sequence from file
seq_len = len(sequence)  # take its length

codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]  # make a list of codons in sequence
codon_pairs = sequence_preparations.count_codon_pair(sequence)  # make a list of codon pairs and its quatity

codon_freq = sequence_preparations.calc_codon_freq(sequence)  # calculate codon frequency
nucl_freq = sequence_preparations.calc_nucl_freq(sequence)  # calculate nucleotide frequency

norm_codon_freq = sequence_preparations.calc_norm_codon_freq(sequence)  # calculate normalized codon frequency
norm_nucl_freq = sequence_preparations.calc_norm_nucl_freq(sequence)  # calculate normalized nucleotide frequency

positional_nucl_freq = sequence_preparations.calc_positional_nucl_freq(sequence)  # calculate nucleotide frequency depending on its position in codon
norm_positional_nucl_freq = sequence_preparations.calc_norm_positional_nucl_freq(sequence)  # calculate normalized nucleotide frequency depending on its position in codon

# PRINT THE LIST OF NAMES OF ALL AVAILABLE TABLES
# print ('Available tables:', pct.available_codon_tables_names)

# LOAD ONE TABLE BY NAME
codon_usage_table = pct.get_codons_table("h_sapiens_9606")  # codon usage table

aa_table = {  # list of amino acids and codons which codes them 
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


# its all in the func name, dude 
def translate_codon(codon):
    genetic_code = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    return genetic_code[codon]

# get codon frequency from codon usage table
def get_codon_frequency(codon: str) -> float:
    for amino_acid, codon_dict in codon_usage_table.items():
        if codon in codon_dict:
            return codon_dict[codon]
    return None

# get normalized nucleotide positional frequency referenced from codon usage table
def calculate_nucleotide_frequency(nucleotide: str, position: int) -> float:
    frequency = 0.0
    counter = 0
    for amino_acid, codon_dict in codon_usage_table.items():
        for codon, codon_frequency in codon_dict.items():
            if codon[position-1] == nucleotide:
                frequency += codon_frequency
                counter += 1
    return frequency / counter

# Relative Codon Adaptation index (RCA)
def RCA_calc(sequence: str, seq_len: int) -> float:
    RCA = 1
    L = (seq_len / 3)  # Calculate the sequence length in codons
    for l in range(0, seq_len, 3):
        codon = sequence[l:l+3]
        f = get_codon_frequency(codon)
        f1 = calculate_nucleotide_frequency(codon[0], 1)
        f2 = calculate_nucleotide_frequency(codon[1], 2)
        f3 = calculate_nucleotide_frequency(codon[2], 3)
        RCAxyz = f / (f1 * f2 * f3)
        RCA += log(RCAxyz) 
    # log П x_i^n = /Sigma log x_i^n = n * /Sigma log x_i
    # /Sigma == RCA, n == 1 / L
    return exp(RCA * (1 / L))


# Relative Codon Bias Strength index (RCBS)
def RCBS_calc(sequence: str, seq_len: int) -> float:
    RCBS = 0
    L = seq_len // 3  # Calculate the sequence length in codons
    
    for l in range(0, seq_len, 3):
        codon = sequence[l:l+3]
        f = norm_codon_freq[codon]
        f1 = norm_positional_nucl_freq[0][codon[0]]
        f2 = norm_positional_nucl_freq[1][codon[1]]
        f3 = norm_positional_nucl_freq[2][codon[2]]
        RCBxyz = f / (f1 * f2 * f3)
        RCBS += log(RCBxyz)
    # log П x_i^n = /Sigma log x_i^n = n * /Sigma log x_i
    # /Sigma == RCA, n == 1 / L
    return exp((RCBS * (1 / L))) - 1


# Codon Pair Strength/Codon Pair Bias (CPS/CPB)
CPS = codonpair.CodonPair.from_named_reference('E. coli')

def codon_pair_freq(AB: str) -> float:
    for codonpair, freq in codon_pairs.items():
        if AB == codonpair:
            return freq / len(codon_pairs)


# Codon Pair Strength for Homo sapience (CPS)
def CPS_human_calc(codon: str) -> float:
    with open('human_CPS.txt', 'r') as f:
        next(f)
        for line in f:
            line = line.strip()
            fields = line.split()
            if len(fields) > 2 and fields[1] == codon:
                return float(fields[5])
    return 0
    


# Codon Pair Bias for Homo sapience (CPB)
def CPB_human_calc(sequence: str) -> float:
    CPS = 0
    for codon_pair in range(0, seq_len-6, 3):
        # print(sequence[codon_pair:codon_pair+6])
        CPS += CPS_human_calc(sequence[codon_pair:codon_pair+6])
    CPB = CPS / ((seq_len / 3) - 1)
    return CPS, CPB
    


# GC frequency HANDWRITTEN
def GC_freq_calc() -> float:
    G_freq = nucl_freq['G'] 
    C_freq = nucl_freq['C'] 
    A_freq = nucl_freq['A']
    T_freq = nucl_freq['T']
    GC_freq = (G_freq + C_freq) / (G_freq + C_freq + A_freq + T_freq)
    return GC_freq * 100


# Relative Synonimouse Codon Usage (RSCU)
def calculate_RSCU(codon_count):
    # Calculate total counts of each amino acid
    amino_acid_counts = {}
    for codon, count in codon_count.items():
        amino_acid = translate_codon(codon)
        amino_acid_counts[amino_acid] = amino_acid_counts.get(amino_acid, 0) + count

    # Calculate RSCU for each codon
    rscu_values = {}
    for codon, count in codon_count.items():
        amino_acid = translate_codon(codon)
        if amino_acid == '*':
            expected_count = 0  # Set count to 0 for stop codons
        else:
            synonymous_codons = aa_table[amino_acid]
            synonymous_count = sum(codon_count.get(c, 0) for c in synonymous_codons)
            if synonymous_count == 0:
                expected_count = 0  # Set count to 0 for synonymous codons with zero count
            else:
                expected_count = amino_acid_counts[amino_acid] / synonymous_count
        
        if expected_count == 0:
            expected_count = 0.1
            
        rscu = count / expected_count

        rscu_values[codon] = rscu

    return rscu_values

# Average ration of Relative Synonymouse Codon Usage (ARSCU)
def ARSCU_calc(codons: list) -> float:
    codon_count = Counter(codons)
    RSCU = calculate_RSCU(codon_count)
    GC_codons = {key: value for key, value in RSCU.items() if key[-1] in ('G', 'C')}
    AT_codons = {key: value for key, value in RSCU.items() if key[-1] in ('A', 'T')}
    
    GC_ending_RSCU_sum = sum(value for value in GC_codons.values())
    AT_ending_RSCU_sum = sum(value for value in AT_codons.values())
    ARSCU = GC_ending_RSCU_sum / (AT_ending_RSCU_sum / 18)
    return ARSCU


# Effective Number of Codons (ENC)
def ENC_calc() -> float:
    G_3_position_frec = norm_positional_nucl_freq[2]['G']
    C_3_position_frec = norm_positional_nucl_freq[2]['C']
    s = G_3_position_frec + C_3_position_frec
    ENC = 2 + s + (29 / (s*s + (1 - s) ** 2))
    return ENC


# print all the metrics results 
print('RCA metric:', RCA_calc(sequence, seq_len))
print('RCBS metric:', RCBS_calc(sequence, seq_len))
# print('CPS/CPB in E.coli metrics:', CPS.cpb(sequence))
print('CPS/CPB in Homo sapience metrics:', CPB_human_calc(sequence))
print('GC frequency metric:', gc_fraction(sequence) * 100)
print('ARSCU metric:', ARSCU_calc(codons))
print('ENC frequency metric:', ENC_calc())
print(len(sequence))

# print(codon_usage_table)
# print(codon_pairs)