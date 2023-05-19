from collections import Counter
import sequence_preparations
import codonpair

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
    tmp = 1
    L = seq_len // 3  # Calculate the sequence length in codons
    
    for l in range(0, seq_len, 3):
        codon = sequence[l:l+3]
        f = norm_codon_freq[codon]
        f1 = norm_positional_nucl_freq[0][codon[0]]
        f2 = norm_positional_nucl_freq[1][codon[1]]
        f3 = norm_positional_nucl_freq[2][codon[2]]
        RCBxyz = f / (f1 * f2 * f3)
        tmp = RCBxyz ** (1 / L)
        RCBS *= tmp
    
    return RCBS - 1


# Codon Pair Strength/Codon Pair Bias (CPS/CPB)
CPS = codonpair.CodonPair.from_named_reference('E. coli')


# GC frequency 
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


# print all the metrics results 
print('RCA metric:', RCA_calc(sequence, seq_len))
print('RCBS metric:', RCBS_calc(sequence, seq_len))
print('CPS/CPB metrics:', CPS.cpb(sequence))
print('GC frequency metric:', GC_freq_calc())
print('ARSCU metric:', ARSCU_calc(codons))
print(len(sequence))
