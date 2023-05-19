# read the nucleotide sequence a.k.a. gene 
with open("D:/Диплом/Sequences/Short/human il21 opt1.fa", 'r') as file:
    file.readline()
    lines = file.readlines()

# concatenate all strings in the list
sequence = "".join(lines)

# remove any newline characters from the sequence
sequence = sequence.replace("\n", "")

# make letters BIG
sequence = sequence.upper()

# calculate the codon frecuency 
def calc_codon_freq(sequence: str) -> dict: 
    codon_frequency = {}
    for i in range(0, len(sequence), 3):  
        codon = sequence[i:i+3]
        if codon in codon_frequency:
            codon_frequency[codon] += 1
        else:
            codon_frequency[codon] = 1
    
    return codon_frequency
    
    
# calculate nucleotide frequency
def calc_nucl_freq(sequence: str) -> dict:
    nucleotide_freq = {}
    for nucleotide in sequence:
        if nucleotide in nucleotide_freq:
            nucleotide_freq[nucleotide] += 1
        else:
            nucleotide_freq[nucleotide] = 1
    
    return nucleotide_freq


# calculate NORMALIZED codon frecuency 
def calc_norm_codon_freq(sequence: str) -> dict:
    norm_codon_frequency = {}
    total_codons = len(sequence) // 3  # Calculate the total number of codons in the sequence
    
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        if codon in norm_codon_frequency:
            norm_codon_frequency[codon] += 1
        else:
            norm_codon_frequency[codon] = 1
    
    # Calculate the frequency of each codon and divide by the total number of codons
    for codon in norm_codon_frequency:
        norm_codon_frequency[codon] /= total_codons
    
    return norm_codon_frequency
    
    
# calculate NORMALIZED nucleotide frequency
def calc_norm_nucl_freq(sequence: str) -> dict:
    norm_nucleotide_freq = {}
    for nucleotide in sequence:
        if nucleotide in norm_nucleotide_freq:
            norm_nucleotide_freq[nucleotide] += 1
        else:
            norm_nucleotide_freq[nucleotide] = 1
            
    # calculate the total number of codons in sequence 
    total_nucleotides = sum(norm_nucleotide_freq.values())

    # calculate the frequency of each codon and store it in the dictionary
    for codon in norm_nucleotide_freq:
        norm_nucleotide_freq[codon] /= total_nucleotides
    
    return norm_nucleotide_freq


# calculate nucleotide frequency DEPENDING ON THE POSITION IN CODON
def calc_positional_nucl_freq(sequence: str) -> dict:
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]

    freq_dict = {1: {'A': 0, 'T': 0, 'G': 0, 'C': 0},
                2: {'A': 0, 'T': 0, 'G': 0, 'C': 0},
                3: {'A': 0, 'T': 0, 'G': 0, 'C': 0}}

    for codon in codons:
        for i in range(3):
            freq_dict[i+1][codon[i]] += 1
            
    return freq_dict


# calculate NORMALIZED nucleotide frequency DEPENDING ON THE POSITION IN CODON
def calc_norm_positional_nucl_freq(sequence: str) -> dict:
    codon_positions = [0, 1, 2]
    nucleotides = ['A', 'T', 'G', 'C']
    codon_freqs = {position: {nt: 0 for nt in nucleotides} for position in codon_positions}
    total_freqs = {position: 0 for position in codon_positions}

    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        for j, nt in enumerate(codon):
            codon_freqs[j][nt] += 1
            total_freqs[j] += 1

    # calculate normalized frequencies
    normalized_freqs = {position: {nt: freq / total_freqs[position] for nt, freq in codon_freqs[position].items()} 
                        for position in codon_positions}
            
    return normalized_freqs


# codon pair count
def count_codon_pair(sequence: str) -> dict:
    codon_pair_count = {}
    for i in range(0, len(sequence)-6, 3):
        codon_pair = sequence[i:i+6]
        if codon_pair in codon_pair_count:
            codon_pair_count[codon_pair] += 1
        else:
            codon_pair_count[codon_pair] = 1
    return codon_pair_count
