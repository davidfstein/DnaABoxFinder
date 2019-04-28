import sys
from argparse import parseargs

def minimum_skew(genome):
    minimum_index = []
    minimum = 0
    skew = 0
    skews = []
    for nuc in genome:
        if nuc == "C":
            skew -= 1
        if nuc == "G":
            skew += 1
        if skew < minimum:
            minimum = skew
        skews.append(skew)
    for i in range(0, len(skews)):
        if skews[i] == minimum:
            minimum_index.append(i + 1)
    return minimum_index

def frequent_words_with_mistmatches_and_reverse_complements(text, k, d):
    frequent_patterns = set()
    frequency_array = computing_frequencies_with_mismatches(text, k, d)
    max_count = max(frequency_array)
    for i in range(0, 4**k - 1):
        if frequency_array[i] == max_count:
            pattern = number_to_pattern(i, k)
            frequent_patterns.add(pattern)
    return frequent_patterns, max_count

def computing_frequencies_with_mismatches(text, k, d):
    freq_arr = initialize_frequency_array(k)
    for i in range(0, len(text) - k + 1):
        pattern = text[i:i+k]
        rc_pattern = reverse_complement(pattern)
        rc_neighborhood = neighbors(rc_pattern, d)
        neighborhood = neighbors(pattern, d)
        for neighbor in neighborhood:
            index = convert_sequence_to_num(neighbor)
            freq_arr[index] += 1
        for neighbor in rc_neighborhood:
            index = convert_sequence_to_num(neighbor)
            freq_arr[index] += 1
    return freq_arr
        
def reverse_complement(pattern):
    complements = {'A':'T', 'G':'C', 'C':'G', 'T':'A'}
    complement = [complements[nuc] for nuc in pattern]
    complement.reverse()
    return ('').join(complement)

def convert_sequence_to_num(sequence):
    nuc_vals = {"A":0, "C":1, "G":2, "T":3}
    seq_hex = 0
    power = len(sequence) - 1
    for nuc in sequence:
        hex_val = nuc_vals[nuc] * 4**power
        seq_hex += hex_val
        power -= 1
    return seq_hex

def number_to_pattern(index, k):
    num_to_nuc = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    if k == 1:
        return num_to_nuc[index]
    quotient = index // 4
    remainder = index % 4
    symbol = num_to_nuc[remainder]
    pattern = number_to_pattern(quotient, k - 1)
    return pattern + symbol
    
def initialize_frequency_array(k):
    return [0 for i in range(0, 4**k)]

def neighbors(pattern, d):
    nucs = {'A', 'G', 'T', 'C'}
    if d == 0:
        return pattern
    if len(pattern) == 1:
        return nucs
    neighborhood = set()
    suffix_neighbors = neighbors(pattern[1:], d)
    for suffix in suffix_neighbors:
        if hamming_distance(pattern[1:], suffix) < d:
            for nuc in nucs:
                neighborhood.add(nuc + suffix)
        else:
            neighborhood.add(pattern[0] + suffix)
    return neighborhood

def hamming_distance(p, q):
    hamming_distance = 0
    p_len = len(p)
    q_len = len(q)
    for i in range(0, p_len):
        if i < p_len and i < q_len and p[i] != q[i]:
            hamming_distance += 1
        elif not i < p_len:
            hamming_distance = hamming_distance + q_len - 1
        elif not i < q_len:
            hamming_distance = hamming_distance + p_len - 1
    return hamming_distance

def find_dnaa_box(genome, window, num_kmers, dnaa_box_length, hamming_distance):
    potential_dnaa_boxes = set()
    potential_ori_index = minimum_skew(genome)[0]

    genome_window_ori_at_center = genome[potential_ori_index - (window // 2):potential_ori_index + (window // 2)]
    frequent_words_ori_at_center, count_at_center = frequent_words_with_mistmatches_and_reverse_complements(genome_window_ori_at_center, 
                                                                                                            dnaa_box_length, 
                                                                                                            hamming_distance)

    genome_window_ori_at_start = genome[potential_ori_index:potential_ori_index + window]
    frequent_words_ori_at_start, count_at_start = frequent_words_with_mistmatches_and_reverse_complements(genome_window_ori_at_start, 
                                                                                                          dnaa_box_length, 
                                                                                                          hamming_distance)

    genome_window_ori_at_end = genome[potential_ori_index - window:potential_ori_index]
    frequent_words_ori_at_end, count_at_end = frequent_words_with_mistmatches_and_reverse_complements(genome_window_ori_at_end, 
                                                                                                      dnaa_box_length, 
                                                                                                      hamming_distance)
                                                                                            
    if count_at_center >= num_kmers:
        for kmer in frequent_words_ori_at_center:
            potential_dnaa_boxes.add(kmer)
    if count_at_start >= num_kmers:
        for kmer in frequent_words_ori_at_start:
            potential_dnaa_boxes.add(kmer)
    if count_at_end >= num_kmers:
        for kmer in frequent_words_ori_at_end:
            potential_dnaa_boxes.add(kmer)
    
    return potential_dnaa_boxes

if __name__ == '__main__':
    genome_file_path, window, num_kmers, dnaa_box_length, hamming_dist = parseargs(sys.argv[1:])
    with open(genome_file_path) as text:
        genome = ''
        for line in text.readlines():
            genome = genome + line.strip().strip('\n')
        print(*find_dnaa_box(genome, window, num_kmers, dnaa_box_length, hamming_dist), sep=' ')

        
        