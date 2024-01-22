def extract_ngrams(rna_sequences, n):
    ngrams_list = []
    for rna_sequence in rna_sequences:
        length = len(rna_sequence)
        ngrams_list.append([rna_sequence[i:i + n] for i in range(length - n + 1)])
    return ngrams_list

