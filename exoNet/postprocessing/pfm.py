import numpy as np

def compute_pfm(data, motif_scores, motif_len=10, n_filters=32):
    position_frequency_matrix = np.zeros(shape=(32, motif_len, 4))
    for i in range(n_filters):
        for j in range(motif_scores.shape[0]):
            k = np.argmax(motif_scores[j, :, i]) 
            subsequence = data[j, k:k + motif_len, :]
            position_frequency_matrix[i] += subsequence
    return position_frequency_matrix