import numpy as np
from keras_preprocessing.sequence import pad_sequences


def seq_encoder(seq_data, char_encoder, max_len, unknown_char=True, unknown_value=0.25):
    one_hot_encoder = {}
    if unknown_char:
        one_hot_encoder[-1] = [0 for i in range(len(char_encoder.keys()) - 1)]
        for i in range(len(char_encoder.keys())):
            one_hot_encoder[i] = [0 for _ in range(i)] + [1] + [0 for __ in range(len(char_encoder.keys()) - 2 - i)]
    else:
        one_hot_encoder[-1] = [0 for i in range(len(char_encoder.keys()))]
        for i in range(len(char_encoder.keys())):
            one_hot_encoder[i] = [0 for _ in range(i)] + [1] + [0 for __ in range(len(char_encoder.keys()) - 1 - i)]
    encoded_sequences = []
    for sequence in seq_data:
        encoded_sequence = [char_encoder[char] for char in sequence]
        encoded_sequences.append(encoded_sequence)

    encoded_sequences = pad_sequences(encoded_sequences, maxlen=max_len, padding='post', truncating='post', value=-1)
    onehot_sequences = []
    for encoded_sequence in encoded_sequences.tolist():
        onehot_sequence = [one_hot_encoder[enc] for enc in encoded_sequence]
        onehot_sequences.append(onehot_sequence)

    onehot_sequences = np.array(onehot_sequences)
    return onehot_sequences
