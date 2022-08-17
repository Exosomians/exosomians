import os
import exoNet
import numpy as np
import pandas as pd
import argparse
import warnings

parser = argparse.ArgumentParser(description='Mutation Map Generator')
arguments_group = parser.add_argument_group("Parameters")
arguments_group.add_argument('-d', '--data_path', type=str, required=True,
                                help='data path')
arguments_group.add_argument('-m', '--model_path', type=str, required=True,
                                help='model path')
arguments_group.add_argument('-o', '--output_path', type=str, required=True,
                                help='')
arguments_group.add_argument('-f', '--output_file', type=str, required=True,
                                help='')

args = vars(parser.parse_args())

data_path = args['data_path']
model_path = args['model_path']
output_path = args['output_path']
output_filename = args['output_file']

sequences_df = pd.read_csv(data_path)

identifiers = sequences_df['id'].values
sequences = sequences_df['seq'].values

char_encoder = {
    'N': -1,
    'A': 0,
    'C': 1,
    'G': 2,
    'U': 3
}

max_len = 50
sequences_encoded = exoNet.pep.seq_encoder(sequences, char_encoder, max_len, unknown_char=True)


network = exoNet.models.ExoCNN(seq_len=sequences_encoded.shape[1],
                               n_channels=sequences_encoded.shape[2],
                               n_classes=2,
                               padding="same",
                               use_batchnorm=True,
                               lr=0.0001,
                               model_path=f"./models/ExoCNN/{model_path}/",
                               dropout_rate=0.25,
                               )

network.restore_model_weights()

predictions = network.model.predict(sequences_encoded)

results = pd.DataFrame({'id': identifiers, 'seq': sequences, 'secretion_prob': predictions[:, 1]})

# ev_results = results[results['id'].str.contains('ev')]
# ev_extreme_results = ev_results[ev_results['secretion_prob'] >= 0.95]

# ic_results = results[results['id'].str.contains('ic')]
# ic_extreme_results = ic_results[ic_results['secretion_prob'] <= 0.05]

os.makedirs(os.path.join(output_path), exist_ok=True)

results.to_csv(os.path.join(output_path, f"ExoCNN.{output_filename}.seqs.tsv"), index=None, sep='\t')
# ev_extreme_results.to_csv(os.path.join(output_path, f"ExoCNN.{output_filename}.ev.extreme.95.refined.seqs.tsv"), sep='\t', index=None)
# ic_extreme_results.to_csv(os.path.join(output_path, f"ExoCNN.{output_filename}.ic.extreme.95.refined.seqs.tsv"), sep='\t', index=None)