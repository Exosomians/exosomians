import os
import exopy
import numpy as np
import pandas as pd
import argparse
import warnings
warnings.simplefilter("ignore")
warnings.filterwarnings("ignore")
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 

parser = argparse.ArgumentParser(description='Mutation Map Generator')
arguments_group = parser.add_argument_group("Parameters")
arguments_group.add_argument('-d', '--data_path', type=str, required=True,
                                help='data path')
arguments_group.add_argument('-n', '--n_generated', type=int, required=False, default=3,
                                help='')
arguments_group.add_argument('-o', '--output_path', type=str, required=True,
                                help='')
arguments_group.add_argument('-f', '--output_file', type=str, required=True,
                                help='')
arguments_group.add_argument('-m', '--model_name', type=str, required=True,
                                help='')


args = vars(parser.parse_args())

data_path = args['data_path']
n_generated = args['n_generated']
output_path = args['output_path']
output_filename = args['output_file']
model_name = args['model_name']

data_name = data_path.split("/")[-1]

sequences = pd.read_csv(data_path, delimiter='\t')['seq']

char_encoder = {
    'N': -1,
    'A': 0,
    'C': 1,
    'G': 2,
    'U': 3
}
char_reverse_encoder = {
    0: 'A',
    1: 'C',
    2: 'G', 
    3: 'U'
}

one_hot_encoder = {
    'A': [1, 0, 0, 0],
    'C': [0, 1, 0, 0],
    'G': [0, 0, 1, 0],
    'U': [0, 0, 0, 1]
}

one_hot_reverse_encoder = {
    0 : [1, 0, 0, 0],
    1 : [0, 1, 0, 0],
    2 : [0, 0, 1, 0],
    3 : [0, 0, 0, 1]
}

random_sequences = exopy.gen.generate_random_sequence(sequences.values, NUM=n_generated)

max_len = 50
sequences_encoded = exopy.pep.seq_encoder(random_sequences, char_encoder, max_len, unknown_char=True)

network = exopy.models.ExoCNN(seq_len=sequences_encoded.shape[1],
                              n_channels=sequences_encoded.shape[2],
                              n_classes=2,
                              padding="same",
                              use_batchnorm=True,
                              lr=0.0001,
                              model_path=f"./models/ExoCNN/{model_name}/",
                              dropout_rate=0.25,
                              )
network.restore_model_weights()

secretion_probs = network.model.predict(sequences_encoded)[:, 1]

os.makedirs(output_path, exist_ok=True)
random_sequences_df = pd.DataFrame({"id": [f'random_{n_generated}_{i}' for i in range(len(random_sequences))], 
                                    "seq": random_sequences, "secretion_prob": secretion_probs})
random_sequences_df.to_csv(os.path.join(output_path, f"{output_filename}.random.sequences.{n_generated}.tsv"), index=None, sep='\t')
