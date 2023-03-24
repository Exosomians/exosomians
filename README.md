# Exosomians

## Getting Started

## Installation

### Installation with pip

To install the latest version from PyPI, simply use the following bash script:

```bash
pip install exopy
```

or you can clone this repository and install via setup.py file:

```bash
git clone https://github.com/Exosomians/exosomians
cd exosomians
python setup.py -q
``` 

## Examples

### Inference

You can use pre-trained models to make predictions on your own datasets

```python
import exopy as exo

# Load the pre-trained model
model = exo.ml.ExoGRU.load('./saved_models/ExoGRU/exogru_best-v2.ckpt')

# Prepare dataset
data = model.prepare_data('/path/to/fasta/data.fasta', seq_key='seq')

# Get the predictions 
df_results = model.predict(data, batch_size=128)  # Results will be stored in a pandas dataframe 


```

### Train from scratch

```python
import exopy as exo

exo.ml.ExoGRU.setup_dataset(path='/path/to/dataset/design.mat.csv',
                            seq_key='seq',
                            target_key='label',
                            fraction=1.0)

config = {
    'activation_fn': 'relu',
    'batch_size': 32,
    'bidirectional': False,
    'dropout_rate': 0.1,
    'lr': 0.00011342016019358544,
    'n_head_hidden': 512, 'n_head_layers': 2,
    'n_hidden': 1024,
    'n_layers': 1,
    'network': 'exogru',
    'use_batch_norm': True,
    'use_layer_norm': False
}

model = exo.ml.ExoGRU(**config)

model.fit(max_epochs=2000,
          train_size=0.8,
          batch_size=128,
          early_stopping_patience=5,
          check_val_every_n_epoch=3,
          save_path='./saved_models/ExoGRU/',
          )

```

### Sample Notebooks

| Model   | Path                                                                                                                                |
|---------|-------------------------------------------------------------------------------------------------------------------------------------|
| ExoGRU  | [notebooks/training_demos/ExoGRU.ipynb](https://github.com/Exosomians/exosomians/blob/main/notebooks/training_demos/ExoGRU.ipynb)   |
| ExoCNN  | [notebooks/training_demos/ExoCNN.ipynb](https://github.com/Exosomians/exosomians/blob/main/notebooks/training_demos/ExoCNN.ipynb)   |
| ExoLSTM | [notebooks/training_demos/ExoLSTM.ipynb](https://github.com/Exosomians/exosomians/blob/main/notebooks/training_demos/ExoLSTM.ipynb) |

