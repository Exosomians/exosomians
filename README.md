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
git clone https://github.com/saberi1/Exosomians
cd Exosomians
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

exo.ml.ExoGRU.setup_dataset(path='/home/mohsen/projects/exosomians-v2/data/design.mat.csv',
                            seq_key='seq',
                            target_key='label',
                            fraction=0.5)

config = {
    'n_conv_blocks': 2,
    'n_conv_layers': 3,
    'kernel_size': 7,
    'n_filters': 32,
    'lr': 1e-3,
    'pooling': 'avg',
    'pooling_size': 2,
    'activation_fn': 'tanh',
    'use_batch_norm': True,
    'use_layer_norm': False,
    'dropout_rate': 0.1,
    'n_head_layers': 2,
    'n_head_hidden': 128,
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

