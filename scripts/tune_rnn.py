import sys
import os
import time

import numpy as np

os.chdir('/home/mohsen/projects/exosomians-v2/')

from ray import tune, air
from ray.tune import CLIReporter
from ray.tune.schedulers import ASHAScheduler
from ray.tune.integration.pytorch_lightning import TuneReportCallback

import exopy as exo

np.random.seed(int(time.time()) % 1000)


def train_fn(config):
    batch_size = config.get('batch_size')

    exo.ml.ExoNet.setup_dataset(path='/home/mohsen/projects/exosomians-v2/data/design.mat.csv',
                                seq_key='seq',
                                target_key='label',
                                fraction=0.5)

    model = exo.ml.ExoNet(**config)

    callback = TuneReportCallback(metrics={"loss": "loss",
                                           "val_loss": "val_loss",
                                           "val_auroc": "val_auroc",
                                           "val_precision": "val_precision",
                                           "val_recall": "val_recall",
                                           "val_specificity": "val_specificity",
                                           },
                                  on="validation_end")

    model.fit(max_epochs=2000,
              train_size=0.8,
              batch_size=batch_size,
              early_stopping_patience=5,
              check_val_every_n_epoch=3,
              callbacks=[callback],
              save_path=None,
              )

config = {
    'n_layers': tune.choice([1, 2, 3]),
    'n_hidden': tune.choice([64, 128, 256, 512, 1024]),
    'network': tune.choice(['exogru', 'exolstm']),
    'bidirectional': tune.choice([True, False]),
    'lr': tune.loguniform(1e-6, 1e-2),
    'activation_fn': tune.choice(['relu', 'leaky_relu', 'tanh']),
    'use_batch_norm': tune.choice([True, False]),
    'use_layer_norm': tune.sample_from(
        lambda spec: False if spec.config.use_batch_norm else np.random.choice([True, False])),
    'dropout_rate': tune.choice([0.0, 0.1, 0.2, 0.25]),
    'n_head_layers': tune.choice([1, 2, 3]),
    'n_head_hidden': tune.choice([512, 256, 128, 64, 32]),
    'batch_size': tune.choice([32, 64, 128, 256, 512, 1024]),
}

scheduler = ASHAScheduler(
    max_t=500000,
    grace_period=5,
    reduction_factor=4)

reporter = CLIReporter(
    parameter_columns=list(config.keys()),
    metric_columns=['loss', 'val_loss', 'val_auroc', 'val_precision', 'val_recall', 'val_specificity'])

tuner = tune.Tuner(
    tune.with_resources(
        train_fn,
        resources={
            'cpu': 5,
            'gpu': 0.1,
        }
    ),
    tune_config=tune.TuneConfig(
        metric="val_loss",
        mode="min",
        scheduler=scheduler,
        num_samples=50000,
    ),
    run_config=air.RunConfig(
        name="tune_ExoRNN",
        progress_reporter=reporter,
        log_to_file=True,
        local_dir='/data/mohsen/ray_results/',
    ),
    param_space=config,
)
results = tuner.fit()

print("Best hyperparameters found were: ", results.get_best_result().config)
