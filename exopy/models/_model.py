import os
import json
import pickle
import random
from collections import defaultdict
from typing import Optional, List, Union

import numpy as np
import torch
import torch.nn as nn
from pytorch_lightning.callbacks import EarlyStopping, ModelCheckpoint, TQDMProgressBar
from torch.autograd import Variable
from torch.optim import lr_scheduler
from torch.utils.data import DataLoader
from tqdm import tqdm

from ._data import RNASeqDataset, SpecialTokens
from ._base import RNAClassifierBase
import pytorch_lightning as pl
import torch.nn.functional as F

import scanpy as sc

from ._module import ExoNetModule


class ExoNet(RNAClassifierBase):
    token_encoder: dict = {}
    module: ExoNetModule = None
    data: RNASeqDataset = None

    @classmethod
    def load(cls, path: str):
        cls.module = ExoNetModule.load_from_checkpoint(path)
        cls.token_encoder = cls.module.config.token_encoder
        return cls()

    def __init__(self,
                 network: str = 'exonet',
                 **kwargs):
        super().__init__()
        if ExoNet.module is not None:
            self.module = ExoNet.module
        else:
            pos_weight = self.data.weights[1]
            ExoNet.module = self.module = ExoNetModule(network=network,
                                                       token_encoder=self.token_encoder,
                                                       pos_weight=pos_weight,
                                                       **kwargs)

    @classmethod
    def setup_dataset(cls, path, seq_key, target_key, categorical_keys=[], continuous_keys=[], fraction=1.0):
        dataset = RNASeqDataset(path=path,
                                seq_key=seq_key,
                                target_key=target_key,
                                categorical_keys=categorical_keys,
                                continuous_keys=continuous_keys,
                                fraction=fraction,
                                )
        cls.n_tokens = dataset.n_tokens
        cls.token_encoder = dataset.c2i
        cls.data = dataset

    def get_embeddings(self, batch_size: int = 128, **kwargs) -> sc.AnnData:
        self.module.eval()
        dataset = self.data

        dataloader = DataLoader(dataset,
                                batch_size=batch_size,
                                pin_memory=False,
                                num_workers=0,
                                shuffle=False,
                                collate_fn=dataset.collate_fn)

        embeddings = []
        probs = []
        for batch in tqdm(dataloader):
            # for key, value in batch.items():
            #     if isinstance(batch[key], torch.Tensor):
            #         batch[key].to(self.module.device)
            #     else:
            #         batch[key] = torch.tensor(batch[key]).to(self.module.device)

            with torch.no_grad():
                logits, batch_embeddings = self.module.forward(batch, return_embeddings=True)
                batch_probs = torch.sigmoid(logits)

            embeddings.append(batch_embeddings.detach().cpu().numpy())
            probs.append(batch_probs.detach().cpu().numpy())

        adata = sc.AnnData(X=np.concatenate(embeddings, axis=0))
        adata.obs = self.data.df.copy()
        adata.obs['prob'] = np.concatenate(probs, axis=0).reshape(-1,)

        return adata

    def fit(self, max_epochs: int = 500,
            batch_size: int = 128,
            early_stopping_patience: int = 5,
            train_size: Optional[float] = 0.9,
            check_val_every_n_epoch: Optional[Union[int, float]] = 1,
            save_path: Optional[str] = None, **kwargs):
        train_dataloader, valid_dataloader = self.data.split(train_size=train_size, batch_size=batch_size)

        es_callback = EarlyStopping(monitor='val_loss', patience=early_stopping_patience, mode='min')
        if 'callbacks' not in kwargs.keys():
            kwargs['callbacks'] = []

        kwargs['callbacks'] += [es_callback]

        if save_path is not None:
            checkpoint_callback = ModelCheckpoint(
                dirpath=save_path,
                filename=f'{self.module.network}_best',
                save_last=True,
                verbose=False,
                monitor='val_loss',
                mode='min',
            )
            kwargs['callbacks'] += [checkpoint_callback]

        kwargs['callbacks'] += [TQDMProgressBar(refresh_rate=10)]

        self.trainer = pl.Trainer(max_epochs=max_epochs,
                                  accelerator='gpu' if torch.cuda.is_available() else 'cpu',
                                  check_val_every_n_epoch=check_val_every_n_epoch,
                                  enable_progress_bar=True,
                                  default_root_dir=save_path, **kwargs)

        self.trainer.fit(self.module, train_dataloader, valid_dataloader)

    def save(self, path: str) -> None:
        pass
