from collections import defaultdict
from typing import Optional, Union

import numpy as np
import pandas as pd
import pytorch_lightning as pl
import scanpy as sc
import torch
from pytorch_lightning.callbacks import EarlyStopping, ModelCheckpoint, TQDMProgressBar
from torch.utils.data import DataLoader
from tqdm import tqdm

from ._base import RNAClassifierBase
from ._data import RNASeqCSVDataset
from ._module import ExoNetModule
from ..utils import load_fasta


class ExoNet(RNAClassifierBase):
    token_encoder: dict = {}
    module: ExoNetModule = None
    data: RNASeqCSVDataset = None

    @classmethod
    def load(cls, path: str, **kwargs):
        ExoNet.module = ExoNetModule.load_from_checkpoint(path, **kwargs)
        ExoNet.token_encoder = ExoNet.module.config.token_encoder
        return cls()

    def __init__(self, network: str = "exonet", **kwargs):
        super().__init__()
        if ExoNet.module is not None:
            print("Model successfully restored!")
            self.module = ExoNet.module
        else:
            pos_weight = self.data.weights[1]
            ExoNet.module = self.module = ExoNetModule(
                network=network,
                token_encoder=self.token_encoder,
                pos_weight=pos_weight,
                **kwargs,
            )

    @classmethod
    def setup_dataset(
        cls,
        path: str,
        seq_key: str,
        target_key: str,
        categorical_keys=[],
        continuous_keys=[],
        fraction=1.0,
    ):
        dataset = cls.prepare_data(
            path, seq_key, target_key, categorical_keys, continuous_keys, fraction
        )

        cls.n_tokens = dataset.n_tokens
        cls.token_encoder = dataset.c2i
        cls.data = dataset

    @staticmethod
    def prepare_data(
        path_or_df: Union[str, pd.DataFrame],
        seq_key: str,
        target_key: Optional[str] = None,
        categorical_keys=[],
        continuous_keys=[],
        fraction: Optional[float] = 1.0,
    ):
        if isinstance(path_or_df, str) and path_or_df.endswith(".fasta"):
            path_or_df = load_fasta(path_or_df, seq_key, target_key)
            categorical_keys = []
            continuous_keys = []

        return RNASeqCSVDataset(
            path_or_df=path_or_df,
            seq_key=seq_key,
            target_key=target_key,
            categorical_keys=categorical_keys,
            continuous_keys=continuous_keys,
            fraction=fraction,
        )

    def fit(
        self,
        max_epochs: int = 500,
        batch_size: int = 128,
        early_stopping_patience: int = 5,
        train_size: Optional[float] = 0.9,
        check_val_every_n_epoch: Optional[Union[int, float]] = 1,
        save_path: Optional[str] = None,
        **kwargs,
    ):
        train_dataloader, valid_dataloader = self.data.split(
            train_size=train_size, batch_size=batch_size
        )

        es_callback = EarlyStopping(
            monitor="val_loss", patience=early_stopping_patience, mode="min"
        )
        if "callbacks" not in kwargs.keys():
            kwargs["callbacks"] = []

        kwargs["callbacks"] += [es_callback]

        if save_path is not None:
            checkpoint_callback = ModelCheckpoint(
                dirpath=save_path,
                filename=f"{self.module.network}_best",
                save_last=True,
                verbose=False,
                monitor="val_loss",
                mode="min",
            )
            kwargs["callbacks"] += [checkpoint_callback]

        kwargs["callbacks"] += [TQDMProgressBar(refresh_rate=10)]

        self.trainer = pl.Trainer(
            max_epochs=max_epochs,
            accelerator="gpu" if torch.cuda.is_available() else "cpu",
            check_val_every_n_epoch=check_val_every_n_epoch,
            enable_progress_bar=True,
            default_root_dir=save_path,
            **kwargs,
        )

        self.trainer.fit(
            self.module,
            train_dataloaders=train_dataloader,
            val_dataloaders=valid_dataloader,
        )

    @torch.no_grad()
    def get_embeddings(
        self, data: Optional[RNASeqCSVDataset] = None, batch_size: int = 128, **kwargs
    ) -> sc.AnnData:
        self.module.eval()
        dataset = self.data if data is None else data

        dataloader = DataLoader(
            dataset,
            batch_size=batch_size,
            pin_memory=False,
            num_workers=0,
            shuffle=False,
            collate_fn=dataset.collate_fn,
        )

        embeddings = []
        probs = []
        for batch in tqdm(dataloader):
            for key in batch.keys():
                if isinstance(batch[key], torch.Tensor):
                    batch[key] = batch[key].to(self.module.device)
                elif isinstance(batch[key], list):
                    batch[key] = [x.to(self.module.device) for x in batch[key]]
                else:
                    print('WTF!')
            logits, batch_embeddings = self.module.forward(
                batch, return_embeddings=True
            )
            if logits.shape[1] == 1:
                batch_probs = torch.sigmoid(logits)
            else:
                batch_probs = torch.softmax(logits, dim=1)  # (batch_size, n_output)
                batch_probs = batch_probs[torch.arange(batch_probs.shape[0]), batch['label']]

            embeddings.append(batch_embeddings.detach().cpu().numpy())
            probs.append(batch_probs.detach().cpu().numpy())

        adata = sc.AnnData(X=np.concatenate(embeddings, axis=0))
        adata.obs = self.data.df.copy()
        adata.obs["prob"] = np.concatenate(probs, axis=0).reshape(
            -1,
        )

        return adata

    @torch.no_grad()
    def predict(
        self, data: Optional[RNASeqCSVDataset] = None, batch_size: Optional[int] = 32
    ) -> pd.DataFrame:
        self.module.eval()
        dataset = self.data if data is None else data

        dataloader = DataLoader(
            dataset,
            batch_size=batch_size,
            pin_memory=False,
            num_workers=0,
            shuffle=False,
            collate_fn=dataset.collate_fn,
        )

        df = defaultdict(list)
        df["sequence"] = list(dataset.sequences)
        probs = []
        for batch in tqdm(dataloader):
            logits = self.module.forward(batch, return_embeddings=False)
            batch_probs = torch.sigmoid(logits)

            probs.append(batch_probs.detach().cpu().numpy())

        df["prob"] = list(
            np.concatenate(probs, axis=0).reshape(
                -1,
            )
        )

        return pd.DataFrame(df)

    def save(self, path: str) -> None:
        pass
