# Inspired by MOSES Github repository (https://github.com/molecularsets/moses/blob/master/moses/utils.py)
from typing import List, Union, Optional

import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset, DataLoader
from torch.utils.data import SubsetRandomSampler

import logging
import multiprocessing

from tqdm import tqdm

from exopy.models._constants import ExoNetCONSTANTS
from exopy.models._utils import categorical_encoder

logger = logging.getLogger()
logger.setLevel(logging.INFO)


class SpecialTokens:
    bos = '<bos>'
    eos = '<eos>'
    pad = '<pad>'
    unk = '<unk>'


class RNASeqCSVDataset(Dataset):
    def __init__(self,
                 path_or_df: Union[str, pd.DataFrame],
                 seq_key: str,
                 target_key: Optional[str] = None,
                 categorical_keys: Optional[Union[List[str]]] = [],
                 continuous_keys: Optional[Union[List[str]]] = [],
                 fraction: float = 1.0,
                 **kwargs):

        if isinstance(categorical_keys, str):
            categorical_keys = [categorical_keys]

        if isinstance(continuous_keys, str):
            continuous_keys = [continuous_keys]

        self.seq_key = seq_key
        self.target_key = target_key
        self.cat_keys = categorical_keys
        self.cont_keys = continuous_keys

        ExoNetCONSTANTS.SEQ_KEY = seq_key
        ExoNetCONSTANTS.Y_KEY = target_key
        ExoNetCONSTANTS.CAT_COVS_KEYS = categorical_keys
        ExoNetCONSTANTS.CONT_COVS_KEYS = continuous_keys

        if isinstance(path_or_df, str):
            self.df = pd.read_csv(path_or_df, **kwargs)
        else:
            self.df = path_or_df

        if 0.0 < fraction < 1.0:
            indices = np.random.choice(np.arange(len(self.df)), int(fraction * len(self.df)))
            print(len(indices))
            self.df = self.df.iloc[indices, :]

        self.sequences = list(self.df[seq_key].values)
        self.target_encoder, self.unique_target_values = {'IC': 0, 'EV': 1}, ['IC', 'EV']

        if target_key is not None:
            self.targets = self.df[target_key].values

            weights = []
            n_total_samples = len(self.df)
            for target in self.unique_target_values:
                n_samples = len(self.df[self.df[target_key] == target])
                weights.append((n_total_samples / n_samples) - 1.)
            self.weights = weights
        else:
            self.targets = None
            self.weights = None

        self.cat_info = {}
        for cat_key in categorical_keys:
            cat_encoder, unique_cat_values = categorical_encoder(self.df[cat_key].values)
            self.cat_info[cat_key] = {
                'encoder': cat_encoder,
                'unique_values': unique_cat_values
            }

        self.cont_info = {}
        for cont_key in continuous_keys:
            self.cont_info[cont_key] = {}  # :D

        self.nucs = list('ACGU')
        self.special_tokens = [SpecialTokens.bos, SpecialTokens.eos, SpecialTokens.pad, SpecialTokens.unk]
        self.c2i = {}
        for i, token in enumerate(self.nucs + self.special_tokens):
            self.c2i[token] = i

        self.i2c = {v: k for k, v in self.c2i.items()}

        self.n_tokens = len(self.c2i)

        print(f'Successfully loaded the dataset with {len(self.df)} sequences!')

    @property
    def bos(self):
        return self.c2i[SpecialTokens.bos]

    @property
    def eos(self):
        return self.c2i[SpecialTokens.eos]

    @property
    def pad(self):
        return self.c2i[SpecialTokens.pad]

    @property
    def unk(self):
        return self.c2i[SpecialTokens.unk]

    def seq2ids(self, seq: str, add_bos=False, add_eos=False):
        encoded_drug = []
        if add_bos:
            encoded_drug.append(self.bos)

        encoded_drug += [self.c2i.get(nuc, self.unk) for nuc in seq]

        if add_eos:
            encoded_drug.append(self.eos)

        return encoded_drug

    def ids2seq(self, ids, rem_bos=True, rem_eos=True):
        if len(ids) == 0:
            return ''
        if rem_bos and ids[0] == self.bos:
            ids = ids[1:]
        if rem_eos and ids[-1] == self.eos:
            ids = ids[:-1]

        seq = ''.join([self.i2c.get(idx) for idx in ids])

        return seq

    def __len__(self):
        return len(self.df)

    def collate_fn(self, data):
        data.sort(key=lambda x: len(x[self.seq_key]), reverse=True)
        collated_data = {
            self.seq_key: [x[self.seq_key] for x in data],
        }

        if self.target_key is not None:
            collated_data[self.target_key] = torch.tensor([x[self.target_key] for x in data], dtype=torch.long)

        for key in self.cat_keys:
            collated_data[key] = torch.tensor([x[key] for x in data], dtype=torch.long)

        for key in self.cont_keys:
            collated_data[key] = torch.tensor([x[key] for x in data], dtype=torch.float)

        return collated_data

    def __getitem__(self, index):
        encoded_seq = torch.tensor(self.seq2ids(self.sequences[index], add_bos=True, add_eos=True),
                                   dtype=torch.long)
        item = {self.seq_key: encoded_seq}
        if self.target_key is not None:
            item[self.target_key] = self.target_encoder[self.targets[index]]

        for key in self.cat_keys:
            item[key] = self.cat_info[key]['encoder'][self.df[key].values[index]]

        for key in self.cont_keys:
            item[key] = self.df[key].values[index]

        return item

    def split(self, train_size: float, batch_size: int):
        n_samples = len(self.sequences)

        indices = np.arange(n_samples)
        np.random.shuffle(indices)

        train_indices = indices[:int(train_size * n_samples)]
        valid_indices = indices[int(train_size * n_samples):]

        train_dataloader = DataLoader(self, sampler=SubsetRandomSampler(train_indices),
                                      num_workers=max(1, min(6, multiprocessing.cpu_count() // 2)), pin_memory=False,
                                      batch_size=batch_size, collate_fn=self.collate_fn)
        valid_dataloader = DataLoader(self, sampler=SubsetRandomSampler(valid_indices),
                                      num_workers=max(1, min(6, multiprocessing.cpu_count() // 2)), pin_memory=False,
                                      batch_size=batch_size, collate_fn=self.collate_fn)

        return train_dataloader, valid_dataloader
