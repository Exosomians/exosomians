from typing import List, Tuple

import torch
import torch.nn as nn
from torch.autograd import Variable


class MLP(nn.Module):
    def __init__(self, n_input: int,
                 n_hidden: int,
                 n_layers: int,
                 n_output: int,
                 use_batch_norm: bool,
                 use_layer_norm: bool,
                 dropout_rate: float,
                 activation_fn: nn.Module,
                 ):
        super().__init__()
        self.n_input = n_input
        self.n_hidden = n_hidden
        self.n_layers = n_layers
        self.use_batch_norm = use_batch_norm
        self.use_layer_norm = use_layer_norm
        self.dropout_rate = dropout_rate
        self.activation_fn = activation_fn

        layers = [n_input] + [n_hidden for _ in range(n_layers)]

        self.network = nn.ModuleList()
        for n_in, n_out in zip(layers[:-1], layers[1:]):
            self.network.append(
                nn.Linear(n_in, n_out, bias=True)
            )

            if self.use_batch_norm:
                self.network.append(
                    nn.BatchNorm1d(n_out)
                )

            if self.use_layer_norm:
                self.network.append(
                    nn.LayerNorm(n_out)
                )

            self.network.append(
                self.activation_fn()
            )

            if 1.0 > self.dropout_rate > 0.0:
                self.network.append(
                    nn.Dropout(p=dropout_rate)
                )

        self.classifier = nn.Linear(n_hidden, n_output)

        self.network = nn.Sequential(*self.network)

    def forward(self, x, return_embeddings=False):
        """
            x: (batch_size, n_input)
        """
        x = self.network(x)
        y = self.classifier(x)
        if return_embeddings:
            return y, x
        return y


class RNNCell(nn.Module):
    def __init__(self, n_input: int,
                 n_layers: int,
                 n_hidden: int,
                 rnn_cell: str = 'GRU',
                 dropout_rate: float = 0.0,
                 bidirectional: bool = False):
        assert rnn_cell.upper() in ['GRU', 'LSTM']
        super().__init__()

        self.n_input = n_input
        self.n_layers = n_layers
        self.n_hidden = n_hidden
        self.dropout_rate = dropout_rate
        self.bidirectional = bidirectional
        self.rnn_cell = rnn_cell.upper()

        if self.rnn_cell == 'GRU':
            self.rnn = nn.GRU(input_size=self.n_input,
                              hidden_size=self.n_hidden,
                              num_layers=self.n_layers,
                              batch_first=True,
                              bidirectional=self.bidirectional,
                              dropout=self.dropout_rate)
        else:
            self.rnn = nn.LSTM(input_size=self.n_input,
                               hidden_size=self.n_hidden,
                               num_layers=self.n_layers,
                               batch_first=True,
                               bidirectional=self.bidirectional,
                               dropout=self.dropout_rate)

    def forward(self, x, hidden=None):
        """
            x: (batch_size, seq_len, n_input)

        """
        return self.rnn(x, hidden)


class ConvBlock(nn.Module):
    def __init__(self,
                 in_channels: int,
                 n_filters: int,
                 n_layers: int = 3,
                 kernel_size: int = 3,
                 use_batch_norm: bool = False,
                 use_layer_norm: bool = False,
                 dropout_rate: float = 0.0,
                 pooling: str = 'max',
                 pooling_size: int = 2,
                 activation_fn: nn.Module = nn.ReLU,
                 ):
        assert pooling.lower() in ['max', 'avg']
        super().__init__()
        self.in_channels = in_channels
        self.n_filters = n_filters
        self.n_layers = n_layers
        self.kernel_size = kernel_size
        self.use_batch_norm = use_batch_norm
        self.use_layer_norm = use_layer_norm
        self.dropout_rate = dropout_rate
        self.pooling = nn.MaxPool1d if pooling == 'max' else nn.AvgPool1d
        self.pooling_size = pooling_size
        self.activation_fn = activation_fn

        self.network = nn.ModuleList()
        for i in range(self.n_layers):
            if i == 0:
                in_channels = self.in_channels
            else:
                in_channels = self.n_filters

            self.network.append(
                nn.Conv1d(in_channels=in_channels,
                          out_channels=self.n_filters,
                          kernel_size=kernel_size))

            if self.use_batch_norm:
                self.network.append(nn.BatchNorm1d(self.n_filters))

            if self.use_layer_norm:
                self.network.append(nn.GroupNorm(1, self.n_filters))

            self.network.append(self.activation_fn())

            if 1.0 > self.dropout_rate > 0.0:
                self.network.append(nn.Dropout(p=self.dropout_rate))

        self.network.append(self.pooling(self.pooling_size))

        self.network = nn.Sequential(*self.network)

    def forward(self, x):
        """
            x: (batch_size, seq_len, in_channels)
        """

        return self.network(x)
