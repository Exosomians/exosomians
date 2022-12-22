import numpy as np
import pytorch_lightning as pl
import torch
from torch import nn
from torch.optim import lr_scheduler
from torchmetrics import AUROC, Precision, Recall, Specificity

from ._base_modules import ConvBlock, MLP, RNNCell
from ._constants import ExoNetCONSTANTS
from ._data import SpecialTokens
from ._utils import AttributeDict, _get_cnn_output_shape

activations = {
    'relu': nn.ReLU,
    'leaky_relu': nn.LeakyReLU,
    'tanh': nn.Tanh,
}


class ExoNetModule(pl.LightningModule):
    def __init__(self, network, **network_params):
        assert network in ['exocnn', 'exogru', 'exolstm', 'deepbind', 'exonet']
        super().__init__()

        self.network = network
        self.config = AttributeDict()

        self.set_config(network, **network_params)
        self.use_cnn = network in ['exocnn', 'deepbind']
        self.use_lstm = network in ['exolstm'] or network in ['exonet'] and self.config.rnn_cell in ['LSTM']

        self.pos_weight = network_params.get('pos_weight', 1.0)

        self.save_hyperparameters(network, self.config)

        self.input_embedding = nn.Embedding(self.config.n_tokens, self.config.n_tokens)
        self.input_embedding.weight.data.copy_(torch.eye(self.config.n_tokens))
        self.input_embedding.weight.requires_grad = False

        if network in ['exocnn', 'deepbind']:
            self.cnn = nn.ModuleList()
            if isinstance(self.config.kernel_size, int):
                kernel_sizes = [self.config.kernel_size] * self.config.n_conv_blocks
            else:
                kernel_sizes = self.config.kernel_size
            for i in range(self.config.n_conv_blocks):
                if i == 0:
                    in_channels = self.config.n_tokens
                else:
                    in_channels = self.config.n_filters
                self.cnn.append(
                    ConvBlock(in_channels=in_channels,
                              n_filters=self.config.n_filters,
                              n_layers=self.config.n_conv_layers,
                              pooling=self.config.pooling,
                              kernel_size=kernel_sizes[i],
                              pooling_size=self.config.pooling_size,
                              use_batch_norm=self.config.use_batch_norm,
                              use_layer_norm=self.config.use_layer_norm,
                              dropout_rate=self.config.dropout_rate,
                              activation_fn=activations[self.config.activation_fn],
                              )
                )
            self.cnn = nn.Sequential(*self.cnn)

            self.head = MLP(
                n_input=_get_cnn_output_shape(self.cnn, (self.config.n_tokens, self.config.max_len)),
                n_layers=self.config.n_head_layers,
                n_hidden=self.config.n_head_hidden,
                n_output=self.config.n_output,
                use_batch_norm=self.config.use_batch_norm,
                use_layer_norm=self.config.use_layer_norm,
                dropout_rate=self.config.dropout_rate,
                activation_fn=activations[self.config.activation_fn],
            )

        elif network in ['exolstm', 'exogru']:
            self.rnn = RNNCell(n_input=self.config.n_tokens,
                               n_layers=self.config.n_layers,
                               n_hidden=self.config.n_hidden,
                               rnn_cell=self.config.rnn_cell,
                               bidirectional=self.config.bidirectional
                               )

            self.head = MLP(
                n_input=2 * self.config.n_layers * self.config.n_hidden if self.config.bidirectional else self.config.n_hidden * self.config.n_layers,
                n_hidden=self.config.n_head_hidden,
                n_layers=self.config.n_head_layers,
                n_output=self.config.n_output if self.config.n_output > 2 else 1,
                use_batch_norm=self.config.use_batch_norm,
                use_layer_norm=self.config.use_layer_norm,
                dropout_rate=self.config.dropout_rate,
                activation_fn=activations[self.config.activation_fn],
            )
        elif network in ['exonet']:
            self.cnn = nn.ModuleList()
            if isinstance(self.config.kernel_size, int):
                kernel_sizes = [self.config.kernel_size] * self.config.n_conv_blocks
            else:
                kernel_sizes = self.config.kernel_size
            for i in range(self.config.n_conv_blocks):
                if i == 0:
                    in_channels = self.config.n_tokens
                else:
                    in_channels = self.config.n_filters
                self.cnn.append(
                    ConvBlock(in_channels=in_channels,
                              n_filters=self.config.n_filters,
                              n_layers=self.config.n_conv_layers,
                              pooling=self.config.pooling,
                              kernel_size=kernel_sizes[i],
                              pooling_size=self.config.pooling_size,
                              use_batch_norm=self.config.use_batch_norm,
                              use_layer_norm=self.config.use_layer_norm,
                              dropout_rate=self.config.dropout_rate,
                              activation_fn=activations[self.config.activation_fn],
                              )
                )
            self.cnn = nn.Sequential(*self.cnn)

            self.rnn = RNNCell(
                n_input=self.config.n_filters,
                n_layers=self.config.n_layers,
                n_hidden=self.config.n_hidden,
                rnn_cell=self.config.rnn_cell,
                bidirectional=self.config.bidirectional
            )

            self.head = MLP(
                n_input=2 * self.config.n_layers * self.config.n_hidden if self.config.bidirectional else self.config.n_hidden * self.config.n_layers,
                n_hidden=self.config.n_head_hidden,
                n_layers=self.config.n_head_layers,
                n_output=self.config.n_output if self.config.n_output > 2 else 1,
                use_batch_norm=self.config.use_batch_norm,
                use_layer_norm=self.config.use_layer_norm,
                dropout_rate=self.config.dropout_rate,
                activation_fn=activations[self.config.activation_fn],
            )

        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

        self.loss = nn.CrossEntropyLoss() if self.config.n_output > 2 else nn.BCEWithLogitsLoss(
            pos_weight=torch.tensor(self.pos_weight, device=device))

        self.eval_metrics = {
            'auroc': AUROC(pos_label=1).to(device),
            'precision': Precision().to(device),
            'recall': Recall().to(device),
            'specificity': Specificity().to(device),
        }

    def set_config(self, network, **config):
        if network == 'exocnn':
            self.config = AttributeDict({
                'n_conv_blocks': config.get('n_conv_blocks', 3),
                'n_conv_layers': config.get('n_conv_layers', 3),
                'kernel_size': config.get('kernel_size', [10, 5, 3]),
                'n_filters': config.get('n_filters', 32),
                'pooling': config.get('pooling', 'max'),
                'pooling_size': config.get('pooling_size', 2),
                'activation_fn': config.get('activation_fn', 'relu'),
                'use_batch_norm': config.get('use_batch_norm', False),
                'use_layer_norm': config.get('use_layer_norm', False),
                'dropout_rate': config.get('dropout_rate', 0.1),
                'n_head_layers': config.get('n_head_layers', 3),
                'n_head_hidden': config.get('n_head_hidden', 256),
            })
        elif network == 'deepbind':
            self.config = AttributeDict({
                'n_conv_blocks': config.get('n_conv_blocks', 1),
                'n_conv_layers': config.get('n_conv_layers', 1),
                'kernel_size': config.get('kernel_size', 5),
                'n_filters': config.get('n_filters', 32),
                'pooling': config.get('pooling', 'avg'),
                'pooling_size': config.get('pooling_size', 2),
                'activation_fn': config.get('activation_fn', 'relu'),
                'use_batch_norm': config.get('use_batch_norm', True),
                'use_layer_norm': config.get('use_layer_norm', False),
                'dropout_rate': config.get('dropout_rate', 0.1),
                'n_head_layers': config.get('n_head_layers', 3),
                'n_head_hidden': config.get('n_head_hidden', 256),
            })
        elif network in ['exolstm', 'exogru']:
            self.config = AttributeDict({
                'n_hidden': config.get('n_hidden', 128),
                'n_layers': config.get('n_layers', 3),
                'rnn_cell': config.get('network', self.network)[3:].upper(),
                'bidirectional': config.get('bidirectional', True),
                'activation_fn': config.get('activation_fn', 'relu'),
                'use_batch_norm': config.get('use_batch_norm', True),
                'use_layer_norm': config.get('use_layer_norm', False),
                'dropout_rate': config.get('dropout_rate', 0.1),
                'n_head_layers': config.get('n_head_layers', 3),
                'n_head_hidden': config.get('n_head_hidden', 256),
            })
        elif network in ['exonet']:
            self.config = AttributeDict({
                'n_conv_blocks': config.get('n_conv_blocks', 3),
                'n_conv_layers': config.get('n_conv_layers', 3),
                'kernel_size': config.get('kernel_size', [10, 5, 3]),
                'n_filters': config.get('n_filters', 32),
                'pooling': config.get('pooling', 'max'),
                'pooling_size': config.get('pooling_size', 2),
                'activation_fn': config.get('activation_fn', 'relu'),
                'use_batch_norm': config.get('use_batch_norm', False),
                'use_layer_norm': config.get('use_layer_norm', False),
                'dropout_rate': config.get('dropout_rate', 0.1),
                'n_head_layers': config.get('n_head_layers', 3),
                'n_head_hidden': config.get('n_head_hidden', 256),

                'n_hidden': config.get('n_hidden', 128),
                'n_layers': config.get('n_layers', 3),
                'rnn_cell': config.get('rnn_cell', 'GRU'),
                'bidirectional': config.get('bidirectional', False),
            })

        self.config['token_encoder'] = config.get('token_encoder')
        self.config['n_tokens'] = len(self.config.token_encoder)
        self.config['max_len'] = config.get('max_len', 512)
        self.config['n_output'] = config.get('n_output', 1)
        self.config['lr'] = config.get('lr', 1e-4)

    # def mixup_data(self, batch, alpha: float = 1.0):
    #     """
    #         Returns mixed inputs, pairs of targets, and lambda
    #     """
    #     alpha = max(0.0, alpha)
    #
    #     if alpha == 0.0:
    #         mixup_lambda = 1.0
    #     else:
    #         mixup_lambda = np.random.beta(alpha, alpha)
    #
    #     x = batch[ExoNetCONSTANTS.SEQ_KEY]
    #     y = batch[ExoNetCONSTANTS.Y_KEY]
    #
    #     batch_size = x.size()[0]
    #     index = torch.randperm(batch_size).to(x.device)
    #
    #     mixed_x = mixup_lambda * x + (1. - mixup_lambda) * x[index, :]
    #
    #     batch[ExoNetCONSTANTS.SEQ_KEY] = mixed_x
    #     batch[GeneCPA_REGISTRY_KEYS.X_KEY + '_true'] = x
    #     batch[GeneCPA_REGISTRY_KEYS.GENE_ADV_KEY + '_mixup'] = y_knockouts[index]
    #     batch[GeneCPA_REGISTRY_KEYS.GENE_KO_KEY] = mixed_knockouts
    #     batch[GeneCPA_REGISTRY_KEYS.GENE_KO_KEY + '_true'] = knockouts
    #
    #     for covar, encoder in self.covars_encoder.items():
    #         batch[covar + '_mixup'] = batch[covar][index]
    #
    #     return batch, mixup_lambda

    def forward(self, tensors, return_embeddings=False):
        x = [self.input_embedding(x_i) for x_i in tensors[ExoNetCONSTANTS.SEQ_KEY]]
        x = nn.utils.rnn.pack_sequence(x, enforce_sorted=False)

        if self.network in ['exonet']:
            x, _ = nn.utils.rnn.pad_packed_sequence(x,
                                                    batch_first=True,
                                                    padding_value=self.config.token_encoder[SpecialTokens.pad],
                                                    total_length=self.config.max_len)
            x = x.permute(0, 2, 1)  # (batch_size, in_channels, max_len)
            x = self.cnn(x)  # (batch_size, n_filters, seq_len)
            x = x.permute(0, 2, 1)
            if self.use_lstm:
                _, (h, _) = self.rnn(x)
            else:
                _, h = self.rnn(x)

            bs = h.shape[1]
            h = h.view(bs, -1)
            head_outputs = self.head(h, return_embeddings=return_embeddings)

        elif self.use_cnn:
            x, _ = nn.utils.rnn.pad_packed_sequence(x,
                                                    batch_first=True,
                                                    padding_value=self.config.token_encoder[SpecialTokens.pad],
                                                    total_length=self.config.max_len)
            x = x.permute(0, 2, 1)  # (batch_size, in_channels, max_len)
            x = self.cnn(x)
            x = x.view(x.shape[0], -1)
            head_outputs = self.head(x, return_embeddings=return_embeddings)

        else:
            if self.use_lstm:
                _, (h, _) = self.rnn(x)
            else:
                _, h = self.rnn(x)

            bs = h.shape[1]
            h = h.view(bs, -1)

            head_outputs = self.head(h, return_embeddings=return_embeddings)

        if return_embeddings:
            y, x = head_outputs
            return y, x
        return head_outputs

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(list(filter(lambda p: p.requires_grad, self.parameters())), lr=self.config.lr)
        scheduler = lr_scheduler.ExponentialLR(optimizer, gamma=0.9)
        return [optimizer], [scheduler]

    def training_step(self, batch, batch_idx):
        y_pred = self.forward(batch)
        y_true = batch[ExoNetCONSTANTS.Y_KEY]
        batch_size = y_pred.shape[0]

        loss = self.loss(y_pred.squeeze(), y_true.float())
        results = {
            key: f(y_pred.squeeze(), y_true).item() for key, f in self.eval_metrics.items()
        }
        results['loss'] = loss

        for key in results.keys():
            if key != 'loss':
                self.log(f'{key}', results[key], prog_bar=False, on_step=True, batch_size=batch_size)

        return results

    def validation_step(self, batch, batch_idx):
        y_pred = self.forward(batch)
        y_true = batch[ExoNetCONSTANTS.Y_KEY]

        loss = self.loss(y_pred.squeeze(), y_true.float())
        results = {
            key: f(y_pred.squeeze(), y_true).item() for key, f in self.eval_metrics.items()
        }
        results['loss'] = loss.item()

        return results

    def validation_epoch_end(self, outputs):
        for key in outputs[0].keys():
            epoch_val = float(np.mean([output[key] for output in outputs]))
            self.log(f'val_{key}', epoch_val, prog_bar=True, on_epoch=True)
