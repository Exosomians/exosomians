import numpy as np
import torch
from torch.autograd import Variable


class AttributeDict(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


def categorical_encoder(values):
    unique_values = list(sorted(np.unique(values)))
    encoder = {k: i for i, k in enumerate(unique_values)}

    return encoder, unique_values


def _get_cnn_output_shape(network, input_shape, flatten=True):
    bs = 1
    x = Variable(torch.rand(bs, *input_shape))
    output_feat = network(x)
    if flatten:
        n_size = output_feat.data.view(bs, -1).size(1)
        return n_size
    else:
        return output_feat.shape[1]

