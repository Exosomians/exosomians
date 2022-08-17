import numpy as np
from keras import backend as K
from keras.layers import Layer


class FixedDilatedCNN(Layer):
    """
    This layer takes a tensor in shape (n_samples, rna_length, 4) and returns a
    tensor in shape (n_samples, rna_length, n_structural_features)

    It uses a kind of dilated CNN with no trainable parameters to extract the
    RNA structural features

    Developed by Alireza Omidi
    """

    def __init__(self, dilation_start=2, dilation_end=None, n_bulge=2, **kwargs):
        self.dilation_start = dilation_start
        self.dilation_end = dilation_end
        self.n_bulge = n_bulge
        self.n_dilations = None
        super(FixedDilatedCNN, self).__init__(**kwargs)

    def build(self, input_shape):
        if self.dilation_end is None:
            self.dilation_end = int(input_shape[1] / 2 - 1)
        self.n_dilations = self.dilation_end - self.dilation_start
        super(FixedDilatedCNN, self).build(input_shape)

    def dilate(self, x, space, weights):
        weights = np.array(weights).reshape((1, 1, -1))
        x_padded_left = K.temporal_padding(x, (space+1, 0))
        x_padded_right = K.temporal_padding(x, (0, space+1))
        dilation = x_padded_right[:, space+1:, :] * x_padded_left[:, :-(space+1), ::-1]
        dilation = dilation * weights
        return K.sum(dilation, axis=2)
    
    def left_bulge_dilate(self, x, space, weights, bulge_size):
        x_padded_left = K.temporal_padding(x, (space+1+bulge_size, 0))
        x_padded_right = K.temporal_padding(x, (0, space+1))
        dilation = x_padded_right[:, space+1:, :] * x_padded_left[:, :-(space+1+bulge_size), ::-1]
        weights = np.array(weights).reshape((1, 1, -1))
        dilation = dilation * weights
        return K.sum(dilation, axis=2)
    
    def right_bulge_dilate(self, x, space, weights, bulge_size):
        x_padded_left = K.temporal_padding(x, (space+1, 0))
        x_padded_right = K.temporal_padding(x, (0, space+1+bulge_size))
        dilation = x_padded_right[:, space+1+bulge_size:, :] * x_padded_left[:, :-(space+1), ::-1]
        weights = np.array(weights).reshape((1, 1, -1))
        dilation = dilation * weights
        return K.sum(dilation, axis=2)
        
        
    def call(self, x):
        result = []

        for s in range(self.dilation_start, self.dilation_end):
            AU_CG = self.dilate(x, s, weights=[2, 3, 3, 2])
            GU = self.dilate(x[:, :, 2:], s, weights=[2, 2])
            result.append(AU_CG + GU)

            for d in range(1, self.n_bulge+1):
                # left bulges
                dilation = self.left_bulge_dilate(x, s, [2, 3, 3, 2], d)
                result.append(dilation)
                
            for d in range(1, self.n_bulge+1):
                # right bulges
                dilation = self.right_bulge_dilate(x, s, [2, 3, 3, 2], d)
                result.append(dilation)
 
        result = K.stack(result)
        result = K.permute_dimensions(result, (1, 2, 0))
        return result

    def compute_output_shape(self, input_shape):
        return (input_shape[0], input_shape[1], self.n_dilations * (self.n_bulge * 2 + 1))