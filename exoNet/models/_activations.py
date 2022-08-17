from keras.layers import Activation
from keras.layers.advanced_activations import LeakyReLU


ACTIVATIONS = {
    "relu": Activation("relu"),
    'leaky_relu': LeakyReLU(),
    'linear': Activation("linear"),
    'sigmoid': Activation('sigmoid'),
}
