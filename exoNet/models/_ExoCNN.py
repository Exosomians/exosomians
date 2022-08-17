import math
import os

import keras
import numpy as np
from keras.callbacks import EarlyStopping, ReduceLROnPlateau
from keras.layers import Input, Dense, BatchNormalization, Dropout, Conv1D, MaxPooling1D, Flatten
from keras.models import Model, load_model
from keras.optimizers import Nadam, Adam
from keras.utils import to_categorical
from keras.layers import ReLU
from keras.utils import plot_model

from exoNet.models._losses import METRICS, LOSSES
from exoNet.models._network import Network
from exoNet.utils import label_encoder, train_test_split_data


class ExoCNN(Network):
    def __init__(self, seq_len, n_channels, n_classes=2, **kwargs):
        super().__init__()
        self.seq_len = seq_len
        self.n_channels = n_channels
        self.n_classes = n_classes
        self.model_name = ExoCNN.__name__

        self.lr = kwargs.get("learning_rate", 0.001)
        self.dr_rate = kwargs.get("dropout_rate", 0.2)
        self.model_path = kwargs.get("model_path", "./models/FCN/")
        self.loss_fn = kwargs.get("loss_fn", 'cce')
        self.lambda_l1 = kwargs.get('lambda_l1', 1e-6)
        self.lambda_l2 = kwargs.get('lambda_l2', 1e-6)
        self.use_batchnorm = kwargs.get("use_batchnorm", False)
        self.padding = kwargs.get("padding", "valid")
        self.class_weight = kwargs.get("class_weight", {0: 1.0, 1:1.0})

        self.sequence = Input(shape=(self.seq_len, self.n_channels,), name="data")

        self.network_kwargs = {
            "seq_len": self.seq_len,
            "n_channels": self.n_channels,
            "n_classes": self.n_classes,
            "dropout_rate": self.dr_rate,
            "loss_fn": self.loss_fn,
            "use_batchnorm": self.use_batchnorm,
        }

        self.training_kwargs = {
            "learning_rate": self.lr,
            "model_path": self.model_path,
            "lambda_l1": self.lambda_l1,
            "lambda_l2": self.lambda_l2,
        }

        self.aux_models = {}

        self.init_w = keras.initializers.glorot_normal()
        self.regularizer = keras.regularizers.l1_l2(self.lambda_l1, self.lambda_l2)
        self._create_network()
        self._compile_models()

        print_summary = kwargs.get("print_summary", False)

        if print_summary:
            self.model.summary()

    def _create_network(self):
        conv = Conv1D(filters=32, kernel_size=10, padding=self.padding, kernel_initializer=self.init_w, use_bias=False,
                      kernel_regularizer=self.regularizer)(self.sequence)
        if self.use_batchnorm:
            conv = BatchNormalization(trainable=True)(conv)
        conv = ReLU()(conv)

        conv = Conv1D(filters=32, kernel_size=3, padding=self.padding, kernel_initializer=self.init_w, use_bias=False,
                      kernel_regularizer=self.regularizer)(conv)
        if self.use_batchnorm:
            conv = BatchNormalization(trainable=True)(conv)
        conv = ReLU()(conv)

        max_pool = MaxPooling1D(pool_size=2)(conv)

        flat = Flatten()(max_pool)

        dense = Dense(64, kernel_initializer=self.init_w, kernel_regularizer=self.regularizer, use_bias=False)(flat)
        if self.use_batchnorm:
            dense = BatchNormalization()(dense)
        dense = ReLU()(dense)

        probs = Dense(self.n_classes, activation='softmax', kernel_initializer=self.init_w, kernel_regularizer=self.regularizer)(dense)

        self.model = Model(inputs=self.sequence, outputs=probs)
        self.aux_models['latent'] = Model(inputs=self.sequence, outputs=dense)


    def _compile_models(self):
        self.optimizer = Nadam(lr=self.lr)
        self.model.compile(optimizer=self.optimizer, loss=LOSSES[self.loss_fn],
                           metrics=['acc', METRICS['sensitivity'], METRICS['specificity']])

    def to_latent(self, data):
        latent_model = Model(self.model.input, self.model.layers[-2].output)
        return latent_model.predict(data)

    def predict(self, data):
        return self.label_encoder.inverse_transform(np.argmax(self.model.predict(data), axis=1))

    def save_model(self):
        os.makedirs(self.model_path, exist_ok=True)
        self.model.save(os.path.join(self.model_path, f"{self.model_name}.h5"), overwrite=True)
    
    def save_model_weights(self):
        os.makedirs(self.model_path, exist_ok=True)
        self.model.save_weights(os.path.join(self.model_path, f"weights.h5"), overwrite=True)

    def restore_model(self):
        self.model = load_model(os.path.join(self.model_path, f"{self.model_name}.h5"), compile=False)
        self._compile_models()
    
    def restore_model_weights(self):
        self.model.load_weights(os.path.join(self.model_path, f"weights.h5"))
        self._compile_models()
        print("model has been successfully restored")
    
    def plot_ExoCNN(self):
        os.makedirs(self.model_path, exist_ok=True)
        plot_model(self.model, os.path.join(self.model_path, "architecture.pdf"), show_layer_names=False, show_shapes=True)

    def train(self, seq_data, labels, le=None, n_epochs=500, batch_size=32, early_stopping_kwargs={},
              lr_reducer_kwargs={}, verbose=2, save=False, retrain=True):

        x_train, x_valid, y_train, y_valid = train_test_split_data(seq_data, labels, 0.80, stratify=True)

        y_train, self.label_encoder = label_encoder(y_train, label_encoder=le)
        y_train = to_categorical(y_train, num_classes=self.n_classes)

        y_valid, self.label_encoder = label_encoder(y_valid, label_encoder=le)
        y_valid = to_categorical(y_valid, num_classes=self.n_classes)

        if not retrain and os.path.exists(os.path.join(self.model_path, f"weights.h5")):
            self.restore_model_weights()
            return self.model.evaluate(x_valid, y_valid, verbose=0)

        callbacks = []

        if early_stopping_kwargs != {}:
            callbacks.append(EarlyStopping(**early_stopping_kwargs))

        if lr_reducer_kwargs != {}:
            callbacks.append(ReduceLROnPlateau(**lr_reducer_kwargs))

        self.model.fit(x=x_train,
                       y=y_train,
                       validation_data=(x_valid, y_valid),
                       epochs=n_epochs,
                       class_weight=self.class_weight,
                       batch_size=batch_size,
                       verbose=verbose,
                       callbacks=callbacks,
                       )
        if save:
            self.save_model_weights()

        return self.model.evaluate(x_valid, y_valid, verbose=0)
