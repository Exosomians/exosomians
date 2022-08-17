import os

import keras
import numpy as np
from keras.callbacks import EarlyStopping, ReduceLROnPlateau
from keras.layers import Input, Dense, BatchNormalization, Dropout, Conv1D, MaxPooling1D, Flatten, concatenate
from keras.layers.advanced_activations import LeakyReLU
from keras.models import Model, load_model
from keras.optimizers import Nadam
from keras.utils import to_categorical

from exoNet.models._losses import METRICS, LOSSES
from exoNet.models._network import Network
from exoNet.utils import label_encoder, train_test_split_data


class ExoFCNN(Network):
    def __init__(self, n_features, seq_len, n_channels, n_classes=2, **kwargs):
        super().__init__()
        self.n_features = n_features
        self.seq_len = seq_len
        self.n_channels = n_channels
        self.n_classes = n_classes
        self.model_name = ExoFCNN.__name__

        self.lr = kwargs.get("learning_rate", 0.001)
        self.dr_rate = kwargs.get("dropout_rate", 0.2)
        self.model_path = kwargs.get("model_path", "./models/exoFCNN/")
        self.loss_fn = kwargs.get("loss_fn", 'cce')
        self.lambda_l1 = kwargs.get('lambda_l1', 0.0)
        self.lambda_l2 = kwargs.get('lambda_l2', 0.0)
        self.use_batchnorm = kwargs.get("use_batchnorm", False)

        self.sequence = Input(shape=(self.seq_len, self.n_channels,), name="sequences")
        self.mlp_x = Input(shape=(self.n_features,), name="other_features")

        self.network_kwargs = {
            "n_features": self.n_features,
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

        self.init_w = keras.initializers.glorot_normal()
        self.regularizer = keras.regularizers.l1_l2(self.lambda_l1, self.lambda_l2)
        self._create_network()
        self._compile_models()

        print_summary = kwargs.get("print_summary", True)

        if print_summary:
            self.model.summary()

    def _create_network(self):
        def cnn():
            h = Conv1D(filters=32, kernel_size=8, activation='linear', padding='same',
                       kernel_initializer=self.init_w, kernel_regularizer=self.regularizer)(self.sequence)
            if self.use_batchnorm:
                h = BatchNormalization(trainable=True)(h)
            h = LeakyReLU()(h)
            h = MaxPooling1D(pool_size=2)(h)
            if self.dr_rate > 0:
                h = Dropout(self.dr_rate)(h)
            h = Flatten()(h)
            return h

        def fcn():
            h = Dense(64, activation='linear', kernel_initializer=self.init_w,
                      kernel_regularizer=self.regularizer)(self.mlp_x)
            if self.use_batchnorm:
                h = BatchNormalization(trainable=True)(h)
            h = LeakyReLU()(h)
            if self.dr_rate > 0:
                h = Dropout(self.dr_rate)(h)
            return h

        h = concatenate([cnn(), fcn()], axis=1)

        h = Dense(32, activation='linear', kernel_initializer=self.init_w,
                  kernel_regularizer=self.regularizer)(h)
        if self.use_batchnorm:
            h = BatchNormalization()(h)
        h = LeakyReLU()(h)
        if self.dr_rate > 0:
            h = Dropout(self.dr_rate)(h)

        output = Dense(self.n_classes, kernel_initializer=self.init_w, kernel_regularizer=self.regularizer,
                       activation='softmax')(h)

        self.model = Model(inputs=[self.sequence, self.mlp_x], outputs=output, name='exoFCNN')

    def _compile_models(self):
        self.optimizer = Nadam(lr=self.lr)
        self.model.compile(optimizer=self.optimizer, loss=LOSSES[self.loss_fn],
                           metrics=['acc', METRICS['sensitivity'], METRICS['specificity']])

    def to_latent(self):
        pass

    def predict(self, data):
        return self.label_encoder.inverse_transform(np.argmax(self.model.predict(data), axis=1))

    def save_model(self):
        self.model.save(os.path.join(self.model_path, f"{self.model_name}.h5"), overwrite=True)

    def restore_model(self):
        self.model = load_model(os.path.join(self.model_path, f"{self.model_name}.h5"), compile=False)
        self._compile_models()

    def train(self, seq_data, fcn_data, labels, le=None, n_epochs=500, batch_size=32, early_stopping_kwargs={},
              lr_reducer_kwargs={}, verbose=2):
        train_labels, self.label_encoder = label_encoder(labels, label_encoder=le)
        train_labels = to_categorical(train_labels, num_classes=self.n_classes)

        callbacks = []

        if early_stopping_kwargs != {}:
            callbacks.append(EarlyStopping(**early_stopping_kwargs))

        if lr_reducer_kwargs != {}:
            callbacks.append(ReduceLROnPlateau(**lr_reducer_kwargs))

        x_train = [np.reshape(seq_data, (-1, self.seq_len, self.n_channels)), fcn_data]
        y_train = train_labels

        self.model.fit(x=x_train,
                       y=y_train,
                       validation_split=0.2,
                       epochs=n_epochs,
                       batch_size=batch_size,
                       verbose=verbose,
                       callbacks=callbacks,
                       )
