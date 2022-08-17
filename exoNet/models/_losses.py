from keras import backend as K


def sensitivity(y_true, y_pred):
    y_true = K.argmax(y_true, axis=1)
    y_pred = K.argmax(y_pred, axis=1)
    true_positives = K.cast(K.sum(K.round(K.clip(y_true * y_pred, 0, 1))), 'float32')
    possible_positives = K.cast(K.sum(K.round(K.clip(y_true, 0, 1))), 'float32')
    return true_positives / (possible_positives + K.constant(K.epsilon()))


def specificity(y_true, y_pred):
    y_true = K.argmax(y_true, axis=1)
    y_pred = K.argmax(y_pred, axis=1)
    true_negatives = K.cast(K.sum(K.round(K.clip((1 - y_true) * (1 - y_pred), 0, 1))), 'float32')
    possible_negatives = K.cast(K.sum(K.round(K.clip(1 - y_true, 0, 1))), 'float32')
    return true_negatives / (possible_negatives + K.constant(K.epsilon()))


LOSSES = {
    "cce": 'categorical_crossentropy',
}

METRICS = {
    "sensitivity": sensitivity,
    "specificity": specificity,
}
