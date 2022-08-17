import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder


def train_test_split_data(data, labels=None, train_frac=0.85, stratify=True):
    if labels is not None:
        x_train, x_valid, y_train, y_valid = train_test_split(data, labels, test_size=1. - train_frac,
                                                              stratify=labels)
        return x_train, x_valid, y_train, y_valid
    else:
        x_train, x_valid = train_test_split(data.X, test_size=1. - train_frac)

        return x_train, x_valid


def label_encoder(labels, label_encoder):
    labels_encoded = np.zeros(labels.shape[0])
    if isinstance(label_encoder, dict):
        for condition, label in label_encoder.items():
            labels_encoded[labels == condition] = label
    elif isinstance(label_encoder, LabelEncoder):
        labels_encoded = label_encoder.transform(labels)
    else:
        label_encoder = LabelEncoder()
        labels_encoded = label_encoder.fit_transform(labels)
    return labels_encoded.reshape(-1, 1), label_encoder
