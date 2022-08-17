from sklearn.preprocessing import normalize


def normalize_count_data(data, lengths):
    data = data / lengths
    return data


def normalize_data(data, norm='l2', axis=0):
    return normalize(data, norm=norm, axis=axis)
