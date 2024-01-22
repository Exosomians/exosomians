from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.feature_extraction.text import CountVectorizer
from typing import Iterable, Optional

def extract_kmers(rna_sequences: Iterable[str], k: int):
    kmers_list = []
    for rna_sequence in rna_sequences:
        length = len(rna_sequence)
        kmers_list.append([rna_sequence[i:i + k] for i in range(length - k + 1)])
    return kmers_list

def encode_kmers(kmer_representations: Iterable[Iterable[str]], 
                 max_features: Optional[int] = 5000):
    # Step 2: Encode the n-grams as features (bag-of-words representation)
    vectorizer = CountVectorizer(analyzer=lambda x: x, max_features=max_features, lowercase=False)
    X = vectorizer.fit_transform(kmer_representations)

    return X

def preprocess_kmers(X, method: Optional[str] = 'tf-idf'):
    if method == 'tf-idf':
        from sklearn.feature_extraction.text import TfidfTransformer
        tfidf_transformer = TfidfTransformer()
        X = tfidf_transformer.fit_transform(X)
    elif method == 'binary':
        X[X > 0] = 1
    elif method == 'standard-scaling':
        from sklearn.preprocessing import StandardScaler
        scaler = StandardScaler(with_mean=False)
        X = scaler.fit_transform(X)
    else:
        raise ValueError('Invalid method: {}'.format(method))

    return X