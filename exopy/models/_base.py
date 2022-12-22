from abc import abstractmethod


class RNAClassifierBase:
    @classmethod
    @abstractmethod
    def load(cls, *args, **kwargs):
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def setup_dataset(cls, *args, **kwargs):
        raise NotImplementedError

    @abstractmethod
    def get_embeddings(self, *args, **kwargs):
        raise NotImplementedError

    @abstractmethod
    def fit(self, *args, **kwargs):
        raise NotImplementedError

    @abstractmethod
    def predict(self, *args, **kwargs):
        raise NotImplementedError

    @abstractmethod
    def save(self, *args, **kwargs) -> None:
        raise NotImplementedError
