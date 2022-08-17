from abc import abstractmethod


class Network(object):
    def __init__(self):
        pass

    @abstractmethod
    def _create_network(self):
        pass

    @abstractmethod
    def _compile_models(self):
        pass

    @abstractmethod
    def to_latent(self):
        pass

    @abstractmethod
    def predict(self):
        pass

    @abstractmethod
    def train(self):
        pass
