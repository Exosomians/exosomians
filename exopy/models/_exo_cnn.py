from ._model import ExoNet


class ExoCNN(ExoNet):
    def __init__(self, **kwargs):
        kwargs.pop('network', None)
        super().__init__(network='exocnn', **kwargs)
