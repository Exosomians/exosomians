from ._model import ExoNet


class ExoFCN(ExoNet):
    def __int__(self, **kwargs):
        kwargs.pop('network', None)
        super().__init__(network='exofcn', **kwargs)
