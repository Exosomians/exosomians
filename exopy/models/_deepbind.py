from ._model import ExoNet


class DeepBind(ExoNet):
    def __init__(self, **kwargs):
        kwargs.pop('network', 'deepbind')
        super().__init__(network='deepbind', **kwargs)
