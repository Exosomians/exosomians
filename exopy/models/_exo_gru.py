from ._model import ExoNet


class ExoGRU(ExoNet):
    def __init__(self, **kwargs):
        kwargs['network'] = 'exogru'
        super().__init__(**kwargs)
