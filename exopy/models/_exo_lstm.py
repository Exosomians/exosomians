from ._model import ExoNet


class ExoLSTM(ExoNet):
    def __init__(self, **kwargs):
        kwargs['network'] = 'exolstm'
        super().__init__(**kwargs)
