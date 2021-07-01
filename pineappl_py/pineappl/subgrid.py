try:
    from .pineappl import PySubgridParams
except:
    import warnings

    warnings.warn("binary files missing")

from .utils import PyWrapper


class SubgridEnum(PyWrapper):
    def __init__(self):
        self._raw = PySubgridEnum()


class SubgridParams(PyWrapper):
    def __init__(self):
        self._raw = PySubgridParams()
