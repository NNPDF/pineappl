try:
    from .pineappl import PyExtraSubgridParams, PySubgridParams
except:
    import warnings

    warnings.warn("binary files missing")

from .utils import PyWrapper


class SubgridParams(PyWrapper):
    def __init__(self):
        self._raw = PySubgridParams()


class ExtraSubgridParams(PyWrapper):
    def __init__(self):
        self._raw = PyExtraSubgridParams()
