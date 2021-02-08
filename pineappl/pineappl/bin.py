try:
    from .pineappl import PyBinRemapper
except:
    import warnings

    warnings.warn("binary files missing")

from .utils import PyWrapper


class BinRemapper(PyWrapper):
    def __init__(self, normalization, limits):
        self._raw = PyBinRemapper(normalization, limits)
