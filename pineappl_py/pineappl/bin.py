try:
    from .pineappl import PyBinRemapper
except:
    import warnings

    warnings.warn("binary files missing")


class BinRemapper:
    def __init__(self, normalization, limits):
        self._br = PyBinRemapper(normalization, limits)
