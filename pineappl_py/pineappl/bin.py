from .pineappl import PyBinRemapper
from .utils import PyWrapper


class BinRemapper(PyWrapper):
    """
    Python wrapper object for :class:`~pineappl.pineappl.PyBinRemapper`.

    Parameters
    ----------
        normalization : list(float)
            list with normalizations
        limits : list(tuple(float,float))
            all bin limits as a flat list
    """

    def __init__(self, normalization, limits):
        self._raw = PyBinRemapper(normalization, limits)
