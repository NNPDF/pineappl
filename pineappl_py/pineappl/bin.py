import numpy as np

from .pineappl import PyBinRemapper
from .utils import PyWrapper


class BinRemapper(PyWrapper):
    """
    Python wrapper object for :class:`~pineappl.pineappl.PyBinRemapper`.

    Parameters
    ----------
        normalizations : sequence(float)
            list with normalizations
        limits : list(tuple(float,float))
            all bin limits as a flat list
    """

    def __init__(self, normalizations, limits):
        self._raw = PyBinRemapper(np.array(normalizations), limits)
