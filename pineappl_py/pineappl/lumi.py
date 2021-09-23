try:
    from .pineappl import PyLumiEntry
except:
    import warnings

    warnings.warn("binary files missing")

from .utils import PyWrapper


class LumiEntry(PyWrapper):
    """
    Python wrapper object to :class:`~pineappl.pineappl.PyLumiEntry`.

    Parameters
    ----------
        lumis : list(tuple(int,int,float))
            sequence describing a luminosity function.
    """

    def __init__(self, lumis):
        self._raw = PyLumiEntry(lumis)
