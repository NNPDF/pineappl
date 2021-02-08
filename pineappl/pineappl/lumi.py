try:
    from .pineappl import PyLumiEntry
except:
    import warnings

    warnings.warn("binary files missing")

from .utils import PyWrapper


class LumiEntry(PyWrapper):
    """Luminosity function.

    Each entry consists of a tuple, which contains, in the following order:

    1. the PDG id of the first incoming parton
    2. the PDG id of the second parton
    3. a numerical factor that will multiply the result for this specific
    combination.

    Parameters
    ----------
    lumis : Sequence(tuple)
        A sequence describing a luminosity function.

    """

    def __init__(self, lumis):
        self._raw = PyLumiEntry(lumis)
