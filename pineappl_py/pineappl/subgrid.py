try:
    from .pineappl import PySubgridParams
except:
    import warnings

    warnings.warn("binary files missing")

from .utils import PyWrapper


class SubgridParams(PyWrapper):
    """
    Python wrapper object to :class:`~pineappl.pineappl.PySubgridParams`.
    """
    def __init__(self):
        self._raw = PySubgridParams()
