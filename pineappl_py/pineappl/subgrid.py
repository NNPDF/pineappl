from .pineappl import PySubgridParams, PyMu2
from .utils import PyWrapper


class SubgridParams(PyWrapper):
    """
    Python wrapper object to :class:`~pineappl.pineappl.PySubgridParams`.
    """

    def __init__(self):
        self._raw = PySubgridParams()
    
class Mu2(PyWrapper):

    def __init__(self, ren, fac):
        self._raw = PyMu2(ren, fac)
