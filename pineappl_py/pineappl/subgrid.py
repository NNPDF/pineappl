from .pineappl import PySubgridParams
from .utils import PyWrapper


class SubgridParams(PyWrapper):
    """
    Python wrapper object to :class:`~pineappl.pineappl.PySubgridParams`.
    """

    def __init__(self):
        self._raw = PySubgridParams()
