from .pineappl import PyFkTable
from .utils import PyWrapper


class FkTable(PyWrapper):
    """
    Python wrapper object to interface :class:`~pineappl.pineappl.PyFkTable`.

    Parameters
    ----------
        grid : PyFkTable
            raw wrapper object
    """

    def __init__(self, grid):
        self._raw = PyFkTable(grid._raw)
