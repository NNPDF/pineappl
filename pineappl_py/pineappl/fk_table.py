from .pineappl import PyFkTable, PyFkAssumptions
from .utils import PyWrapper


class FkTable(PyWrapper):
    """
    Python wrapper object to interface :class:`~pineappl.pineappl.PyFkTable`.

    Parameters
    ----------
        pyfktable : PyFkTable
            raw wrapper object
    """

    def __init__(self, pyfktable):
        self._raw = pyfktable

    @classmethod
    def read(cls, path):
        """
        Load an existing grid from file.

        Convenience wrapper for :meth:`pineappl.pineappl.PyFkTable.read()`.

        Parameters
        ----------
            path : pathlike
                file path

        Returns
        -------
            FkTable
                grid object
        """
        return cls(PyFkTable.read(path))


class FkAssumptions(PyWrapper):
    """
    Python wrapper object to interface
    :class:`~pineappl.pineappl.PyFkAssumptions`.

    Parameters
    ----------
        assumption : str
            assumption identifier
    """

    def __init__(self, assumption):
        self._raw = PyFkAssumptions(assumption)
