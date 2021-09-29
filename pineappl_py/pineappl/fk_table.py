from .pineappl import PyFkTable
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
            path : str
                file path

        Returns
        -------
            FkTable
                grid object
        """
        return cls(PyFkTable.read(path))
