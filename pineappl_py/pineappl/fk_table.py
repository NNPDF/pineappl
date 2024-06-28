from .pineappl import PyFkTable, PyFkAssumptions
from .utils import PyWrapper


class FkTable(PyWrapper):
    """Python wrapper object to interface
    :class:`~pineappl.pineappl.PyFkTable`.

    Parameters
    ----------
        pyfktable : PyFkTable
            raw wrapper object
    """

    def __init__(self, pyfktable):
        self._raw = pyfktable

    @classmethod
    def from_grid(cls, grid):
        return cls(PyFkTable(grid.raw))

    @classmethod
    def read(cls, path):
        """Load an existing grid from file.

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

    def optimize(self, assumptions="Nf6Ind"):
        """Optimize FK table storage.

        In order to perform any relevant optimization, assumptions are needed, and they are passed
        as parameters to the function.

        Parameters
        ----------
        assumptions : FkAssumptions or str
            assumptions about the FkTable properties, declared by the user, deciding which
            optimizations are possible
        """
        if not isinstance(assumptions, FkAssumptions):
            assumptions = FkAssumptions(assumptions)
        return self._raw.optimize(assumptions._raw)


class FkAssumptions(PyWrapper):
    """Python wrapper object to interface
    :class:`~pineappl.pineappl.PyFkAssumptions`.

    Parameters
    ----------
        assumption : str
            assumption identifier
    """

    def __init__(self, assumption):
        self._raw = PyFkAssumptions(assumption)
