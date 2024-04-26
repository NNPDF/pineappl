import numpy as np

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

    def convolute_with_one(
        self,
        pdg_id,
        xfx,
        bin_indices=np.array([], dtype=np.uint64),
        lumi_mask=np.array([], dtype=bool),
    ):
        r"""Convolute FkTable with a pdf.

        Parameters
        ----------
            pdg_id : int
                PDG Monte Carlo ID of the hadronic particle `xfx` is the PDF for
            xfx : callable
                lhapdf like callable with arguments `pid, x, Q2` returning x*pdf for :math:`x`-grid
            bin_indices : sequence(int)
                A list with the indices of the corresponding bins that should be calculated. An
                empty list means that all orders should be calculated.
            lumi_mask : sequence(bool)
                Mask for selecting specific luminosity channels. The value `True` means the
                corresponding channel is included. An empty list corresponds to all channels being
                enabled.

        Returns
        -------
            list(float) :
                cross sections for all bins, for each scale-variation tuple (first all bins, then
                the scale variation)
        """
        return self.raw.convolute_with_one(
            pdg_id,
            xfx,
            np.array(bin_indices),
            np.array(lumi_mask),
        )


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
