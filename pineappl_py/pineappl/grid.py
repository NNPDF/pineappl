import numpy as np

from .pineappl import PyGrid, PyOrder
from .utils import PyWrapper


class Order(PyWrapper):
    r"""
    Python wrapper object to interface :class:`~pineappl.pineappl.PyOrder`.

    Parameters
    ----------
        alphas : int
            power of :math:`\alpha_s`
        alpha : int
            power of :math:`\alpha`
        logxir : int
            power of :math:`\log(\xi_r)`
        logxif : int
            power of :math:`\log(\xi_f)`
    """

    def __init__(self, alphas, alpha, logxir, logxif):
        self._raw = PyOrder(alphas, alpha, logxir, logxif)


class Grid(PyWrapper):
    r"""
    Python wrapper object to interface :class:`~pineappl.pineappl.PyGrid`.

    To create an object, you should call either :meth:`create`
    or :meth:`read`.

    Parameters
    ----------
        pygrid : PyGrid
            raw wrapper object
    """

    def __init__(self, pygrid):
        self._raw = pygrid

    @classmethod
    def create(cls, lumi, orders, bin_limits, subgrid_params):
        """
        Create a grid object from its ingredients

        Parameters
        ---------
            lumi : list(LumiEntry)
                List of active luminosities
            orders: list(Order)
                List of available orders
            bin_limits: BinRemapper
                bins
            subgrid_params : SubgridParams
                subgrid parameters
        """
        lumi = [l.raw for l in lumi]
        orders = [o.raw for o in orders]
        return cls(PyGrid(lumi, orders, bin_limits, subgrid_params.raw))

    def set_subgrid(self, order, bin_, lumi, subgrid):
        """
        Set the subgrid at the given position.

        Convenience wrapper for :meth:`pineappl.pineappl.PyGrid.set_subgrid()`.

        Parameters
        ----------
            order : int
                index of order
            bin_ : int
                index of bin
            lumi : int
                index of luminosity
            subgrid : ImportOnlySubgridV1
                subgrid content
        """
        self.raw.set_subgrid(order, bin_, lumi, subgrid.into())

    def set_remapper(self, remapper):
        """
        Set the normalizations.

        Convenience wrapper for :meth:`pineappl.pineappl.PyGrid.set_remapper()`.

        Parameters
        ----------
            remapper: BinRemapper
                Remapper object
        """
        self.raw.set_remapper(remapper.raw)

    def convolute_eko(self, operators):
        """
        Create an FKTable with the EKO.

        Convenience wrapper for :meth:`pineappl.pineappl.PyGrid.convolute_eko()`.

        Parameters
        ----------
            operators : dict
                EKO Output

        Returns
        ------
            PyFkTable :
                raw grid as an FKTable
        """
        operator_grid = np.array(
            [op["operators"] for op in operators["Q2grid"].values()]
        )
        q2grid = list(operators["Q2grid"].keys())
        alphas_values = [op["alphas"] for op in operators["Q2grid"].values()]
        return self.raw.convolute_eko(
            operators["q2_ref"],
            alphas_values,
            operators["targetpids"],
            operators["interpolation_xgrid"],
            q2grid,
            operator_grid.flatten(),
            operator_grid.shape,
        )

    @classmethod
    def read(cls, path):
        """
        Load an existing grid from file.

        Convenience wrapper for :meth:`pineappl.pineappl.PyGrid.read()`.

        Parameters
        ----------
            path : str
                file path

        Returns
        -------
            Grid
                grid object
        """
        return cls(PyGrid.read(path))
