try:
    from .pineappl import PyGrid, PyOrder
except:
    import warnings

    warnings.warn("binary files missing")

import numpy as np

from .utils import PyWrapper


class Order(PyWrapper):
    r"""
    Python wrapper object to interface `Order`.

    Parameters
    ----------
        alphas : int
            power of :math:`\alphas_s`
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
    Python wrapper object to interface `Grid`.

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

        Parameters
        ----------
            remapper: BinRemapper
                Remapper object
        """
        self.raw.set_remapper(remapper.raw)

    def set_key_value(self, key, value):
        """
        Set a metadata key-value pair in the grid.

        Parameters
        ----------
            key : str
                key
            value : str
                value
        """
        self.raw.set_key_value(key, value)

    def convolute(self, xfx1, xfx2, alphas):
        """
        Convolute grid with given PDF sets.

        Parameters
        ----------
            xfx1 : callable
                LHAPDF like callable for x1
            xfx2 : callable
                LHAPDF like callable for x2
            alphas : float
                value of strong coupling

        Returns
        -------
            list(float) :
                predictions for all bins
        """
        return self.raw.convolute(xfx1, xfx2, alphas)

    def eko_info(self):
        """
        Extract the necessary informations for EKO.

        Returns
        -------
            x_grid: list(float)
                interpolation grid
            muf2_grid : list(float)
                factorization scale list
        """
        return self.raw.eko_info()

    def convolute_eko(self, operators):
        """
        Create an FKTable with the EKO.

        Parameters
        ----------
            operators : dict
                EKO Output

        Returns
        ------
            PyGrid
                raw Grid as an FKTable
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

    def write(self, path):
        """
        Write the grid to file.

        Parameters
        ----------
            path : str
                file path
        """
        self.raw.write(path)

    def optimize(self):
        """
        Optimize grid content.
        """
        return self.raw.optimize()

    def bin_dimensions(self):
        """
        Extract the number of dimensions for bins.

        E.g.: two differential cross-sections will return 2.

        Returns
        -------
            int :
                bin dimension
        """
        return self.raw.bin_dimensions()

    def bin_left(self, dimension):
        """
        Extract the left edges of a specific bin dimension.

        Parameters
        ----------
            dimension : int
                bin dimension

        Returns
        -------
            list(float) :
                left edges of bins
        """
        return self.raw.bin_left(dimension)

    def bin_right(self, dimension):
        """
        Extract the right edges of a specific bin dimension.

        Parameters
        ----------
            dimension : int
                bin dimension

        Returns
        -------
            list(float) :
                right edges of bins
        """
        return self.raw.bin_right(dimension)
