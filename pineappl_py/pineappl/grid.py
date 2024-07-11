import numpy as np

from .fk_table import FkTable
from .pineappl import PyGrid, PyOrder, PyOperatorSliceInfo, PyPidBasis
from .utils import PyWrapper


class Order(PyWrapper):
    r"""Python wrapper object to interface :class:`~pineappl.pineappl.PyOrder`.

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

    @staticmethod
    def create_mask(orders, max_as, max_al, logs):
        r"""
        Return a mask suitable to pass as the `order_mask` parameter of
        :meth:`Grid.convolve`.

        Parameters
        ----------
        orders : list(Order)
            list of available orders
        max_as : int
            maximum power of :math:`\alpha_s`
        max_al : int
            maximum power of :math:`\alpha`
        logs : bool
            whether to include log grids or not

        Returns
        -------
        list(bool)
            boolean mask

        """
        return PyOrder.create_mask([o._raw for o in orders], max_as, max_al, logs)


class Grid(PyWrapper):
    r"""Python wrapper object to interface :class:`~pineappl.pineappl.PyGrid`.

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
        """Create a grid object from its ingredients.

        Parameters
        ---------
            lumi : list(LumiEntry)
                List of active luminosities
            orders: list(Order)
                List of available orders
            bin_limits: sequence(float)
                Bin limits
            subgrid_params : SubgridParams
                subgrid parameters
        """
        lumi = [lentry.raw for lentry in lumi]
        orders = [o.raw for o in orders]
        return cls(PyGrid(lumi, orders, np.array(bin_limits), subgrid_params.raw))

    def subgrid(self, order, bin_, lumi):
        """Retrieve the subgrid at the given position.

        Convenience wrapper for :meth:`pineappl.pineappl.PyGrid.set_subgrid()`.

        Parameters
        ----------
            order : int
                index of order
            bin_ : int
                index of bin
            lumi : int
                index of luminosity

        Returns
        -------
            subgrid : Subgrid
                subgrid content
        """
        return self.raw.subgrid(order, bin_, lumi)

    def __getitem__(self, key):
        """Retrieve the subgrid at the given position.

        Syntactic sugar for :meth:`subgrid`

        Parameters
        ----------
            key : (int, int, int)
                a 3-element integers tuple, consisting in `(order, bin, lumi)`

        Returns
        -------
            subgrid : Subgrid
                subgrid content
        """
        if len(key) != 3:
            raise ValueError("A tuple with `(order, bin, lumi)` is required as key.")

        return self.subgrid(*key)

    def set_subgrid(self, order, bin_, lumi, subgrid):
        """Set the subgrid at the given position.

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

    def __setitem__(self, key, subgrid):
        """Set the subgrid at the given position.

        Syntactic sugar for :meth:`set_subgrid`

        Parameters
        ----------
            key : (int, int, int)
                a 3-element integers tuple, consisting in `(order, bin, lumi)`
            subgrid : ImportOnlySubgridV1
                subgrid content
        """
        if len(key) != 3:
            raise ValueError("A tuple with `(order, bin, lumi)` is required as key.")

        self.set_subgrid(*key, subgrid)

    def set_remapper(self, remapper):
        """Set the normalizations.

        Convenience wrapper for :meth:`pineappl.pineappl.PyGrid.set_remapper()`.

        Parameters
        ----------
            remapper: BinRemapper
                Remapper object
        """
        self.raw.set_remapper(remapper.raw)

    def orders(self):
        """Extract the available perturbative orders and scale variations.

        Convenience wrapper for :meth:`pineappl.pineappl.PyGrid.orders()`.

        Parameters
        ----------
            list(Order) :
                list with perturbative orders and scale variations
        """
        return [Order(*pyorder.as_tuple()) for pyorder in self.raw.orders()]

    def convolve_with_one(
        self,
        pdg_id,
        xfx,
        alphas,
        order_mask=np.array([], dtype=bool),
        bin_indices=np.array([], dtype=np.uint64),
        lumi_mask=np.array([], dtype=bool),
        xi=((1.0, 1.0),),
    ):
        r"""Convolute grid with pdf.

        Parameters
        ----------
            pdg_id : int
                PDG Monte Carlo ID of the hadronic particle `xfx` is the PDF for
            xfx : callable
                lhapdf like callable with arguments `pid, x, Q2` returning x*pdf for :math:`x`-grid
            alphas : callable
                lhapdf like callable with arguments `Q2` returning :math:`\alpha_s`
            order_mask : sequence(bool)
                Mask for selecting specific orders. The value `True` means the corresponding order
                is included. An empty list corresponds to all orders being enabled.
            bin_indices : sequence(int)
                A list with the indices of the corresponding bins that should be calculated. An
                empty list means that all orders should be calculated.
            lumi_mask : sequence(bool)
                Mask for selecting specific luminosity channels. The value `True` means the
                corresponding channel is included. An empty list corresponds to all channels being
                enabled.
            xi : list((float, float))
                A list with the scale variation factors that should be used to calculate
                scale-varied results. The first entry of a tuple corresponds to the variation of
                the renormalization scale, the second entry to the variation of the factorization
                scale. If only results for the central scale are need the list should contain
                `(1.0, 1.0)`.

        Returns
        -------
            list(float) :
                cross sections for all bins, for each scale-variation tuple (first all bins, then
                the scale variation)
        """
        return self.raw.convolve_with_one(
            pdg_id,
            xfx,
            alphas,
            np.array(order_mask),
            np.array(bin_indices),
            np.array(lumi_mask),
            xi,
        )

    def convolve_with_two(
        self,
        pdg_id1,
        xfx1,
        pdg_id2,
        xfx2,
        alphas,
        order_mask=np.array([], dtype=bool),
        bin_indices=np.array([], dtype=np.uint64),
        lumi_mask=np.array([], dtype=bool),
        xi=((1.0, 1.0),),
    ):
        r"""Convolute grid with two pdfs.

        Parameters
        ----------
            pdg_id1 : int
                PDG Monte Carlo ID of the hadronic particle `xfx1` is the PDF for
            xfx1 : callable
                lhapdf like callable with arguments `pid, x, Q2` returning x*pdf for :math:`x`-grid
            pdg_id2 : int
                PDG Monte Carlo ID of the hadronic particle `xfx2` is the PDF for
            xfx2 : callable
                lhapdf like callable with arguments `pid, x, Q2` returning x*pdf for :math:`x`-grid
            alphas : callable
                lhapdf like callable with arguments `Q2` returning :math:`\alpha_s`
            order_mask : sequence(bool)
                Mask for selecting specific orders. The value `True` means the corresponding order
                is included. An empty list corresponds to all orders being enabled.
            bin_indices : sequence(int)
                A list with the indices of the corresponding bins that should be calculated. An
                empty list means that all orders should be calculated.
            lumi_mask : sequence(bool)
                Mask for selecting specific luminosity channels. The value `True` means the
                corresponding channel is included. An empty list corresponds to all channels being
                enabled.
            xi : list((float, float))
                A list with the scale variation factors that should be used to calculate
                scale-varied results. The first entry of a tuple corresponds to the variation of
                the renormalization scale, the second entry to the variation of the factorization
                scale. If only results for the central scale are need the list should contain
                `(1.0, 1.0)`.

        Returns
        -------
            list(float) :
                cross sections for all bins, for each scale-variation tuple (first all bins, then
                the scale variation)
        """
        return self.raw.convolve_with_two(
            pdg_id1,
            xfx1,
            pdg_id2,
            xfx2,
            alphas,
            np.array(order_mask),
            np.array(bin_indices),
            np.array(lumi_mask),
            xi,
        )

    def evolve(
        self,
        operators,
        mur2_grid,
        alphas_values,
        lumi_id_types="pdg_mc_ids",
        order_mask=(),
        xi=(1.0, 1.0),
    ):
        """Create an FKTable with the EKO.

        Convenience wrapper for :meth:`pineappl.pineappl.PyGrid.evolve()`.

        Parameters
        ----------
        operators : dict
            EKO Output
        mur2_grid : list[float]
            renormalization scales
        alphas_values : list[float]
            alpha_s values associated to the renormalization scales
        lumi_id_types : str
            kind of lumi types (e.g. "pdg_mc_ids" for flavor basis, "evol"
            for evolution basis)
        order_mask : list(bool)
            Mask for selecting specific orders. The value `True` means the corresponding order
            is included. An empty list corresponds to all orders being enabled.
        xi : (float, float)
            A tuple with the scale variation factors that should be used.
            The first entry of a tuple corresponds to the variation of
            the renormalization scale, the second entry to the variation of the factorization
            scale. If only results for the central scale are need the tuple should be
            `(1.0, 1.0)`.

        Returns
        ------
        PyFkTable :
            raw grid as an FKTable
        """
        operator_grid = np.array(
            [op["operators"] for op in operators["Q2grid"].values()]
        )
        q2grid = list(operators["Q2grid"].keys())
        return FkTable(
            self.raw.evolve(
                np.array(operator_grid),
                operators["q2_ref"],
                np.array(operators["inputpids"], dtype=np.int32),
                np.array(operators["inputgrid"]),
                np.array(q2grid, dtype=np.float64),
                np.array(operators["targetpids"], dtype=np.int32),
                np.array(operators["targetgrid"]),
                np.array(mur2_grid, dtype=np.float64),
                np.array(alphas_values, dtype=np.float64),
                xi,
                lumi_id_types,
                np.array(order_mask, dtype=bool),
            )
        )

    @classmethod
    def read(cls, path):
        """Load an existing grid from file.

        Convenience wrapper for :meth:`pineappl.pineappl.PyGrid.read()`.

        Parameters
        ----------
            path : pathlike
                file path

        Returns
        -------
            Grid
                grid object
        """
        return cls(PyGrid.read(path))

    def merge(self, other: "Grid"):
        """Merge a second grid in the current one."""
        self.raw.merge(other.raw)

    def delete_bins(self, bin_indices):
        """Delete bins.

        Repeated bins and those exceeding length are ignored.

        Parameters
        ----------
        bin_indices : sequence(int)
            list of indices of bins to removed
        """
        self.raw.delete_bins(np.array(bin_indices, dtype=np.uint))
