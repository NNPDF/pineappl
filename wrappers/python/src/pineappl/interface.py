"""
PineAPPL interface using ctypes
"""
import ctypes
from pineappl.loader import pineappl_lib, AVAILABLE_KEYS
from pineappl.loader import pineappl_lumi, pineappl_keyval, pineappl_grid


class lumi:
    """Luminosity object, creates a new luminosity function.

    Example:
        ::

            import pineappl

            # create a luminosity object
            lumi = pineappl.lumi()

            # adds linear combination of initial states
            lumi.add(pdg_ids, ckm_factors)

            # (optional) force cleanup / class destructor
            del lumi
    """

    def __init__(self):
        self._lumi = ctypes.byref(
            pineappl_lumi.from_address(pineappl_lib.pineappl_lumi_new())
        )

    def __del__(self):
        """The luminosity destructor method."""
        if self._lumi is not None:
            pineappl_lib.pineappl_lumi_delete(self._lumi)

    def add(self, pdg_id_pairs, factors):
        """Adds a linear combination of initial states to the
        luminosity function.

        Args:
            pdg_id_pair (list): list with pairs of parton distribution
                functions ids following the pdg convention.
            factors (list): list of ckm factors.
        """
        combinations = len(factors)
        if len(pdg_id_pairs) != 2 * combinations:
            raise RuntimeError("pdg_id_pairs and factors size mismatch")
        type1 = ctypes.c_int32 * len(pdg_id_pairs)
        type2 = ctypes.c_double * len(factors)
        pineappl_lib.pineappl_lumi_add(
            self._lumi, combinations, type1(*pdg_id_pairs), type2(*factors))


class grid:
    """The PineAPPL grid object. Creates a new and empty grid.
       The creation requires four different sets of parameters.

    Args:
        luminosity (pineappl.interface.lumi): the luminosity function that
            specifies how the cross section should be reconstructed.
        order_params (list): number of different perturbative orders.
        bin_limits (list): entries denoting the left and right limit
            for each bin.
        key_vals (pineappl.keyval): a key-value storage object.

    Example:
        ::

            import pineappl

            # create luminosity and key val objects
            lumi = pineappl.lumi()
            lumi.add(pdg_ids, factors)

            # create grid
            grid = pineappl.grid(lumi, orders, bins)

            # fill grid
            grid.fill(x1, x2, q2, orders, observable, lumi, weight)

            # write grid to file
            grid.write('my_grid.pineappl')
    """

    def __init__(self, luminosity, order_params, bin_limits, key_vals=None):
        # initialize
        self._keyval = None
        self._grid = None

        # basic checks
        if len(bin_limits) < 2:
            raise RuntimeError("bin_limits lenght must be larger than 2.")
        if len(order_params) % 4 != 0:
            raise RuntimeError(
                "order_params lenght must be a multiple of four.")
        if key_vals is not None:
            if not isinstance(key_vals, dict):
                raise ValueError("key_vals must be a dictionnary.")

        # allocating keyvals
        self._keyval = ctypes.byref(
            pineappl_keyval.from_address(pineappl_lib.pineappl_keyval_new())
        )

        if key_vals is not None:
            for key, value in key_vals.items():
                if key in AVAILABLE_KEYS:
                    if isinstance(value, str):
                        value = value.encode('utf-8')
                    AVAILABLE_KEYS.get(key)(self._keyval,
                                            key.encode('utf-8'),
                                            value)

        orders = len(order_params) // 4
        bins = len(bin_limits) - 1
        type1 = ctypes.c_uint32 * len(order_params)
        type2 = ctypes.c_double * len(bin_limits)
        self._grid = ctypes.byref(
            pineappl_grid.from_address(
                pineappl_lib.pineappl_grid_new(
                    luminosity._lumi, orders, type1(*order_params),
                    bins, type2(*bin_limits),
                    self._keyval)
            )
        )

    def __del__(self):
        """The grid destructor method."""
        if self._keyval is not None:
            pineappl_lib.pineappl_keyval_delete(self._keyval)
        if self._grid is not None:
            pineappl_lib.pineappl_grid_delete(self._grid)

    def fill(self, x1, x2, q2, order, observable, lumi, weight):
        """Fills the grid for a specific kinematics.

        Args:
            x1 (float): the partonic momentum fraction.
            x2 (float): the partonic momentum fraction.
            q2 (float) the event energy scale.
            order (int): the order value.
            observable (float): the observable value.
            lumi (int): the luminosity channel value.
            weight (float): the event weight.
        """
        pineappl_lib.pineappl_grid_fill(self._grid,
                                    x1, x2, q2,
                                    order, observable,
                                    lumi, weight)

    def write(self, filename):
        """Write grid to disk.

        Args:
            filename (str): the filename for the grid storage.
        """
        pineappl_lib.pineappl_grid_write(self._grid, filename.encode('utf-8'))
