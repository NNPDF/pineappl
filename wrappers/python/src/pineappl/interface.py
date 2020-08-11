"""
PineAPPL interface using ctypes
"""
import ctypes
import pkgconfig


# Loading pineppl library.
if not pkgconfig.exists('pineappl_capi'):
    raise RuntimeError('Cannot find the PineAPPL C-API, please make' \
                       'sure pkgconfig is able to access the pineappl')

paths = pkgconfig.libs('pineappl_capi').split(' ')
libdir = f'{paths[0][2:]}/lib{paths[1][2:]}.so'
pineappl = ctypes.CDLL(libdir)

# Mirror C structures in python.
class pineappl_lumi(ctypes.Structure):
    pass


class pineappl_keyval(ctypes.Structure):
    pass


class pineappl_grid(ctypes.Structure):
    pass


# Specify the return and argument types.
pineappl.pineappl_lumi_new.restype = ctypes.c_void_p
pineappl.pineappl_keyval_new.restype = ctypes.c_void_p
pineappl.pineappl_grid_new.restype = ctypes.c_void_p
pineappl.pineappl_grid_fill.argtypes = [
        ctypes.POINTER(pineappl_grid),
        ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_uint, ctypes.c_double, ctypes.c_uint, ctypes.c_double
    ]


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
            pineappl_lumi.from_address(pineappl.pineappl_lumi_new())
        )

    def __del__(self):
        """The luminosity destructor method."""
        pineappl.pineappl_lumi_delete(self._lumi)

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
        pineappl.pineappl_lumi_add(self._lumi, combinations, type1(*pdg_id_pairs), type2(*factors))


class keyval:
    """Auxiliary key values object."""

    def __init__(self):
        self._keyval = ctypes.byref(
            pineappl_keyval.from_address(pineappl.pineappl_keyval_new())
        )

    def __del__(self):
        """The key vals destructor method."""
        pineappl.pineappl_keyval_delete(self._keyval)


class grid:
    """The PineAPPL grid object. Creates a new and empty grid.
       The creation requires four different sets of parameters.

    Args:
        luminosity (pineappl.interface.lumi): the luminosity function that
            specifies how the cross section should be reconstructed.
        order_params (list): number of different perturbative orders.
        bin_limits (list): entries denoting the left and right limit
            for each bin.
        key_vals (pineappl.interface.keyval): a key-value storage object.

    Example:
        ::

            import pineappl

            # create luminosity and key val objects
            lumi = pineappl.lumi()
            lumi.add(pdg_ids, factors)
            keyval = pineappl.keyval()

            # create grid
            grid = pineappl.grid(lumi, orders, bins, keyval)

            # fill grid
            grid.fill(x1, x2, q2, orders, observable, lumi, weight)

            # write grid to file
            grid.write('my_grid.pineappl')
    """

    def __init__(self, luminosity, order_params, bin_limits, key_vals):
        if len(bin_limits) < 2:
            raise RuntimeError("bin_limits lenght must be larger than 2.")
        if len(order_params) % 4 != 0:
            raise RuntimeError("order_params lenght must be a multiple of four.")
        orders = len(order_params) // 4
        bins = len(bin_limits) - 1
        type1 = ctypes.c_uint32 * len(order_params)
        type2 = ctypes.c_double * len(bin_limits)
        self._grid = ctypes.byref(
            pineappl_grid.from_address(
                pineappl.pineappl_grid_new(
                    luminosity._lumi, orders, type1(*order_params),
                    bins, type2(*bin_limits),
                    key_vals._keyval)
            )
        )

    def __del__(self):
        """The grid destructor method."""
        pineappl.pineappl_grid_delete(self._grid)

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
        pineappl.pineappl_grid_fill(self._grid,
                                    x1, x2, q2,
                                    order, observable,
                                    lumi, weight)

    def write(self, filename):
        """Write grid to disk.

        Args:
            filename (str): the filename for the grid storage.
        """
        pineappl.pineappl_grid_write(self._grid, filename.encode('utf-8'))
