"""
PineAPPL interface using ctypes
"""
import ctypes
import numpy as np
from pineappl.loader import pineappl_lib, AVAILABLE_KEYS
from pineappl.loader import pineappl_lumi, pineappl_keyval, pineappl_grid
from pineappl.loader import xfx_callback_prototype, as_callback_prototype


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
       You can load the pineappl grid from file using the `grid.read` method.

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

            # load pdf for testing
            pdf = lhapdf.mkPDF('NNPDF31_nlo_as_0118_luxqed', 0)

            def xfx(id, x, q2, p):
                return pdf.xfxQ2(id, x, q2)

            def alphas(q2, p):
                return pdf.alphasQ2(q2)

            # perform convolution
            result = grid.convolute(xfx, xfx, alphas, None, None, None, 1.0, 1.0)

            # write grid to file
            grid.write('my_grid.pineappl')
    """

    def __init__(self, luminosity=None, order_params=None, bin_limits=None, key_vals=None):
        # initialize
        self._keyval = None
        self._grid = None

        if None not in [luminosity, order_params, bin_limits]:
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

    def convolute(self, xfx1, xfx2, alphas, order_mask, lumi_mask, xi_ren, xi_fac):
        """Convolutes the grid with the PDFs and strong coupling.

        Args:
            xfx1 (function): xfxQ2 function of (id, x, q2, p) for parton 1.
            xfx2 (function): xfxQ2 function of (id, x, q2, p) for parton 2.
            alphas (function): alphasQ2 function of (q2, p).
            order_mask (array): must be provided as long as there are perturbative orders.
            lumi_mask (array): specify the luminosity mask.
            xi_ren (float): renomarlization scale.
            xi_fac (floag): factorization scale.

        Return:
            A numpy array with the convolution output for each bin.
        """
        cxfx1 = xfx_callback_prototype(xfx1)
        cxfx2 = xfx_callback_prototype(xfx2)
        calphas = as_callback_prototype(alphas)

        if order_mask is not None:
            type1 = ctypes.c_bool * len(order_mask)
            order_mask = type1(*order_mask)
        if lumi_mask is not None:
            type2 = ctypes.c_bool * len(lumi_mask)
            lumi_mask = type2(*lumi_mask)

        result = np.ones(self.bin_count(), dtype=np.float64)
        pineappl_lib.pineappl_grid_convolute(
            self._grid, cxfx1, cxfx2, calphas, None,
            order_mask, lumi_mask,
            ctypes.c_double(xi_ren), ctypes.c_double(xi_fac),
            ctypes.c_void_p(result.ctypes.data)
        )
        return result

    def write(self, filename):
        """Write grid to disk.

        Args:
            filename (str): the filename for the grid storage.
        """
        pineappl_lib.pineappl_grid_write(self._grid, filename.encode('utf-8'))

    def bin_count(self):
        """Return bin count."""
        return pineappl_lib.pineappl_grid_bin_count(self._grid)

    @classmethod
    def read(cls, filename):
        """Read grid from file.

        Args:
            filename (str): the filename for the pineappl grid.

        Returns:
            A pineappl grid object.
        """
        obj = cls()
        obj._grid = ctypes.byref(
            pineappl_grid.from_address(pineappl_lib.pineappl_grid_read(filename.encode('utf-8')))
            )
        return obj