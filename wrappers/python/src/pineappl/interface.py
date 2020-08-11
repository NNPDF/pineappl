"""
PineAPPL interface using ctypes
"""
import ctypes
import pkgconfig


# loading pineppl library
if not pkgconfig.exists('pineappl_capi'):
    raise RuntimeError('Cannot find the PineAPPL C-API, please make' \
                       'sure pkgconfig is able to access the pineappl')

paths = pkgconfig.libs('pineappl_capi').split(' ')
libdir = f'{paths[0][2:]}/lib{paths[1][2:]}.so'
pineappl = ctypes.CDLL(libdir)


class pineappl_lumi(ctypes.Structure):
    pass


class pineappl_keyval(ctypes.Structure):
    pass


class pineappl_grid(ctypes.Structure):
    pass


pineappl.pineappl_lumi_new.restype = ctypes.c_void_p
pineappl.pineappl_keyval_new.restype = ctypes.c_void_p
pineappl.pineappl_grid_new.restype = ctypes.c_void_p
pineappl.pineappl_grid_fill.argtypes = [
        ctypes.POINTER(pineappl_grid),
        ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_uint, ctypes.c_double, ctypes.c_uint, ctypes.c_double
    ]


class lumi:

    def __init__(self):
        self._lumi = ctypes.byref(
            pineappl_lumi.from_address(pineappl.pineappl_lumi_new())
        )

    def __del__(self):
        pineappl.pineappl_lumi_delete(self._lumi)

    def add(self, combinations, pdg_id_pairs, factors):
        type1 = ctypes.c_int32 * len(pdg_id_pairs)
        type2 = ctypes.c_double * len(factors)
        pineappl.pineappl_lumi_add(self._lumi, 2, type1(*pdg_id_pairs), type2(*factors))


class keyval:

    def __init__(self):
        self._keyval = ctypes.byref(
            pineappl_keyval.from_address(pineappl.pineappl_keyval_new())
        )

    def __del__(self):
        pineappl.pineappl_keyval_delete(self._keyval)


class grid:

    def __init__(self, luminosity, orders, order_params,
                 bins, bin_limits, key_vals):
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
        pineappl.pineappl_grid_delete(self._grid)

    def fill(self, x1, x2, q2, order, observable, lumi, weight):
        pineappl.pineappl_grid_fill(self._grid,
                                    x1, x2, q2,
                                    order, observable,
                                    lumi, weight)

    def write(self, filename):
        pineappl.pineappl_grid_write(self._grid, filename.encode('utf-8'))