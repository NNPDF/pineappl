"""
PineAPPL interface using ctypes
"""
import ctypes
import pkgconfig


# loading pineppl library
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


# refers to CAPI documentation
def lumi_new():
    return pineappl_lumi.from_address(pineappl.pineappl_lumi_new())


def lumi_add(lumi, combinations, pdg_id_pairs, factors):
    type1 = ctypes.c_int32 * len(pdg_id_pairs)
    type2 = ctypes.c_double * len(factors)
    pineappl.pineappl_lumi_add(ctypes.byref(lumi), 2, type1(*pdg_id_pairs), type2(*factors))


def keyval_new():
    return pineappl_keyval.from_address(pineappl.pineappl_keyval_new())


def grid_new(lumi, orders, order_params, bins, bin_limits, key_vals):
    type1 = ctypes.c_uint32 * len(order_params)
    type2 = ctypes.c_double * len(bin_limits)
    return pineappl_grid.from_address(pineappl.pineappl_grid_new(
        ctypes.byref(lumi), orders, type1(*order_params), bins, type2(*bin_limits),
        ctypes.byref(key_vals)
    ))


def grid_fill(grid, x1, x2, q2, order, observable, lumi, weight):
    pineappl.pineappl_grid_fill(ctypes.byref(grid),
                                 x1, x2, q2, order, observable, lumi, weight)


def grid_write(grid, filename):
    pineappl.pineappl_grid_write(ctypes.byref(grid), filename.encode('utf-8'))


def keyval_delete(keyval):
    pineappl.pineappl_keyval_delete(ctypes.byref(keyval))


def lumi_delete(lumi):
    pineappl.pineappl_lumi_delete(ctypes.byref(lumi))


def grid_delete(grid):
    pineappl.pineappl_grid_delete(ctypes.byref(grid))
