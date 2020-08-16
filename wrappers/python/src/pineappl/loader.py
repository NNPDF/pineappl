import ctypes
import pkgconfig


# Loading pineppl library.
if not pkgconfig.exists('pineappl_capi'):
    raise RuntimeError('Cannot find the PineAPPL C-API, please make sure ' \
                       'the PineAPPL Rust library is properly installed and ' \
                       'that pkgconfig is able to access the pineappl.pc file.')

paths = pkgconfig.libs('pineappl_capi').split(' ')
libdir = f'{paths[0][2:]}/lib{paths[1][2:]}.so'
pineappl_lib = ctypes.CDLL(libdir)

# Mirror C structures in python.
class pineappl_lumi(ctypes.Structure):
    pass


class pineappl_keyval(ctypes.Structure):
    pass


class pineappl_grid(ctypes.Structure):
    pass


# Specify the return and argument types.
pineappl_lib.pineappl_lumi_new.restype = ctypes.c_void_p
pineappl_lib.pineappl_keyval_new.restype = ctypes.c_void_p
pineappl_lib.pineappl_grid_new.restype = ctypes.c_void_p
pineappl_lib.pineappl_grid_read.restype = ctypes.c_void_p
pineappl_lib.pineappl_grid_fill.argtypes = [
    ctypes.POINTER(pineappl_grid),
    ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_uint, ctypes.c_double, ctypes.c_uint, ctypes.c_double
]
pineappl_lib.pineappl_keyval_set_double.argtypes = [
    ctypes.POINTER(pineappl_keyval),
    ctypes.c_char_p,
    ctypes.c_double
]
pineappl_lib.pineappl_keyval_set_int.argtypes = [
    ctypes.POINTER(pineappl_keyval),
    ctypes.c_char_p,
    ctypes.c_int
]
pineappl_lib.pineappl_keyval_set_bool.argtypes= [
    ctypes.POINTER(pineappl_keyval),
    ctypes.c_char_p,
    ctypes.c_bool
]
pineappl_lib.pineappl_keyval_set_string.argtypes = [
    ctypes.POINTER(pineappl_keyval),
    ctypes.c_char_p,
    ctypes.c_char_p
]
xfx_callback_prototype = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_int32,
                                          ctypes.c_double, ctypes.c_double,
                                          ctypes.c_void_p)
as_callback_prototype = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double,
                                         ctypes.c_void_p)
AVAILABLE_KEYS = {
    'q2_bins': pineappl_lib.pineappl_keyval_set_int,
    'q2_max': pineappl_lib.pineappl_keyval_set_double,
    'q2_min': pineappl_lib.pineappl_keyval_set_double,
    'q2_order': pineappl_lib.pineappl_keyval_set_int,
    'reweight': pineappl_lib.pineappl_keyval_set_bool,
    'x_bins': pineappl_lib.pineappl_keyval_set_int,
    'x_max': pineappl_lib.pineappl_keyval_set_double,
    'x_min': pineappl_lib.pineappl_keyval_set_double,
    'x_order': pineappl_lib.pineappl_keyval_set_int,
    'subgrid_type': pineappl_lib.pineappl_keyval_set_string,
}
