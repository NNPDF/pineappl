try:
    from .pineappl import PyGrid, PyOrder
except:
    import warnings

    warnings.warn("binary files missing")

from .utils import PyWrapper


class Order(PyWrapper):
    def __init__(self, alphas, alpha, logxir, logxif):
        self._raw = PyOrder(alphas, alpha, logxir, logxif)


class Grid(PyWrapper):
    def __init__(self, lumi, orders, bin_limits, subgrid_params):
        lumi = [l.raw for l in lumi]
        orders = [o.raw for o in orders]
        self._raw = PyGrid(lumi, orders, bin_limits, subgrid_params.raw)
