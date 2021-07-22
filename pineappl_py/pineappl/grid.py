try:
    from .pineappl import PyGrid, PyOrder
except:
    import warnings

    warnings.warn("binary files missing")

from .utils import PyWrapper

import numpy as np


class Order(PyWrapper):
    def __init__(self, alphas, alpha, logxir, logxif):
        self._raw = PyOrder(alphas, alpha, logxir, logxif)


class Grid(PyWrapper):
    def __init__(self, pygrid):
        self._raw = pygrid

    @classmethod
    def create(cls, lumi, orders, bin_limits, subgrid_params):
        lumi = [l.raw for l in lumi]
        orders = [o.raw for o in orders]
        return cls(PyGrid(lumi, orders, bin_limits, subgrid_params.raw))

    def set_subgrid(self, order, bin_, lumi, subgrid):
        self.raw.set_subgrid(order, bin_, lumi, subgrid.into())

    def set_remapper(self, remapper):
        self.raw.set_remapper(remapper.raw)

    def convolute_eko(self, alphas, operators):
        operator_grid = np.array(
            [op["operators"] for op in operators["Q2grid"].values()]
        )
        q2grid = list(operators["Q2grid"].keys())
        alphas_values = [alphas(q2) for q2 in q2grid]
        return self.raw.convolute_eko(
            operators["q2_ref"],
            alphas_values,
            operators["pids"],
            operators["interpolation_xgrid"],
            q2grid,
            operator_grid.flatten(),
            operator_grid.shape,
        )

    @classmethod
    def read(cls, path):
        return cls(PyGrid.read(path))
