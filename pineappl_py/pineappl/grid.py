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

    def subgrid_q2s(self, order, bin_, lumi):
        return np.array(self.raw.subgrid_q2s(order, bin_, lumi))

    def set_subgrid(self, order, bin_, lumi, subgrid):
        self.raw.set_subgrid(order, bin_, lumi, subgrid.raw)

    def set_remapper(self, remapper):
        self.raw.set_remapper(remapper.raw)

    @classmethod
    def read(cls, path):
        return cls(PyGrid.read(path))
