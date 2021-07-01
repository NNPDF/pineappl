
import numpy as np

from . import grid

class FKTable:
    def __init__(self, fktable, xgrid, basis):
        self.fktable = fktable
        self.xgrid = xgrid
        self.basis = basis

    @property
    def ndata(self):
        return self.fktable.shape[0]

    @property
    def nx(self):
        return self.xgrid.shape[0]

    @property
    def nbasis(self):
        return self.basis.shape[0]

    @classmethod
    def load(cls, filename):
        g = grid.Grid.read(filename)
        obs = g.bins()
        xgrid = np.array([.1,.2])
        basis = np.array([21,200])
        fktable = np.random.rand((obs,len(basis), len(basis), len(xgrid), len(xgrid)))
        return cls(fktable, xgrid, basis)
