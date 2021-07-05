
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
        xgrid = g.subgrid(0,0,0).x1_grid()
        # load lumis - they actually have to be single-valued
        basis = []
        for e in g.lumi():
            le = e.into_array()
            if len(le) != 1:
                raise ValueError("Lumis in FKTables have to be a single element")
            le = le[0]
            if not np.isclose(le[-1],1.):
                raise ValueError("All Lumis in a FKTables have to have weight 1.0")
            basis.append(le[:2])
        basis = np.array(basis)

        fktable = np.random.rand(obs, len(basis), 1, len(xgrid))
        return cls(fktable, xgrid, basis)
