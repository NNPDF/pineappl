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
        # check xgrid
        x1_grid = g.subgrid(0, 0, 0).x1_grid()
        x2_grid = g.subgrid(0, 0, 0).x2_grid()
        if np.isclose(x1_grid, x2_grid).all():
            xgrid = x1_grid
        else:
            if len(x2_grid) == 1:
                xgrid = x1_grid
            elif len(x1_grid) == 1:
                xgrid = x2_grid
            else:
                raise ValueError(
                    "x1_grid and x2_grid have to be either the same or one of the two is [1.0]"
                )
        # load lumis - they actually have to be single-valued
        basis = []
        for e in g.lumi():
            le = e.into_array()
            if len(le) != 1:
                raise ValueError("Lumis in FKTables have to be a single element")
            le = le[0]
            if not np.isclose(le[-1], 1.0):
                raise ValueError("All Lumis in a FKTables have to have weight 1.0")
            basis.append(le[:2])
        basis = np.array(basis)
        # fill numpy representation
        fktable = np.zeros((obs, len(basis), len(x1_grid), len(x2_grid)))
        for bin in range(obs):
            for lumi in range(basis.shape[0]):
                try:
                    fktable[bin, lumi] = g.subgrid(0, bin, lumi).fk_subgrid_array()
                except ValueError:
                    raise ValueError(f"Invalid FK subgrid: (0, {bin}, {lumi})")

        return cls(fktable, xgrid, basis)
