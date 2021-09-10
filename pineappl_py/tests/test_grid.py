import pineappl

import numpy as np


class TestOrder:
    def test_init(self):
        o = pineappl.grid.Order(2, 1, 0, 1)

        assert isinstance(o, pineappl.grid.Order)
        assert isinstance(o.raw, pineappl.pineappl.PyOrder)


class TestGrid:
    def fake_grid(self):
        lumis = [pineappl.lumi.LumiEntry([(1, 21, 0.1)])]
        orders = [pineappl.grid.Order(3, 0, 1, 0)]
        bin_limits = [1e-7, 1e-3, 1]
        subgrid_params = pineappl.subgrid.SubgridParams()
        g = pineappl.grid.Grid.create(lumis, orders, bin_limits, subgrid_params)
        return g

    def test_init(self):
        g = self.fake_grid()
        assert isinstance(g, pineappl.grid.Grid)
        assert isinstance(g.raw, pineappl.pineappl.PyGrid)

    def test_set_subgrid(self):
        g = self.fake_grid()

        # DIS grid
        xs = np.linspace(0.1, 1.0, 5)
        vs = np.random.rand(len(xs))
        subgrid = pineappl.import_only_subgrid.ImportOnlySubgridV1(
            vs[np.newaxis, :, np.newaxis],
            [90],
            xs,
            [1.0],
        )
        g.set_subgrid(0, 0, 0, subgrid)

        # let's mix it for fun with an hadronic one
        x1s = np.linspace(0.1, 1, 2)
        x2s = np.linspace(0.5, 1, 2)
        Q2s = np.linspace(10, 20, 2)
        subgrid = pineappl.import_only_subgrid.ImportOnlySubgridV1(
            np.random.rand(len(Q2s), len(x1s), len(x2s)), Q2s, x1s, x2s
        )
        g.set_subgrid(0, 1, 0, subgrid)

    def test_set_key_value(self):
        g = self.fake_grid()
        g.set_key_value("bla", "blub")
        g.set_key_value("\"", "'")
        g.set_key_value("äöü", "ß\\")

    def test_set_remapper(self):
        g = self.fake_grid()
        normalizations = [1.0] * 2
        limits = [(1,1),(2,2)]
        remapper = pineappl.bin.BinRemapper(normalizations, limits)
        g.set_remapper(remapper)
