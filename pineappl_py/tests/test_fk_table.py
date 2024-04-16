import numpy as np

import pineappl


class TestFkTable:
    def fake_grid(self, bins=None):
        lumis = [pineappl.lumi.LumiEntry([(1, 21, 0.1)])]
        orders = [pineappl.grid.Order(3, 0, 0, 0)]
        bin_limits = np.array([1e-7, 1e-3, 1] if bins is None else bins, dtype=float)
        subgrid_params = pineappl.subgrid.SubgridParams()
        g = pineappl.grid.Grid.create(lumis, orders, bin_limits, subgrid_params)
        return g

    def test_convolute_with_one(self):
        g = self.fake_grid()

        # DIS grid
        xs = np.linspace(0.5, 1.0, 5)
        vs = xs.copy()
        subgrid = pineappl.import_only_subgrid.ImportOnlySubgridV1(
            vs[np.newaxis, :, np.newaxis],
            np.array([90.0]),
            xs,
            np.array([1.0]),
        )
        g.set_subgrid(0, 0, 0, subgrid)
        fk = pineappl.fk_table.FkTable(g)
        np.testing.assert_allclose(
            fk.convolute_with_one(2212, lambda pid, x, q2: 0.0, lambda q2: 0.0),
            [0.0] * 2,
        )
        np.testing.assert_allclose(
            fk.convolute_with_one(2212, lambda pid, x, q2: 1, lambda q2: 1.0),
            [5e6 / 9999, 0.0],
        )
        np.testing.assert_allclose(
            fk.convolute_with_one(2212, lambda pid, x, q2: 1, lambda q2: 2.0),
            [2**3 * 5e6 / 9999, 0.0],
        )
