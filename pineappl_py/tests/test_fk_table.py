import numpy as np

import pineappl


class TestFkTable:
    def fake_grid(self, bins=None):
        lumis = [pineappl.lumi.LumiEntry([(1, 21, 1.0)])]
        orders = [pineappl.grid.Order(0, 0, 0, 0)]
        bin_limits = np.array([1e-7, 1e-3, 1] if bins is None else bins, dtype=float)
        subgrid_params = pineappl.subgrid.SubgridParams()
        g = pineappl.grid.Grid.create(lumis, orders, bin_limits, subgrid_params)
        g.set_key_value("lumi_id_types", "pdg_mc_ids")
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
        fk = pineappl.fk_table.FkTable.from_grid(g)
        np.testing.assert_allclose(
            fk.convolute_with_one(2212, lambda pid, x, q2: 0.0),
            [0.0] * 2,
        )
        np.testing.assert_allclose(
            fk.convolute_with_one(2212, lambda pid, x, q2: 1),
            [5e7 / 9999, 0.0],
        )
