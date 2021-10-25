import pineappl

import numpy as np


class TestOrder:
    def test_init(self):
        args = (2, 1, 0, 1)
        o = pineappl.grid.Order(*args)

        assert isinstance(o, pineappl.grid.Order)
        assert isinstance(o.raw, pineappl.pineappl.PyOrder)
        assert o.as_tuple() == args


class TestGrid:
    def fake_grid(self):
        lumis = [pineappl.lumi.LumiEntry([(1, 21, 0.1)])]
        orders = [pineappl.grid.Order(3, 0, 0, 0)]
        bin_limits = [1e-7, 1e-3, 1]
        subgrid_params = pineappl.subgrid.SubgridParams()
        g = pineappl.grid.Grid.create(lumis, orders, bin_limits, subgrid_params)
        return g

    def test_init(self):
        g = self.fake_grid()
        assert isinstance(g, pineappl.grid.Grid)
        assert isinstance(g.raw, pineappl.pineappl.PyGrid)
        # orders
        assert len(g.orders()) == 1
        assert g.orders()[0].as_tuple() == (3,0,0,0)

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
        g.optimize()

    def test_set_key_value(self):
        g = self.fake_grid()
        g.set_key_value("bla", "blub")
        g.set_key_value("\"", "'")
        g.set_key_value("äöü", "ß\\")

    def test_bins(self):
        g = self.fake_grid()
        # 1D
        normalizations = [1.0] * 2
        limits = [(1,1),(2,2)]
        remapper = pineappl.bin.BinRemapper(normalizations, limits)
        g.set_remapper(remapper)
        assert g.bin_dimensions() == 1
        np.testing.assert_allclose(g.bin_left(0), [1,2])
        np.testing.assert_allclose(g.bin_right(0), [1,2])
        # 2D
        limits = [(1,2),(2,3),(2,4),(3,5)]
        remapper = pineappl.bin.BinRemapper(normalizations, limits)
        g.set_remapper(remapper)
        assert g.bin_dimensions() == 2
        np.testing.assert_allclose(g.bin_left(0), [1,2])
        np.testing.assert_allclose(g.bin_right(0), [2,4])
        np.testing.assert_allclose(g.bin_left(1), [2,3])
        np.testing.assert_allclose(g.bin_right(1), [3,5])

    def test_convolute_with_one(self):
        g = self.fake_grid()

        # DIS grid
        xs = np.linspace(0.5, 1.0, 5)
        vs = xs.copy()
        subgrid = pineappl.import_only_subgrid.ImportOnlySubgridV1(
            vs[np.newaxis, :, np.newaxis],
            [90],
            xs,
            [1.0],
        )
        g.set_subgrid(0, 0, 0, subgrid)
        np.testing.assert_allclose(g.convolute_with_one(2212, lambda pid,x,q2: 0., lambda q2: 0.), [0.]*2)
        np.testing.assert_allclose(g.convolute_with_one(2212, lambda pid,x,q2: 1, lambda q2: 1.), [5e6/9999,0.])
        np.testing.assert_allclose(g.convolute_with_one(2212, lambda pid,x,q2: 1, lambda q2: 2.), [2**3 * 5e6/9999,0.])

    def test_axes(self):
        g = self.fake_grid()

        # add 2 DIS grids
        xs = np.linspace(0.5, 1.0, 5)
        vs = np.random.rand(len(xs))
        subgrid = pineappl.import_only_subgrid.ImportOnlySubgridV1(
            vs[np.newaxis, :, np.newaxis],
            [90.],
            xs,
            [1.0],
        )
        g.set_subgrid(0, 0, 0, subgrid)
        vs2 = np.random.rand(len(xs))
        subgrid = pineappl.import_only_subgrid.ImportOnlySubgridV1(
            vs2[np.newaxis, :, np.newaxis],
            [100.],
            xs,
            [1.0],
        )
        g.set_subgrid(0, 1, 0, subgrid)
        # now get the thing
        ei = g.axes()

        np.testing.assert_allclose(ei[0], xs)
        np.testing.assert_allclose(ei[1], [])
        np.testing.assert_allclose(ei[2], [90., 100.])

    def test_io(self, tmp_path):
        g = self.fake_grid()
        p = tmp_path / "test.pineappl"
        p.write_text("")
        g.write(str(p))
        gg = pineappl.grid.Grid.read(str(p))
        assert isinstance(gg, pineappl.grid.Grid)

    def test_convolute_eko(self):
        g = self.fake_grid()
        fake_eko = {
            "q2_ref": 1.,
            "targetpids": [1],
            "targetgrid": [.1,1.],
            "inputpids": [1],
            "inputgrid": [.1,1.],
            "interpolation_xgrid": [.1,1.],
            "Q2grid": {
                90: {
                    "operators": np.random.rand(1,2,1,2),
                    "alphas": 1.
                }
            }
        }
        g.set_key_value("lumi_id_types", "pdg_mc_ids")
        fk = g.convolute_eko(fake_eko)
        assert isinstance(fk.raw, pineappl.pineappl.PyFkTable)
