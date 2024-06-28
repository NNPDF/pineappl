import numpy as np
import pytest

import pineappl


class TestOrder:
    def test_init(self):
        args = (2, 1, 0, 1)
        o = pineappl.grid.Order(*args)

        assert isinstance(o, pineappl.grid.Order)
        assert isinstance(o.raw, pineappl.pineappl.PyOrder)
        assert o.as_tuple() == args


class TestGrid:
    def fake_grid(self, bins=None):
        lumis = [pineappl.lumi.LumiEntry([(1, 21, 0.1)])]
        orders = [pineappl.grid.Order(3, 0, 0, 0)]
        bin_limits = np.array([1e-7, 1e-3, 1] if bins is None else bins, dtype=float)
        subgrid_params = pineappl.subgrid.SubgridParams()
        g = pineappl.grid.Grid.create(lumis, orders, bin_limits, subgrid_params)
        return g

    def test_init(self):
        g = self.fake_grid()
        assert isinstance(g, pineappl.grid.Grid)
        assert isinstance(g.raw, pineappl.pineappl.PyGrid)
        # orders
        assert len(g.orders()) == 1
        assert g.orders()[0].as_tuple() == (3, 0, 0, 0)

    def test_set_subgrid(self):
        g = self.fake_grid()

        # DIS grid
        xs = np.linspace(0.1, 1.0, 5)
        vs = np.random.rand(len(xs))
        subgrid = pineappl.import_only_subgrid.ImportOnlySubgridV1(
            vs[np.newaxis, :, np.newaxis],
            np.array([90.0]),
            np.array(xs),
            np.array([1.0]),
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
        g.set_key_value('"', "'")
        g.set_key_value("äöü", "ß\\")

    def test_bins(self):
        g = self.fake_grid()
        # 1D
        normalizations = [1.0] * 2
        limits = [(1, 1), (2, 2)]
        remapper = pineappl.bin.BinRemapper(normalizations, limits)
        g.set_remapper(remapper)
        assert g.bin_dimensions() == 1
        np.testing.assert_allclose(g.bin_left(0), [1, 2])
        np.testing.assert_allclose(g.bin_right(0), [1, 2])
        # 2D
        limits = [(1, 2), (2, 3), (2, 4), (3, 5)]
        remapper = pineappl.bin.BinRemapper(normalizations, limits)
        g.set_remapper(remapper)
        assert g.bin_dimensions() == 2
        np.testing.assert_allclose(g.bin_left(0), [1, 2])
        np.testing.assert_allclose(g.bin_right(0), [2, 4])
        np.testing.assert_allclose(g.bin_left(1), [2, 3])
        np.testing.assert_allclose(g.bin_right(1), [3, 5])

    def test_convolve_with_one(self):
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
        np.testing.assert_allclose(
            g.convolve_with_one(2212, lambda pid, x, q2: 0.0, lambda q2: 0.0),
            [0.0] * 2,
        )
        np.testing.assert_allclose(
            g.convolve_with_one(2212, lambda pid, x, q2: 1, lambda q2: 1.0),
            [5e6 / 9999, 0.0],
        )
        np.testing.assert_allclose(
            g.convolve_with_one(2212, lambda pid, x, q2: 1, lambda q2: 2.0),
            [2**3 * 5e6 / 9999, 0.0],
        )

    def test_io(self, tmp_path):
        g = self.fake_grid()
        p = tmp_path / "test.pineappl"
        p.write_text("")
        g.write(str(p))
        gg = pineappl.grid.Grid.read(p)
        assert isinstance(gg, pineappl.grid.Grid)
        _ = pineappl.grid.Grid.read(str(p))

    def test_fill(self):
        g = self.fake_grid()
        g.fill(0.5, 0.5, 10.0, 0, 0.01, 0, 10.0)
        res = g.convolve_with_one(2212, lambda pid, x, q2: x, lambda q2: 1.0)
        pytest.approx(res) == 0.0

    def test_fill_array(self):
        g = self.fake_grid()
        g.fill_array(
            np.array([0.5, 1.0]),
            np.array([0.5, 1.0]),
            np.array([0.5, 1.0]),
            0,
            np.array([1e-3, 1e-2]),
            0,
            np.array([10.0, 100.0]),
        )
        res = g.convolve_with_one(2212, lambda pid, x, q2: x, lambda q2: 1.0)
        pytest.approx(res) == 0.0

    def test_fill_all(self):
        g = self.fake_grid()
        g.fill_all(1.0, 1.0, 1.0, 0, 1e-2, np.array([10.0]))
        res = g.convolve_with_one(2212, lambda pid, x, q2: x, lambda q2: 1.0)
        pytest.approx(res) == 0.0

    def test_merge(self):
        g = self.fake_grid([1, 2, 3])
        g1 = self.fake_grid([3, 4, 5])
        assert g.bins() == 2
        assert g1.bins() == 2

        g.merge(g1)
        assert g.bins() == 4

        g2 = self.fake_grid([1, 2, 3])
        g3 = self.fake_grid([1, 2, 3])
        assert g2.bins() == 2
        assert g3.bins() == 2

        g2.merge(g3)
        assert g2.bins() == 2

        g4 = self.fake_grid([2, 3, 4])
        g5 = self.fake_grid([4, 5, 6])
        assert g4.bins() == 2
        assert g5.bins() == 2

        with pytest.raises(ValueError, match="NonConsecutiveBins"):
            g2.merge(g4)

        with pytest.raises(ValueError, match="NonConsecutiveBins"):
            g2.merge(g5)
