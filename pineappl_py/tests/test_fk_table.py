import numpy as np

import pineappl


class TestFkTable:
    def fake_grid(self, bins=None):
        channels = [pineappl.boc.Channel([(1, 21, 1.0)])]
        orders = [pineappl.grid.Order(0, 0, 0, 0)]
        bin_limits = np.array([1e-7, 1e-3, 1] if bins is None else bins, dtype=float)
        subgrid_params = pineappl.subgrid.SubgridParams()
        g = pineappl.grid.Grid(channels, orders, bin_limits, subgrid_params)
        return g

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
        g.set_subgrid(0, 0, 0, subgrid.into())
        fk = pineappl.fk_table.FkTable(g)
        np.testing.assert_allclose(
            fk.convolve_with_one(2212, lambda pid, x, q2: 0.0),
            [0.0] * 2,
        )
        np.testing.assert_allclose(
            fk.convolve_with_one(2212, lambda pid, x, q2: 1),
            [5e7 / 9999, 0.0],
        )

        info = pineappl.evolution.OperatorSliceInfo(
            1.0, [], [], 1.0, [], [], pineappl.pids.PidBasis.Pdg
        )

        # TODO: write a better test
        try:
            g.evolve_with_slice_iter(
                iter(
                    [
                        (info, np.ndarray([0, 0, 0, 0])),
                        (info, np.ndarray([0, 0, 0, 0])),
                    ]
                ),
                np.array([], dtype=bool),
                (1.0, 1.0),
                [],
                [],
            )

            assert False
        except:  # noqa: E722
            assert True

        # TODO: write a better test
        try:
            g.evolve_with_slice_iter2(
                iter(
                    [
                        (info, np.ndarray([0, 0, 0, 0])),
                        (info, np.ndarray([0, 0, 0, 0])),
                    ]
                ),
                iter(
                    [
                        (info, np.ndarray([0, 0, 0, 0])),
                        (info, np.ndarray([0, 0, 0, 0])),
                    ]
                ),
                np.array([], dtype=bool),
                (1.0, 1.0),
                [],
                [],
            )

            assert False
        except:  # noqa: E722
            assert True
