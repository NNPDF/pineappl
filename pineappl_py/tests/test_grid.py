import numpy as np
import pytest
from pineappl.bin import BinRemapper
from pineappl.boc import Channel, Kinematics
from pineappl.convolutions import Conv, ConvType
from pineappl.grid import Grid, Order
from pineappl.import_subgrid import ImportSubgridV1
from pineappl.interpolation import Interp
from pineappl.pids import PidBasis

# Construct the type of convolutions and the convolution object
# We assume unpolarized proton PDF
TYPECONV = ConvType(polarized=False, time_like=False)
CONVOBJECT = Conv(conv_type=TYPECONV, pid=2212)

# Construct the Channel and Order objetcs
UP_ANTIUP_CHANNEL = [([2, -2], 0.1)]
CHANNELS = [Channel(UP_ANTIUP_CHANNEL)]
ORDERS = [Order(3, 0, 0, 0, 0)]


class TestGrid:
    def fake_grid(
        self,
        channels: list[Channel] = CHANNELS,
        orders: list[Order] = ORDERS,
        bins: list[float] = [1e-7, 1e-3, 1],
    ) -> Grid:
        # We assume symmetrical proton-proton in the initial state
        convolutions = [CONVOBJECT, CONVOBJECT]
        # Define the kinematics. Kinematics are defined as a list of items.
        # 1st item: factorization and renormalization scale
        # 2nd item: parton momentum fraction of the 1st convolution
        # 3rd tiem: parton momentum fraction of the 2nd convolution
        kinematics = [
            Kinematics(0),  # Scale
            Kinematics(1),  # x1 momentum fraction
            Kinematics(2),  # x2 momentum fraction
        ]
        # Define the interpolation specs for each item of the Kinematics
        interpolations = [
            Interp(
                min=1e2,
                max=1e8,
                nodes=40,
                order=3,
            ),  # Interpolation on the Scale
            Interp(
                min=2e-7,
                max=1.0,
                nodes=50,
                order=3,
            ),  # Interpolation on x1 momentum fraction
            Interp(
                min=2e-7,
                max=1.0,
                nodes=50,
                order=3,
            ),  # Interpolation on x2 momentum fraction
        ]
        bin_limits = np.array(bins)
        return Grid(
            pid_basis=PidBasis.Evol,
            channels=channels,
            orders=orders,
            bin_limits=bin_limits,
            convolutions=convolutions,
            interpolations=interpolations,
            kinematics=kinematics,
        )

    def test_init(self):
        g = self.fake_grid()
        assert isinstance(g, Grid)
        assert len(g.orders()) == 1
        assert g.orders()[0].as_tuple() == (3, 0, 0, 0, 0)

    def test_channels(self):
        g = self.fake_grid()
        assert len(g.channels()) == 1
        assert g.channels()[0] == UP_ANTIUP_CHANNEL

    def test_set_subgrid(self):
        g = self.fake_grid()

        # DIS grid
        xs = np.linspace(0.1, 1.0, 5)
        vs = np.random.rand(len(xs))
        subgrid = ImportSubgridV1(
            vs[np.newaxis, :, np.newaxis],
            np.array([90.0]),
            xs,
            np.array([1.0]),
        )
        g.set_subgrid(0, 0, 0, subgrid.into())

        # let's mix it for fun with an hadronic one
        x1s = np.linspace(0.1, 1, 2)
        x2s = np.linspace(0.5, 1, 2)
        Q2s = np.linspace(10, 20, 2)
        subgrid = ImportSubgridV1(
            np.random.rand(len(Q2s), len(x1s), len(x2s)), Q2s, x1s, x2s
        )
        g.set_subgrid(0, 1, 0, subgrid.into())
        g.optimize()

    def test_bins(self):
        g = self.fake_grid()
        # 1D
        normalizations = np.array([1.0, 1.0])
        limits = [(1, 1), (2, 2)]
        remapper = BinRemapper(normalizations, limits)
        g.set_remapper(remapper)
        assert g.bin_dimensions() == 1
        np.testing.assert_allclose(g.bin_left(0), [1, 2])
        np.testing.assert_allclose(g.bin_right(0), [1, 2])
        # 2D
        limits = [(1, 2), (2, 3), (2, 4), (3, 5)]
        remapper = BinRemapper(normalizations, limits)
        g.set_remapper(remapper)
        assert g.bin_dimensions() == 2
        np.testing.assert_allclose(g.bin_left(0), [1, 2])
        np.testing.assert_allclose(g.bin_right(0), [2, 4])
        np.testing.assert_allclose(g.bin_left(1), [2, 3])
        np.testing.assert_allclose(g.bin_right(1), [3, 5])

    def test_convolve_with_two(self):
        g = self.fake_grid()

        # Fill the subgrid-part of the GRID object
        xs = np.linspace(0.5, 1.0, 5)
        vs = xs.copy()
        subgrid = ImportSubgridV1(
            array=vs[np.newaxis, :, np.newaxis],
            scales=np.array([90.0]),
            x1_grid=xs,
            x2_grid=np.array([1.0]),
        )
        g.set_subgrid(0, 0, 0, subgrid.into())

        # Check the convolutions of the GRID
        np.testing.assert_allclose(
            g.convolve_with_two(
                pdg_conv1=CONVOBJECT,
                xfx1=lambda pid, x, q2: 0.0,
                pdg_conv2=CONVOBJECT,
                xfx2=lambda pid, x, q2: 0.0,
                alphas=lambda q2: 0.0,
            ),
            [0.0] * 2,
        )
        np.testing.assert_allclose(
            g.convolve_with_two(
                pdg_conv1=CONVOBJECT,
                xfx1=lambda pid, x, q2: 1.0,
                pdg_conv2=CONVOBJECT,
                xfx2=lambda pid, x, q2: 1.0,
                alphas=lambda q2: 1.0,
            ),
            [5e6 / 9999, 0.0],
        )
        np.testing.assert_allclose(
            g.convolve_with_two(
                pdg_conv1=CONVOBJECT,
                xfx1=lambda pid, x, q2: 1.0,
                pdg_conv2=CONVOBJECT,
                xfx2=lambda pid, x, q2: 1.0,
                alphas=lambda q2: 2.0,
            ),
            [2**3 * 5e6 / 9999, 0.0],
        )

    def test_io(self, tmp_path):
        g = self.fake_grid()
        p = tmp_path / "test.pineappl"
        p.write_text("")
        g.write(str(p))
        gg = Grid.read(p)
        assert isinstance(gg, Grid)
        _ = Grid.read(str(p))

    def test_fill(self):
        g = self.fake_grid()
        # Fill the Grid with some values
        n_tuple = [0.5, 0.5, 10.0]
        g.fill(
            order=0,
            observable=0.01,
            channel=0,
            ntuple=n_tuple,
            weight=10,
        )
        # Peform convolutions using Toy LHPDF & AlphasQ2 functions
        res = g.convolve_with_two(
            pdg_conv1=CONVOBJECT,
            xfx1=lambda pid, x, q2: x,
            pdg_conv2=CONVOBJECT,
            xfx2=lambda pid, x, q2: x,
            alphas=lambda q2: 1.0,
        )
        pytest.approx(res) == 0.0

    def test_merge(self):
        g = self.fake_grid(bins=[1, 2, 3])
        g1 = self.fake_grid(bins=[3, 4, 5])
        assert g.bins() == 2
        assert g1.bins() == 2

        g.merge(g1)
        assert g.bins() == 4

        g2 = self.fake_grid(bins=[1, 2, 3])
        g3 = self.fake_grid(bins=[1, 2, 3])
        assert g2.bins() == 2
        assert g3.bins() == 2

        g2.merge(g3)
        assert g2.bins() == 2

        g4 = self.fake_grid(bins=[2, 3, 4])
        g5 = self.fake_grid(bins=[4, 5, 6])
        assert g4.bins() == 2
        assert g5.bins() == 2

        with pytest.raises(ValueError, match="NonConsecutiveBins"):
            g2.merge(g4)

        with pytest.raises(ValueError, match="NonConsecutiveBins"):
            g2.merge(g5)
