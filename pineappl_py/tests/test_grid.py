import numpy as np
import pytest
import tempfile

from pineappl.bin import BinRemapper
from pineappl.boc import Channel
from pineappl.convolutions import Conv, ConvType
from pineappl.evolution import OperatorSliceInfo
from pineappl.fk_table import FkTable
from pineappl.grid import Grid, Order
from pineappl.import_subgrid import ImportSubgridV1
from pineappl.pids import PidBasis

# Construct the type of convolutions and the convolution object
# We assume unpolarized protons in the initial state
TYPECONV = ConvType(polarized=False, time_like=False)
CONVOBJECT = Conv(conv_type=TYPECONV, pid=2212)

# Construct the Channel and Order objetcs
UP_ANTIUP_CHANNEL = [([2, -2], 0.1)]
CHANNELS = [Channel(UP_ANTIUP_CHANNEL)]
ORDERS = [Order(3, 0, 0, 0, 0)]

# Testing specs for Convolution checks. Each element of the list is
# a tuple with two elements where the first element is a dictionary
# whose keys are the arguments of the `convolve` function and the
# second element is the expected results.
REF_VALUE = 5e6 / 9999
TESTING_SPECS = [
    (
        {
            "pdg_convs": [CONVOBJECT, CONVOBJECT],
            "xfxs": [lambda pid, x, q2: 0.0, lambda pid, x, q2: 0.0],
            "alphas": lambda q2: 0.0,
        },
        [0.0] * 2,
    ),  # fixed alphas(Q2) == 0.0
    (
        {
            "pdg_convs": [CONVOBJECT, CONVOBJECT],
            "xfxs": [lambda pid, x, q2: 1.0, lambda pid, x, q2: 1.0],
            "alphas": lambda q2: 1.0,
        },
        [REF_VALUE, 0.0],
    ),  # fixed alphas(Q2) == 1.0
    (
        {
            "pdg_convs": [CONVOBJECT, CONVOBJECT],
            "xfxs": [lambda pid, x, q2: 1.0, lambda pid, x, q2: 1.0],
            "alphas": lambda q2: 2.0,
        },
        [2**3 * REF_VALUE, 0.0],
    ),  # fixed alphas(Q2) == 2.0
    (
        {
            "pdg_convs": [CONVOBJECT, CONVOBJECT],
            "xfxs": [lambda pid, x, q2: 1.0, lambda pid, x, q2: 1.0],
            "alphas": lambda q2: 1.0,
            "bin_indices": [0],
        },
        [REF_VALUE],
    ),  # block first Bin without argument
    (
        {
            "pdg_convs": [CONVOBJECT, CONVOBJECT],
            "xfxs": [lambda pid, x, q2: 1.0, lambda pid, x, q2: 1.0],
            "alphas": lambda q2: 1.0,
            "bin_indices": [0],
            "order_mask": [False],
        },
        [0.0],
    ),  # block first Bin with order_mask
    (
        {
            "pdg_convs": [CONVOBJECT, CONVOBJECT],
            "xfxs": [lambda pid, x, q2: 1.0, lambda pid, x, q2: 1.0],
            "alphas": lambda q2: 1.0,
            "bin_indices": [0],
            "channel_mask": [False],
        },
        [0.0],
    ),  # block first Bin with channel_mask
    (
        {
            "pdg_convs": [CONVOBJECT, CONVOBJECT],
            "xfxs": [lambda pid, x, q2: 1.0, lambda pid, x, q2: 1.0],
            "alphas": lambda q2: 1.0,
            "bin_indices": [1],
        },
        [0.0],
    ),  # second Bin is empty
]

# Define the raw and target PIDS for testing the Evolution
EVOL_BASIS_PIDS = (
    22,
    100,
    21,
    200,
    203,
    208,
    215,
    224,
    235,
    103,
    108,
    115,
    124,
    135,
)

TARGET_PIDS = [22, -6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6]


class TestGrid:
    def test_init(self, fake_grids):
        g = fake_grids.grid_with_two_convolutions(
            channels=CHANNELS,
            orders=ORDERS,
            convolutions=[CONVOBJECT, CONVOBJECT],
        )
        assert isinstance(g, Grid)
        assert len(g.orders()) == 1
        assert g.orders()[0].as_tuple() == (3, 0, 0, 0, 0)

    def test_channels(self, fake_grids):
        g = fake_grids.grid_with_two_convolutions(
            channels=CHANNELS,
            orders=ORDERS,
            convolutions=[CONVOBJECT, CONVOBJECT],
        )
        assert len(g.channels()) == 1
        assert g.channels()[0] == UP_ANTIUP_CHANNEL

    def test_write(self, fake_grids):
        g = fake_grids.grid_with_two_convolutions(
            channels=CHANNELS,
            orders=ORDERS,
            convolutions=[CONVOBJECT, CONVOBJECT],
        )

        # Test writing/dumping the FK table into disk
        with tempfile.TemporaryDirectory() as tmpdir:
            g.write(f"{tmpdir}/toy_grid.pineappl")
            g.write_lz4(f"{tmpdir}/toy_grid.pineappl.lz4")

    def test_set_subgrid(self, fake_grids):
        g = fake_grids.grid_with_two_convolutions(
            channels=CHANNELS,
            orders=ORDERS,
            convolutions=[CONVOBJECT, CONVOBJECT],
        )

        # DIS grid
        xs = np.linspace(0.1, 1.0, 5)
        vs = np.random.rand(len(xs))
        subgrid = ImportSubgridV1(
            array=vs[np.newaxis, :, np.newaxis],
            node_values=[np.array([90.0]), xs, np.array([1.0])],
        )
        g.set_subgrid(0, 0, 0, subgrid.into())

        # let's mix it for fun with an hadronic one
        x1s = np.linspace(0.1, 1, 2)
        x2s = np.linspace(0.5, 1, 2)
        Q2s = np.linspace(10, 20, 2)
        subgrid = ImportSubgridV1(
            array=np.random.rand(len(Q2s), len(x1s), len(x2s)),
            node_values=[Q2s, x1s, x2s],
        )
        g.set_subgrid(0, 1, 0, subgrid.into())
        g.optimize()

    def test_bins(self, fake_grids):
        g = fake_grids.grid_with_two_convolutions(
            channels=CHANNELS,
            orders=ORDERS,
            convolutions=[CONVOBJECT, CONVOBJECT],
        )
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

    def test_grid(
        self,
        download_objects,
        gridname: str = "GRID_STAR_WMWP_510GEV_WP-AL-POL.pineappl.lz4",
    ):
        grid = download_objects(f"{gridname}")
        g = Grid.read(grid)

        # Get the types of convolutions for this grid
        for conv in g.convolutions():
            assert isinstance(conv, Conv)

        # Check that the scalings work, ie run without error
        # TODO: implement method to check the actual values
        g.scale(factor=10.0)
        g.scale_by_bin(factors=[10.0, 20.0])
        g.delete_bins(bin_indices=[0, 1, 2])

    def test_incosistent_convs(
        self,
        pdf,
        download_objects,
        gridname: str = "GRID_DYE906R_D_bin_1.pineappl.lz4",
    ):
        """Check that if the passed convolution types do not match the
        information in the grid the fail with `PanicException`.
        """
        grid = download_objects(f"{gridname}")
        g = Grid.read(grid)

        # The following grid has UNPOLARIZED proton, ie should be
        # `polarized=False`.
        h = ConvType(polarized=True, time_like=False)
        h_conv = Conv(conv_type=h, pid=2212)

        with pytest.raises(BaseException) as err_func:
            g.convolve(
                pdg_convs=[h_conv],  # Requires ONE single convolutions
                xfxs=[pdf.polarized_pdf],  # Requires ONE single PDF
                alphas=pdf.alphasQ,
            )
        assert "called `Option::unwrap()` on a `None` value" == str(
            err_func.value
        )

    @pytest.mark.parametrize("params,expected", TESTING_SPECS)
    def test_toy_convolution(self, fake_grids, params, expected):
        g = fake_grids.grid_with_two_convolutions(
            channels=CHANNELS,
            orders=ORDERS,
            convolutions=[CONVOBJECT, CONVOBJECT],
        )

        # Fill the subgrid-part of the GRID object
        xs = np.linspace(0.5, 1.0, 5)
        vs = xs.copy()
        subgrid = ImportSubgridV1(
            array=vs[np.newaxis, :, np.newaxis],
            node_values=[np.array([90.0]), xs, np.array([1.0])],
        )
        g.set_subgrid(0, 0, 0, subgrid.into())

        # Check the convolutions of the GRID
        np.testing.assert_allclose(g.convolve(**params), expected)

    def test_unpolarized_convolution(
        self,
        pdf,
        download_objects,
        gridname: str = "GRID_DYE906R_D_bin_1.pineappl.lz4",
    ):
        """Tes convolution with an actual Grid. In the following example,
        it is a DIS grid that involves a single unique hadron/proton.
        """
        expected_results = [
            +3.71019208e4,
            +3.71019208e4,
            +2.13727492e4,
            -1.83941398e3,
            +3.22728612e3,
            +5.45646897e4,
        ]

        grid = download_objects(f"{gridname}")
        g = Grid.read(grid)

        # Convolution object of the Unpolarized proton
        h = ConvType(polarized=False, time_like=False)
        h_conv = Conv(conv_type=h, pid=2212)

        np.testing.assert_allclose(
            g.convolve(
                pdg_convs=[h_conv],  # Requires ONE single convolutions
                xfxs=[pdf.polarized_pdf],  # Requires ONE single PDF
                alphas=pdf.alphasQ,
            ),
            expected_results,
        )

    def test_polarized_convolution(
        self,
        pdf,
        download_objects,
        gridname: str = "GRID_STAR_WMWP_510GEV_WP-AL-POL.pineappl.lz4",
    ):
        expected_results = [
            +5.50006832e6,
            +1.68117895e6,
            +3.08224445e5,
            -2.65602464e5,
            -1.04664085e6,
            -5.19002089e6,
        ]

        grid = download_objects(f"{gridname}")
        g = Grid.read(grid)

        # Convolution object of the 1st hadron - Polarized
        h1 = ConvType(polarized=True, time_like=False)
        h1_conv = Conv(conv_type=h1, pid=2212)

        # Convolution object of the 2nd hadron - Unpolarized
        h2 = ConvType(polarized=False, time_like=False)
        h2_conv = Conv(conv_type=h2, pid=2212)

        np.testing.assert_allclose(
            g.convolve(
                pdg_convs=[h1_conv, h2_conv],
                xfxs=[pdf.polarized_pdf, pdf.unpolarized_pdf],
                alphas=pdf.alphasQ,
            ),
            expected_results,
        )

    def test_evolve_with_two_ekos(
        self,
        pdf,
        download_objects,
        gridname: str = "GRID_STAR_WMWP_510GEV_WP-AL-POL.pineappl.lz4",
    ):
        """Test the evolution on a grid that contains two different convolutions,
        ie. requires two different EKOs.
        """
        grid = download_objects(f"{gridname}")
        g = Grid.read(grid)

        # Extract oder mask
        order_mask = Order.create_mask(g.orders(), 2, 2, True)
        evinfo = g.evolve_info(order_mask=order_mask)

        # Define the convolution types objects
        h1 = ConvType(polarized=True, time_like=False)
        h2 = ConvType(polarized=False, time_like=False)
        conv_type = [h1, h2]

        input_xgrid = np.geomspace(2e-7, 1, num=50)
        slices = []
        for conv_id, cvtype in enumerate(conv_type):
            sub_slices = []
            for q2 in evinfo.fac1:
                info = OperatorSliceInfo(
                    fac0=1.0,
                    fac1=q2,
                    x0=input_xgrid,
                    x1=evinfo.x1,
                    pids0=EVOL_BASIS_PIDS,
                    pids1=TARGET_PIDS,
                    pid_basis=PidBasis.Evol,
                    conv_type=cvtype,
                )
                op = np.random.uniform(
                    low=1,
                    high=10,
                    size=(
                        len(TARGET_PIDS),
                        evinfo.x1.size,
                        len(EVOL_BASIS_PIDS),
                        input_xgrid.size,
                    ),
                )
                sub_slices.append((info, op))
            slices.append(sub_slices)

        fktable = g.evolve(
            slices=slices,
            order_mask=order_mask,
            xi=(1.0, 1.0, 1.0),
            ren1=evinfo.fac1,
            alphas=[0.12029247510152144],
        )
        assert isinstance(fktable, FkTable)

    def test_io(self, tmp_path, fake_grids):
        g = fake_grids.grid_with_two_convolutions(
            channels=CHANNELS,
            orders=ORDERS,
            convolutions=[CONVOBJECT, CONVOBJECT],
        )
        p = tmp_path / "test.pineappl"
        p.write_text("")
        g.write(str(p))
        gg = Grid.read(p)
        assert isinstance(gg, Grid)
        _ = Grid.read(str(p))

    def test_set_key_value(self, fake_grids):
        g = fake_grids.grid_with_two_convolutions(
            channels=CHANNELS,
            orders=ORDERS,
            convolutions=[CONVOBJECT, CONVOBJECT],
        )
        g.set_key_value("bla", "blub")
        g.set_key_value('"', "'")
        g.set_key_value("äöü", "ß\\")

    def test_fill(self, fake_grids):
        g = fake_grids.grid_with_two_convolutions(
            channels=CHANNELS,
            orders=ORDERS,
            convolutions=[CONVOBJECT, CONVOBJECT],
        )
        # Fill the Grid with some values
        n_tuple = [10.0, 0.5, 0.5]
        g.fill(
            order=0,
            observable=0.01,
            channel=0,
            ntuple=n_tuple,
            weight=10,
        )
        # Peform convolutions using Toy LHPDF & AlphasQ2 functions
        res = g.convolve(
            pdg_convs=[CONVOBJECT, CONVOBJECT],
            xfxs=[lambda pid, x, q2: x, lambda pid, x, q2: x],
            alphas=lambda q2: 1.0,
        )
        np.testing.assert_allclose(res, [0.0, 0.0])

    def test_fill_array(self, fake_grids):
        g = fake_grids.grid_with_two_convolutions(
            channels=CHANNELS,
            orders=ORDERS,
            convolutions=[CONVOBJECT, CONVOBJECT],
        )
        g.fill_array(
            order=0,
            observables=np.array([1e-3, 1e-2]),
            channel=0,
            ntuples=[
                np.array([15.0, 1.0, 1.0]),
                np.array([15.0, 1.0, 1.0]),
                np.array([15.0, 1.0, 1.0]),
            ],
            weights=np.array([10.0, 100.0]),
        )

        # Convolution of two symmetrical hadrons
        h = ConvType(polarized=False, time_like=False)
        h_conv = Conv(conv_type=h, pid=2212)
        res = g.convolve(
            pdg_convs=[h_conv, h_conv],
            xfxs=[lambda pid, x, q2: x, lambda pid, x, q2: x],
            alphas=lambda q2: 1.0,
        )
        np.testing.assert_allclose(res, [0.0, 0.0])

    def test_fill_all(self, fake_grids):
        g = fake_grids.grid_with_two_convolutions(
            channels=CHANNELS,
            orders=ORDERS,
            convolutions=[CONVOBJECT, CONVOBJECT],
        )
        g.fill_all(
            order=0,
            observable=1e-2,
            ntuple=[1.0, 1.0, 1.0],
            weights=np.array([10.0]),
        )

        # Convolution of two symmetrical hadrons
        h = ConvType(polarized=False, time_like=False)
        h_conv = Conv(conv_type=h, pid=2212)
        res = g.convolve(
            pdg_convs=[h_conv, h_conv],
            xfxs=[lambda pid, x, q2: x, lambda pid, x, q2: x],
            alphas=lambda q2: 1.0,
        )
        np.testing.assert_allclose(res, [0.0, 0.0])

    def test_merge(self, fake_grids):
        g0 = fake_grids.grid_with_two_convolutions(
            channels=CHANNELS,
            orders=ORDERS,
            convolutions=[CONVOBJECT, CONVOBJECT],
            bins=[1, 2, 3],
        )
        g1 = fake_grids.grid_with_two_convolutions(
            channels=CHANNELS,
            orders=ORDERS,
            convolutions=[CONVOBJECT, CONVOBJECT],
            bins=[3, 4, 5],
        )
        assert g0.bins() == 2
        assert g1.bins() == 2
        g0.merge(g1)
        assert g0.bins() == 4

        g2 = fake_grids.grid_with_two_convolutions(
            channels=CHANNELS,
            orders=ORDERS,
            convolutions=[CONVOBJECT, CONVOBJECT],
            bins=[1, 2, 3],
        )
        g3 = fake_grids.grid_with_two_convolutions(
            channels=CHANNELS,
            orders=ORDERS,
            convolutions=[CONVOBJECT, CONVOBJECT],
            bins=[1, 2, 3],
        )
        assert g2.bins() == 2
        assert g3.bins() == 2

        g2.merge(g3)
        assert g2.bins() == 2

        g4 = fake_grids.grid_with_two_convolutions(
            channels=CHANNELS,
            orders=ORDERS,
            convolutions=[CONVOBJECT, CONVOBJECT],
            bins=[2, 3, 4],
        )
        g5 = fake_grids.grid_with_two_convolutions(
            channels=CHANNELS,
            orders=ORDERS,
            convolutions=[CONVOBJECT, CONVOBJECT],
            bins=[4, 5, 6],
        )
        assert g4.bins() == 2
        assert g5.bins() == 2

        with pytest.raises(ValueError, match="NonConsecutiveBins"):
            g2.merge(g4)

        with pytest.raises(ValueError, match="NonConsecutiveBins"):
            g2.merge(g5)

    def test_evolveinfo(
        self,
        download_objects,
        gridname: str = "GRID_STAR_WMWP_510GEV_WP-AL-POL.pineappl.lz4",
    ):
        grid = download_objects(f"{gridname}")
        g = Grid.read(grid)
        g_evinfo = g.evolve_info(order_mask=[True, False, False, False])

        np.testing.assert_allclose(g_evinfo.fac1, [6463.838404])
        np.testing.assert_allclose(g_evinfo.ren1, [6463.838404])
        np.testing.assert_allclose(g_evinfo.pids1, [-5, -3, -1, 2, 4])
        assert g_evinfo.x1.size == 23
        np.testing.assert_allclose(g_evinfo.x1[0], 0.01437507)
