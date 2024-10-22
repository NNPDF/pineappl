"""Test module for the interface of the `fk_table`.

It checks the cases in which we have one, two, and
three (general) convolutions.
"""

import numpy as np
import tempfile
from pineappl.boc import Channel, Kinematics
from pineappl.convolutions import Conv, ConvType
from pineappl.fk_table import FkTable, FkAssumptions
from pineappl.grid import Grid, Order
from pineappl.import_subgrid import ImportSubgridV1
from pineappl.interpolation import Interp
from pineappl.pids import PidBasis
from typing import List


class TestFkTable:
    def fake_grid(
        self,
        channels: List[Channel],
        orders: List[Order],
        convolutions: List[Conv],
        bins: List[float] = [1e-7, 1e-3, 1],
    ) -> Grid:
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
                reweight_meth="noreweight",
                map="applgrid_h0",
            ),  # Interpolation on the Scale
            Interp(
                min=2e-7,
                max=1.0,
                nodes=50,
                order=3,
                reweight_meth="applgrid",
                map="applgrid_f2",
            ),  # Interpolation on x1 momentum fraction
            Interp(
                min=2e-7,
                max=1.0,
                nodes=50,
                order=3,
                reweight_meth="applgrid",
                map="applgrid_f2",
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

    def test_convolve(self):
        # Define convolution types and the initial state hadrons
        # We consider an initial state Polarized Proton
        h = ConvType(polarized=True, time_like=False)
        h_conv = Conv(conv_type=h, pid=2212)
        # The length of the convolutions has to match the nb of hadrons
        convolutions = [h_conv]
        # We define the PIDs of the partons out of the Proton
        down_channel = [([1], 1.0)]  # DIS-case
        up_channel = [([2], 1.0)]  # DIS-case
        channels = [Channel(down_channel), Channel(up_channel)]
        # Now we define the perturbative orders
        orders = [Order(0, 0, 0, 0, 0)]
        g = self.fake_grid(channels, orders, convolutions)

        # DIS grid
        xs = np.linspace(0.5, 1.0, 5)
        vs = xs.copy()
        subgrid = ImportSubgridV1(
            vs[np.newaxis, :],  # DIS shape: (len(q2), len(x_grid))
            np.array([90.0]),
            xs,
        )
        g.set_subgrid(0, 0, 0, subgrid.into())

        # Convert the Grid -> FkTable
        fk = FkTable(g)

        # Test a simple convolution of the FK table
        np.testing.assert_allclose(
            fk.convolve(
                pdg_convs=[h_conv],
                xfxs=[lambda pid, x, q2: 0.0],
            ),
            [0.0] * 2,
        )
        np.testing.assert_allclose(
            fk.convolve(
                pdg_convs=[h_conv],
                xfxs=[lambda pid, x, q2: 1.0],
            ),
            [5e7 / 9999, 0.0],
        )

        # Test writing/dumping the FK table into disk
        with tempfile.TemporaryDirectory() as tmpdir:
            fk.write(f"{tmpdir}/toy_fktable.pineappl")
            fk.write_lz4(f"{tmpdir}/toy_fktable.pineappl.lz4")

    def test_fktable(
        self,
        download_objects,
        fkname: str = "FKTABLE_CMSTTBARTOT8TEV-TOPDIFF8TEVTOT.pineappl.lz4",
    ):
        fk_table = download_objects(f"{fkname}")
        fk = FkTable.read(fk_table)

        assert fk.table().shape == (1, 51, 34, 34)
        np.testing.assert_allclose(fk.muf2(), 2.7224999999999997)

        # Check the various aspects of the Bins
        assert fk.bins() == 1
        np.testing.assert_allclose(fk.bin_normalizations(), [1.0])
        np.testing.assert_allclose(fk.bin_right(dimension=0), [1.0])
        np.testing.assert_allclose(fk.bin_left(dimension=0), [0.0])

        # Check the various aspects of the Channels
        channels = fk.channels()
        assert len(channels) == 51
        assert [21, 200] in channels

        # Check the contents of the x-grid
        x_grid = fk.x_grid()
        assert x_grid.size == 34
        np.testing.assert_allclose(x_grid[0], 1.57456056e-04)

        # Test FK optimization
        assumption = FkAssumptions("Nf6Sym")
        fk.optimize(assumption)

    def test_unpolarized_convolution(
        self,
        pdf,
        download_objects,
        fkname: str = "FKTABLE_CMSTTBARTOT8TEV-TOPDIFF8TEVTOT.pineappl.lz4",
    ):
        """Check the convolution of an actual FK table that involves two
        symmetrical unpolarized protons:
        """
        expected_results = [3.72524538e04]
        fk_table = download_objects(f"{fkname}")
        fk = FkTable.read(fk_table)

        # Convolution object of the 1st hadron - Polarized
        h = ConvType(polarized=False, time_like=False)
        h_conv = Conv(conv_type=h, pid=2212)

        np.testing.assert_allclose(
            fk.convolve(
                pdg_convs=[h_conv, h_conv],
                xfxs=[pdf.unpolarized_pdf, pdf.unpolarized_pdf],
            ),
            expected_results,
        )

    def test_polarized_convolution(
        self,
        pdf,
        download_objects,
        fkname: str = "FKTABLE_STAR_WMWP_510GEV_WM-AL-POL.pineappl.lz4",
    ):
        """Check the convolution of an actual FK table that involves two
        different initial states:
            - 1st hadron: polarized proton
            - 2nd hadron: unpolarized proton
        """
        expected_results = [
            -1.00885071e6,
            -2.40862657e5,
            -1.66407218e5,
            -2.96098362e5,
            -5.67594297e5,
            +6.59245015e4,
        ]
        fk_table = download_objects(f"{fkname}")
        fk = FkTable.read(fk_table)

        # Convolution object of the 1st hadron - Polarized
        h1 = ConvType(polarized=True, time_like=False)
        h1_conv = Conv(conv_type=h1, pid=2212)

        # Convolution object of the 2nd hadron - Unpolarized
        h2 = ConvType(polarized=False, time_like=False)
        h2_conv = Conv(conv_type=h2, pid=2212)

        np.testing.assert_allclose(
            fk.convolve(
                pdg_convs=[h1_conv, h2_conv],
                xfxs=[pdf.polarized_pdf, pdf.unpolarized_pdf],
            ),
            expected_results,
        )
