"""Test module for the interface of the `fk_table`.

It checks the cases in which we have one, two, and
three (general) convolutions.
"""

import numpy as np
import pytest
import subprocess
from pineappl.boc import Channel, Kinematics
from pineappl.convolutions import Conv, ConvType
from pineappl.fk_table import FkTable
from pineappl.grid import Grid, Order
from pineappl.import_subgrid import ImportSubgridV1
from pineappl.interpolation import Interp
from pineappl.pids import PidBasis
from typing import List


@pytest.fixture
def download_fktable(tmp_path_factory):
    def _download_fk(fkname: str) -> None:
        download_dir = tmp_path_factory.mktemp("data")
        file_path = download_dir / f"{fkname}"
        args = [
            "wget",
            "--no-verbose",
            "--no-clobber",
            "-P",
            f"{download_dir}",
            f"https://data.nnpdf.science/pineappl/test-data/{fkname}",
        ]

        try:
            _ = subprocess.run(
                args,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            return file_path
        except OSError as error:
            msg = f"Failed to execute the command {args}."
            raise EnvironmentError(msg) from error

    return _download_fk


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
            vs[np.newaxis, :, np.newaxis],
            np.array([90.0]),
            xs,
            np.array([1.0]),
        )
        g.set_subgrid(0, 0, 0, subgrid.into())
        fk = FkTable(g)  # Convert Grid -> FkTable
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

    def test_unpolarized_convolution(
        self,
        download_fktable,
        fkname: str = "CMSTTBARTOT8TEV-TOPDIFF8TEVTOT.pineappl.lz4",
    ):
        """Check the convolution of an actual FK table that involves two
        symmetrical unpolarized protons:
        """
        expected_results = [3.72524538e04]
        fk_table = download_fktable(f"{fkname}")
        fk = FkTable.read(fk_table)

        # Convolution object of the 1st hadron - Polarized
        h = ConvType(polarized=False, time_like=False)
        h_conv = Conv(conv_type=h, pid=2212)

        # Define the Toy Unpolarized PDF set
        def _unpolarized_pdf(pid, x, q2):
            return 1.0

        np.testing.assert_allclose(
            fk.convolve(
                pdg_convs=[h_conv, h_conv],
                xfxs=[_unpolarized_pdf, _unpolarized_pdf],
            ),
            expected_results,
        )

    def test_polarized_convolution(
        self,
        download_fktable,
        fkname: str = "GRID_STAR_WMWP_510GEV_WM-AL-POL.pineappl.lz4",
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
        fk_table = download_fktable(f"{fkname}")
        fk = FkTable.read(fk_table)

        # Convolution object of the 1st hadron - Polarized
        h1 = ConvType(polarized=True, time_like=False)
        h1_conv = Conv(conv_type=h1, pid=2212)

        # Convolution object of the 2nd hadron - Unpolarized
        h2 = ConvType(polarized=False, time_like=False)
        h2_conv = Conv(conv_type=h2, pid=2212)

        # Define the Toy Polarized PDF set
        def _polarized_pdf(pid, x, q2):
            return 2.0

        # Define the Toy Unpolarized PDF set
        def _unpolarized_pdf(pid, x, q2):
            return 1.0

        np.testing.assert_allclose(
            fk.convolve(
                pdg_convs=[h1_conv, h2_conv],
                xfxs=[_polarized_pdf, _unpolarized_pdf],
            ),
            expected_results,
        )
