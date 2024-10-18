"""Test module for the interface of the `fk_table`.

It checks the cases in which we have one, two, and
three (general) convolutions.
"""

import numpy as np
from pineappl.boc import Channel, Kinematics
from pineappl.convolutions import Conv, ConvType
from pineappl.fk_table import FkTable
from pineappl.grid import Grid, Order
from pineappl.import_subgrid import ImportSubgridV1
from pineappl.interpolation import Interp
from pineappl.pids import PidBasis


class TestFkTable:
    def fake_grid(
        self,
        channels: list[Channel],
        orders: list[Order],
        convolutions: list[Conv],
        bins: list[float] = [1e-7, 1e-3, 1],
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

    def test_convolve_with_one(self):
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
            fk.convolve_with_one(
                pdg_conv=h_conv,
                xfx=lambda pid, x, q2: 0.0,
            ),
            [0.0] * 2,
        )
        np.testing.assert_allclose(
            fk.convolve_with_one(
                pdg_conv=h_conv,
                xfx=lambda pid, x, q2: 1.0,
            ),
            [5e7 / 9999, 0.0],
        )

    def test_convolve_with_two(self):
        # TODO
        pass

    def test_convolve_with_many(self):
        # TODO
        pass
