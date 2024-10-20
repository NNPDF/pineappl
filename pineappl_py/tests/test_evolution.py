"""Test module for the interface of the `evolution`.

It checks the cases in which we have evolve with one,
two, and three (general) EKOs.
"""

import itertools

import numpy as np
from pineappl.boc import Channel, Kinematics
from pineappl.convolutions import Conv, ConvType
from pineappl.evolution import OperatorSliceInfo, EvolveInfo
from pineappl.grid import Grid, Order
from pineappl.interpolation import Interp
from pineappl.pids import PidBasis
from typing import List


class TestEvolution:
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
                min=1.0,
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

    def test_evolveinfo(self):
        evinfo = EvolveInfo(
            fac1=[0.5, 1.0, 2.0],
            pids1=[-2, 0, 2],
            x1=[1e-3, 0.5, 1],
            ren1=[0.5, 1.0, 2.0],
        )
        np.testing.assert_array_equal(evinfo.fac1, [0.5, 1.0, 2.0])
        np.testing.assert_array_equal(evinfo.pids1, [-2, 0, 2])
        np.testing.assert_array_equal(evinfo.x1, [1e-3, 0.5, 1.0])
        np.testing.assert_array_equal(evinfo.fac1, [0.5, 1.0, 2.0])

    def test_with_one_eko(self):
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

        # Construct the Grid and fill with some values
        grid = self.fake_grid(channels, orders, convolutions)

        x1g = np.linspace(0.5, 1.0, 5)
        x2g = x1g.copy()
        q2g = np.array([10, 90, 100])

        for x1, x2, q2 in itertools.product(x1g, x2g, q2g):
            grid.fill(
                order=0,
                observable=0.01,
                channel=0,
                ntuple=[x1, x2, q2],
                weight=10,
            )

        # Check the Evolution of the Grid
        info = OperatorSliceInfo(
            fac0=1.0,
            pids0=[],
            x0=[],
            fac1=1.0,
            pids1=[],
            x1=[],
            pid_basis=PidBasis.Pdg,
            conv_type=h,
        )

        # TODO: check a Toy evolution
        assert isinstance(info, OperatorSliceInfo)
