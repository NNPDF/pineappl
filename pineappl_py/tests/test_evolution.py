"""Test module for the interface of the `evolution`.

It checks the cases in which we have evolve with one,
two, and three (general) EKOs.
"""

import numpy as np
from pineappl.convolutions import ConvType
from pineappl.evolution import EvolveInfo, OperatorSliceInfo
from pineappl.pids import PidBasis


class TestEvolution:
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

    def test_init_operatorsliceinfo(self):
        info = OperatorSliceInfo(
            fac0=1.0,
            pids0=[],
            x0=[],
            fac1=1.0,
            pids1=[],
            x1=[],
            pid_basis=PidBasis.Pdg,
            conv_type=ConvType(polarized=False, time_like=False),
        )

        assert isinstance(info, OperatorSliceInfo)
