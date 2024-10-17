import pytest
import numpy as np

from pineappl.pids import PidBasis
from pineappl.boc import Channel, Order, Kinematics
from pineappl.convolutions import Conv, ConvType
from pineappl.grid import Grid
from pineappl.subgrid import SubgridParams
from pineappl.interpolation import Interp


# See `test_grid.py` for more detailed information
TYPECONV = ConvType(polarized=False, time_like=False)
CONVOBJECT = Conv(conv_type=TYPECONV, pid=2212)


class TestSubgridParams:
    def test_init(self):
        sp = SubgridParams()
        assert isinstance(sp, SubgridParams)


def test_issue_164(pdf):
    # https://github.com/NNPDF/pineappl/issues/164
    channels = [Channel([([1, 2], 1.0)])]
    orders = [Order(0, 0, 0, 0, 0)]
    params = SubgridParams()

    def convolve_grid():
        bin_limits = np.array([0.0, 1.0])
        # See `test_grid.py` for more detailed information
        # on the meaning of the following parameters
        convolutions = [CONVOBJECT]  # Consider DIS-case
        kinematics = [
            Kinematics(0),  # Scale
            Kinematics(1),  # x1 momentum fraction
            Kinematics(2),  # x2 momentum fraction
        ]
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
        grid = Grid(
            pid_basis=PidBasis.Evol,
            channels=channels,
            orders=orders,
            bin_limits=bin_limits,
            convolutions=convolutions,
            interpolations=interpolations,
            kinematics=kinematics,
        )
        grid.fill(
            order=0,
            observable=0.5,
            channel=0,
            ntuple=[0.2, 0.2, 10],
            weight=0.5,
        )
        return grid.convolve_with_one(
            pdg_conv=CONVOBJECT,
            xfx=pdf.xfxQ,
            alphas=pdf.alphasQ,
        )

    # default minimum is q2=100
    res = convolve_grid()
    assert res == 0.0

    # lower minimum to q2=1
    params.set_q2_min(1.0)
    res = convolve_grid()
    assert res == 0.0
