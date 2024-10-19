import numpy as np
import pytest
from pineappl.boc import Channel, Kinematics, Order
from pineappl.convolutions import Conv, ConvType
from pineappl.grid import Grid
from pineappl.import_subgrid import ImportSubgridV1
from pineappl.interpolation import Interp
from pineappl.pids import PidBasis
from pineappl.subgrid import SubgridEnum

# Define some default for the minimum value of `Q2`
Q2_MIN = 1e2

# See `test_grid.py` for more detailed information
TYPECONV = ConvType(polarized=False, time_like=False)
CONVOBJECT = Conv(conv_type=TYPECONV, pid=2212)


def fake_dis_grid(
    orders: Order, channels: Channel, q2_min: float = Q2_MIN
) -> Grid:
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
            min=q2_min,
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
    return Grid(
        pid_basis=PidBasis.Evol,
        channels=channels,
        orders=orders,
        bin_limits=bin_limits,
        convolutions=convolutions,
        interpolations=interpolations,
        kinematics=kinematics,
    )


def test_issue_164(pdf):
    # https://github.com/NNPDF/pineappl/issues/164
    # DIS-like convolution now ONLY requires one entry of `PID`
    channels = [Channel([([2], 1.0)])]  # DIS-case
    orders = [Order(0, 0, 0, 0, 0)]

    def convolve_grid(q2_min: float = Q2_MIN) -> Grid:
        grid = fake_dis_grid(orders, channels, q2_min)
        # Fill the Grid with some values
        grid.fill(
            order=0,
            observable=0.5,
            channel=0,
            ntuple=[0.2, 0.2, 10],
            weight=0.5,
        )
        return grid.convolve(
            pdg_convs=[CONVOBJECT],
            xfxs=[pdf.xfxQ],
            alphas=pdf.alphasQ,
        )

    # Using default minimum
    res = convolve_grid()
    assert res == 0.0
    # lower minimum to q2=1
    res = convolve_grid(q2_min=1.0)
    assert res == 0.0


class TestSubgrid:
    def fake_grid(self):
        channels = [Channel([([2], 1.0)]), Channel([([3], 0.5)])]
        orders = [Order(0, 0, 0, 0, 0)]
        return fake_dis_grid(orders, channels)

    def fake_importonlysubgrid(self) -> tuple:
        x1s = np.linspace(0.1, 1, 2)
        x2s = np.linspace(0.5, 1, 2)
        Q2s = np.linspace(10, 20, 2)
        scale = [q2 for q2 in Q2s]
        array = np.random.rand(len(Q2s), len(x1s), len(x2s))
        subgrid = ImportSubgridV1(
            array=array,
            scales=scale,
            x1_grid=x1s,
            x2_grid=x2s,
        )
        return subgrid, [x1s, x2s, scale, array]

    def test_subgrid_methods(self):
        # TODO: extract the values of the scales and x grids
        grid = self.fake_grid()
        test_subgrid, infos = self.fake_importonlysubgrid()
        x1s, x2s, mu2s, _ = (obj for obj in infos)
        grid.set_subgrid(0, 0, 0, test_subgrid.into())
        extr_subgrid = grid.subgrid(0, 0, 0)
        assert isinstance(extr_subgrid, SubgridEnum)

    @pytest.mark.skip(reason="No implementation of Array3 for subgrid.")
    def test_to_array3(self):
        # TODO: extract and check the dense array of the subgrid
        # requires `impl From<&SubgridEnum> for Array3<f64>`
        grid = self.fake_grid()
        test_subgrid, infos = self.fake_importonlysubgrid()
        _, _, _, array = (obj for obj in infos)
        grid.set_subgrid(0, 0, 0, test_subgrid.into())
        extr_subgrid = grid.subgrid(0, 0, 0)
        test_array = extr_subgrid.to_array3()
        print(test_array)
        print(array)
        np.testing.assert_allclose(test_array, array)
