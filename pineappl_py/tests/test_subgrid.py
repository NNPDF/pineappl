import numpy as np
import pytest

from dataclasses import dataclass
from typing import List, Tuple

from pineappl.boc import Channel, Order
from pineappl.convolutions import Conv, ConvType
from pineappl.grid import Grid
from pineappl.import_subgrid import ImportSubgridV1
from pineappl.subgrid import SubgridEnum

# Define some default for the minimum value of `Q2`
Q2_MIN = 1e2

# See `test_grid.py` for more detailed information
TYPECONV = ConvType(polarized=False, time_like=False)
CONVOBJECT = Conv(convolution_types=TYPECONV, pid=2212)


@dataclass
class OperatorInfo:
    x_grids: List[np.ndarray]
    scale: List[float]
    array: np.ndarray


def test_issue_164(pdf, fake_grids):
    # https://github.com/NNPDF/pineappl/issues/164
    # DIS-like convolution now ONLY requires one entry of `PID`
    channels = [Channel([([2], 1.0)])]  # DIS-case
    orders = [Order(0, 0, 0, 0, 0)]

    def convolve_grid(q2_min: float = Q2_MIN) -> np.ndarray:
        grid = fake_grids.grid_with_generic_convolution(
            nb_convolutions=1,
            orders=orders,
            channels=channels,
            convolutions=[CONVOBJECT],
            q2min=q2_min,
            bins=np.array([0.0, 0.1]),
        )
        # Fill the Grid with some values
        grid.fill(
            order=0,
            observable=0.5,
            channel=0,
            ntuple=[10.0, 0.2, 0.2],
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
    def fake_grid(self, fake_grids) -> Grid:
        channels = [Channel([([2], 1.0)]), Channel([([3], 0.5)])]
        orders = [Order(0, 0, 0, 0, 0)]
        return fake_grids.grid_with_generic_convolution(
            nb_convolutions=1,
            orders=orders,
            channels=channels,
            convolutions=[CONVOBJECT],
        )

    def fake_importonlysubgrid(
        self, nb_xdim: int = 1
    ) -> Tuple[ImportSubgridV1, OperatorInfo]:
        x_grids = [np.linspace(0.1, 1, 2) for _ in range(nb_xdim)]
        xgrid_size = [x.size for x in x_grids]
        Q2s = np.linspace(10, 20, 2)
        scale = [q2 for q2 in Q2s]  # One single scale Q2
        array = np.random.rand(len(Q2s), *xgrid_size)
        infos = OperatorInfo(x_grids, scale, array)
        subgrid = ImportSubgridV1(array=array, node_values=[scale, *x_grids])
        return subgrid, infos

    def test_subgrid_methods(self, fake_grids):
        grid = self.fake_grid(fake_grids)
        test_subgrid, infos = self.fake_importonlysubgrid()
        grid.set_subgrid(0, 0, 0, test_subgrid.into())
        extr_subgrid = grid.subgrid(0, 0, 0)
        assert isinstance(extr_subgrid, SubgridEnum)

        # Check that the subgrid can be scaled
        extr_subgrid.scale(factor=100)
        assert isinstance(extr_subgrid.into(), SubgridEnum)

    @pytest.mark.parametrize("nb_xdim", [1, 2, 3, 4])
    def test_subgrid_arrays(self, nb_xdim: int):
        """This simply checks that the commands run without raising any
        errors and that the objects have been succesfully instantiated.
        """
        subgrid, info = self.fake_importonlysubgrid(nb_xdim=nb_xdim)
        assert isinstance(subgrid, ImportSubgridV1)

    def test_to_array(self, fake_grids):
        grid = self.fake_grid(fake_grids)
        test_subgrid, infos = self.fake_importonlysubgrid()
        grid.set_subgrid(0, 0, 0, test_subgrid.into())
        extr_subgrid = grid.subgrid(0, 0, 0)

        # Check that the shape of the subgrid matches specs
        extr_subgrid_shape = extr_subgrid.shape
        assert tuple(extr_subgrid_shape) == infos.array.shape

        # Check that the `node_values` correspond with the Kinematics
        node_values = extr_subgrid.node_values
        np.testing.assert_allclose(node_values, [infos.scale, *infos.x_grids])

        test_array = extr_subgrid.to_array(shape=extr_subgrid_shape)
        np.testing.assert_allclose(test_array, infos.array)
