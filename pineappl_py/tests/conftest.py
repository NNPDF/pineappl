import numpy as np
import pytest
import subprocess
from typing import List

from pineappl.boc import Channel, Kinematics, ScaleFuncForm, Scales
from pineappl.convolutions import Conv
from pineappl.grid import Grid, Order
from pineappl.interpolation import (
    Interp,
    InterpolationMethod,
    MappingMethod,
    ReweightingMethod,
)
from pineappl.pids import PidBasis


class PDF:
    """PDF class whose attributes are some toy PDF functions."""

    def xfxQ2(self, pid, x, q2):
        if pid in range(-6, 6):
            return x * (1 - x)
        else:
            return 0.0

    def xfxQ(self, pid, x, q):
        return self.xfxQ2(pid, x, q**2)

    def alphasQ(self, q):
        return 1.0

    # Define the Toy Polarized PDF set
    def polarized_pdf(self, pid, x, q2):
        return 2.0

    # Define the Toy Unpolarized PDF set
    def unpolarized_pdf(self, pid, x, q2):
        return 1.0


class FakeGrid:
    """Class that mocks a PineAPPL grid. This should contain functions
    that return all the possible number of convolutions.

    TODO: combine the different convolutions
    """

    def grid_with_one_convolution(
        self,
        channels: List[Channel],
        orders: List[Order],
        convolutions: List[Conv],
        bins: List[float] = [1e-7, 1e-3, 1],
        q2min: float = 1e2,
        q2max: float = 1e8,
        q2nodes: int = 40,
        xmin: float = 2e-7,
        xmax: float = 1,
        xnodes: int = 50,
    ) -> Grid:
        kinematics = [
            Kinematics.Scale(0),  # Scale
            Kinematics.X(0),  # momentum fraction x
        ]
        # Define the interpolation specs for each item of the Kinematics
        interpolations = [
            Interp(
                min=q2min,
                max=q2max,
                nodes=q2nodes,
                order=3,
                reweight_meth=ReweightingMethod.NoReweight,
                map=MappingMethod.ApplGridH0,
                interpolation_meth=InterpolationMethod.Lagrange,
            ),  # Interpolation on the Scale
            Interp(
                min=xmin,
                max=xmax,
                nodes=xnodes,
                order=3,
                reweight_meth=ReweightingMethod.ApplGridX,
                map=MappingMethod.ApplGridF2,
                interpolation_meth=InterpolationMethod.Lagrange,
            ),  # Interpolation on momentum fraction x
        ]
        # Construct the `Scales` object
        scale_funcs = Scales(
            ren=ScaleFuncForm.Scale(0),
            fac=ScaleFuncForm.Scale(0),
            frg=ScaleFuncForm.NoScale(0),
        )

        return Grid(
            pid_basis=PidBasis.Evol,
            channels=channels,
            orders=orders,
            bin_limits=np.array(bins),
            convolutions=convolutions,
            interpolations=interpolations,
            kinematics=kinematics,
            scale_funcs=scale_funcs,
        )

    def grid_with_two_convolutions(
        self,
        channels: List[Channel],
        orders: List[Order],
        convolutions: List[Conv],
        bins: List[float] = [1e-7, 1e-3, 1],
        q2min: float = 1e2,
        q2max: float = 1e8,
        q2nodes: int = 40,
        xmin: float = 2e-7,
        xmax: float = 1,
        xnodes: int = 50,
    ) -> Grid:
        kinematics = [
            Kinematics.Scale(0),  # Scale
            Kinematics.X(0),  # x1 momentum fraction
            Kinematics.X(1),  # x2 momentum fraction
        ]
        # Define the interpolation specs for each item of the Kinematics
        interpolations = [
            Interp(
                min=q2min,
                max=q2max,
                nodes=q2nodes,
                order=3,
                reweight_meth=ReweightingMethod.NoReweight,
                map=MappingMethod.ApplGridH0,
                interpolation_meth=InterpolationMethod.Lagrange,
            ),  # Interpolation on the Scale
            Interp(
                min=xmin,
                max=xmax,
                nodes=xnodes,
                order=3,
                reweight_meth=ReweightingMethod.ApplGridX,
                map=MappingMethod.ApplGridF2,
                interpolation_meth=InterpolationMethod.Lagrange,
            ),  # Interpolation on x1 momentum fraction
            Interp(
                min=xmin,
                max=xmax,
                nodes=xnodes,
                order=3,
                reweight_meth=ReweightingMethod.ApplGridX,
                map=MappingMethod.ApplGridF2,
                interpolation_meth=InterpolationMethod.Lagrange,
            ),  # Interpolation on x2 momentum fraction
        ]
        # Construct the `Scales` object
        scale_funcs = Scales(
            ren=ScaleFuncForm.Scale(0),
            fac=ScaleFuncForm.Scale(0),
            frg=ScaleFuncForm.NoScale(0),
        )

        return Grid(
            pid_basis=PidBasis.Evol,
            channels=channels,
            orders=orders,
            bin_limits=np.array(bins),
            convolutions=convolutions,
            interpolations=interpolations,
            kinematics=kinematics,
            scale_funcs=scale_funcs,
        )


@pytest.fixture
def pdf():
    return PDF()


@pytest.fixture
def fake_grids():
    return FakeGrid()


@pytest.fixture
def download_objects(tmp_path_factory):
    def _download_fk(objname: str) -> None:
        download_dir = tmp_path_factory.mktemp("data")
        file_path = download_dir / f"{objname}"
        args = [
            "wget",
            "--no-verbose",
            "--no-clobber",
            "-P",
            f"{download_dir}",
            f"https://data.nnpdf.science/pineappl/test-data/{objname}",
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
