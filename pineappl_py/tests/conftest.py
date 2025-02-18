import numpy as np
import pytest
import subprocess
from typing import List

from pineappl.boc import (
    BinsWithFillLimits,
    Channel,
    Kinematics,
    ScaleFuncForm,
    Scales,
    Order,
)
from pineappl.convolutions import Conv
from pineappl.grid import Grid
from pineappl.interpolation import (
    Interp,
    InterpolationMethod,
    MappingMethod,
    ReweightingMethod,
)
from pineappl.pids import PidBasis


class CustomGrid(Grid):
    def _get_bins_array(self) -> np.ndarray:
        """Get the bin values as a combined array.

        Returns
        -------
        np.ndarray:
            an array of shape (n_bins, bin_dim, 2)
        """
        return np.array(self.bin_limits())

    def bin_left(self, bin_dim_index: int) -> np.ndarray:
        """Get the left-hand bin for a given dimension of bin.

        Parameters
        ----------
        bin_index: int
            the index-th dimension of the bins

        Returns
        -------
        np.ndarray:
            the array of bins with shape (n_bins,)
        """
        return self._get_bins_array()[:, bin_dim_index, 0]

    def bin_right(self, bin_dim_index: int) -> np.ndarray:
        """Get the right-hand bin for a given dimension of bin.

        Parameters
        ----------
        bin_index: int
            the index-th dimension of the bins

        Returns
        -------
        np.ndarray:
            the array of bins with shape (n_bins,)
        """
        return self._get_bins_array()[:, bin_dim_index, 1]


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

    # Define the Toy Unpolarized PDF set
    def unpolarized_pdf(self, pid, x, q2):
        return 1.0

    # Define the Toy Polarized PDF set
    def polarized_pdf(self, pid, x, q2):
        return 2.0

    # Define the Toy Fragmentation set
    def ff_set(self, pid, x, q2):
        return 3.0


class FakeGrid:
    """Class that mocks a PineAPPL grid. This should contain functions
    that return all the possible number of convolutions.

    TODO: Expose the index that defines the `ScaleFuncForm`.
    """

    def grid_with_generic_convolution(
        self,
        nb_convolutions: int,
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
    ) -> CustomGrid:
        """A function to generate fake GRIDs that can take any number of convolutions.
        Note that the `nb_convolutions` can be different from the number of convolution
        types passed to `convolutions`. Indeed, if all the types of convolutions are
        the same, then only one single element can be passed to `convolutions`.
        """
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

        # Extend the Kinematics and Interpolations
        if nb_convolutions > 1:
            for i in range(1, nb_convolutions):
                kinematics.append(Kinematics.X(i))
                interpolations.append(
                    Interp(
                        min=xmin,
                        max=xmax,
                        nodes=xnodes,
                        order=3,
                        reweight_meth=ReweightingMethod.ApplGridX,
                        map=MappingMethod.ApplGridF2,
                        interpolation_meth=InterpolationMethod.Lagrange,
                    )
                )

        # Construct the `Scales` object
        fragmentation_scale = (
            ScaleFuncForm.Scale(0) if nb_convolutions >= 3 else ScaleFuncForm.NoScale(0)
        )
        scale_funcs = Scales(
            ren=ScaleFuncForm.Scale(0),
            fac=ScaleFuncForm.Scale(0),
            frg=fragmentation_scale,
        )

        #  Construct the bin object
        bin_limits = BinsWithFillLimits.from_fill_limits(fill_limits=bins)

        return CustomGrid(
            pid_basis=PidBasis.Evol,
            channels=channels,
            orders=orders,
            bins=bin_limits,
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
