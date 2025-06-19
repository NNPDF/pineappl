"""Example used in the PineAPPL v1 paper."""

import math
import numpy as np

from dataclasses import dataclass
from numpy.random import Generator, PCG64
from pineappl.boc import (
    BinsWithFillLimits,
    Channel,
    Kinematics,
    ScaleFuncForm,
    Scales,
    Order,
)
from pineappl.convolutions import Conv, ConvType
from pineappl.grid import Grid
from pineappl.interpolation import (
    Interp,
    InterpolationMethod,
    MappingMethod,
    ReweightingMethod,
)
from pineappl.pids import PidBasis


RNG = Generator(PCG64())


@dataclass
class Psp2to2Hadron:
    s: np.ndarray  # Mandelstam variable s
    t: np.ndarray  # Mandelstam variable s
    u: np.ndarray  # Mandelstam variable s
    x1: np.ndarray  # momemtum fraction of the 1st initial-state hadron
    x2: np.ndarray  # momemtum fraction of the 2nd initial-state hadron
    z: np.ndarray  # momentum fraction of the final-state hadron
    pth: np.ndarray  # transverse momentum of the final-state hadron
    yh: np.ndarray  # rapidity of the final-state hadron
    jc: np.ndarray  # jacobian factor


def create_grid(pt_bins: np.ndarray) -> Grid:
    """Instantiate the Grid object for the single-inclusive pion production in pp."""
    # define the sub-partonic channels for gg -> qqbar
    sub_channels = [([21, 21, pid], 1.0) for pid in range(-3, 4) if pid != 0]
    channels = [Channel(sub_channels)]  # construct the channel object

    # define the perturbative order that will be filled into the grid
    # orders specifies the power of the tuple `orders = (αs, α, lR, lF, lD)`
    # in this example, we only fill the LO QCD
    orders = [Order(1, 0, 0, 0, 0)]

    # define the convolution object: 2 initial-state protons + 1 final-state pion
    proton_conv = Conv(
        convolution_types=ConvType(polarized=False, time_like=False), pid=2212
    )
    pion_conv = Conv(
        convolution_types=ConvType(polarized=False, time_like=True), pid=211
    )

    # define the kinematics object: (mu, x1, x2, z)
    kinematics = [
        Kinematics.Scale(0),
        Kinematics.X(0),
        Kinematics.X(1),
        Kinematics.X(2),
    ]

    # define the specifities of the interpolations corresponding to the kinematics
    scale_interp = Interp(
        min=1e2,
        max=1e8,
        nodes=40,
        order=3,
        reweight_meth=ReweightingMethod.NoReweight,
        map=MappingMethod.ApplGridH0,
        interpolation_meth=InterpolationMethod.Lagrange,
    )  # interpolation on the scale
    x_interp = Interp(
        min=2e-7,
        max=1.0,
        nodes=50,
        order=3,
        reweight_meth=ReweightingMethod.ApplGridX,
        map=MappingMethod.ApplGridF2,
        interpolation_meth=InterpolationMethod.Lagrange,
    )  # interpolation on the momentum fraction x
    interpolations = [scale_interp, x_interp, x_interp, x_interp]

    # Construc the scale object for the (muR, muF, muD) scales
    scale_funcs = Scales(
        ren=ScaleFuncForm.Scale(0),
        fac=ScaleFuncForm.Scale(0),
        frg=ScaleFuncForm.Scale(0),
    )

    # construct the bin object
    bin_limits = BinsWithFillLimits.from_fill_limits(fill_limits=pt_bins)

    return Grid(
        pid_basis=PidBasis.Evol,  # use the Evolution basis to represent the grid
        channels=channels,
        orders=orders,
        bins=bin_limits,
        convolutions=[proton_conv, proton_conv, pion_conv],
        interpolations=interpolations,
        kinematics=kinematics,
        scale_funcs=scale_funcs,
    )


def me_gg2qqbar(_s: np.ndarray, t: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Calculate the matrix element squared for the gg -> qqbar channel."""
    as2 = 0.118 * 0.118
    PI2 = math.pi * math.pi
    # TODO: double-check
    return (16 * PI2 * as2 / 6.0) * (np.power(u, 2) + np.power(t, 2)) / (u * t)


def psgen_pp2hadron(
    n_calls: int,
    mmin: float,
    mmax: float,
    pt_min: float,
    pt_max: float,
    abs_y_max: float,
) -> Psp2to2Hadron:
    """Generate the kinematics for the entire phase-space."""

    smin = mmin * mmin
    smax = mmax * mmax

    r1: np.ndarray = RNG.uniform(0, 1, size=n_calls)
    r2: np.ndarray = RNG.uniform(0, 1, size=n_calls)
    r3: np.ndarray = RNG.uniform(0, 1, size=n_calls)
    r4: np.ndarray = RNG.uniform(0, 1, size=n_calls)
    r5: np.ndarray = RNG.uniform(0, 1, size=n_calls)
    r6: np.ndarray = RNG.uniform(0, 1, size=n_calls)

    tau0 = smin / smax
    tau = np.power(tau0, r1)
    y = np.power(tau, 1.0 - r2)
    x1 = y
    x2 = tau / y
    s = tau * smax

    jacobian = tau * np.log(tau0) * np.log(tau0) * r1

    # `theta` integration
    cos_theta = 2.0 * r3 - 1.0
    jacobian *= 2.0

    t = -0.5 * s * (1.0 - cos_theta)
    u = -0.5 * s * (1.0 + cos_theta)

    # `phi` integration
    jacobian *= 2.0 * math.acos(-1.0)

    # sample hadron `pT` uniformly in log scale
    log_pt_min = math.log(pt_min)
    log_pt_max = math.log(pt_max)

    pt_hadron = np.exp(log_pt_min + (log_pt_max - log_pt_min) * r4)
    jacobian *= pt_hadron * (log_pt_max - log_pt_min)

    # sample hadron rapidity uniformly
    y_hadron = 2.0 * abs_y_max * r5 - abs_y_max
    jacobian *= 2.0 * abs_y_max

    # define the momentum fracion `z`
    z_min = pt_hadron * np.exp(-y_hadron) / np.sqrt(s)
    z_max_kin = pt_hadron * np.exp(y_hadron) / np.sqrt(s)
    z_max = np.minimum(z_max_kin, 1.0)

    # ensure that `z` is physical - set unphysical values to zero
    check_kin1: np.ndarray = z_min >= 1
    check_kin2: np.ndarray = z_min >= z_max
    z_min[check_kin1] = 0.0
    z_min[check_kin2] = 0.0
    jacobian[check_kin1] = 0.0
    jacobian[check_kin2] = 0.0

    # sample `z` uniformly between the kinematic limits
    z = z_min + (z_max - z_min) * r6
    jacobian *= z_max - z_min

    return Psp2to2Hadron(
        s=s, t=t, u=u, x1=x1, x2=x2, z=z, pth=pt_hadron, yh=y_hadron, jc=jacobian
    )


def fill_grid(grid: Grid, n_calls: int, pt_bins: list) -> Grid:
    """Fill the Grid with the phase-space points."""
    hbarc2 = 389379372.1
    abs_y_max = 2.4  # rapidity cut

    # compute the phase-space kinematics
    ps = psgen_pp2hadron(n_calls, 3e3, 14e3, pt_bins[0], pt_bins[-1], abs_y_max)

    # calculate the partonic cross-section
    jacobian = ps.jc * (hbarc2 / n_calls)
    weight = jacobian * me_gg2qqbar(ps.s, ps.t, ps.u)

    # pass the kinematics as a list of n-tuples
    scale = ps.pth * ps.pth
    kinematics_ntuples = [
        np.array([mu, x1, x2, z]) for mu, x1, x2, z in zip(scale, ps.x1, ps.x2, ps.z)
    ]

    # fill the grid by passing **at once** the weights and observables for a given order and channel
    for pto_idx in range(len(grid.orders())):
        for channel_idx in range(len(grid.channels())):
            grid.fill_array(
                order=pto_idx,
                observables=ps.pth,
                channel=channel_idx,
                ntuples=kinematics_ntuples,
                weights=weight,
            )

    return grid


def main() -> None:
    pt_bins = [
        5.0,
        7.0,
        10.0,
        15.0,
        20.0,
        25.0,
        30.0,
        35.0,
        40.0,
        45.0,
        50.0,
        60.0,
        70.0,
        80.0,
        90.0,
        100.0,
    ]
    n_calls = 100000
    grid_instance = create_grid(pt_bins=np.array(pt_bins))

    # fill the grid
    grid = fill_grid(grid=grid_instance, n_calls=n_calls, pt_bins=pt_bins)

    # add some metadata and write the Grid onto disk
    grid.set_metadata("x1_label", "pT")
    grid.set_metadata("y_label", "dsig/dpT")
    grid.set_metadata("x1_unit", "GeV")
    grid.set_metadata("y_unit", "pb/GeV")

    # write the Grid onto disk
    grid.write_lz4("pp2hadron-pt-py.pineappl.lz4")

    return


if __name__ == "__main__":
    main()
