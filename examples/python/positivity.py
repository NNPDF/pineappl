#!/usr/bin/env python

import numpy as np
import pineappl


def main(filename, Q2):
    # setup data
    pid = 4
    xgrid = np.geomspace(5e-5, 0.7, 10)
    bins_length = len(xgrid)
    bin_limits = [float(i) for i in range(0, bins_length + 1)]

    # Instantiate the objecs required to construct a new Grid
    channels = [pineappl.boc.Channel([([pid], 1.0)])]
    orders = [pineappl.boc.Order(0, 0, 0, 0, 0)]
    convolution_types = pineappl.convolutions.ConvType(polarized=False, time_like=False)
    convolutions = [
        pineappl.convolutions.Conv(convolution_types=convolution_types, pid=2212)
    ]
    kinematics = [pineappl.boc.Kinematics.Scale(0), pineappl.boc.Kinematics.X(0)]
    scale_funcs = pineappl.boc.Scales(
        ren=pineappl.boc.ScaleFuncForm.Scale(0),
        fac=pineappl.boc.ScaleFuncForm.Scale(0),
        frg=pineappl.boc.ScaleFuncForm.NoScale(0),
    )
    bin_limits = pineappl.boc.BinsWithFillLimits.from_fill_limits(
        fill_limits=bin_limits
    )
    interpolations = [
        pineappl.interpolation.Interp(
            min=1e2,
            max=1e3,
            nodes=50,
            order=3,
            reweight_meth=pineappl.interpolation.ReweightingMethod.NoReweight,
            map=pineappl.interpolation.MappingMethod.ApplGridH0,
            interpolation_meth=pineappl.interpolation.InterpolationMethod.Lagrange,
        ),  # Interpolation on the Scale
        pineappl.interpolation.Interp(
            min=1e-5,
            max=1,
            nodes=40,
            order=3,
            reweight_meth=pineappl.interpolation.ReweightingMethod.ApplGridX,
            map=pineappl.interpolation.MappingMethod.ApplGridF2,
            interpolation_meth=pineappl.interpolation.InterpolationMethod.Lagrange,
        ),  # Interpolation on momentum fraction x
    ]

    grid = pineappl.grid.Grid(
        pid_basis=pineappl.pids.PidBasis.Evol,
        channels=channels,
        orders=orders,
        bins=bin_limits,
        convolutions=convolutions,
        interpolations=interpolations,
        kinematics=kinematics,
        scale_funcs=scale_funcs,
    )

    limits = []
    # add each point as a bin
    for bin_, x in enumerate(xgrid):
        # keep DIS bins
        limits.append([(Q2, Q2), (x, x)])
        # Fill the subgrid with delta functions
        array_subgrid = np.zeros((1, xgrid.size))
        array_subgrid[0][bin_] = 1
        # create and set the subgrid
        subgrid = pineappl.subgrid.ImportSubgridV1(
            array=array_subgrid,
            node_values=[[Q2], xgrid],
        )
        grid.set_subgrid(0, bin_, 0, subgrid.into())
    # set the correct observables
    normalizations = [1.0] * bins_length
    bin_configs = pineappl.boc.BinsWithFillLimits.from_limits_and_normalizations(
        limits=limits,
        normalizations=normalizations,
    )
    grid.set_bwfl(bin_configs)

    # set the initial state PDF ids for the grid
    grid.set_metadata(
        "runcard",
        f"positivity constraint for quark {pid}",
    )

    # dump file
    grid.optimize()
    grid.write(filename)


if __name__ == "__main__":
    main("charm-positivity.pineappl", 1.5**2)
