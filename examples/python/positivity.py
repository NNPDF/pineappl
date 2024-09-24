#!/usr/bin/env python

import numpy as np
import pineappl


def main(filename, Q2):
    # setup data
    xgrid = np.geomspace(5e-5, 0.7, 10)
    lepton_pid = 11
    pid = 4

    # init pineappl objects
    lumi_entries = [pineappl.boc.Channel([(pid, lepton_pid, 1.0)])]
    orders = [pineappl.grid.Order(0, 0, 0, 0)]
    bins = len(xgrid)
    # NOTE: `bin_limits` have to be `np.ndarray`
    bin_limits = np.array([float(i) for i in range(0, bins + 1)])
    # subgrid params - default is just sufficient
    params = pineappl.subgrid.SubgridParams()
    # inti grid
    grid = pineappl.grid.Grid(lumi_entries, orders, bin_limits, params)
    limits = []
    # add each point as a bin
    for bin_, x in enumerate(xgrid):
        # keep DIS bins
        limits.append((Q2, Q2))
        limits.append((x, x))
        # delta function
        array = np.zeros(len(xgrid))
        array[bin_] = 1
        # create and set
        subgrid = pineappl.import_only_subgrid.ImportOnlySubgridV1(
            array[np.newaxis, :, np.newaxis],
            np.array([Q2]),  # `q2_grid` has to be `np.ndarrary`
            np.array(xgrid),  # `x_grid` has to be `np.ndarrary`
            np.array([1.0]),  # `x_grid` has to be `np.ndarrary`
        )
        grid.set_subgrid(0, bin_, 0, subgrid.into())
    # set the correct observables
    normalizations = np.array(
        [1.0] * bins
    )  # `normalizations` has to be `np.ndarray`
    remapper = pineappl.bin.BinRemapper(normalizations, limits)
    grid.set_remapper(remapper)

    # set the initial state PDF ids for the grid
    grid.set_key_value("initial_state_1", "2212")
    grid.set_key_value("initial_state_2", str(lepton_pid))
    grid.set_key_value(
        "runcard",
        f"positivity constraint for quark {pid}",
    )

    # dump file
    grid.optimize()
    grid.write(filename)


if __name__ == "__main__":
    main("charm-positivity.pineappl", 1.5**2)
