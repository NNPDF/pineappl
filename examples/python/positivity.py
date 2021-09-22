#!/usr/bin/env python

import pineappl
import numpy as np


def main(filename):
    # setup data
    xgrid = np.geomspace(5e-5, 0.7, 10)
    Q2 = 1.5 ** 2
    lepton_pid = 11
    pid = 4

    # init pineappl objects
    lumi_entries = [pineappl.lumi.LumiEntry([(pid, lepton_pid, 1.0)])]
    orders = [pineappl.grid.Order(0, 0, 0, 0)]
    bins = len(xgrid)
    bin_limits = list(map(float, range(0, bins + 1)))
    # subgrid params - default is just sufficient
    params = pineappl.subgrid.SubgridParams()
    # inti grid
    grid = pineappl.grid.Grid.create(lumi_entries, orders, bin_limits, params)
    limits = []
    # add each point as a bin
    for bin_, x in enumerate(xgrid):
        x = x

        limits.append((Q2, Q2))
        limits.append((x, x))

        array = np.zeros(len(xgrid))
        array[bin_] = 1

        subgrid = pineappl.import_only_subgrid.ImportOnlySubgridV1(
            array[np.newaxis, :, np.newaxis],
            [Q2],
            xgrid,
            [1.0],
        )
        grid.set_subgrid(0, bin_, 0, subgrid)
    # set the correct observables
    normalizations = [1.0] * bins
    remapper = pineappl.bin.BinRemapper(normalizations, limits)
    grid.set_remapper(remapper)

    # set the initial state PDF ids for the grid
    grid.set_key_value("initial_state_1", "2212")
    grid.set_key_value("initial_state_2", str(lepton_pid))
    grid.set_key_value(
        "runcard",
        f"positivity constraint for quark {pid}",
    )
    grid.set_key_value("lumi_id_types", "pdg_mc_ids")

    # dump file
    grid.optimize()
    grid.write(filename)


if __name__ == "__main__":
    main("charm-positivity.pineappl")
