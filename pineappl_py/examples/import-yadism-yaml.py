import argparse
import pathlib

import yaml
import numpy as np

import pineappl


def make_pineappl(input_yaml, output_pineappl):
    with open(input_yaml) as f:
        yad_out = yaml.safe_load(f)

    # extract output objects
    interpolation_is_log = yad_out["interpolation_is_log"]
    interpolation_polynomial_degree = yad_out["interpolation_polynomial_degree"]
    interpolation_xgrid = yad_out["interpolation_xgrid"]
    pids = yad_out["pids"]

    # hardcoded options
    lepton_pid = 11

    # init pineappl objects
    lumi = [pineappl.lumi.LumiEntry([(pid, lepton_pid, 1.0)]) for pid in pids]
    orders = [pineappl.grid.Order(0, 0, 0, 0)]
    bins = len(yad_out["F2total"])
    bin_limits = list(map(float, range(0, bins + 1)))
    # subgrid params
    params = pineappl.subgrid.SubgridParams()
    params.set_x_bins(len(interpolation_xgrid))
    params.set_x_max(interpolation_xgrid[-1])
    params.set_x_min(interpolation_xgrid[0])
    params.set_x_order(interpolation_polynomial_degree)

    extra = pineappl.subgrid.ExtraSubgridParams()
    extra.set_reweight2(False)
    extra.set_x2_bins(1)
    extra.set_x2_max(1.0)
    extra.set_x2_min(1.0)
    extra.set_x2_order(0)

    grid = pineappl.grid.Grid(
        lumi, orders, bin_limits, pineappl.subgrid.SubgridParams()
    )
    limits = []

    for bin_, obs in enumerate(yad_out["F2total"]):
        q2 = obs["Q2"]
        x = obs["x"]

        limits.append((q2, q2))
        limits.append((x, x))

        params.set_q2_bins(1)
        params.set_q2_max(q2)
        params.set_q2_min(q2)
        params.set_q2_order(0)

        order = 0

        for lumi, values in enumerate(obs["values"]):
            values = list(reversed(values))

            assert len(values) == params.x_bins()

            if any(np.array(values) != 0):
                subgrid = pineappl.lagrange_subgrid.LagrangeSubgridV2(params, extra)
                subgrid.write_q2_slice(0, values)
                grid.set_subgrid(order, bin_, lumi, subgrid)

    # set the correct observables
    normalizations = [1.0] * bins
    remapper = pineappl.bin.BinRemapper(normalizations, limits)
    grid.set_remapper(remapper)

    # set the initial state PDF ids for the grid
    grid.set_key_value("initial_state_1", "2212")
    grid.set_key_value("initial_state_2", str(lepton_pid))

    # TODO: find a way to open file in python
    # with open(output_pineappl, "wb") as f:
    grid.write(output_pineappl)

    print(f"\nI wrote the file '{pathlib.Path(output_pineappl).absolute()}'")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="produce a *pineappl* grid from *yadism* yaml output"
    )

    parser.add_argument("input_yaml", help="input yaml file, from yadism output")
    parser.add_argument("output_pineappl", help="name for output pineappl grid")

    args = parser.parse_args()
    print(args)
    make_pineappl(args.input_yaml, args.output_pineappl)
