import argparse

import yaml

import pineappl


def make_pineappl(input_yaml, output_pineappl):
    with open(input_yaml) as f:
        yad_out = yaml.safe_load(f)
    # print(yad_out)

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

    print("\nI'm doing something")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="produce a *pineappl* grid from *yadism* yaml output"
    )

    parser.add_argument("input_yaml", help="input yaml file, from yadism output")
    parser.add_argument("output_pineappl", help="name for output pineappl grid")

    args = parser.parse_args()
    print(args)
    make_pineappl(args.input_yaml, args.output_pineappl)
