#!/usr/bin/env python3

import argparse
import eko
import pathlib
import pineappl

from eko.io import manipulate


def main():
    parser = argparse.ArgumentParser(
        description="Re-interpolate EKO to match PineAPPL Grid x-nodes."
    )
    parser.add_argument("grid_path", help="Path to the PineAPPL Grid file")
    parser.add_argument("eko_path", help="Path to the EKO tar file")
    parser.add_argument(
        "--max-as", type=int, default=3, help="Maximum as order (default: 3)"
    )
    parser.add_argument(
        "--max-al", type=int, default=0, help="Maximum al order (default: 0)"
    )

    args = parser.parse_args()

    grid = pineappl.grid.Grid.read(args.grid_path)
    mask = pineappl.boc.Order.create_mask(grid.orders(), args.max_as, args.max_al, True)
    evinfo = grid.evolve_info(mask)
    x_grid = evinfo.x1

    with eko.EKO.edit(pathlib.Path(args.eko_path)) as operator:
        for (q2, _), op in operator.items():
            op_slice = manipulate.xgrid_reshape(
                op,
                eko.interpolation.XGrid(x_grid),
                operator.operator_card.configs.interpolation_polynomial_degree,
                targetgrid=eko.interpolation.XGrid(x_grid),
            )
            operator[(q2, _)] = op_slice
            operator.xgrid = eko.interpolation.XGrid(x_grid)


if __name__ == "__main__":
    main()
