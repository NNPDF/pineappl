#!/usr/bin/env python3

import eko
import numpy
import pineappl
import sys

#grid = pineappl.grid.Grid.read('{}.pineappl.lz4'.format(sys.argv[1]))
grid = pineappl.grid.Grid.read('test-data/ATLASPHT15_Et_1bin_last_three_bins.pineappl.lz4')
mask = pineappl.grid.Order.create_mask(grid.orders(), 2, 0, True)
x_grid = grid.evolve_info(mask).x1

#operator = eko.EKO.edit('{}.tar'.format(sys.argv[1]))
operator = eko.EKO.edit('test-data/ATLASPHT15-ATLASPHT15_Et_1bin.tar')

eko.io.manipulate.xgrid_reshape(operator, targetgrid=eko.interpolation.XGrid(x_grid))
operator.close()
