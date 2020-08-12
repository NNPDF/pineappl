Basic examples
==============


How to create and fill a PineAPPL grid?
---------------------------------------

Here is an example of Drell-Yan photon induced contributions:

.. code-block::  python

    # create a new luminosity function for the $\gamma\gamma$ initial state
    lumi = pineappl.lumi()
    pdg_ids = [22, 22]
    ckm_factors = [1.0]
    lumi.add(pdg_ids, ckm_factors)

    # only LO $\alpha_\mathrm{s}^0 \alpha^2 \log^0(\xi_\mathrm{R}) \log^0(\xi_\mathrm{F})$
    orders = [0, 2, 0, 0]
    bins = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
            1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4]
    grid = pineappl.grid(lumi, orders, bins)

    # fill the grid with phase-space points
    print('Generating events, please wait...')
    fill_grid(grid, 100000)

    # write the grid to disk
    filename = 'DY-LO-AA.pineappl'
    print(f'Writing PineAPPL grid to disk: {filename}')
    grid.write(filename)

