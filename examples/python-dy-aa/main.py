#!/usr/bin/env python
import pineappl
import numpy as np


def int_photo(s, t, u):
    alpha0 = 1.0 / 137.03599911
    return alpha0 * alpha0 / 2.0 / s * (t / u + u / t)


def hadronic_pspgen(mmin, mmax):
    smin = mmin * mmin
    smax = mmax * mmax

    r1 = np.random.uniform()
    r2 = np.random.uniform()
    r3 = np.random.uniform()

    tau0 = smin / smax
    tau = pow(tau0, r1)
    y = pow(tau, 1.0 - r2)
    x1 = y
    x2 = tau / y
    s = tau * smax

    jacobian = tau * np.log(tau0) * np.log(tau0) * r1

    # theta integration (in the CMS)
    cos_theta = 2.0 * r3 - 1.0
    jacobian *= 2.0

    t = -0.5 * s * (1.0 - cos_theta)
    u = -0.5 * s * (1.0 + cos_theta)

    # phi integration
    jacobian *= 2.0 * np.math.acos(-1.0);

    return [s, t, u, x1, x2, jacobian]


def fill_grid(grid, calls):

    # in GeV^2 pbarn
    hbarc2 = 389379372.1

    for i in range(calls):
        s, t, u, x1, x2, jacobian = hadronic_pspgen(10.0, 7000.0)

        ptl = np.sqrt((t * u / s))
        mll = np.sqrt(s)
        yll = 0.5 * np.log(x1 / x2)
        ylp = np.abs(yll + np.math.acosh(0.5 * mll / ptl))
        ylm = np.abs(yll - np.math.acosh(0.5 * mll / ptl))

        jacobian *= hbarc2 / calls;

        # cuts for LO for the invariant-mass slice containing the
        # Z-peak from CMSDY2D11
        if ptl < 14.0 or np.abs(yll) > 2.4 or ylp > 2.4 \
            or ylm > 2.4 or mll < 60.0 or mll > 120.0:
            continue

        weight = jacobian * int_photo(s, u, t)
        q2 = 90.0 * 90.0

        pineappl.grid_fill(grid, x1, x2, q2, 0, np.abs(yll), 0, weight)


def main():
    # create a new luminosity function for the $\gamma\gamma$ initial state
    lumi = pineappl.lumi_new()
    pdg_ids = [22, 22]
    ckm_factors = [1.0]
    pineappl.lumi_add(lumi, 1, pdg_ids, ckm_factors)

    # only LO $\alpha_\mathrm{s}^0 \alpha^2 \log^0(\xi_\mathrm{R}) \log^0(\xi_\mathrm{F})$
    orders = [0, 2, 0, 0]
    bins = [ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
             1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4]
    keyval = pineappl.keyval_new()
    grid = pineappl.grid_new(lumi, 1, orders, 24, bins, keyval)

    # now we no longer need keyval and lumi
    pineappl.keyval_delete(keyval)
    pineappl.lumi_delete(lumi)

    # fill the grid with phase-space points
    fill_grid(grid, 1000)

    # write the grid to disk
    pineappl.grid_write(grid, "DY-LO-AA.pineappl")
    pineappl.grid_delete(grid)


if __name__ == '__main__':
    main()
