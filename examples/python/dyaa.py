#!/usr/bin/env python
import argparse

import numpy as np
import pineappl


def int_photo(s, t, u):
    """
    Photon-photon matrix element.

    Parameters
    ----------
        s : float
            Mandelstam s
        t : float
            Mandelstam t
        u : float
            Mandelstam u

    Returns
    -------
        float :
            matrix element
    """
    alpha0 = 1.0 / 137.03599911
    return alpha0 * alpha0 / 2.0 / s * (t / u + u / t)


def hadronic_pspgen(mmin, mmax):
    r"""
    Hadronic phase space generator.

    Parameters
    ----------
        mmin : float
            minimal energy :math:`\sqrt{s_{min}}`
        mmax : float
            maximal energy :math:`\sqrt{s_{max}}`

    Returns
    -------
        s : float
            Mandelstam s
        t : float
            Mandelstam t
        u : float
            Mandelstam u
        x1 : float
            first momentum fraction
        x2 : float
            second momentum fraction
        jacobian : float
            jacobian from the uniform generation
    """
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
    jacobian *= 2.0 * np.math.acos(-1.0)

    return [s, t, u, x1, x2, jacobian]


def fill_grid(grid, calls):
    """
    Fill grid with points.

    Parameters
    ----------
        grid : pineappl.Grid
            grid to fill
        calls : int
            number of events
    """

    # in GeV^2 pbarn
    hbarc2 = 389379372.1

    for _ in range(calls):
        # compute phase space
        s, t, u, x1, x2, jacobian = hadronic_pspgen(10.0, 7000.0)
        # build observables
        ptl = np.sqrt((t * u / s))
        mll = np.sqrt(s)
        yll = 0.5 * np.log(x1 / x2)
        ylp = np.abs(yll + np.math.acosh(0.5 * mll / ptl))
        ylm = np.abs(yll - np.math.acosh(0.5 * mll / ptl))

        jacobian *= hbarc2 / calls

        # cuts for LO for the invariant-mass slice containing the
        # Z-peak from CMSDY2D11
        if (
            ptl < 14.0
            or np.abs(yll) > 2.4
            or ylp > 2.4
            or ylm > 2.4
            or mll < 60.0
            or mll > 120.0
        ):
            continue

        # build event
        weight = jacobian * int_photo(s, u, t)
        q2 = 90.0 * 90.0
        # put
        grid.fill(x1, x2, q2, 0, np.abs(yll), 0, weight)


def main(calls, pdfname, filename):
    """
    Generate the grid.

    Parameters
    ----------
        calls : int
            number of events
        pdfname : str
            if given, write the predictions
        filename : str
            if given, write the grid to this path
    """
    # create a new luminosity function for the $\gamma\gamma$ initial state
    lumi_entries = [pineappl.lumi.LumiEntry([(22, 22, 1.0)])]
    # only LO $\alpha_\mathrm{s}^0 \alpha^2 \log^0(\xi_\mathrm{R}) \log^0(\xi_\mathrm{F})$
    orders = [pineappl.grid.Order(0, 2, 0, 0)]
    bins = np.arange(0, 2.4, 0.1)
    params = pineappl.subgrid.SubgridParams()
    grid = pineappl.grid.Grid.create(lumi_entries, orders, bins, params)

    # fill the grid with phase-space points
    print(f"Generating {calls} events, please wait...")
    fill_grid(grid, calls)

    # load pdf for testing
    if pdfname:
        import lhapdf  # pylint: disable=import-outside-toplevel

        pdf = lhapdf.mkPDF(pdfname, 0)
        pdg_id = int(pdf.set().get_entry("Particle"))
        # perform convolution
        dxsec = grid.convolve_with_one(pdg_id, pdf.xfxQ2, pdf.alphasQ2)
        for i in range(len(dxsec)):
            print(f"{bins[i]:.1f} {bins[i + 1]:.1f} {dxsec[i]:.3e}")

    # write the grid to disk
    if filename:
        print(f"Writing PineAPPL grid to disk: {filename}")
        grid.write(filename)


if __name__ == "__main__":
    # expose the arguments as script
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--calls", default=100000, type=int, help="Number of events"
    )
    ap.add_argument(
        "--pdf",
        default="NNPDF31_nlo_as_0118_luxqed",
        type=str,
        help="plot with PDF set",
    )
    ap.add_argument("--filename", type=str, help="Output path")
    args = ap.parse_args()
    main(args.calls, args.pdf, args.filename)
