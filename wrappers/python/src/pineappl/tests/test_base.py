"""
Test imports and basic functionality that is indepedent of calculation backend.
"""
import pathlib
import pineappl
import numpy as np


REGRESSION_FOLDER = pathlib.Path(__file__).with_name('regressions')


def assert_regression_fixture(array, filename):
    """Check array matches data inside filename.

    Args:
        array: numpy array/
        filename: fixture filename

    If filename does not exists, this function
    creates the missing file otherwise it loads
    from file and compare.
    """
    def load(filename):
        return np.loadtxt(filename)
    try:
        array_fixture = load(filename)
    except: # pragma: no cover
        np.savetxt(filename, array)
        array_fixture = load(filename)
    np.testing.assert_allclose(array, array_fixture, rtol=1e-12)


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
    jacobian *= 2.0 * np.math.acos(-1.0)

    return [s, t, u, x1, x2, jacobian]


def fill_grid(grid, calls):

    # in GeV^2 pbarn
    hbarc2 = 389379372.1

    for _ in range(calls):
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

        grid.fill(x1, x2, q2, 0, np.abs(yll), 0, weight)


def simulate():
    # create a new luminosity function for the $\gamma\gamma$ initial state
    np.random.seed(0)

    lumi = pineappl.lumi()
    pdg_ids = [22, 22]
    ckm_factors = [1.0]
    lumi.add(pdg_ids, ckm_factors)

    # only LO
    orders = [0, 2, 0, 0]
    bins = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
            1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4]
    grid = pineappl.grid(lumi, orders, bins)

    # fill the grid with phase-space points
    fill_grid(grid, 100000)

    return grid


def convolute(grid):
    def xfx(id, x, q2, p):
        return 1
    def alphas(q2, p):
        return 1
    return grid.convolute(xfx, xfx, alphas, None, None, 1.0, 1.0)


def test_sanity():
    """Check that pineappl is consistent."""
    grid = simulate()

    # perform convolution
    dxsec = convolute(grid)
    assert_regression_fixture(dxsec, REGRESSION_FOLDER/'dxsec.out')


def test_read_and_write(filename='DY-LO-AA.pineappl'):
    """Check that read/write works"""
    original_grid = simulate()

    # write the grid to disk
    original_grid.write(filename)
    original_dxsec = convolute(original_grid)

    # read and check consistency
    read_grid = pineappl.grid.read(filename)
    read_dxsec = convolute(read_grid)
    np.testing.assert_allclose(original_dxsec, read_dxsec)
