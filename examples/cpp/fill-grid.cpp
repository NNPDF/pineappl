#include <pineappl_capi.h>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <random>

double int_photo(double s, double t, double u) {
    double alpha0 = 1.0 / 137.03599911;
    return alpha0 * alpha0 / 2.0 / s * (t / u + u / t);
}

struct Psp2to2 {
    double s;
    double t;
    double u;
    double x1;
    double x2;
    double jacobian;
};

Psp2to2 hadronic_pspgen(std::mt19937& rng, double mmin, double mmax) {
    using std::acos;
    using std::log;
    using std::pow;

    double smin = mmin * mmin;
    double smax = mmax * mmax;

    double r1 = std::generate_canonical<double, 53>(rng);
    double r2 = std::generate_canonical<double, 53>(rng);
    double r3 = std::generate_canonical<double, 53>(rng);

    double tau0 = smin / smax;
    double tau = pow(tau0, r1);
    double y = pow(tau, 1.0 - r2);
    double x1 = y;
    double x2 = tau / y;
    double s = tau * smax;

    double jacobian = tau * log(tau0) * log(tau0) * r1;

    // theta integration (in the CMS)
    double  cos_theta = 2.0 * r3 - 1.0;
    jacobian *= 2.0;

    double  t = -0.5 * s * (1.0 - cos_theta);
    double  u = -0.5 * s * (1.0 + cos_theta);

    // phi integration
    jacobian *= 2.0 * acos(-1.0);

    return { s, t, u, x1, x2, jacobian };
}

void fill_grid(pineappl_grid* grid, std::size_t calls) {
    using std::acosh;
    using std::fabs;
    using std::log;
    using std::sqrt;

    auto rng = std::mt19937();

    // in GeV^2 pbarn
    double  hbarc2 = 389379372.1;

    for (std::size_t i = 0; i != calls; ++i) {
        // generate a phase-space point
        auto tmp = hadronic_pspgen(rng, 10.0, 7000.0);
        auto s = tmp.s;
        auto t = tmp.t;
        auto u = tmp.u;
        auto x1 = tmp.x1;
        auto x2 = tmp.x2;
        auto jacobian = tmp.jacobian;

        double ptl = sqrt((t * u / s));
        double mll = sqrt(s);
        double yll = 0.5 * log(x1 / x2);
        double ylp = fabs(yll + acosh(0.5 * mll / ptl));
        double ylm = fabs(yll - acosh(0.5 * mll / ptl));

        jacobian *= hbarc2 / calls;

        // cuts for LO for the invariant-mass slice containing the
        // Z-peak from CMSDY2D11
        if ((ptl < 14.0) || (fabs(yll) > 2.4) || (ylp > 2.4)
            || (ylm > 2.4) || (mll < 60.0) || (mll > 120.0))
        {
            continue;
        }

        auto weight = jacobian * int_photo(s, t, u);
        double q2 = 90.0 * 90.0;
        std::size_t order = 0;
        std::size_t channel = 0;

        // fill the LO `weight` into `grid` for parton fractions `x1` and `x2`, and the (squared)
        // renormalization/factorization scale `q2`. The parameters `order` and `channel` are
        // indices defined from the arrays `orders` and `channel` used in creating the grid. In this
        // case they are both `0` and denote the order #0 (leading order) and the channel #0
        // (photon-photon channel), respectively
        pineappl_grid_fill(grid, x1, x2, q2, order, fabs(yll), channel, weight);
    }
}

int main() {
    // ---
    // Create all channels

    // this object will contain all channels (initial states) that we define
    auto* channels = pineappl_lumi_new();

    // photon-photon initial state, where `22` is the photon (PDG MC ids)
    int32_t pids1[] = { 22, 22 };

    // factor that each channel is multiplied with when convoluting with PDFs
    double factors1[] = { 1.0 };

    // define the channel #0
    pineappl_lumi_add(channels, 1, pids1, factors1);

    // create another channel, which we won't fill, however

    // this channel is the down-type-antidown-type quark channel; here we combine down-antidown,
    // strange-antistrange and bottom-antibottom into a single channel, which is often done if the
    // CKM matrix is taken to be diagonal
    int32_t pids2[] = { 1, -1, 3, -3, 5, -5 };

    // for each pair of particle ids we need to give a factor; in case of a non-diagonal CKM matrix
    // we could factor out the CKM matrix elements here
    double factors2[] = { 1.0, 1.0, 1.0 };

    // define the channel #1
    pineappl_lumi_add(channels, 3, pids2, factors2);

    // ---
    // Specify the perturbative orders that will be filled into the grid

    // in this example we only fill the LO, which has the exponents
    // - 0 of alphas,
    // - 2 of alpha (electroweak coupling),
    // - 0 of log (xiR^2) (renormalization scale logarithm) and
    // - 0 of log (xiF^2) (factorization scale logarithm)
    uint32_t orders[] = {
        0, 2, 0, 0, // order #0: LO
        1, 2, 0, 0, // order #1: NLO QCD
        1, 2, 0, 1  // order #2: NLO QCD factorization log
    };

    // ---
    // Specify the bin limits

    // Similar to many Monte Carlo integrators PineAPPL supports only one-dimensional differential
    // distributions, and only one distribution for each grid. However, one can generate multiple
    // grids to support multiple distributions, and since every n-dimensional distribution can be
    // written as a one-dimensional one (by using the bin index as a new binning variable, for
    // instance), this isn't a limitation.

    // we bin the rapidity of the final-state lepton pair from 0 to 2.4 in steps of 0.1
    double bins[] = {
        0.0,
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
        1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4
    };

    // ---
    // Create the grid using the previously set information about orders, bins and channels

    // create the PineAPPL grid with default interpolation and binning parameters
    auto* keyval = pineappl_keyval_new();
    auto* grid = pineappl_grid_new(channels, 1, orders, 24, bins, keyval);

#ifdef USE_CUSTOM_GRID_PARAMETERS
    // set custom grid parameters. If left out, the standard values will be used, which are the ones
    // used below
    pineappl_keyval_set_int(keyval, "q2_bins", 40);
    pineappl_keyval_set_double(keyval, "q2_max", 1e8);
    pineappl_keyval_set_double(keyval, "q2_min", 1e2);
    pineappl_keyval_set_int(keyval, "q2_order", 3);
    pineappl_keyval_set_bool(keyval, "reweight", true);

    // Settings for all x-values (x1 and x2)
    pineappl_keyval_set_int(keyval, "x_bins", 50);
    pineappl_keyval_set_double(keyval, "x_max", 1.0);
    pineappl_keyval_set_double(keyval, "x_min", 2e-7);
    pineappl_keyval_set_int(keyval, "x_order", 3);

    // these parameters can be used to override the values specifically for x1
    pineappl_keyval_set_int(keyval, "x1_bins", 50);
    pineappl_keyval_set_double(keyval, "x1_max", 1.0);
    pineappl_keyval_set_double(keyval, "x1_min", 2e-7);
    pineappl_keyval_set_int(keyval, "x1_order", 3);

    // these parameters can be used to override the values specifically for x2
    pineappl_keyval_set_int(keyval, "x2_bins", 50);
    pineappl_keyval_set_double(keyval, "x2_max", 1.0);
    pineappl_keyval_set_double(keyval, "x2_min", 2e-7);
    pineappl_keyval_set_int(keyval, "x2_order", 3);
#endif

    // now we no longer need `keyval` and `lumi`
    pineappl_keyval_delete(keyval);
    pineappl_lumi_delete(channels);

    // ---
    // Fill the grid with phase-space points
    fill_grid(grid, 10000000);

    // ---
    // Write the grid to disk - with `.lz4` suffix the grid is automatically LZ4 compressed
    char const* filename =
#ifdef USE_CUSTOM_GRID_PARAMETERS
        "drell-yan-rap-ll-custom-grid.pineappl.lz4";
#else
        "drell-yan-rap-ll.pineappl.lz4";
#endif
    pineappl_grid_write(grid, filename);

    // destroy the object
    pineappl_grid_delete(grid);

    std::cout << "Generated " << filename << " containing a a -> l+ l-.\n\n"
        "Try running (PDF sets must contain non-zero photon PDF):\n"
        "  - pineappl convolute " << filename << " NNPDF31_nnlo_as_0118_luxqed\n"
        "  - pineappl --silence-lhapdf plot " << filename
        << " NNPDF31_nnlo_as_0118_luxqed MSHT20qed_nnlo > plot_script.py\n"
        "  - pineappl --help\n";
}
