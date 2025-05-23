#include <cstdint>
#include <pineappl_capi.h>

#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <random>
#include <string>
#include <vector>

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

        // Values of the kinematic variables
        std::vector<double> ntuples = {q2, x1, x2};

        // fill the LO `weight` into `grid` for parton fractions `x1` and `x2`, and the (squared)
        // renormalization/factorization scale `q2`. The parameters `order` and `channel` are
        // indices defined from the arrays `orders` and `channel` used in creating the grid. In this
        // case they are both `0` and denote the order #0 (leading order) and the channel #0
        // (photon-photon channel), respectively
        pineappl_grid_fill2(grid, order, fabs(yll), channel, ntuples.data(), weight);
    }
}

int main() {
    // ---
    // Create all channels

    // this object will contain all channels (for two initial states) that we define
    std::size_t nb_convolutions = 2;
    auto* channels = pineappl_channels_new(nb_convolutions);

    // photon-photon initial state, where `22` is the photon (PDG MC ids)
    int32_t pids1[] = { 22, 22 };

    // factor that each channel is multiplied with when convoluting with PDFs
    double factors1[] = { 1.0 };

    // define the channel #0
    pineappl_channels_add(channels, 1, pids1, factors1);

    // create another channel, which we won't fill, however

    // this channel is the down-type-antidown-type quark channel; here we combine down-antidown,
    // strange-antistrange and bottom-antibottom into a single channel, which is often done if the
    // CKM matrix is taken to be diagonal
    int32_t pids2[] = { 1, -1, 3, -3, 5, -5 };

    // for each pair of particle ids we need to give a factor; in case of a non-diagonal CKM matrix
    // we could factor out the CKM matrix elements in this array and still treat the down-type
    // contributions in a single channel. In this case, however, all factors are `1.0`, for which we
    // can also pass `nullptr`

    // define the channel #1
    pineappl_channels_add(channels, 3, pids2, nullptr);

    // ---
    // Specify the perturbative orders that will be filled into the grid

    // in this example we only fill the LO, which has the exponents
    // - 0 of alphas,
    // - 2 of alpha (electroweak coupling),
    // - 0 of log (xiR^2) (renormalization scale logarithm) and
    // - 0 of log (xiF^2) (factorization scale logarithm)
    std::vector<uint8_t> orders = {
        0, 2, 0, 0, 0, // order #0: LO
        1, 2, 0, 0, 0, // order #1: NLO QCD
        1, 2, 0, 1, 0  // order #2: NLO QCD factorization log
    };

    // ---
    // Specify the bin limits

    // Similar to many Monte Carlo integrators PineAPPL supports only one-dimensional differential
    // distributions, and only one distribution for each grid. However, one can generate multiple
    // grids to support multiple distributions, and since every n-dimensional distribution can be
    // written as a one-dimensional one (by using the bin index as a new binning variable, for
    // instance), this isn't a limitation.

    // we bin the rapidity of the final-state lepton pair from 0 to 2.4 in steps of 0.1
    std::vector<double> bins = {
        0.0,
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
        1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4
    };

    // ---
    // Construct the objects that are needed to fill the Grid

    // First we define the types of convolutions required by the involved initial-/final-state
    // hadrons. Then we add the corresponding PID of each of the hadrons, and finally define the
    // Basis onto which the partons are mapped.
    pineappl_pid_basis pid_basis = PINEAPPL_PID_BASIS_EVOL;
    pineappl_conv convs[] = {
        { PINEAPPL_CONV_TYPE_UNPOL_PDF, 2212 },
        { PINEAPPL_CONV_TYPE_UNPOL_PDF, 2212 },
    };

    // Define the kinematics required for this process. In the following example we have ONE
    // single scale and two momentum fractions (corresponding to the two initial-state hadrons).
    // The format of the kinematics is: { type, value }.
    pineappl_kinematics scales = { PINEAPPL_KINEMATICS_SCALE, 0 };
    pineappl_kinematics x1 = { PINEAPPL_KINEMATICS_X, 0 };
    pineappl_kinematics x2 = { PINEAPPL_KINEMATICS_X, 1 };
    pineappl_kinematics kinematics[3] = { scales, x1, x2 };

    // Define the specificities of the interpolations for each of the kinematic variables.
    pineappl_reweight_meth scales_reweight = PINEAPPL_REWEIGHT_METH_NO_REWEIGHT; // Reweighting method
    pineappl_reweight_meth moment_reweight = PINEAPPL_REWEIGHT_METH_APPL_GRID_X;
    pineappl_map scales_mapping = PINEAPPL_MAP_APPL_GRID_H0; // Mapping method
    pineappl_map moment_mapping = PINEAPPL_MAP_APPL_GRID_F2;
    pineappl_interp_meth interpolation_meth = PINEAPPL_INTERP_METH_LAGRANGE;
    pineappl_interp interpolations[3] = {
        { 1e2, 1e8, 40, 3, scales_reweight, scales_mapping, interpolation_meth },  // Interpolation fo `scales`
        { 2e-7, 1.0, 50, 3, moment_reweight, moment_mapping, interpolation_meth }, // Interpolation fo `x1`
        { 2e-7, 1.0, 50, 3, moment_reweight, moment_mapping, interpolation_meth }, // Interpolation fo `x2`
    };

    // Define the unphysical scale objecs
    pineappl_scale_func_form scale_mu = { PINEAPPL_SCALE_FUNC_FORM_SCALE, 0 };
    pineappl_scale_func_form no_scale_mu = { PINEAPPL_SCALE_FUNC_FORM_NO_SCALE, 0 }; // Here `.scale=0` is dummy value
    pineappl_scale_func_form mu_scales[3] = { scale_mu, scale_mu, no_scale_mu };

    // ---
    // Create the grid using the previously set information about orders, bins and channels

    // create a new grid with the previously defined channels, 3 perturbative orders defined by the
    // exponents in `orders`, 24 bins given as the 25 limits in `bins` and potential extra
    // parameters in `keyval`.
    auto* grid = pineappl_grid_new2(bins.size() - 1, bins.data(), orders.size() / 5, orders.data(),
        channels, pid_basis, convs, 3, interpolations, kinematics, mu_scales);

    // now we no longer need `keyval` and `channels`
    pineappl_channels_delete(channels);

    // ---
    // Fill the grid with phase-space points
    fill_grid(grid, 10000000);

    std::string filename = "drell-yan-rap-ll.pineappl";

    // ---
    // Write the grid to disk - the filename can be anything ...
    pineappl_grid_write(grid, filename.c_str());

    // but if it has an `.lz4` suffix ...
    filename.append(".lz4");
    // the grid is automatically LZ4 compressed
    pineappl_grid_write(grid, filename.c_str());

    // destroy the object
    pineappl_grid_delete(grid);

    std::cout << "Generated " << filename << " containing a a -> l+ l-.\n\n"
        "Try running (PDF sets must contain non-zero photon PDF):\n"
        "  - pineappl convolve " << filename << " NNPDF31_nnlo_as_0118_luxqed\n"
        "  - pineappl --silence-lhapdf plot " << filename
        << " NNPDF31_nnlo_as_0118_luxqed MSHT20qed_nnlo > plot_script.py\n"
        "  - pineappl --help\n";
}
