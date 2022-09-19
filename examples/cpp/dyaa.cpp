#include <LHAPDF/LHAPDF.h>

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <random>
#include <vector>
#include <list>

#include "PineAPPL.hpp"

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

void fill_grid(PineAPPL::Grid* grid, std::size_t calls) {
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

        grid->fill(x1, x2, q2, 0, fabs(yll), 0, weight);
    }
}

int main() {
    // create a new luminosity function for the $\gamma\gamma$ initial state
    PineAPPL::Lumi lumi;
    lumi.add({PineAPPL::LumiEntry {22,22,1.0}});

    // only LO $\alpha_\mathrm{s}^0 \alpha^2 \log^0(\xi_\mathrm{R}) \log^0(\xi_\mathrm{F})$
    std::vector<PineAPPL::Order> orders = {PineAPPL::Order {0,2,0,0}};

    // we bin in rapidity from 0 to 2.4 in steps of 0.1
    std::vector<double> bins;
    double b = 0.;
    for (size_t j = 0; j < 24+1; ++j) {
        bins.push_back(b);
        b += 0.1;
    }

    // create the PineAPPL grid with default interpolation and binning parameters
    PineAPPL::KeyVal kv;
    PineAPPL::Grid grid(lumi, orders, bins, kv);

    // fill the grid with phase-space points
    fill_grid(&grid, 10000000);

    // perform a convolution of the grid with PDFs
    auto* pdf = LHAPDF::mkPDF("NNPDF31_nlo_as_0118_luxqed", 0);
    std::vector<double> dxsec = grid.convolute_with_one(2212, pdf);
    
    // print the results
    for (std::size_t j = 0; j != 24; ++j) {
        std::printf("%02d %.1f %.1f %.3e\n", j, bins[j], bins[j + 1], dxsec[j]);
    }

    // store some metadata in the grid
    grid.set_key_value("events", "10000000");

    // read out the stored value and print it on stdout
    const auto value = grid.get_key_value("events");
    std::printf("Finished running %s events.\n", value.c_str());

    // write the grid to disk - with `.lz4` suffix the grid is automatically LZ4 compressed
    std::string filename = "DY-LO-AA.pineappl.lz4";
    grid.write(filename);

    std::printf("Generated %s containing a a -> l+ l-.\n\n"
        "Try running (PDF sets must contain non-zero photon PDF):\n"
        "  - pineappl convolute %s NNPDF31_nnlo_as_0118_luxqed\n"
        "  - pineappl --silence-lhapdf plot %s NNPDF31_nnlo_as_0118_luxqed MSHT20qed_nnlo > plot_script.py\n"
        "  - pineappl --help\n", filename.c_str(), filename.c_str(), filename.c_str());
}
