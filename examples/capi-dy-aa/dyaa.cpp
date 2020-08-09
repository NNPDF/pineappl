#include <pineappl_capi.h>
#include <LHAPDF/LHAPDF.h>

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <random>
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

        auto weight = jacobian * int_photo(s, u, t);
        double  q2 = 90.0 * 90.0;

        pineappl_grid_fill(grid, x1, x2, q2, 0, fabs(yll), 0, weight);
    }
}

int main() {
    // create a new luminosity function for the $\gamma\gamma$ initial state
    auto* lumi = pineappl_lumi_new();
    int32_t pdg_ids[] = { 22, 22 };
    double ckm_factors[] = { 1.0 };
    pineappl_lumi_add(lumi, 1, pdg_ids, ckm_factors);

    // only LO $\alpha_\mathrm{s}^0 \alpha^2 \log^0(\xi_\mathrm{R}) \log^0(\xi_\mathrm{F})$
    uint32_t orders[] = { 0, 2, 0, 0 };

    // we bin in rapidity from 0 to 2.4 in steps of 0.1
    double bins[] = {
        0.0,
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
        1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4
    };

    // create the PineAPPL grid with default interpolation and binning parameters
    auto* keyval = pineappl_keyval_new();
    auto* grid = pineappl_grid_new(lumi, 1, orders, 24, bins, keyval);

    // now we no longer need `keyval` and `lumi`
    pineappl_keyval_delete(keyval);
    pineappl_lumi_delete(lumi);

    // fill the grid with phase-space points
    fill_grid(grid, 10000000);

    // perform a convolution of the grid with PDFs
    auto* pdf = LHAPDF::mkPDF("NNPDF31_nlo_as_0118_luxqed", 0);
    auto xfx = [](int32_t id, double x, double q2, void* pdf) {
        return static_cast <LHAPDF::PDF*> (pdf)->xfxQ2(id, x, q2);
    };
    auto alphas = [](double q2, void* pdf) {
        return static_cast <LHAPDF::PDF*> (pdf)->alphasQ2(q2);
    };

    std::vector<double> dxsec(24);
    pineappl_grid_convolute(grid, xfx, xfx, alphas, pdf, nullptr,
        nullptr, 1.0, 1.0, dxsec.data());

    // print the results
    for (std::size_t i = 0; i != 24; ++i) {
        std::printf("%.1f %.1f %.3e\n", bins[i], bins[i + 1], dxsec[i]);
    }

    // write the grid to disk
    pineappl_grid_write(grid, "DY-LO-AA.pineappl");

    // destroy the object
    pineappl_grid_delete(grid);
}
