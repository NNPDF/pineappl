#include <LHAPDF/LHAPDF.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

#include "PineAPPL.hpp"
#include "pineappl_capi.h"

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
  double cos_theta = 2.0 * r3 - 1.0;
  jacobian *= 2.0;

  double t = -0.5 * s * (1.0 - cos_theta);
  double u = -0.5 * s * (1.0 + cos_theta);

  // phi integration
  jacobian *= 2.0 * acos(-1.0);

  return {s, t, u, x1, x2, jacobian};
}

void fill_grid(PineAPPL::Grid& grid, std::size_t calls) {
  using std::acosh;
  using std::fabs;
  using std::log;
  using std::sqrt;

  auto rng = std::mt19937();

  // in GeV^2 pbarn
  double hbarc2 = 389379372.1;

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
    if ((ptl < 14.0) || (fabs(yll) > 2.4) || (ylp > 2.4) || (ylm > 2.4) ||
        (mll < 60.0) || (mll > 120.0)) {
      continue;
    }

    auto weight = jacobian * int_photo(s, t, u);
    double q2 = 90.0 * 90.0;

    std::vector<double> ntuples = {q2, x1, x2};
    grid.fill(0, fabs(yll), 0, ntuples, weight);
  }
}

int main() {
  // --- create a new `Channels` function for the $\gamma\gamma$ initial state
  PineAPPL::Channels channels;
  PineAPPL::SubChannelEntry subchannels;
  subchannels.entry.push_back({{22, 22}, 1.0});
  PineAPPL::ChannelsEntry channels_entry;
  channels_entry.channels_entry.push_back(subchannels);
  channels.add(channels_entry);

  // --- Instatiate the Order object
  // only LO $\alpha_\mathrm{s}^0 \alpha^2 \log^0(\xi_\mathrm{R})
  // \log^0(\xi_\mathrm{F}) \log^0(\xi_\mathrm{A})$
  std::vector<PineAPPL::Order> orders = {PineAPPL::Order{0, 2, 0, 0, 0}};

  // --- Define the binning
  // we bin in rapidity from 0 to 2.4 in steps of 0.1
  std::vector<double> bins = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
                              0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
                              1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4};

  // --- Construct the PineAPPL grid
  pineappl_pid_basis pid_basis = Evol;
  std::vector<std::int32_t> pids = {2212, 2212};

  // Define the types of convolutions
  pineappl_conv_type h1 = UnpolPDF;
  pineappl_conv_type h2 = UnpolPDF;
  std::vector<pineappl_conv_type> convolution_types = {h1, h2};

  // Define the Kinematics
  pineappl_kinematics scales = {Scale, 0};
  pineappl_kinematics x1 = {X, 0};
  pineappl_kinematics x2 = {X, 1};
  std::vector<pineappl_kinematics> kinematics = {scales, x1, x2};

  // Define the interpolation configurations
  pineappl_reweight_meth scales_reweight = NoReweight;  // Reweighting method
  pineappl_reweight_meth moment_reweight = ApplGridX;
  pineappl_map scales_mapping = ApplGridH0;  // Mapping method
  pineappl_map moment_mapping = ApplGridF2;
  pineappl_interp_meth interpolation_meth = Lagrange;
  std::vector<pineappl_interp_tuples> interpolations = {
      {1e2, 1e8, 40, 3, scales_reweight, scales_mapping,
       interpolation_meth},  // Interpolation fo `scales`
      {2e-7, 1.0, 50, 3, moment_reweight, moment_mapping,
       interpolation_meth},  // Interpolation fo `x1`
      {2e-7, 1.0, 50, 3, moment_reweight, moment_mapping,
       interpolation_meth},  // Interpolation fo `x2`
  };

  // Define the μ scale
  std::vector<std::size_t> mu_scales = {1, 1, 1};

  PineAPPL::Grid grid(orders, channels, pid_basis, pids, convolution_types,
                        kinematics, interpolations, bins, mu_scales);

  // fill the grid with phase-space points
  fill_grid(grid, 10000000);
  grid.optimize();

  // write the grid to disk - with `.lz4` suffix the grid is automatically LZ4
  // compressed
  const std::string filename = "DY-LO-AA.pineappl.lz4";
  grid.write(filename);

  std::printf(
      "Generated %s containing a a -> l+ l-.\n\n"
      "Try running (PDF sets must contain non-zero photon PDF):\n"
      "  - pineappl convolve %s NNPDF31_nnlo_as_0118_luxqed\n"
      "  - pineappl --silence-lhapdf plot %s NNPDF31_nnlo_as_0118_luxqed "
      "MSHT20qed_nnlo > plot_script.py\n"
      "  - pineappl --help\n",
      filename.c_str(), filename.c_str(), filename.c_str());
}
