/// Example used in the PineAPPL v1 paper ///
#include <cstdint>
#include <pineappl_capi.h>

#include <cassert>
#include <cmath>
#include <cstddef>
#include <random>
#include <string>
#include <vector>

struct Psp2to2Hadron {
    double s;
    double t;
    double u;
    double x1;
    double x2;
    double z;
    double pt_hadron;
    double y_hadron;
    double jacobian;
};

double me_gg2qqbar(double s, double t, double u) {
    (void) s; // ignore dummy variable
    double as2 = 0.118 * 0.118;
    double PI2 = M_PI * M_PI;
    // TODO: double-check
    return (16 * PI2 * as2 / 6.0) * (std::pow(u, 2) + std::pow(t, 2)) / (u * t);
}

Psp2to2Hadron pspgen_pp2hadron(std::mt19937& rng, double mmin,
    double mmax, double pt_min, double pt_max, double abs_y_max) {
    using std::acos;
    using std::log;
    using std::pow;
    using std::exp;
    using std::sqrt;

    double smin = mmin * mmin;
    double smax = mmax * mmax;

    double r1 = std::generate_canonical<double, 53>(rng);
    double r2 = std::generate_canonical<double, 53>(rng);
    double r3 = std::generate_canonical<double, 53>(rng);
    double r4 = std::generate_canonical<double, 53>(rng);
    double r5 = std::generate_canonical<double, 53>(rng);
    double r6 = std::generate_canonical<double, 53>(rng);

    double tau0 = smin / smax;
    double tau = pow(tau0, r1);
    double y = pow(tau, 1.0 - r2);
    double x1 = y;
    double x2 = tau / y;
    double s = tau * smax;

    double jacobian = tau * log(tau0) * log(tau0) * r1;

    // `theta` integration
    double cos_theta = 2.0 * r3 - 1.0;
    jacobian *= 2.0;

    double t = -0.5 * s * (1.0 - cos_theta);
    double u = -0.5 * s * (1.0 + cos_theta);

    // `phi` integration
    jacobian *= 2.0 * acos(-1.0);

    // sample hadron `pT` uniformly in log scale
    double log_pt_min = log(pt_min);
    double log_pt_max = log(pt_max);

    double pt_hadron = exp(log_pt_min + (log_pt_max - log_pt_min) * r4);
    jacobian *= pt_hadron * (log_pt_max - log_pt_min);

    // sample hadron rapidity uniformly
    double y_hadron = 2.0 * abs_y_max * r5 - abs_y_max;
    jacobian *= 2.0 * abs_y_max;

    // define the momentum fracion `z`
    double z_min = pt_hadron * exp(-y_hadron) / sqrt(s);
    double z_max_kin = pt_hadron * exp(y_hadron) / sqrt(s);
    double z_max = std::min(1.0, z_max_kin);

    // ensure that `z` is physical
    if ((z_min >= 1) || (z_min >= z_max)) {
        return {s, t, u, x1, x2, 0.0, pt_hadron, y_hadron, 0.0};
    }

    // sample `z` uniformly between the kinematic limits
    double z = z_min + (z_max - z_min) * r6;
    jacobian *= (z_max - z_min);

    return {s, t, u, x1, x2, z, pt_hadron, y_hadron, jacobian};
}

void fill_grid(pineappl_grid* grid, std::size_t calls) {
    using std::acosh;
    using std::fabs;
    using std::log;
    using std::sqrt;

    auto rng = std::mt19937();
    double hbarc2 = 389379372.1;

    // define hadron kinematic ranges
    double pt_min = 5.0;   // GeV
    double pt_max = 100.0; // GeV
    double abs_y_max = 2.4;    // rapidity range

    for (std::size_t i = 0; i != calls; ++i) {
        auto tmp = pspgen_pp2hadron(rng, 3000.0, 14000.0, pt_min, pt_max, abs_y_max);
        auto s = tmp.s;
        auto t = tmp.t;
        auto u = tmp.u;
        auto x1 = tmp.x1;
        auto x2 = tmp.x2;
        auto z = tmp.z;
        auto pt_hadron = tmp.pt_hadron;
        auto y_hadron = tmp.y_hadron;
        auto jacobian = tmp.jacobian;

        // skip if kinematically forbidden
        if (jacobian == 0.0 || z <= 0.0) {
            continue;
        }

        // apply cuts on hadron kinematics
        if ((pt_hadron < pt_min) || (pt_hadron > pt_max) || (fabs(y_hadron) > abs_y_max)) {
            continue;
        }

        jacobian *= hbarc2 / calls;

        // calculate the partonic cross-section
        auto weight = jacobian * me_gg2qqbar(s, t, u);

        double q2 = pt_hadron * pt_hadron;
        std::size_t order = 0;
        std::size_t channel = 0;

        // define the tuple of kinematic variables `ntuples = (q2, x1, x2, z)`
        std::vector<double> ntuples = {q2, x1, x2, z};

        // Fill the grid using hadron `pT` as the observable
        pineappl_grid_fill2(grid, order, pt_hadron, channel, ntuples.data(), weight);
    }
}

int main() {
    // ---
    // Define the partonic channels and orders that will be filled into the grid

    // specify the number of convolutions: 2 for initial-state PDFs + 1 for FFs
    std::size_t nb_convolutions = 3;

    // instantiate the channel object
    auto* channels = pineappl_channels_new(nb_convolutions);

    // specify the contributing channel(s) and the corresponding factor(s)
    // for the process `gg -> qqbar` we need to sum over the light quarks
    std::vector<int32_t> pids;
    std::vector<double> factors;
    for (int i = -3; i <= 3; ++i) {
        if (i == 0) continue;
        pids.insert(pids.end(), {21, 21, i});
        factors.push_back(1.0);
    }
    pineappl_channels_add(channels, pids.size() / nb_convolutions, pids.data(),
        factors.data());

    // specify the perturbative orders that will be filled into the grid
    // orders specifies the power of the tuple `orders = (αs, α, lR, lF, lD)`
    // in this example, we only fill the LO QCD
    std::vector<uint8_t> orders = {1, 0, 0, 0, 0};

    // bin limits of the final-state hadron transverse momentum
    std::vector<double> bins = {
        5.0, 7.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0,
        60.0, 70.0, 80.0, 90.0, 100.0
    };

    // ---
    // Construct the objects that are needed to fill the grid

    // choose the Evolution Basis to represent the grid
    pineappl_pid_basis pid_basis = PINEAPPL_PID_BASIS_EVOL;
    // define the types of hadrons and set them to be Unpolarised
    pineappl_conv convs[] = {
        {PINEAPPL_CONV_TYPE_UNPOL_PDF, 2212}, // proton
        {PINEAPPL_CONV_TYPE_UNPOL_PDF, 2212}, // proton
        {PINEAPPL_CONV_TYPE_UNPOL_FF, 211},   // pion
    };

    // define the kinematics object `kinematics = (μ, x1, x2, x)`
    pineappl_kinematics scales = {PINEAPPL_KINEMATICS_SCALE, 0};
    pineappl_kinematics x1 = {PINEAPPL_KINEMATICS_X, 0};
    pineappl_kinematics x2 = {PINEAPPL_KINEMATICS_X, 1};
    pineappl_kinematics z = {PINEAPPL_KINEMATICS_X, 2};
    pineappl_kinematics kinematics[4] = {scales, x1, x2, z};

    // define the specificities of the interpolations `interpolations = (μ, x1, x2, z)`
    pineappl_reweight_meth scales_reweight = PINEAPPL_REWEIGHT_METH_NO_REWEIGHT;
    pineappl_reweight_meth moment_reweight = PINEAPPL_REWEIGHT_METH_APPL_GRID_X;
    pineappl_map scales_mapping = PINEAPPL_MAP_APPL_GRID_H0;
    pineappl_map moment_mapping = PINEAPPL_MAP_APPL_GRID_F2;
    pineappl_interp_meth interpolation_meth = PINEAPPL_INTERP_METH_LAGRANGE;
    pineappl_interp interpolations[4] = {
        {1e2, 1e8, 40, 3, scales_reweight, scales_mapping, interpolation_meth},  // μ
        {2e-7, 1.0, 50, 3, moment_reweight, moment_mapping, interpolation_meth}, // x1
        {2e-7, 1.0, 50, 3, moment_reweight, moment_mapping, interpolation_meth}, // x2
        {2e-7, 1.0, 50, 3, moment_reweight, moment_mapping, interpolation_meth},  // z
    };

    // define the values of the unphysical scales `mu_scales = (μR, μF, μD)`
    // where here we do not consider the fragmentation scale μD
    pineappl_scale_func_form scale_mu = {PINEAPPL_SCALE_FUNC_FORM_SCALE, 0};
    pineappl_scale_func_form mu_scales[3] = {scale_mu, scale_mu, scale_mu};

    // ---
    // Create the grid, fill it with Monte Carlo weights, and dump into disk

    auto* grid = pineappl_grid_new2(bins.size() - 1, bins.data(), orders.size() / 5,
        orders.data(), channels, pid_basis, convs, nb_convolutions + 1,
        interpolations, kinematics, mu_scales);

    // delete no longer needed channel object
    pineappl_channels_delete(channels);

    // fill the grid with phase-space points
    fill_grid(grid, 100000);

    // add some metadata to the grid
    pineappl_grid_set_key_value(grid, "x1_label", "pT");
    pineappl_grid_set_key_value(grid, "y_label", "dsig/dpT");
    pineappl_grid_set_key_value(grid, "x1_unit", "GeV");
    pineappl_grid_set_key_value(grid, "y_unit", "pb/GeV");

    // write the grid into disk
    std::string filename = "pp2hadron-pt.pineappl.lz4";
    pineappl_grid_write(grid, filename.c_str());

    // remove grid object from memory
    pineappl_grid_delete(grid);
}
