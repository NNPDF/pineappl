#include <cstdint>
#include <pineappl_capi.h>

#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <random>
#include <string>
#include <tuple>
#include <vector>

using KinematicsTuple = std::tuple<double, double, double>;

struct KinInterpolation {
    std::vector<double> x_interp;
    std::vector<double> z_interp;
};

template<typename T>
std::vector<T> geomspace(T start, T stop, int num, bool endpoint = false) {
    std::vector<T> result(num);

    if (num == 1) {
        result[0] = start;
        return result;
    }

    T log_start = std::log(start);
    T log_stop = std::log(stop);
    T step = (log_stop - log_start) / (endpoint ? (num - 1) : num);

    for (int i = 0; i < num; ++i) {
        result[i] = std::exp(log_start + i * step);
    }

    return result;
}

KinInterpolation kinematics_interpolation_points() {
    std::vector<double> x = geomspace(1e-5, 1.0, 50);
    std::vector<double> z = geomspace(1e-5, 1.0, 50);

    return {x, x};
}

std::vector<double> generate_subgrid_arrays(
    std::mt19937& rng,
    std::vector<double> x,
    std::vector<double> z
) {
    std::size_t kin_length = x.size() * z.size();
    std::vector<double> subgrid(kin_length);

    // NOTE: `subgrid` is a flatten matrix whose layout was `[q2=1][x][z]`.
    // The order of the kinematics shoud match the kinematics declaration.
    for (std::size_t i = 0; i != kin_length; i++) {
        subgrid[i] = std::generate_canonical<double, 53>(rng);
    }

    return subgrid;
}

void fill_grid(pineappl_grid* grid, std::vector<KinematicsTuple> kin_tuples) {
    auto rng = std::mt19937();

    auto* channels = pineappl_grid_channels(grid);
    std::size_t n_bins = pineappl_grid_bin_count(grid);
    std::size_t n_orders = pineappl_grid_order_count(grid);
    std::size_t n_channels = pineappl_channels_count(channels);

    // Get the kinematics
    KinInterpolation kins = kinematics_interpolation_points();
    std::vector<double> x_interp = kins.x_interp;
    std::vector<double> z_interp = kins.z_interp;

    // Extract the shape of the subgrid - Q2 always passed as an array of ONE element
    std::vector<std::size_t> subgrid_shape = {1, x_interp.size(), z_interp.size()};

    for (std::size_t b = 0; b != n_bins; b++) {
        for (std::size_t o = 0; o != n_orders; o++) {
            for (std::size_t c = 0; c != n_channels; c++) {
                // Construct the node values of {Q2, x_inter, z_interp}
                // NOTE: Pay attention to the order, it should match the kinematics declaration
                // and how the subgrid was constructed (see `generate_subgrid_array`).
                std::vector<double> node_values = { std::get<0>(kin_tuples[b]) };
                node_values.insert(node_values.end(), x_interp.begin(), x_interp.end());
                node_values.insert(node_values.end(), z_interp.begin(), z_interp.end());

                // Mock the subgrids for a given bin, order, and channel
                std::vector<double> subgrid_arrays = generate_subgrid_arrays(rng, x_interp, z_interp);

                // set the subgrids
                pineappl_grid_set_subgrid(
                    grid,
                    b, o, c, // kinematics index
                    node_values.data(),
                    subgrid_arrays.data(),
                    subgrid_shape.data()
                );
            }
        }
    }
}

int main() {
    // ---
    // Create all channels

    std::size_t nb_convolutions = 2;
    auto* channels = pineappl_channels_new(nb_convolutions);

    int32_t pids1[] = { 21, 21 };
    double factors1[] = { 1.0 };
    // define the channel #0
    pineappl_channels_add(channels, 1, pids1, factors1);

    // create another channel this channel is the down-type-antidown-type quark channel; here we
    int32_t pids2[] = { 1, -1, 3, -3, 5, -5 };
    // define the channel #1
    pineappl_channels_add(channels, 3, pids2, nullptr);

    // ---
    // Specify the perturbative orders that will be filled into the grid
    std::vector<uint8_t> orders = {
        1, 0, 0, 0, 0, // order #0: LO QCD
        2, 0, 0, 0, 0, // order #1: NLO QCD
    };

    // ---
    // Specify the bin limits

    // In SIDIS, a bin is defined as a tuple (Q2, x, z) values (3D).
    std::vector<KinematicsTuple> kin_obs = {std::make_tuple(1e3, 1e-5, 1e-2), std::make_tuple(1e4, 1e-2, 1e-3)};
    // We are going to define some placeholder 1D bins that we'll overwrite later.
    std::vector<double> bins;
    for (std::size_t i = 0; i < kin_obs.size() + 1; ++i) {
        bins.push_back(static_cast<float>(i));
    }

    // ---
    // Construct the objects that are needed to fill the Grid

    pineappl_pid_basis pid_basis = PINEAPPL_PID_BASIS_EVOL;
    pineappl_conv convs[] = {
        { PINEAPPL_CONV_TYPE_UNPOL_PDF, 2212 },
        { PINEAPPL_CONV_TYPE_UNPOL_FF, 211 }, // Assumes Pion
    };

    // Define the kinematics required for this process. In the following example we have ONE single
    // scale and two momentum fractions (corresponding to the two initial- and final-state hadrons).
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

    // Define the unphysical scale objects
    pineappl_scale_func_form scale_mu = { PINEAPPL_SCALE_FUNC_FORM_SCALE, 0 };
    pineappl_scale_func_form mu_scales[3] = { scale_mu, scale_mu, scale_mu };

    // ---
    // Create the grid using the previously set information about orders, bins and channels

    auto* grid = pineappl_grid_new2(bins.size() - 1, bins.data(), orders.size() / 5, orders.data(),
        channels, pid_basis, convs, 3, interpolations, kinematics, mu_scales);

    pineappl_channels_delete(channels);

    // ---
    // Fill the grid with phase-space points
    fill_grid(grid, kin_obs);

    // ---
    // NOTE: We now have to remap the bins to 3D

    // We need to flatten the array of `KinematicsTuple` and define the normalizations
    std::vector<double> flat_kin_obs;
    for (const auto& kin : kin_obs) {
        flat_kin_obs.push_back(std::get<0>(kin));
        flat_kin_obs.push_back(std::get<1>(kin));
        flat_kin_obs.push_back(std::get<2>(kin));
    }
    std::vector<double> normalizations(bins.size() - 1, 1.0);
    pineappl_grid_set_bwfl(
        grid,
        flat_kin_obs.data(), // lower bin limits
        flat_kin_obs.data(), // upper bin limits
        bins.size() - 1, // number of bins
        3, // Dimension of the bins (Q2, x, z)
        normalizations.data()
    );
    pineappl_grid_optimize(grid);

    // ---
    // Write the grid to disk - the filename can be anything ...
    std::string filename = "sidis-toygrid.pineappl.lz4";
    pineappl_grid_write(grid, filename.c_str());

    // destroy the object
    pineappl_grid_delete(grid);

    std::cout << "Generated " << filename << " containing a toy SIDIS.\n\n"
        "Try running the following command to check the bins:\n"
        "  - pineappl read --bins " << filename << "\n";
}
