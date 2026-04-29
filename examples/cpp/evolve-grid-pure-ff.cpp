#include <LHAPDF/PDF.h>
#include <pineappl_capi.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

/** @brief This struct can contain arbitrary parameters that need to be passed to Evolution
 * Operator Callback (`generate_identity_eko`).
 */
struct OperatorParams {
    std::vector<pineappl_conv_type> conv_types;
};

std::vector<std::size_t> unravel_index(std::size_t flat_index, const std::vector<std::size_t>& shape) {
    std::size_t ndim = shape.size();
    std::vector<std::size_t> coords(ndim);

    for (int i = static_cast<int>(ndim) - 1; i >= 0; --i) {
        coords[i] = flat_index % shape[static_cast<std::size_t>(i)];
        flat_index /= shape[static_cast<std::size_t>(i)];
    }

    return coords;
}

extern "C" void generate_identity_eko(
    std::size_t op_index,
    double /*fac1*/,
    const int* /*pids_in*/,
    const double* /*x_in*/,
    const int* /*pids_out*/,
    const double* /*x_out*/,
    const std::size_t* eko_shape,
    double* eko_buffer,
    void* params_state
) {
    OperatorParams* op_params = static_cast<OperatorParams*>(params_state);
    pineappl_conv_type conv_type = op_params->conv_types[op_index];
    assert(conv_type == PINEAPPL_CONV_TYPE_UNPOL_FF);

    std::vector<std::size_t> shape(eko_shape, eko_shape + 4);
    std::size_t flat_len = std::accumulate(shape.begin(), shape.end(), 1ULL, std::multiplies<std::size_t>());

    for (std::size_t i = 0; i != flat_len; i++) {
        std::vector<std::size_t> coords = unravel_index(i, shape);
        double delta_ik = (coords[0] == coords[2]) ? 1.0 : 0.0;
        double delta_jl = (coords[1] == coords[3]) ? 1.0 : 0.0;
        eko_buffer[i] = delta_ik * delta_jl;
    }
}

void print_results(const std::vector<double>& dxsec_grid, const std::vector<double>& dxsec_fktable) {
    const int idx_width = 6;
    const int num_width = 15;
    const int dif_width = 15;

    std::cout << std::setw(idx_width) << "Bin"
              << std::setw(num_width) << "Grid"
              << std::setw(num_width) << "FkTable"
              << std::setw(dif_width) << "reldiff" << std::endl;

    std::cout << std::setw(idx_width) << std::string(idx_width - 2, '-')
              << std::setw(num_width) << std::string(num_width - 2, '-')
              << std::setw(num_width) << std::string(num_width - 2, '-')
              << std::setw(dif_width) << std::string(dif_width - 2, '-') << std::endl;

    std::cout << std::scientific << std::setprecision(6);
    for (std::size_t i = 0; i < dxsec_grid.size(); ++i) {
        double reldiff = (dxsec_fktable[i] - dxsec_grid[i]) / dxsec_grid[i];
        std::cout << std::setw(idx_width) << i
                  << std::setw(num_width) << dxsec_grid[i]
                  << std::setw(num_width) << dxsec_fktable[i]
                  << std::setw(dif_width) << reldiff
                  << std::endl;
    }
}

pineappl_grid* create_and_fill_pure_ff_grid() {
    std::size_t nb_convolutions = 1;
    auto* channels = pineappl_channels_new(nb_convolutions);

    int32_t pids1[] = { 21 };
    double factors1[] = { 1.0 };
    pineappl_channels_add(channels, 1, pids1, factors1);

    // (alphas^0, alpha^0, logxir^0, logxif^0, logxia^0)
    std::vector<uint8_t> orders = { 0, 0, 0, 0, 0 };

    std::vector<double> bins = { 0.0, 0.5, 1.0 };

    pineappl_pid_basis pid_basis = PINEAPPL_PID_BASIS_PDG;
    pineappl_conv convs[] = { { PINEAPPL_CONV_TYPE_UNPOL_FF, 211 } }; // hadron PID: pi+

    pineappl_kinematics scales = { PINEAPPL_KINEMATICS_SCALE, 0 };
    pineappl_kinematics z = { PINEAPPL_KINEMATICS_X, 0 };
    pineappl_kinematics kinematics[2] = { scales, z };

    pineappl_reweight_meth scales_reweight = PINEAPPL_REWEIGHT_METH_NO_REWEIGHT;
    pineappl_reweight_meth moment_reweight = PINEAPPL_REWEIGHT_METH_APPL_GRID_X;
    pineappl_map scales_mapping = PINEAPPL_MAP_APPL_GRID_H0;
    pineappl_map moment_mapping = PINEAPPL_MAP_APPL_GRID_F2;
    pineappl_interp_meth interpolation_meth = PINEAPPL_INTERP_METH_LAGRANGE;
    pineappl_interp interpolations[2] = {
        { 1e2, 1e2 * (1.0 + 1e-12), 2, 1, scales_reweight, scales_mapping, interpolation_meth },
        { 2e-7, 1.0, 50, 3, moment_reweight, moment_mapping, interpolation_meth },
    };

    pineappl_scale_func_form scale_mu = { PINEAPPL_SCALE_FUNC_FORM_SCALE, 0 };
    pineappl_scale_func_form no_scale_mu = { PINEAPPL_SCALE_FUNC_FORM_NO_SCALE, 0 };
    pineappl_scale_func_form mu_scales[3] = {
        scale_mu,    // ren
        no_scale_mu, // fac
        scale_mu,    // frg
    };

    auto* grid = pineappl_grid_new2(
        bins.size() - 1,
        bins.data(),
        orders.size() / 5,
        orders.data(),
        channels,
        pid_basis,
        convs,
        2,
        interpolations,
        kinematics,
        mu_scales
    );

    pineappl_channels_delete(channels);

    std::size_t order = 0;
    std::size_t channel = 0;
    double q2 = 100.0; // GeV^2
    for (int i = 0; i != 10; ++i) {
        double zi = 0.1 + 0.08 * i;
        double obs = zi; // bin observable
        double weight = 1.0 + 0.1 * i;
        std::vector<double> ntuples = { q2, zi };
        pineappl_grid_fill2(grid, order, obs, channel, ntuples.data(), weight);
    }

    return grid;
}

int main() {
    LHAPDF::setVerbosity(0);

    std::string ffset = "MAPFF10NLOPIsum";
    auto ff = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(ffset, 0));

    auto xfx = [](int32_t id, double x, double q2, void* ff) {
        return static_cast<LHAPDF::PDF*>(ff)->xfxQ2(id, x, q2);
    };
    auto alphas = [](double q2, void* ff) {
        return static_cast<LHAPDF::PDF*>(ff)->alphasQ2(q2);
    };

    auto* grid = create_and_fill_pure_ff_grid();

    std::vector<std::size_t> evinfo_shape(5);
    pineappl_grid_evolve_info_shape(grid, nullptr, evinfo_shape.data());

    std::vector<double> fac1(evinfo_shape[0]);
    std::vector<double> frg1(evinfo_shape[1]);
    std::vector<int> pids_in(evinfo_shape[2]);
    std::vector<double> x_in(evinfo_shape[3]);
    std::vector<double> ren1(evinfo_shape[4]);

    pineappl_grid_evolve_info(grid, nullptr, fac1.data(), frg1.data(), pids_in.data(), x_in.data(), ren1.data());

    const std::vector<double>& proc1_all = fac1.empty() ? frg1 : fac1;
    assert(!proc1_all.empty());

    std::size_t n_convs = pineappl_grid_convolutions_len(grid);
    assert(n_convs == 1);
    std::vector<pineappl_conv_type> conv_types(n_convs);
    pineappl_grid_conv_types(grid, conv_types.data());

    std::vector<pineappl_conv_type> unique_convs;
    for (std::size_t i = 0; i != n_convs; i++) {
        pineappl_conv_type conv = conv_types[i];
        if (std::find(unique_convs.begin(), unique_convs.end(), conv) == unique_convs.end()) {
            unique_convs.push_back(conv);
        }
    }
    assert(unique_convs.size() == 1);
    assert(unique_convs[0] == PINEAPPL_CONV_TYPE_UNPOL_FF);

    pineappl_pid_basis pid_basis = pineappl_grid_pid_basis(grid);
    const double fac0 = proc1_all.front();
    std::vector<pineappl_operator_info> opinfo_slices(unique_convs.size() * proc1_all.size());
    for (std::size_t i = 0; i != unique_convs.size(); i++) {
        for (std::size_t j = 0; j != proc1_all.size(); j++) {
            pineappl_operator_info opinfo = {
                fac0,           // fac0: common starting scale
                proc1_all[j],   // fac1: process scale (here: fragmentation scale)
                pid_basis,
                unique_convs[i],
            };
            opinfo_slices[i * proc1_all.size() + j] = opinfo;
        }
    }

    std::vector<double> alphas_table;
    alphas_table.reserve(ren1.size());
    for (double q2 : ren1) {
        alphas_table.push_back(alphas(q2, ff.get()));
    }

    OperatorParams op_params;
    op_params.conv_types = unique_convs;

    std::vector<double> xi = { 1.0, 1.0, 1.0 };
    std::vector<int> pids_out = pids_in;
    std::vector<std::size_t> tensor_shape = { pids_in.size(), x_in.size(), pids_out.size(), x_in.size() };

    pineappl_grid* fktable = pineappl_grid_evolve(
        grid,
        unique_convs.size(),
        generate_identity_eko,
        opinfo_slices.data(),
        pids_in.data(),
        x_in.data(),
        pids_out.data(),
        x_in.data(),
        tensor_shape.data(),
        &op_params,
        nullptr,
        xi.data(),
        ren1.data(),
        alphas_table.data()
    );

    // Compare grid vs evolved FK by convolution
    std::size_t bins = pineappl_grid_bin_count(grid);
    std::vector<double> mu_scales = { 1.0, 1.0, 1.0 };

    std::vector<LHAPDF::PDF*> ffs = { ff.get() };
    void** ff_states = reinterpret_cast<void**>(ffs.data());

    std::vector<double> dxsec_grid(bins);
    pineappl_grid_convolve(
        grid,
        xfx,
        alphas,
        ff_states,
        ff.get(),
        nullptr,
        nullptr,
        nullptr,
        1,
        mu_scales.data(),
        dxsec_grid.data()
    );

    std::vector<double> dxsec_fktable(bins);
    auto as_one = [](double /*q2*/, void* /*state*/) { return 1.0; };
    pineappl_grid_convolve(
        fktable,
        xfx,
        as_one,
        ff_states,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        1,
        mu_scales.data(),
        dxsec_fktable.data()
    );

    print_results(dxsec_grid, dxsec_fktable);

    pineappl_grid_write(grid, "pure-ff-grid.pineappl.lz4");
    pineappl_grid_write(fktable, "pure-ff-fktable.pineappl.lz4");

    pineappl_grid_delete(grid);
    pineappl_grid_delete(fktable);
}
