#include <pineappl_capi.h>
#include <cstdint>
#include <cstdlib>
#include <cassert>
#include <cstddef>
#include <string>
#include <algorithm>
#include <vector>
#include <random>

double FAC0 = 1.65;

std::vector<double> generate_fake_ekos(
    double q2,
    std::vector<int> pids0,
    std::vector<double> x0,
    std::vector<int> pids1,
    std::vector<double> x1
) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distrib(q2 / 1000, q2 / 100);

    std::size_t flat_len = x0.size() * x1.size() * pids0.size() * pids1.size();
    std::vector<double> ops(flat_len);

    for (std::size_t i = 0; i != flat_len; i++) {
        ops[i] = distrib(gen);
    }

    return ops;
}

int main() {
    // TODO: How to get a Grid that can be evolved??
    std::string filename = "LHCB_WP_7TEV_opt.pineappl.lz4";

    // read the grid from a file
    auto* grid = pineappl_grid_read(filename.c_str());

    // Get the PID basis representation
    pineappl_pid_basis pid_basis = pineappl_grid_pid_basis(grid);

    // Get the number of convolutions
    std::size_t n_convs = pineappl_grid_convolutions_len(grid);

    // Fill the vector of unique convolution types. If the EKOs required for the Grid
    // are the same, then it suffices to only pass ONE single EKO.
    std::vector<pineappl_conv_type> conv_types;
    for (std::size_t i = 0; i != n_convs; i++) {
        pineappl_conv_type conv = pineappl_grid_conv_type(grid, i);
        if (std::find(conv_types.begin(), conv_types.end(), conv) == conv_types.end()) {
            conv_types.push_back(conv);
        }
    }

    // Get the shape of the evolve info objects
    std::vector<std::size_t> evinfo_shape(5);
    std::vector<uint8_t> order_mask = {3, 0};
    pineappl_grid_evolve_info_shape(grid, order_mask.data(), evinfo_shape.data());

    // Get the values of the evolve info parameters. These contain, for example, the
    // information on the `x`-grid and `PID` used to interpolate the Grid.
    // NOTE: These are used to construct the Evolution Operator
    std::vector<double> fac1(evinfo_shape[0]);
    std::vector<double> frg1(evinfo_shape[1]);
    std::vector<int> pids1(evinfo_shape[2]);
    std::vector<double> x1(evinfo_shape[3]);
    std::vector<double> ren1(evinfo_shape[4]);
    pineappl_grid_evolve_info(grid, order_mask.data(), fac1.data(),
        frg1.data(), pids1.data(), x1.data(), ren1.data());

    // ------------------ Construct the Operator Info ------------------
    // The Operator Info is a vector with length `N_conv * N_Q2_slices` whose
    // elements are `OperatorInfo` objects.
    std::vector<pineappl_operator_info> opinfo_slices(conv_types.size() * fac1.size());
    for (std::size_t i = 0; i != conv_types.size(); i++) {
        for (std::size_t j = 0; j != fac1.size(); j++) {
            pineappl_operator_info opinfo = {
                FAC0, // fac0
                fac1[j], // fac1
                pid_basis,
                conv_types[i],
            };
            opinfo_slices[i * fac1.size() + j] = opinfo;
        }
    }

    // ------------------ Construct the Evolution Operator ------------------
    // The Evolution Operator is a vector with length `N_conv * N_Q2_slices * Σ product(OP shape)`
    std::vector<double> op_slices;
    std::size_t flat_len = x1.size() * x1.size() * pids1.size() * pids1.size();
    for (std::size_t i = 0; i != conv_types.size(); i++) {
        for (std::size_t j = 0; j != fac1.size(); j++) {
            std::vector<double> eko = generate_fake_ekos(fac1[j], pids1, x1, pids1, x1);
            for (std::size_t k = 0; k != flat_len; k++) {
                op_slices.push_back(eko[k]);
            }
        }
    }

    // Construct the values of alphas
    std::vector<double> alphas_table(ren1.size(), 0.118);
    std::vector<double> xi = {1.0, 1.0, 1.0};
    std::vector<std::size_t> tensor_shape = {pids1.size(), x1.size(), pids1.size(), x1.size()};

    pineappl_fk_table* fktable = pineappl_grid_evolve(grid, opinfo_slices.data(), order_mask.data(),
        op_slices.data(), x1.data(), x1.data(), pids1.data(), pids1.data(),
        tensor_shape.data(), xi.data(), ren1.data(), alphas_table.data());

    pineappl_fktable_write(fktable, "evolved-grid.pineappl.lz4");

    pineappl_grid_delete(grid);
    pineappl_fk_table_delete(fktable);
}
