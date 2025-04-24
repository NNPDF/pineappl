#include <pineappl_capi.h>

#include <cassert>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <memory>
#include <algorithm>
#include <vector>

double FAC0 = 1.65;

std::vector<OperatorInfo> get_operator_info(
    std::vector<double> fac1,
    std::vector<int> pids1,
    std::vector<double> x1,
    pineappl_pid_basis pid_basis,
    pineappl_conv_type conv_type
) {
    std::vector<OperatorInfo> opinfo_slices(fac1.size());
    std::vector<std::size_t> tensor_shape = {pids1.size(), x1.size(), pids1.size(), x1.size()};

    for (std::size_t q2 = 0; q2 != fac1.size(); q2++) {
        OperatorInfo opinfo = {
            FAC0, // fac0
            fac1[q2], // fac1
            x1.data(), // x0
            x1.data(), // x1
            pids1.data(), // pids0
            pids1.data(), // pids1
            pid_basis,
            conv_type,
            tensor_shape.data()
        };

        opinfo_slices[q2] = opinfo;
    }

    return opinfo_slices;
}

int main() {
    std::string filename = "drell-yan-rap-ll.pineappl.lz4";

    // read the grid from a file
    auto* grid = pineappl_grid_read(filename.c_str());

    // Get the number of perturbative orders
    std::size_t orders = pineappl_grid_order_count(grid);

    // Get the PID basis representation
    pineappl_pid_basis pid_basis = pineappl_grid_pid_basis(grid);

    // Get the number of convolutions
    std::size_t n_convs = pineappl_grid_convolutions_len(grid);

    // Fill the vector of unique convolution types
    std::vector<pineappl_conv_type> conv_types(n_convs);
    for (std::size_t i = 0; i != n_convs; i++) {
        pineappl_conv_type conv = pineappl_grid_conv_type(grid, i);
        if (std::find(conv_types.begin(), conv_types.end(), conv) == conv_types.end()) {
            conv_types.push_back(conv);
        }
    }

    // Get the shape of the evolve info objects
    std::vector<std::size_t> evinfo_shape(5);
    std::unique_ptr<bool[]> order_mask(new bool[orders]());
    order_mask[orders - 1] = true; // Choose the highest order
    pineappl_grid_evolve_info_shape(grid, order_mask.get(), evinfo_shape.data());

    // Get the values of the evolve info parameters. These contain, for example, the
    // information on the `x`-grid and `PID` used to interpolate the Grid.
    // NOTE: These are used to construct the Evolution Operator
    std::vector<double> fac1(evinfo_shape[0]);
    std::vector<double> frg1(evinfo_shape[1]);
    std::vector<int> pids1(evinfo_shape[2]);
    std::vector<double> x1(evinfo_shape[3]);
    std::vector<double> ren1(evinfo_shape[4]);
    pineappl_grid_evolve_info(grid, order_mask.get(), fac1.data(),
        frg1.data(), pids1.data(), x1.data(), ren1.data());

    // Construct the Operator Info
    std::vector<OperatorInfo> opinfo_slices(n_convs * fac1.size());
    for (std::size_t i = 0; i != n_convs; i++) {
        std::vector<OperatorInfo> opinfo = get_operator_info(fac1, pids1, x1, pid_basis, conv_types[i]);
        for (std::size_t j = 0; j != fac1.size(); j++) {
            opinfo_slices.push_back(opinfo[j]);
        }
    }

    pineappl_grid_delete(grid);
}
