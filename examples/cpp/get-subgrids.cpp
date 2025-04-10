#include <pineappl_capi.h>

#include <cstddef>
#include <iomanip>
#include <ios>
#include <iostream>
#include <string>
#include <vector>
#include <numeric>

int main() {
    std::string filename = "drell-yan-rap-ll.pineappl.lz4";

    // read the grid from a file
    auto* grid = pineappl_grid_read(filename.c_str());

    // Determine the number of bins and the index of order and channel
    std::size_t n_bins = pineappl_grid_bin_count(grid);
    std::size_t order = 0;
    std::size_t channel = 0;

    // Get the dimension of the subgrids
    std::size_t subgrid_dim = pineappl_convolutions_len(grid) + 1;

    std::cout << std::right << std::setw(10) << "bin" << std::setw(10) << "sg idx"
        << std::setw(16) << "sg value" << "\n";
    std::cout << std::right << std::setw(10) << "---" << std::setw(10) << "------"
        << std::setw(16) << "------------" << "\n";

    for (std::size_t b = 0; b < n_bins; ++b) {
        std::vector<std::size_t> subgrid_shape(subgrid_dim);
        pineappl_subgrid_shape(grid, b, order, channel, subgrid_shape.data());

        // Check if the subgrid is not empty
        if (subgrid_shape[0] != 0) {
            std::size_t flat_shape = std::accumulate(subgrid_shape.begin(),
                subgrid_shape.end(), 1, std::multiplies<std::size_t>());
            std::vector<double> subgrid_array(flat_shape);

            pineappl_subgrid_array(grid, b, order, channel, subgrid_array.data());

            // TODO: Unravel of the index & check against some reference values
            for (std::size_t value = 0; value < subgrid_array.size(); ++value) {
                if (subgrid_array[value] != 0) {
                    std::cout << std::right << std::setw(10) << b << std::setw(10)
                        << value << std::setw(16) << subgrid_array[value] << "\n";
                    break;
                }
            }
        }
    }

    pineappl_grid_delete(grid);
}
