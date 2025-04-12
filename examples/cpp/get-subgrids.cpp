#include <pineappl_capi.h>

#include <cassert>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>


std::vector<std::size_t> unravel_index(std::size_t flat_index, const std::vector<std::size_t>& shape) {
    std::size_t ndim = shape.size();
    std::vector<std::size_t> coords(ndim);

    for (int i = ndim - 1; i >= 0; --i) {
        coords[i] = flat_index % shape[i];
        flat_index /= shape[i];
    }

    return coords;
}

std::string coords_to_string(const std::vector<std::size_t>& coords) {
    std::ostringstream osstream;

    osstream << "(";
    for (std::size_t i = 0; i < coords.size(); ++i) {
        osstream << coords[i];
        if (i != coords.size() - 1) {
            osstream << ", ";
        }
    }
    osstream << ")";

    return osstream.str();
}


int main() {
    std::string filename = "drell-yan-rap-ll.pineappl.lz4";

    // read the grid from a file
    auto* grid = pineappl_grid_read(filename.c_str());

    // Determine the number of bins and the index of order and channel
    std::size_t n_bins = pineappl_grid_bin_count(grid);
    std::size_t order = 0;
    std::size_t channel = 0;

    // Get the dimension of the subgrids
    std::size_t subgrid_dim = pineappl_grid_kinematics_len(grid);

    std::cout << std::right << std::setw(10) << "bin" << std::setw(10) << "sg idx"
        << std::setw(6 * subgrid_dim) << "sg coordinates" << std::setw(16)
        << "sg value" << "\n";
    std::cout << std::right << std::setw(10) << "---" << std::setw(10) << "------"
        << std::setw(6 * subgrid_dim) << "--------------" << std::setw(16)
        << "------------" << "\n";

    for (std::size_t b = 0; b < n_bins; ++b) {
        std::vector<std::size_t> subgrid_shape(subgrid_dim);
        pineappl_grid_subgrid_shape(grid, b, order, channel, subgrid_shape.data());

        // Check if the subgrid is not empty
        if (subgrid_shape[0] != 0) {
            std::size_t flat_shape = std::accumulate(subgrid_shape.begin(),
                subgrid_shape.end(), 1, std::multiplies<std::size_t>());
            std::vector<double> subgrid_array(flat_shape);

            pineappl_grid_subgrid_array(grid, b, order, channel, subgrid_array.data());

            for (std::size_t index = 0; index < subgrid_array.size(); ++index) {
                if (subgrid_array[index] != 0) {
                    // Unravel the index to recover the standard coordinates
                    std::vector<std::size_t> coords = unravel_index(index, subgrid_shape);

                    std::cout << std::right << std::setw(10) << b << std::setw(10)
                        << index << std::setw(6 * subgrid_dim) << coords_to_string(coords)
                        << std::setw(16) << subgrid_array[index] << "\n";

                    // Compare to some reference value
                    if (b==0 && index==41020) {
                        assert(subgrid_array[index] == -4.936156925096015e-07);
                    }
                    break;
                }
            }
        }
    }

    pineappl_grid_delete(grid);
}
