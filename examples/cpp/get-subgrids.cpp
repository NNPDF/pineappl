#include <pineappl_capi.h>

#include <cassert>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

std::vector<std::size_t> unravel_index(std::size_t flat_index, const std::vector<std::size_t>& shape) {
    std::size_t ndim = shape.size();
    std::vector<std::size_t> coords(ndim);

    for (int i = ndim - 1; i >= 0; --i) {
        coords[i] = flat_index % shape[i];
        flat_index /= shape[i];
    }

    return coords;
}

template <typename T>
std::string vector_to_string(const std::vector<T>& coords) {
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

std::vector<double> get_subgrid_array(
    const pineappl_grid* grid,
    std::vector<std::size_t> subgrid_shape,
    std::size_t bin,
    std::size_t order,
    std::size_t channel
) {
    // Compute the length of the flattened shape by multiplying the entries
    std::size_t flat_shape = std::accumulate(subgrid_shape.begin(),
        subgrid_shape.end(), 1, std::multiplies<std::size_t>());

    // Extract the flattened subgrid values/weights
    std::vector<double> subgrid_array(flat_shape);
    pineappl_grid_subgrid_array(grid, bin, order, channel, subgrid_array.data());

    return subgrid_array;
}

std::vector<double> get_node_values(
    const pineappl_grid* grid,
    std::vector<std::size_t> subgrid_shape,
    std::size_t bin,
    std::size_t order,
    std::size_t channel
) {
    // Compute the length of the flattend nodes
    std::size_t nodes_size = std::accumulate(subgrid_shape.begin(),
        subgrid_shape.end(), 0.0);

    // Extract the values of the nodes as a flattened array
    std::vector<double> node_values(nodes_size);
    pineappl_grid_subgrid_node_values(grid, bin, order, channel, node_values.data());

    return node_values;
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
              << std::setw(6 * subgrid_dim) << "sg coordinates"
              << std::setw(12 * subgrid_dim) << "node values" << std::setw(16)
              << "weight value" << "\n";
    std::cout << std::right << std::setw(10) << "---" << std::setw(10) << "------"
              << std::setw(6 * subgrid_dim) << "--------------"
              << std::setw(12 * subgrid_dim) << "--------------------------------"
              << std::setw(16) << "------------" << "\n";

    for (std::size_t b = 0; b < n_bins; ++b) {
        // Extract the shape of the subgrids
        std::vector<std::size_t> subgrid_shape(subgrid_dim);
        pineappl_grid_subgrid_shape(grid, b, order, channel, subgrid_shape.data());

        // Check if the subgrid is not empty
        if (subgrid_shape[0] != 0) {
            std::vector<double> subgrid_array = get_subgrid_array(grid, subgrid_shape, b, order, channel);
            std::vector<double> node_values = get_node_values(grid, subgrid_shape, b, order, channel);

            for (std::size_t index = 0; index < subgrid_array.size(); ++index) {
                if (subgrid_array[index] != 0) {
                    // Unravel the index to recover the standard coordinates
                    std::vector<std::size_t> coords = unravel_index(index, subgrid_shape);

                    // Store the values of the nodes in a vector. The vector therefore
                    // contains as elements {Scale, x1, x2, ..., xn}
                    std::vector<double> node_values_index(coords.size());
                    std::size_t start_index = 0;
                    for (std::size_t nd = 0; nd < coords.size(); ++nd) {
                        if (nd != 0) { start_index += subgrid_shape[nd - 1]; }
                        node_values_index[nd] = node_values[start_index + coords[nd]];
                    }

                    std::cout << std::right << std::setw(10) << b << std::setw(10)
                              << index << std::setw(6 * subgrid_dim)
                              << vector_to_string(coords) << std::setw(12 * subgrid_dim)
                              << vector_to_string(node_values_index) << std::setw(16)
                              << subgrid_array[index] << "\n";

                    // Compare to some reference value.
                    if (b == 0 && index == 41020) {
                        // Check the unravelled index
                        assert(coords[0] == 16);
                        assert(coords[1] == 20);
                        assert(coords[2] == 20);

                        // Check the values of the node entries.
                        assert(node_values_index[0] == 5442.30542919352900); // PyAPI: `subgrid.node_values[0][16]`
                        assert(node_values_index[1] == 0.03052158400782890); // PyAPI: `subgrid.node_values[1][20]`
                        assert(node_values_index[2] == 0.03052158400782890); // PyAPI: `subgrid.node_values[2][20]`

                        // PyAPI: `grid.subgrid(0, 0, 0).to_array(subgrid.shape)[16][20][20]`
                        assert(subgrid_array[index] == -4.936156925096021e-07);
                    }
                    break;
                }
            }
        }
    }

    pineappl_grid_delete(grid);
}
