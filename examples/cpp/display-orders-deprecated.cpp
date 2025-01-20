#include <pineappl_capi.h>

#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
    std::string filename = "drell-yan-rap-ll.pineappl.lz4";

    switch (argc) {
    case 2:
        filename = argv[1];
    case 1:
        break;

    default:
        std::cout << "Usage: " << argv[0] << " [grid]\n";
    }

    // read the grid from a file
    auto* grid = pineappl_grid_read(filename.c_str());

    // how many perturbative orders does this grid contain?
    std::size_t orders = pineappl_grid_order_count(grid);

    std::vector<std::uint32_t> order_params(4 * orders);

    // read out all exponents of the perturbative orders in the grid
    pineappl_grid_order_params(grid, order_params.data());

    for (std::size_t order = 0; order != orders; ++order) {
        std::cout << std::setw(4) << order << ' ';

        // exponent of the strong coupling
        std::uint32_t exp_as = order_params.at(4 * order + 0);
        // exponent of the electromagnetic/electroweak coupling
        std::uint32_t exp_al = order_params.at(4 * order + 1);
        // exponent of the renormalization log
        std::uint32_t exp_lr = order_params.at(4 * order + 2);
        // exponent of the factorization log
        std::uint32_t exp_lf = order_params.at(4 * order + 3);

        std::cout << "O(as^" << exp_as << " a^" << exp_al << " lr^" << exp_lr << " lf^" << exp_lf
            << ")\n";
    }

    // release memory
    pineappl_grid_delete(grid);
}
