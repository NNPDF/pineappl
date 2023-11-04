#include <pineappl_capi.h>

#include <cstddef>
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
        std::cout << "Usage: " << argv[0] << " [grid] [pdf]\n";
    }

    // read the grid from a file
    auto* grid = pineappl_grid_read(filename.c_str());

    // extract all channels
    auto* channels = pineappl_grid_lumi(grid);

    // how many channels are there?
    auto channel_count = pineappl_lumi_count(channels);

    for (std::size_t channel = 0; channel != channel_count; ++channel) {
        // print channel index
        std::cout << std::setw(4) << channel << ' ';

        // how many partonic combinations does this channel have?
        auto combinations = pineappl_lumi_combinations(channels, channel);

        std::vector<double> factors(combinations);
        std::vector<int> pids(2 * combinations);

        pineappl_lumi_entry(channels, channel, pids.data(), factors.data());

        for (std::size_t combination = 0; combination != combinations; ++combination) {
            auto factor = factors.at(combination);
            auto pida = pids.at(2 * combination + 0);
            auto pidb = pids.at(2 * combination + 1);

            if (combination != 0) {
                std::cout << " + ";
            }

            std::cout << factor << " x (" << std::setw(3) << pida << ',' << std::setw(4) << pidb
                << ")";
        }

        std::cout << '\n';
    }

    // release memory
    pineappl_lumi_delete(channels);
    pineappl_grid_delete(grid);
}
