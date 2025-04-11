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
        std::cout << "Usage: " << argv[0] << " [grid]\n";
    }

    // read the grid from a file
    auto* grid = pineappl_grid_read(filename.c_str());

    // How many convolutions are there?
    auto n_conv = pineappl_grid_convolutions_len(grid);

    // extract all channels
    auto* channels = pineappl_grid_channels(grid);

    // how many channels are there?
    auto channel_count = pineappl_channels_count(channels);

    for (std::size_t channel = 0; channel != channel_count; ++channel) {
        // print channel index
        std::cout << std::setw(4) << channel << ' ';

        // how many partonic combinations does this channel have?
        auto combinations = pineappl_channels_combinations(channels, channel);

        std::vector<double> factors(combinations);
        std::vector<int> pids(n_conv * combinations);

        // read out the channel with index given by `channel`, writing the particle identifiers into
        // `pids` and the corresponding factors into `factors`
        pineappl_channels_entry(channels, channel, n_conv, pids.data(), factors.data());

        for (std::size_t combination = 0; combination != combinations; ++combination) {
            auto factor = factors.at(combination);
            auto pida = pids.at(2 * combination + 0);
            auto pidb = pids.at(2 * combination + 1);

            if (combination != 0) {
                std::cout << " + ";
            }

            // print factor and particle ids
            std::cout << factor << " x (" << std::setw(3) << pida << ',' << std::setw(4) << pidb
                << ")";
        }

        std::cout << '\n';
    }

    // release memory
    pineappl_channels_delete(channels);
    pineappl_grid_delete(grid);
}
