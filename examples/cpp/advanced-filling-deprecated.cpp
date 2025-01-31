#include <pineappl_capi.h>

#include <cassert>
#include <cstddef>
#include <string>
#include <vector>

int main() {
    auto* channels = pineappl_lumi_new();
    int32_t pids[] = { 2, -2, 4, -4 };
    double factors[] = { 1.0, 1.0 };
    pineappl_lumi_add(channels, 1, pids, factors);

    std::size_t channel_count = 1;

    std::vector<uint32_t> orders = {
        0, 2, 0, 0,
        1, 2, 0, 0,
        1, 2, 0, 1
    };

    std::vector<double> bins = {
        0.0,
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
        1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4
    };

    auto* keyval = pineappl_keyval_new();

    auto* grid = pineappl_grid_new(channels, orders.size() / 4, orders.data(), bins.size() - 1,
        bins.data(), keyval);

    pineappl_keyval_delete(keyval);
    pineappl_lumi_delete(channels);

    // arbitrary numbers
    double x1 = 0.001;
    double x2 = 0.02;
    double q2 = 10000.0;
    double yll = 1.3;
    std::size_t order = 0;
    std::size_t channel = 0;
    double weight = 1.23e-3;

    // fill a weight for a single order and channel
    pineappl_grid_fill(grid, x1, x2, q2, order, yll, channel, weight);

    // fill weights for a single order and all channels
    std::vector<double> weights(channel_count, weight);
    pineappl_grid_fill_all(grid, x1, x2, q2, order, yll, weights.data());

    // fill multiple events at once
    std::vector<double> weight_array(100, 1.3637e-4);
    std::vector<double> x1_array(weight_array.size(), x1);
    std::vector<double> x2_array(weight_array.size(), x2);
    std::vector<double> q2_array(weight_array.size(), q2);
    std::vector<std::size_t> order_array(weight_array.size(), 0);
    std::vector<double> yll_array(weight_array.size(), yll);
    std::vector<std::size_t> channel_array(weight_array.size(), 0);
    pineappl_grid_fill_array(grid, x1_array.data(), x2_array.data(), q2_array.data(),
        order_array.data(), yll_array.data(), channel_array.data(), weight_array.data(),
        weight_array.size());

    pineappl_grid_write(grid, "advanced-filling-deprecated.pineappl.lz4");

    // release memory
    pineappl_grid_delete(grid);
}
