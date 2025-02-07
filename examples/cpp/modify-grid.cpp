#include <pineappl_capi.h>

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

    // read the grid from file
    auto* grid = pineappl_grid_read(filename.c_str());

    // how many bins does our grid have?
    auto bins = pineappl_grid_bin_count(grid);

    // 1. merge all bins into a single one. This adds the cross sections of all bins together, and
    // the new single bin has the left limit of the first old bin and the right limit of the last
    // old bin
    pineappl_grid_merge_bins(grid, 0, bins);

    // 2a. scale a grid with a number. This multiplies all subgrids with the given number
    pineappl_grid_scale(grid, 1.0);

    // 2b. scale the grid depending on its bins. Bins for which a factor isn't given are not
    // rescaled. If more factors are given than there are bins, these additional factors are ignored
    std::vector<double> factors(bins, 1.0);
    pineappl_grid_scale_by_bin(grid, factors.size(), factors.data());

    // 2c. scale the grid depending on its orders. If a subgrid is quadratic in the strong coupling,
    // it will be scaled by the square of the value `alphas` below. This is useful to convert
    // between differently defined cross sections. In Madgraph5, for instance, the cross sections
    // are filled into the grid factorized in terms of power of gs^2, but PineAPPL requires them to
    // be factorized in terms of alphas. The difference is a factor 4pi, which `alphas` would be set
    // to
    double alphas = 1.0;
    double alpha  = 1.0;
    double logxir = 1.0;
    double logxif = 1.0;
    double global = 1.0;
    pineappl_grid_scale_by_order(grid, alphas, alpha, logxir, logxif, global);

    // 3a. split channels. A grid with multiple initial states in a single channel will then have
    // multiple channels with one initial state
    pineappl_grid_split_channels(grid);

    // 3b. undo the previous operation, detecting equal subgrids by allowing them to differ by up to
    // 64 ULPS
    pineappl_grid_dedup_channels(grid, 64);

    // 4. optimize grid selectively. The following example removes all perturbative orders whose
    // subgrids are empty
    pineappl_grid_optimize_using(grid, PINEAPPL_GOF_STRIP_EMPTY_ORDERS);

    // 5. set a remapper. This function is important if one wants to generate multi-dimensional
    // differential distributions, which first must be generated one-dimensional, because
    // `pineappl_grid_fill` only supports one variable. Afterwards the multi dimensionality can be
    // restored by setting the multi-dimensional limits with this call:
    std::vector<double> normalizations = { 1.0 };
    std::vector<double> limits = {
        60.0, 120.0, // dimension #0 bin #0
         0.0,   2.4, // dimension #1 bin #0
                     // dimension #0 bin #1 - in this example we only have one bin
                     // dimension #1 bin #1
                     // ...
    };
    std::size_t dimensions = limits.size() / 2;
    pineappl_grid_set_remapper(grid, dimensions, normalizations.data(), limits.data());

    // write out the modified grid
    pineappl_grid_write(grid, "modified-grid.pineappl.lz4");

    // release memory
    pineappl_grid_delete(grid);
}
