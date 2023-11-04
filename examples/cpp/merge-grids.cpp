#include <pineappl_capi.h>

#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    std::string filename1 = "drell-yan-rap-ll.pineappl.lz4";
    std::string filename2 = "drell-yan-rap-ll.pineappl.lz4";

    switch (argc) {
    case 3:
        filename1 = argv[2];
        // fall through
    case 2:
        filename1 = argv[1];
    case 1:
        break;

    default:
        std::cout << "Usage: " << argv[0] << " [grid1] [grid2]\n";
    }

    // read the grids from file
    auto* grid1 = pineappl_grid_read(filename1.c_str());
    auto* grid2 = pineappl_grid_read(filename2.c_str());

    // clone `grid2`
    auto* clone = pineappl_grid_clone(grid2);

    // merge `grid2` into `grid1`. This automatically deletes `grid2`, which is why we clone it in
    // the previous line
    pineappl_grid_merge_and_delete(grid1, grid2);

    // this also works and doesn't do anything
    pineappl_grid_merge_and_delete(grid1, nullptr);

    // write out the merged grid
    pineappl_grid_write(grid1, "merged-grids.pineappl.lz4");

    // release memory - `grid2` was already deleted by `pineappl_grid_merge_and_delete`
    pineappl_grid_delete(clone);
    pineappl_grid_delete(grid1);
}

