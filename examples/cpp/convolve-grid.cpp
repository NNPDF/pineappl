#include <LHAPDF/PDF.h>
#include <cstdint>
#include <pineappl_capi.h>

#include <cstddef>
#include <iomanip>
#include <ios>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
    std::string filename = "drell-yan-rap-ll.pineappl.lz4";
    std::string pdfset = "NNPDF31_nlo_as_0118_luxqed";

    switch (argc) {
    case 3:
        pdfset = argv[2];
        // fall through
    case 2:
        filename = argv[1];
    case 1:
        break;

    default:
        std::cout << "Usage: " << argv[0] << " [grid] [pdf]\n";
    }

    // disable LHAPDF banners to guarantee deterministic output
    LHAPDF::setVerbosity(0);

    // read the grid from a file
    auto* grid = pineappl_grid_read(filename.c_str());

    auto* pdf = LHAPDF::mkPDF(pdfset, 0);

    // define callables for the PDFs and alphas
    auto xfx = [](int32_t id, double x, double q2, void* pdf) {
        return static_cast <LHAPDF::PDF*> (pdf)->xfxQ2(id, x, q2);
    };
    auto alphas = [](double q2, void* pdf) {
        return static_cast <LHAPDF::PDF*> (pdf)->alphasQ2(q2);
    };

    // how many bins does this grid have?
    std::size_t bins = pineappl_grid_bin_count(grid);

    // how many dimensions does each bin have?
    std::size_t dims = pineappl_grid_bin_dimensions(grid);

    // allocate a vector holding the left and right bin limits for each dimension
    std::vector<double> bin_limits(2 * bins * dims);

    for (std::size_t dim = 0; dim != dims; ++dim) {
        pineappl_grid_bin_limits_left(grid, dim, &bin_limits.at((2 * dim + 0) * bins));
        pineappl_grid_bin_limits_right(grid, dim, &bin_limits.at((2 * dim + 1) * bins));
    }

    // allocate a vector holding the differential cross sections
    std::vector<double> dxsec(bins);

    auto order_mask = nullptr;
    auto channel_mask = nullptr;
    double xir = 1.0;
    double xif = 1.0;

    // perform the convolution of `grid` with the PDFs given as `xfx` and the strong coupling in
    // `alphas` and the extra parameter `pdf`, which is passed to `xfx` and `alphas` as the last
    // parameter. The integer `2212` is the PDG MC id for a proton and signals and `xfx` is the PDF
    // of a proton. In this case we assume that both initial state hadrons' PDFs can derived from
    // that of a proton. If this isn't the case, for instance for a proton-lead collision, both PDFs
    // must be given separately and the function `pineappl_grid_convolve_with_two` must be used.
    // The parameters `order_mask` and `channel_mask` can be used to select specific orders and
    // channels, respectively. Using `xir` and `xif` the renormalization and factorization scales
    // can be varied around its central values, respectively.
    std::vector<double> mu_scales = { xir, xif, 1.0 };
    using LambdaType = double(*)(int32_t, double, double, void *);
    LambdaType xfxs[] = { xfx, xfx };
    pineappl_grid_convolve(grid, xfxs, alphas, pdf, order_mask, channel_mask, nullptr, 1,
        mu_scales.data(), dxsec.data());

    std::vector<double> normalizations(bins);

    // read out the bin normalizations, which is usually the size of each bin
    pineappl_grid_bin_normalizations(grid, normalizations.data());

    // print table header
    std::cout << "idx";
    for (std::size_t dim = 0; dim != dims; ++dim) {
        std::cout << "   left  right";
    }
    std::cout << "    dsig/dx     dx\n";
    std::cout << "---";
    for (std::size_t dim = 0; dim != dims; ++dim) {
        std::cout << " ----dim #" << dim << "---";
    }
    std::cout << " ------------ ------\n";

    for (std::size_t bin = 0; bin != bins; ++bin) {
        // print the bin index
        std::cout << std::setw(3) << bin << ' ';

        for (std::size_t dim = 0; dim != dims; ++dim) {
            double left_limit = bin_limits.at((2 * dim + 0) * bins + bin);
            double right_limit = bin_limits.at((2 * dim + 1) * bins + bin);

            // print the left and right bin limit for each dimension
            std::cout << std::setw(6) << left_limit << ' ' << std::setw(6) << right_limit << ' ';
        }

        // print the result together with the normalization
        std::cout << std::scientific << dxsec.at(bin) << std::defaultfloat << ' '
            << std::setw(6) << normalizations.at(bin) << '\n';
    }

    pineappl_grid_delete(grid);
}
