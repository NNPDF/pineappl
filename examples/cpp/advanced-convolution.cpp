#include <LHAPDF/PDF.h>
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

    // define callables for the PDFs and alphas, the grid collides two protons
    auto xfx1 = [](int32_t id, double x, double q2, void* pdf) {
        return static_cast <LHAPDF::PDF*> (pdf)->xfxQ2(id, x, q2);
    };
    auto xfx2 = [](int32_t id, double x, double q2, void* pdf) {
        return static_cast <LHAPDF::PDF*> (pdf)->xfxQ2(id, x, q2);
    };
    auto alphas = [](double q2, void* pdf) {
        return static_cast <LHAPDF::PDF*> (pdf)->alphasQ2(q2);
    };

    // how many perturbative orders does the grid contain?
    std::size_t orders = pineappl_grid_order_count(grid);

    // how many bins does this grid have?
    std::size_t bins = pineappl_grid_bin_count(grid);

    // how many channels does the grid have?
    auto* channels = pineappl_grid_channels(grid);
    std::size_t channels_length = pineappl_channels_count(channels);

    pineappl_channels_delete(channels);

    // std::vector<bool> doesn't have `.data()` member
    std::unique_ptr<bool[]> order_mask(new bool[orders]());

    // allocate a vector holding the differential cross sections
    std::vector<double> dxsec1(bins);

    std::unique_ptr<bool[]> channel_mask(new bool[channels_length]());

    // use the variables to select the included orders and channels
    order_mask[0] = true;
    channel_mask[0] = true;

    // use these variables to perform scale variations
    double xir = 1.0;
    double xif = 1.0;

    // with this choice `order_mask` and `channel_mask` we extract the contribution of the first
    // perturbative order and first channel stored in the grid. If the grid contains cross sections
    // of either a proton-proton, proton-antiproton or antiproton-antiproton collision PineAPPL will
    // perform the necessary charge conjugations to yield the correct convolutions.
    // In the case where the convolution requires two different PDFs, it suffices to modify `xfx1`
    // and `xfx2`.
    std::vector<double> mu_scales = { xir, xif, 1.0 };
    using LambdaType = double(*)(int32_t, double, double, void *);
    LambdaType xfxs[] = { xfx1, xfx2 };

    pineappl_grid_convolve(grid, xfxs, alphas, pdf, order_mask.get(), channel_mask.get(), nullptr, 1,
        mu_scales.data(), dxsec1.data());

    // test with both masks set to `nullptr`
    std::vector<double> dxsec2(bins);

    pineappl_grid_convolve(grid, xfxs, alphas, pdf, nullptr, nullptr, nullptr, 1,
        mu_scales.data(), dxsec2.data());

    std::vector<double> normalizations(bins);

    // read out the bin normalizations, which is usually the size of each bin
    pineappl_grid_bin_normalizations(grid, normalizations.data());

    // print table header
    std::cout << "idx  p-p c#0 l#0      p-d       dx\n"
                 "--- ------------ ------------ ------\n";

    for (std::size_t bin = 0; bin != bins; ++bin) {
        // print the bin index
        std::cout << std::setw(3) << bin << ' ';

        // print the result together with the normalization
        std::cout << std::scientific << dxsec1.at(bin) << ' ' << dxsec2.at(bin) << ' '
            << std::defaultfloat << std::setw(6) << normalizations.at(bin) << '\n';
    }

    pineappl_grid_delete(grid);
}
