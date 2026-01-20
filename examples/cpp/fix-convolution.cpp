#include <cassert>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <LHAPDF/LHAPDF.h>
#include <pineappl_capi.h>

void print_results(const std::vector<double>& r1, const std::vector<double>& r2) {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Bin | Original Grid   | Fixed Grid      | Rel. Diff." << std::endl;
    std::cout << "--- | --------------- | --------------- | -----------------" << std::endl;
    for (size_t i = 0; i < r1.size(); ++i) {
        double rel_diff = 0.0;
        if (r1[i] != 0.0) {
            rel_diff = std::abs((r1[i] - r2[i]) / r1[i]);
        }
        std::cout << std::setw(3) << i << " | "
                    << std::setw(15) << r1[i] << " | "
                    << std::setw(15) << r2[i] << " | "
                    << std::setw(17) << rel_diff << std::endl;
        assert(rel_diff < 1e-9);
    }
}

int main() {
    // Grid with 3 convolutions (2x pol PDF, 1x unpol FF)
    const std::string filename = "../../test-data/SIHP-PP-POLARIZED-STAR-NLO.pineappl.lz4";
    const std::string pol_pdf_set = "NNPDFpol11_100";
    const std::string ff_set = "MAPFF10NLOPIsum";

    LHAPDF::setVerbosity(0);
    auto pol_pdf = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(pol_pdf_set, 0));
    auto ff = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(ff_set, 0));

    pineappl_grid* grid = pineappl_grid_read(filename.c_str());
    assert(grid != nullptr);

    auto xfx = [](int32_t id, double x, double q2, void* pdf) {
        return static_cast <LHAPDF::PDF*> (pdf)->xfxQ2(id, x, q2);
    };
    auto alphas = [](double q2, void* pdf) {
        return static_cast <LHAPDF::PDF*> (pdf)->alphasQ2(q2);
    };

    // Convolve original grid
    size_t bins = pineappl_grid_bin_count(grid);
    std::vector<double> results_orig(bins);
    std::vector<LHAPDF::PDF*> pdfs_orig_vec = { pol_pdf.get(), pol_pdf.get(), ff.get() };
    void** pdfs_orig = reinterpret_cast<void**>(pdfs_orig_vec.data());
    std::vector<double> mu_scales = { 1.0, 1.0, 1.0 };
    pineappl_grid_convolve(
        grid,
        xfx,
        alphas,
        pdfs_orig,
        pol_pdf.get(),
        nullptr,
        nullptr,
        nullptr,
        1,
        mu_scales.data(),
        results_orig.data()
    );

    // Fix the third convolution (fragmentation function)
    pineappl_grid* grid_fixed = pineappl_grid_fix_convolution(grid, 2, xfx, ff.get(), 1.0);
    assert(grid_fixed != nullptr);

    // Convolve the new grid
    std::vector<double> results_fixed(bins);
    std::vector<LHAPDF::PDF*> pdfs_fixed_vec = { pol_pdf.get(), pol_pdf.get() };
    void** pdfs_fixed = reinterpret_cast<void**>(pdfs_fixed_vec.data());
    pineappl_grid_convolve(
        grid_fixed,
        xfx,
        alphas,
        pdfs_fixed,
        pol_pdf.get(),
        nullptr,
        nullptr,
        nullptr,
        1,
        mu_scales.data(),
        results_fixed.data()
    );

    print_results(results_orig, results_fixed);

    std::cout << "\nSuccess: results from original and fixed grid match." << std::endl;

    pineappl_grid_delete(grid);
    pineappl_grid_delete(grid_fixed);

    return 0;
}
