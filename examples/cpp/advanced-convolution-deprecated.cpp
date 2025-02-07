#include <LHAPDF/PDF.h>
#include <pineappl_capi.h>

#include <cstddef>
#include <iomanip>
#include <ios>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

// An object to hold the state of the different PDFs
struct PDFState {
    LHAPDF::PDF* pdfs[2];
};

// The callable PDF function using the 1st PDF
double xfx1(int32_t id, double x, double q2, void* state) {
    auto* pdf = static_cast <PDFState*> (state)->pdfs[0];
    return pdf -> xfxQ2(id, x, q2);
};

// The callable PDF function using the 2nd PDF
double xfx2(int32_t id, double x, double q2, void* state) {
    auto* pdf = static_cast <PDFState*> (state)->pdfs[1];
    return pdf -> xfxQ2(id, x, q2);
};

// The callable alphasQ2 function using the 2nd PDF
double alphas(double q2, void* state) {
    auto* pdf = static_cast <PDFState*> (state)->pdfs[1];
    return pdf -> alphasQ2(q2);
};

int main(int argc, char* argv[]) {
    std::string filename = "drell-yan-rap-ll-deprecated.pineappl.lz4";
    std::string pdfset1 = "NNPDF31_nlo_as_0118_luxqed";
    std::string pdfset2 = "MSHT20qed_nnlo";

    switch (argc) {
    case 4:
        pdfset1 = argv[2];
        pdfset2 = argv[3];
        // fall through
    case 3:
        pdfset1 = argv[2];
        pdfset2 = argv[2];
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

    auto pdf1 = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(pdfset1, 0));
    auto pdf2 = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(pdfset2, 0));
    PDFState state = {pdf1.get(), pdf2.get()};

    // how many perturbative orders does the grid contain?
    std::size_t orders = pineappl_grid_order_count(grid);

    // how many bins does this grid have?
    std::size_t bins = pineappl_grid_bin_count(grid);

    auto* lumi = pineappl_grid_lumi(grid);

    // how many channels does the grid have?
    std::size_t channels = pineappl_lumi_count(lumi);

    pineappl_lumi_delete(lumi);

    // std::vector<bool> doesn't have `.data()` member
    std::unique_ptr<bool[]> order_mask(new bool[orders]());

    // allocate a vector holding the differential cross sections
    std::vector<double> dxsec1(bins);

    std::unique_ptr<bool[]> channel_mask(new bool[channels]());

    // use the variables to select the included orders and channels
    order_mask[0] = true;
    channel_mask[0] = true;

    // use these variables to perform scale variations
    double xir = 1.0;
    double xif = 1.0;

    // with this choice `order_mask` and `channel_mask` we extract the contribution of the first
    // perturbative order and first channel stored in the grid. The parameter `2212` declares that
    // `xfx1` and `alphas` are *proton* PDFs. If the grid contains cross sections of either a
    // proton-proton, proton-antiproton or antiproton-antiproton collision PineAPPL will perform the
    // necessary charge conjugations to yield the correct convolutions
    pineappl_grid_convolve_with_one(grid, 2212, xfx1, alphas, &state, order_mask.get(),
        channel_mask.get(), xir, xif, dxsec1.data());

    // how does the grid know which PDFs it must be convolved with? This is determined by the
    // metadata keys `initial_state_1` and `initial_state_2`, which are by default set to `2212`,
    // the PDG MC ID for the proton. Let's change the second value to an antiproton:
    pineappl_grid_set_key_value(grid, "initial_state_1", "2212");
    char* initial_state_1 = pineappl_grid_key_value(grid, "initial_state_1");
    assert( std::string(initial_state_1) == "2212" );
    // don't forget to deallocate!
    pineappl_string_delete(initial_state_1);

    pineappl_grid_set_key_value(grid, "initial_state_2", "-2212");
    char* initial_state_2 = pineappl_grid_key_value(grid, "initial_state_2");
    assert( std::string(initial_state_2) == "-2212" );
    // don't forget to deallocate!
    pineappl_string_delete(initial_state_2);

    std::vector<double> dxsec2(bins);

    // this calculates the corresponding proton-antiproton differential cross sections. Since the
    // grid itself is unchanged, this change effectively means that for the second PDF the charge
    // convolved PDFs are used
    pineappl_grid_convolve_with_one(grid, 2212, xfx1, alphas, &state, order_mask.get(),
        channel_mask.get(), xir, xif, dxsec2.data());

    // what if we have a collision where we actually need two PDFs? Let's simulate the collision of
    // protons with deuterons:
    pineappl_grid_set_key_value(grid, "initial_state_2", "1000010020"); // 1000010020 = deuteron
    char* new_initial_state_2 = pineappl_grid_key_value(grid, "initial_state_2");
    assert( std::string(new_initial_state_2) == "1000010020" );
    // don't forget to deallocate!
    pineappl_string_delete(new_initial_state_2);

    std::vector<double> dxsec3(bins);

    // For proton-deuteron collisions we can't easily relate the PDFs and have to actually pass two
    // different PDFs, each with their ID of the particle they represent. In this example, we use
    // the first PDF to compute alphasQ2.
    pineappl_grid_convolve_with_two(grid, 2212, xfx1, 1000010020, xfx2, alphas, &state,
        order_mask.get(), channel_mask.get(), xir, xif, dxsec3.data());

    std::vector<double> dxsec4(bins);

    // test with both masks set to `nullptr`. And in the following example, we use the second PDF
    // to compue the values of alphasQ2
    pineappl_grid_convolve_with_two(grid, 2212, xfx1, 1000010020, xfx2, alphas, &state, nullptr,
        nullptr, xir, xif, dxsec4.data());

    std::vector<double> normalizations(bins);

    // read out the bin normalizations, which is usually the size of each bin
    pineappl_grid_bin_normalizations(grid, normalizations.data());

    // print table header
    std::cout << "idx  p-p c#0 l#0 p-p~ c#0 l#  p-d c#0 l#0       p-d       dx\n"
                 "--- ------------ ----------- ------------- ------------ ------\n";

    for (std::size_t bin = 0; bin != bins; ++bin) {
        // print the bin index
        std::cout << std::setw(3) << bin << ' ';

        // print the result together with the normalization
        std::cout << std::scientific << dxsec1.at(bin) << ' ' << dxsec2.at(bin) << ' '
            << dxsec3.at(bin) << ' ' << dxsec4.at(bin) << ' ' << std::defaultfloat << std::setw(6)
            << normalizations.at(bin) << '\n';
    }

    pineappl_grid_delete(grid);
}
