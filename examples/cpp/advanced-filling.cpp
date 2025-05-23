#include <LHAPDF/PDF.h>
#include <pineappl_capi.h>

#include <cassert>
#include <cstddef>
#include <iomanip>
#include <memory>
#include <vector>

int main() {
    // Construct the channel object based on the nb of convolutions
    std::size_t nb_convolutions = 2;
    std::size_t nb_channels = 2;
    auto* channels = pineappl_channels_new(nb_convolutions);
    int32_t pids[] = { 2, -2, 4, -4 };
    double factors[] = { 1.0, 1.0 };
    pineappl_channels_add(channels, nb_channels, pids, factors);

    std::size_t channel_count = 1;

    std::vector<uint8_t> orders = {
        0, 2, 0, 0, 0,
        1, 2, 0, 0, 0,
        1, 2, 0, 1, 0
    };

    std::vector<double> bins = {
        0.0,
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
        1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4
    };

    // ---
    // Construct the objects that are needed to fill the Grid

    // First we define the types of convolutions required by the involved initial-/final-state
    // hadrons. Then we add the corresponding PID of each of the hadrons, and finally define the
    // Basis onto which the partons are mapped.
    pineappl_pid_basis pid_basis = PINEAPPL_PID_BASIS_EVOL;
    pineappl_conv convs[] = {
        { PINEAPPL_CONV_TYPE_UNPOL_PDF, 2212 },
        { PINEAPPL_CONV_TYPE_UNPOL_PDF, 2212 },
    };

    // Define the kinematics required for this process. In the following example we have ONE
    // single scale and two momentum fractions (corresponding to the two initial-state hadrons).
    // The format of the kinematics is: { type, value }.
    pineappl_kinematics scales = { PINEAPPL_KINEMATICS_SCALE, 0 };
    pineappl_kinematics x1k = { PINEAPPL_KINEMATICS_X, 0 };
    pineappl_kinematics x2k = { PINEAPPL_KINEMATICS_X, 1 };
    pineappl_kinematics kinematics[3] = { scales, x1k, x2k };

    // Define the specificities of the interpolations for each of the kinematic variables.
    pineappl_reweight_meth scales_reweight = PINEAPPL_REWEIGHT_METH_NO_REWEIGHT; // Reweighting method
    pineappl_reweight_meth moment_reweight = PINEAPPL_REWEIGHT_METH_APPL_GRID_X;
    pineappl_map scales_mapping = PINEAPPL_MAP_APPL_GRID_H0; // Mapping method
    pineappl_map moment_mapping = PINEAPPL_MAP_APPL_GRID_F2;
    pineappl_interp_meth interpolation_meth = PINEAPPL_INTERP_METH_LAGRANGE;
    pineappl_interp interpolations[3] = {
        { 1e2, 1e8, 40, 3, scales_reweight, scales_mapping, interpolation_meth },  // Interpolation fo `scales`
        { 2e-7, 1.0, 50, 3, moment_reweight, moment_mapping, interpolation_meth }, // Interpolation fo `x1`
        { 2e-7, 1.0, 50, 3, moment_reweight, moment_mapping, interpolation_meth }, // Interpolation fo `x2`
    };

    // Define the unphysical scale objecs
    pineappl_scale_func_form scale_mu = { PINEAPPL_SCALE_FUNC_FORM_SCALE, 0 };
    pineappl_scale_func_form no_scale_mu = { PINEAPPL_SCALE_FUNC_FORM_NO_SCALE, 0 }; // Here `.scale=0` is dummy value
    pineappl_scale_func_form mu_scales[3] = { scale_mu, scale_mu, no_scale_mu };

    // ---
    // Create the grid using the previously set information about orders, bins and channels

    // create a new grid with the previously defined channels, 3 perturbative orders defined by the
    // exponents in `orders`, 24 bins given as the 25 limits in `bins` and potential extra
    // parameters in `keyval`.
    auto* grid = pineappl_grid_new2(bins.size() - 1, bins.data(), orders.size() / 5, orders.data(),
        channels, pid_basis, convs, 3, interpolations, kinematics, mu_scales);

    // now we no longer need `channels`
    pineappl_channels_delete(channels);

    // arbitrary numbers
    double x1 = 0.001;
    double x2 = 0.02;
    double q2 = 10000.0;
    double yll = 1.3;
    std::size_t order = 0;
    std::size_t channel = 0;
    double weight = 1.23e-3;

    // Values of the kinematic variables
    std::vector<double> ntuples = { q2, x1, x2 };

    // fill a weight for a single order and channel
    pineappl_grid_fill2(grid, order, yll, channel, ntuples.data(), weight);

    // fill weights for a single order and all channels
    std::vector<double> weights(channel_count, weight);
    pineappl_grid_fill_all2(grid, order, yll, ntuples.data(), weights.data());

    // fill multiple events at once
    std::vector<double> weight_array(100, 1.3637e-4);
    std::vector<double> x1_array(weight_array.size(), x1);
    std::vector<double> x2_array(weight_array.size(), x2);
    std::vector<double> q2_array(weight_array.size(), q2);
    // Construct the `ntuples` kinematic array
    std::vector<double> ntuples_array;
    // Zipping loop: Add elements of q2_array, x1_array, x2_array to `ntuples_array`
    for (size_t i = 0; i < weight_array.size(); ++i) {
        ntuples_array.push_back(q2_array[i]);
        ntuples_array.push_back(x1_array[i]);
        ntuples_array.push_back(x2_array[i]);
    }
    // Construct the array of order, observable, and channel
    std::vector<std::size_t> order_array(weight_array.size(), 0);
    std::vector<double> yll_array(weight_array.size(), yll);
    std::vector<std::size_t> channel_array(weight_array.size(), 0);

    pineappl_grid_fill_array2(grid, order_array.data(), yll_array.data(),
        ntuples_array.data(), channel_array.data(), weight_array.data(), weight_array.size());

    //-------------------- Check Convolution ----------------------//
    std::string pdfset = "NNPDF31_nlo_as_0118_luxqed";
    // disable LHAPDF banners to guarantee deterministic output
    LHAPDF::setVerbosity(0);
    auto pdf = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(pdfset, 0));

    // define callables for the PDFs and alphas
    auto xfx = [](int32_t id, double x, double q2, void* pdf) {
        return static_cast <LHAPDF::PDF*> (pdf)->xfxQ2(id, x, q2);
    };
    auto alphas = [](double q2, void* pdf) {
        return static_cast <LHAPDF::PDF*> (pdf)->alphasQ2(q2);
    };

    auto order_mask = nullptr;
    auto channel_mask = nullptr;
    std::vector<double> mmu_scales = { 1.0, 1.0, 1.0 };

    std::vector<LHAPDF::PDF*> pdfs = {pdf.get(), pdf.get()};
    void** pdf_states = reinterpret_cast<void**>(pdfs.data());

    // allocate a vector holding the differential cross sections
    std::vector<double> dxsec(bins.size() - 1);
    pineappl_grid_convolve(grid, xfx, alphas, pdf_states, pdf.get(), order_mask, channel_mask,
        nullptr, 1, mmu_scales.data(), dxsec.data());

    // Print table header
    std::cout << std::setw(10) << "bin left"
              << std::setw(12) << "bin right"
              << std::setw(15) << "dsig/dx" << std::endl;
    std::cout << std::string(37, '-') << std::endl;

    // Loop through bins and print results
    for (size_t i = 0; i < dxsec.size(); ++i) {
        std::cout << std::setw(10) << bins[i]
                  << std::setw(12) << bins[i + 1]
                  << std::setw(15) << std::scientific << std::setprecision(3) << dxsec[i]
                  << std::endl;
    }
    //-----------------------------------------------------------------------//

    // Check that the Grid contains an empty subgrid at (b, o, c)=(0, 0, 0)
    auto subgrid_dim = pineappl_grid_kinematics_len(grid);
    std::vector<std::size_t> subgrid_shape(subgrid_dim);
    std::vector<std::size_t> zero_vector(subgrid_dim, 0);
    pineappl_grid_subgrid_shape(grid, 0, 0, 0, subgrid_shape.data());
    assert(subgrid_shape == zero_vector);

    // Check that a call to an empty subgrid does not panic
    std::vector<double> subgrid_array(0);
    pineappl_grid_subgrid_array(grid, 0, 0, 0, subgrid_array.data());

    pineappl_grid_write(grid, "advanced-filling.pineappl.lz4");

    // release memory
    pineappl_grid_delete(grid);
}
