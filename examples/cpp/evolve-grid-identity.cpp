#include <LHAPDF/PDF.h>
#include <pineappl_capi.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cassert>
#include <cstddef>
#include <string>
#include <vector>
#include <numeric>

// NOTE: Uses the scale of the Grid as the starting scale such that we can use an IDENTITY EKO.
double FAC0 = 6456.44;

/** @brief This struct can contain arbitrary parameters that need to be passed to Evolution
 * Operator Callback (`generated_fake_ekos`).
 */
struct OperatorParams {
    std::vector<pineappl_conv_type> conv_types;
};

std::vector<std::size_t> unravel_index(std::size_t flat_index, const std::vector<std::size_t>& shape) {
    std::size_t ndim = shape.size();
    std::vector<std::size_t> coords(ndim);

    for (int i = ndim - 1; i >= 0; --i) {
        coords[i] = flat_index % shape[i];
        flat_index /= shape[i];
    }

    return coords;
}

extern "C" void generate_fake_ekos(
    std::size_t op_index,
    double /*fac1*/,
    const int* /*pids_in*/,
    const double* /*x_in*/,
    const int* /*pids_out*/,
    const double* /*x_out*/,
    const std::size_t* eko_shape,
    double* eko_buffer,
    void* params_state
) {
    // select the type of convolution based on the Operator index
    OperatorParams* op_params = static_cast<OperatorParams*>(params_state);
    pineappl_conv_type _ = op_params->conv_types[op_index];

    // NOTE: This has to work because the Evolution Operator is always 4D
    std::vector<std::size_t> shape(eko_shape, eko_shape + 4);
    // Compute the length of the flattened shape by multiplying the entries
    std::size_t flat_len = std::accumulate(shape.begin(),
        shape.end(), 1, std::multiplies<std::size_t>());
    for (std::size_t i = 0; i != flat_len; i++) {
        std::vector<std::size_t> coords = unravel_index(i, shape);

        double delta_ik = (coords[0] == coords[2]) ? 1.0 : 0.0;
        double delta_jl = (coords[1] == coords[3]) ? 1.0 : 0.0;

        eko_buffer[i] = delta_ik * delta_jl;
    }
}

void print_results(std::vector<double> dxsec_grid, std::vector<double> dxsec_fktable) {
    const int idx_width = 6;
    const int num_width = 15;
    const int dif_width = 15;

    // Print headers
    std::cout << std::setw(idx_width) << "Bin"
              << std::setw(num_width) << "Grid"
              << std::setw(num_width) << "FkTable"
              << std::setw(dif_width) << "reldiff" << std::endl;

    // Print dashed lines
    std::cout << std::setw(idx_width) << std::string(idx_width - 2, '-')
              << std::setw(num_width) << std::string(num_width - 2, '-')
              << std::setw(num_width) << std::string(num_width - 2, '-')
              << std::setw(dif_width) << std::string(dif_width - 2, '-') << std::endl;

    // Print the data
    std::cout << std::scientific << std::setprecision(6);
    for (size_t i = 0; i < dxsec_grid.size(); ++i) {
        double reldiff = (dxsec_fktable[i] - dxsec_grid[i]) / dxsec_grid[i];
        std::cout << std::setw(idx_width) << i
                  << std::setw(num_width) << dxsec_grid[i]
                  << std::setw(num_width) << dxsec_fktable[i]
                  << std::setw(dif_width) << reldiff
                  << std::endl;
    }
}

int main() {
    // TODO: How to get a Grid that can be evolved??
    std::string filename = "../../test-data/LHCB_WP_7TEV_opt.pineappl.lz4";

    // disable LHAPDF banners to guarantee deterministic output
    LHAPDF::setVerbosity(0);
    std::string pdfset = "NNPDF31_nlo_as_0118_luxqed";
    auto pdf = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(pdfset, 0));

    auto xfx = [](int32_t id, double x, double q2, void* pdf) {
        return static_cast <LHAPDF::PDF*> (pdf)->xfxQ2(id, x, q2);
    };
    auto alphas = [](double q2, void* pdf) {
        return static_cast <LHAPDF::PDF*> (pdf)->alphasQ2(q2);
    };

    std::vector<LHAPDF::PDF*> pdfs = {pdf.get(), pdf.get()};
    void** pdf_states = reinterpret_cast<void**>(pdfs.data());

    // read the grid from a file
    auto* grid = pineappl_grid_read(filename.c_str());

    // Get the PID basis representation
    pineappl_pid_basis pid_basis = pineappl_grid_pid_basis(grid);
    assert(pid_basis == PINEAPPL_PID_BASIS_PDG);

    // Get the number of convolutions and their types
    std::size_t n_convs = pineappl_grid_convolutions_len(grid);
    std::vector<pineappl_conv_type> conv_types(n_convs);
    pineappl_grid_conv_types(grid, conv_types.data());

    // Fill the vector of unique convolution types. If the Operators required for the Grid
    // are the same, then it suffices to only pass ONE single Operator.
    std::vector<pineappl_conv_type> unique_convs;
    for (std::size_t i = 0; i != n_convs; i++) {
        pineappl_conv_type conv = conv_types[i];
        if (std::find(unique_convs.begin(), unique_convs.end(), conv) == unique_convs.end()) {
            unique_convs.push_back(conv);
        }
    }

    // Get the shape of the evolve info objects
    std::vector<std::size_t> evinfo_shape(5);
    // NOTE: The argument of `pineappl_grid_evolve_info_shape` must follow the following orders:
    //     - `grid`: PineAPPL Grid
    //     - `order_mask`: array of booleans to mask the order(s) to apply the Evolution to,
    //                     `nullptr` selects all the orders
    //     - `evinfo_shape`: placeholder to store the shape of the Evolution Operator
    pineappl_grid_evolve_info_shape(grid, nullptr, evinfo_shape.data());

    // Get the values of the evolve info parameters. These contain, for example, the
    // information on the `x`-grid and `PID` used to interpolate the Grid.
    // NOTE: These are used to construct the Evolution Operator
    std::vector<double> fac1(evinfo_shape[0]);
    std::vector<double> frg1(evinfo_shape[1]);
    std::vector<int> pids_in(evinfo_shape[2]);
    std::vector<double> x_in(evinfo_shape[3]);
    std::vector<double> ren1(evinfo_shape[4]);
    // NOTE: The argument of `pineappl_grid_evolve_info` must follow the following orders:
    //     - `grid`: PineAPPL Grid
    //     - `order_mask`: array of booleans to mask the order(s) to apply the Evolution to,
    //                     `nullptr` selects all the orders
    // The rest of the arguments are placeholders to store data
    pineappl_grid_evolve_info(grid, nullptr, fac1.data(),
        frg1.data(), pids_in.data(), x_in.data(), ren1.data());

    // ------------------ Construct the Operator Info ------------------
    // The Operator Info is a vector with length `N_conv * N_Q2_slices` whose
    // elements are `OperatorInfo` objects.
    std::vector<pineappl_conv_type> convtypes(unique_convs.size());
    std::vector<pineappl_operator_info> opinfo_slices(unique_convs.size() * fac1.size());
    for (std::size_t i = 0; i != unique_convs.size(); i++) {
        for (std::size_t j = 0; j != fac1.size(); j++) {
            pineappl_operator_info opinfo = {
                FAC0, // fac0
                fac1[j], // fac1
                pid_basis,
                unique_convs[i],
            };
            opinfo_slices[i * fac1.size() + j] = opinfo;
        }
        convtypes[i] = unique_convs[i];
    }

    // ------------------ Construct the Evolution Operator ------------------
    // Choose a different PID basis for the FK table
    // std::vector<int> pids_out = {-6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 21, 22};
    std::vector<int> pids_out = pids_in;

    // Construct the values of alphas table
    std::vector<double> alphas_table;
    for (double q2 : ren1) {
        double alpha = alphas(q2, pdf.get());
        alphas_table.push_back(alpha);
    }

    // Construct the Parameters that will get passed to the Callback
    OperatorParams* op_params = new OperatorParams;
    op_params->conv_types = convtypes;
    void* params = static_cast<void*>(op_params);

    std::vector<double> xi = {1.0, 1.0, 1.0};
    // NOTE: The EKO has to have as shape: (pids_in, x_in, pids_out, x_out)
    std::vector<std::size_t> tensor_shape = {pids_in.size(), x_in.size(), pids_out.size(), x_in.size()};

    // NOTE: The arguments of `pineappl_grid_evolve` must follow the following orders:
    //     - `grid`: PineAPPL Grid
    //     - `nb_slices`: the number of convolution(s)/Evolution Operator(s) required
    //     - `slices`: callback that returns the evolution operator(s) in slices
    //     - `operator_info`: operator info
    //     - `pids_in`: PIDs basis representation of the Grid
    //     - `x_in`: x-grid of the Grid
    //     - `pids_out`: PIDs basis representation of the FK table
    //     - `x_out`: x-grid of the FK table
    //     - `state`: parameters that get passed to `operator`
    //     - `order_mask`: array of booleans to mask the order(s) to apply the Evolution to,
    //                     `nullptr` selects all the orders
    //     - `xi`: scale variation
    //     - `ren1`: values of the renormalization scales
    //     - `alphas_table`: values of alphas for each renormalization scales
    //     - `eko_shape`: shape of the evolution operators
    pineappl_grid* fktable = pineappl_grid_evolve(
        grid,                 // `grid`
        unique_convs.size(),  // `nb_slices`
        generate_fake_ekos,   // `slices`
        opinfo_slices.data(), // `operator_info`
        pids_in.data(),       // `pids_in`
        x_in.data(),          // `x_in`
        pids_out.data(),      // `pids_out`
        x_in.data(),          // `x_out`
        tensor_shape.data(),  // `eko_shape`
        params,               // `state`
        nullptr,              // `order_mask`
        xi.data(),            // `xi`
        ren1.data(),          // `ren1`
        alphas_table.data()   // `alphas_table`
    );

    // ------------------ Compare Grid & FK after convolution ------------------
    // how many bins does this grid have?
    std::size_t bins = pineappl_grid_bin_count(grid);

    // [ convolve the Grid ]
    std::vector<double> mu_scales = { 1.0, 1.0, 1.0 };
    std::vector<double> dxsec_grid(bins);
    pineappl_grid_convolve(grid, xfx, alphas, pdf_states, pdf.get(),
        nullptr, nullptr, nullptr, 1,
        mu_scales.data(), dxsec_grid.data());

    // [ convolve the FK Table ]
    std::vector<double> dxsec_fktable(bins);
    auto as_one = [](double /*q2*/, void* /*pdf*/) { return 1.0; };
    pineappl_grid_convolve(fktable, xfx, as_one, pdf_states, nullptr,
        nullptr, nullptr, nullptr, 1,
        mu_scales.data(), dxsec_fktable.data());

    // Print the results
    print_results(dxsec_grid, dxsec_fktable);

    pineappl_grid_write(fktable, "evolved-grid-identity.pineappl.lz4");

    pineappl_grid_delete(grid);
    pineappl_grid_delete(fktable);
}
