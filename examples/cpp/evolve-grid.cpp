#include <LHAPDF/PDF.h>
#include <pineappl_capi.h>

#include <fstream>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cassert>
#include <cstddef>
#include <string>
#include <vector>

// NOTE: Uses the scale of the Grid as the starting scale such that we can use an IDENTITY EKO.
double FAC0 = 2.7224999999999997;

// x-grid nodes definition for the `in` and `out`
std::vector<double> XGRID = {
    2.00000000000000e-07,
    3.03430476586795e-07,
    4.60350147489639e-07,
    6.98420853070036e-07,
    1.05960949591010e-06,
    1.60758549847081e-06,
    2.43894329289168e-06,
    3.70022720698550e-06,
    5.61375771693015e-06,
    8.51680667757335e-06,
    1.29210156907473e-05,
    1.96025050023917e-05,
    2.97384953722449e-05,
    4.51143839496404e-05,
    6.84374491896790e-05,
    1.03811729865769e-04,
    1.57456056008414e-04,
    2.38787829185619e-04,
    3.62054496381397e-04,
    5.48779532367080e-04,
    8.31406883648814e-04,
    1.25867971442728e-03,
    1.90346340228674e-03,
    2.87386758128175e-03,
    4.32850063882081e-03,
    6.49620619463380e-03,
    9.69915957404340e-03,
    1.43750685810901e-02,
    2.10891866837872e-02,
    3.05215840078289e-02,
    4.34149174170227e-02,
    6.04800287544474e-02,
    8.22812212620489e-02,
    1.09143757463307e-01,
    1.41120806444403e-01,
    1.78025660425694e-01,
    2.19504126500389e-01,
    2.65113704158282e-01,
    3.14387400769276e-01,
    3.66875318648224e-01,
    4.22166775358965e-01,
    4.79898902961025e-01,
    5.39757233788045e-01,
    6.01472197967335e-01,
    6.64813948247382e-01,
    7.29586844241431e-01,
    7.95624252292276e-01,
    8.62783932390611e-01,
    9.30944080871754e-01,
    1.00000000000000e+00
};

// Particle PIDs for both `in` and `out`
std::vector<int> PIDS = {- 22 , -6 , -5 , -4 , -3 , -2 , -1 , 21 , 1 , 2 , 3 , 4 , 5 , 6};

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
    const int* /*pids_in*/,
    const double* /*x_in*/,
    const int* /*pids_out*/,
    const double* /*x_out*/,
    double* eko_buffer,
    void* params_state,
    pineappl_conv_type /*conv_type*/,
    double /*fac1*/,
    std::size_t /*pids_in_len*/,
    std::size_t /*x_in_len*/,
    std::size_t /*pids_out_len*/,
    std::size_t /*x_out_len*/
) {
    // Check to get the μ0 from the PDF
    const double _ = static_cast<LHAPDF::PDF*> (params_state)->q2Min();

    std::ifstream input_file("../../test-data/EKO_LHCB_WP_7TEV.txt");
    double weight_value;

    std::size_t count = 0;
    while (input_file >> weight_value) {
        eko_buffer[count++] = weight_value;
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
    std::vector<int> _pids_in(evinfo_shape[2]);
    std::vector<double> _x_in(evinfo_shape[3]);
    std::vector<double> ren1(evinfo_shape[4]);
    // NOTE: The argument of `pineappl_grid_evolve_info` must follow the following orders:
    //     - `grid`: PineAPPL Grid
    //     - `order_mask`: array of booleans to mask the order(s) to apply the Evolution to,
    //                     `nullptr` selects all the orders
    // The rest of the arguments are placeholders to store data
    pineappl_grid_evolve_info(grid, nullptr, fac1.data(),
        frg1.data(), _pids_in.data(), _x_in.data(), ren1.data());

    // ------------------ Construct the Operator Info ------------------
    // The Operator Info is a vector with length `N_conv * N_Q2_slices` whose
    // elements are `OperatorInfo` objects.
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
    }

    // ------------------ Construct the Evolution Operator ------------------
    // Define the same PIDs for both `in` and `out` according to the EKO
    std::vector<int> pids_in = PIDS;
    std::vector<int> pids_out = PIDS;

    // Construct the values of alphas table
    std::vector<double> alphas_table;
    for (double q2 : ren1) {
        double alpha = alphas(q2, pdf.get());
        alphas_table.push_back(alpha);
    }

    std::vector<double> xi = {1.0, 1.0, 1.0};
    // NOTE: The EKO has to have as shape: (pids_in, x_in, pids_out, x_out)
    std::vector<std::size_t> tensor_shape = {pids_in.size(), XGRID.size(), pids_out.size(), XGRID.size()};

    // NOTE: The arguments of `pineappl_grid_evolve` must follow the following orders:
    //     - `grid`: PineAPPL Grid
    //     - `op_info`: operator info
    //     - `operator`: callback that returns an evolution operator
    //     - `order_mask`: array of booleans to mask the order(s) to apply the Evolution to,
    //                     `nullptr` selects all the orders
    //     - `params_state`: parameters that get passed to `operator`
    //     - `nb_convolutions`: the number of convolutions/Evolution Operators required
    //     - `x_in`: x-grid of the Grid
    //     - `x_out`: x-grid of the FK table
    //     - `pids_in`: PIDs basis representation of the Grid
    //     - `pids_out`: PIDs basis representation of the FK table
    //     - `eko_shape`: shape of the evolution operators
    //     - `xi`: scale variation
    //     - `ren1`: values of the renormalization scales
    //     - `alphas_table`: values of alphas for each renormalization scales
    pineappl_grid* fktable = pineappl_grid_evolve(grid, opinfo_slices.data(),
        generate_fake_ekos, nullptr, pdf.get(),
        unique_convs.size(), XGRID.data(), XGRID.data(),
        pids_in.data(), pids_out.data(), tensor_shape.data(),
        xi.data(), ren1.data(), alphas_table.data());

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

    // write the unoptimised FK table into disk
    pineappl_grid_write(fktable, "evolved-grid.pineappl.lz4");

    // optimise the FK table and then write into disk
    pineappl_fktable_optimize(fktable, PINEAPPL_FK_ASSUMPTIONS_NF3_SYM);
    pineappl_grid_write(fktable, "evolved-grid-optimised.pineappl.lz4");

    pineappl_grid_delete(grid);
    pineappl_grid_delete(fktable);
}
