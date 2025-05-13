#include <LHAPDF/PDF.h>
#include <pineappl_capi.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cassert>
#include <cstddef>
#include <string>
#include <vector>

// NOTE: Uses the scale of the Grid as the starting scale such that we can use an IDENTITY EKO.
double FAC0 = 6456.44;

std::vector<std::size_t> unravel_index(std::size_t flat_index, const std::vector<std::size_t>& shape) {
    std::size_t ndim = shape.size();
    std::vector<std::size_t> coords(ndim);

    for (int i = ndim - 1; i >= 0; --i) {
        coords[i] = flat_index % shape[i];
        flat_index /= shape[i];
    }

    return coords;
}

std::vector<double> generate_fake_ekos(
    std::vector<int> pids_in,
    std::vector<double> x_in,
    std::vector<int> pids_out,
    std::vector<double> x_out
) {
    std::size_t flat_len = x_out.size() * x_in.size() * pids_out.size() * pids_in.size();
    std::vector<double> ops(flat_len);

    // NOTE: The EKO has to have as shape: (pids_in, x_in, pids_out, x_out)
    std::vector<std::size_t> shape = {pids_in.size(), x_in.size(), pids_out.size(), x_out.size()};
    for (std::size_t i = 0; i != flat_len; i++) {
        std::vector<std::size_t> coords = unravel_index(i, shape);

        double delta_ik = (coords[0] == coords[2]) ? 1.0 : 0.0;
        double delta_jl = (coords[1] == coords[3]) ? 1.0 : 0.0;

        ops[i] = delta_ik * delta_jl;
    }

    return ops;
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
    std::string filename = "LHCB_WP_7TEV_opt.pineappl.lz4";

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

    // Get the number of convolutions
    std::size_t n_convs = pineappl_grid_convolutions_len(grid);

    // Fill the vector of unique convolution types. If the EKOs required for the Grid
    // are the same, then it suffices to only pass ONE single EKO.
    std::vector<pineappl_conv_type> conv_types;
    for (std::size_t i = 0; i != n_convs; i++) {
        pineappl_conv_type conv = pineappl_grid_conv_type(grid, i);
        if (std::find(conv_types.begin(), conv_types.end(), conv) == conv_types.end()) {
            conv_types.push_back(conv);
        }
    }

    // Get the shape of the evolve info objects
    std::vector<std::size_t> evinfo_shape(5);
    std::vector<uint8_t> max_orders = {2, 3};
    pineappl_grid_evolve_info_shape(grid, max_orders.data(), evinfo_shape.data());

    // Get the values of the evolve info parameters. These contain, for example, the
    // information on the `x`-grid and `PID` used to interpolate the Grid.
    // NOTE: These are used to construct the Evolution Operator
    std::vector<double> fac1(evinfo_shape[0]);
    std::vector<double> frg1(evinfo_shape[1]);
    std::vector<int> pids_in(evinfo_shape[2]);
    std::vector<double> x_in(evinfo_shape[3]);
    std::vector<double> ren1(evinfo_shape[4]);
    pineappl_grid_evolve_info(grid, max_orders.data(), fac1.data(),
        frg1.data(), pids_in.data(), x_in.data(), ren1.data());

    // ------------------ Construct the Operator Info ------------------
    // The Operator Info is a vector with length `N_conv * N_Q2_slices` whose
    // elements are `OperatorInfo` objects.
    std::vector<pineappl_operator_info> opinfo_slices(conv_types.size() * fac1.size());
    for (std::size_t i = 0; i != conv_types.size(); i++) {
        for (std::size_t j = 0; j != fac1.size(); j++) {
            pineappl_operator_info opinfo = {
                FAC0, // fac0
                fac1[j], // fac1
                pid_basis,
                conv_types[i],
            };
            opinfo_slices[i * fac1.size() + j] = opinfo;
        }
    }

    // ------------------ Construct the Evolution Operator ------------------
    // Choose a different PID basis for the FK table
    // std::vector<int> pids_out = {-6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 21, 22};
    std::vector<int> pids_out = pids_in;

    // The Evolution Operator is a vector with length `N_conv * N_Q2_slices * Σ product(OP shape)`
    std::vector<double> op_slices;
    std::size_t flat_len = x_in.size() * x_in.size() * pids_in.size() * pids_out.size();
    for (std::size_t _i = 0; _i != conv_types.size(); _i++) {
        for (std::size_t j = 0; j != fac1.size(); j++) {
            std::vector<double> eko = generate_fake_ekos(pids_in, x_in, pids_out, x_in);
            for (std::size_t k = 0; k != flat_len; k++) {
                op_slices.push_back(eko[k]);
            }
        }
    }

    // Construct the values of alphas table
    std::vector<double> alphas_table;
    for (double q2 : ren1) {
        double alpha = alphas(q2, pdf.get());
        alphas_table.push_back(alpha);
    }

    std::vector<double> xi = {1.0, 1.0, 1.0};
    // NOTE: The EKO has to have as shape: (pids_in, x_in, pids_out, x_out)
    std::vector<std::size_t> tensor_shape = {pids_in.size(), x_in.size(), pids_out.size(), x_in.size()};

    // NOTE: The arguments of `pineappl_grid_evolve` must follow the following orders:
    //     - `grid`: PineAPPL Grid
    //     - `op_info`: operator info
    //     - `max_orders`: max orders to apply the evolution
    //     - `operators`: evolution operator
    //     - `x_in`: x-grid of the Grid
    //     - `x_out`: x-grid of the FK table
    //     - `pids_in`: PIDs basis representation of the Grid
    //     - `pids_out`: PIDs basis representation of the FK table
    //     - `eko_shape`: shape of the evolution operators
    //     - `xi`: scale variation
    //     - `ren1`: values of the renormalization scales
    //     - `alphas_table`: values of alphas for each renormalization scales
    pineappl_fk_table* fktable = pineappl_grid_evolve(grid, opinfo_slices.data(),
        max_orders.data(), op_slices.data(), x_in.data(),
        x_in.data(), pids_in.data(), pids_out.data(),
        tensor_shape.data(), xi.data(), ren1.data(), alphas_table.data());

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
    pineappl_fk_table_convolve(fktable, xfx, pdf_states, nullptr,
        nullptr, dxsec_fktable.data());

    // Print the results
    print_results(dxsec_grid, dxsec_fktable);

    pineappl_fktable_write(fktable, "evolved-grid.pineappl.lz4");

    pineappl_grid_delete(grid);
    pineappl_fk_table_delete(fktable);
}
