#include <cmath>
#include <cstdlib>
#include <string>

#include <fastnlotk/fastNLOLHAPDF.h>
#include <pineappl_capi.h>

std::unique_ptr<LHAPDF::PDF> pdf;

extern "C" double xfx(int32_t pdg_id, double x, double q2, void* state)
{
    return static_cast <LHAPDF::PDF*> (state)->xfxQ2(pdg_id, x, q2);
}

extern "C" double alphas(double q2, void* state)
{
    return static_cast <LHAPDF::PDF*> (state)->alphasQ2(q2);
}

pineappl_grid* convert_coeff_add_fix(
    fastNLOCoeffAddFix* table,
    std::vector<double> const& bin_limits,
    uint32_t alpha
) {
    std::vector<uint32_t> order_params = { table->GetNpow(), alpha, 0, 0 };

    auto* lumi = pineappl_lumi_new();
    auto const& pdf = table->GetPDFCoeff();

    for (auto const& pdf_entries : pdf)
    {
        std::vector<int32_t> combinations;
        std::vector<double> factors;

        for (auto const& entry : pdf_entries)
        {
            combinations.push_back(entry.first);
            combinations.push_back(entry.second);
            factors.push_back(1.0);
        }

        pineappl_lumi_add(lumi, factors.size(), combinations.data(), factors.data());
    }

    auto* key_vals = pineappl_keyval_new();
    auto* pgrid = pineappl_grid_new(lumi, 1, order_params.data(), bin_limits.size() - 1,
        bin_limits.data(), key_vals);
    pineappl_keyval_delete(key_vals);
    pineappl_lumi_delete(lumi);

    //auto descriptions = table->GetContributionDescription();

    //for (auto const& desc : descriptions)
    //{
    //    std::cout << "-- " << desc << '\n';
    //}

    std::size_t n_obs_bin = table->GetNObsBin();
    //std::size_t n_scalevar = table->GetNScalevar();
    //std::size_t n_scale_node = table->GetNScaleNode();
    //std::size_t n_scales = table->GetNScales();
    std::size_t n_subproc = table->GetNSubproc();
    std::size_t total_scalevars = table->GetTotalScalevars();
    std::size_t total_scalenodes = table->GetTotalScalenodes();

    //std::cout << table->GetScaleDescription(0) << '\n';
    //auto const& descr = table->GetScaleDescr();
    //std::cout << descr.at(0).at(0) << std::endl;
    //std::cout << table->GetNevt() << '\n';

    for (std::size_t obs = 0; obs != n_obs_bin; ++obs)
    {
        std::size_t n_xmax = table->GetNxmax(obs);

        auto const& x1_values = table->GetXNodes1(obs);
        auto const& x2_values = table->GetXNodes2(obs);

        for (std::size_t subproc = 0; subproc != n_subproc; ++subproc)
        {
            double const factor = table->GetNevt(obs, subproc);

            for (std::size_t j = 0; j != total_scalevars; ++j)
            {
                auto q2_values = table->GetScaleNodes(obs, j);

                // the values are the unsquared q values, correct that
                for (auto& value : q2_values)
                {
                    value *= value;
                }

                auto* subgrid = pineappl_subgrid_new(q2_values.size(), q2_values.data(),
                    x1_values.size(), x1_values.data(), x2_values.size(), x2_values.data());

                // TODO: figure out what the general case is supposed to be
                assert( j == 0 );

                //std::cout << table->GetScaleFactor(j) << '\n';

                bool non_zero_subgrid = false;

                for (std::size_t k = 0; k != total_scalenodes; ++k)
                {
                    std::vector<double> slice(x1_values.size() * x2_values.size());
                    bool non_zero = false;

                    for (std::size_t ix = 0; ix != n_xmax; ++ix)
                    {
                        auto const value = table->GetSigmaTilde(obs, j, k, ix, subproc);

                        if (value != 0.0)
                        {
                            non_zero = true;

                            auto const ix1 = ix % x1_values.size();
                            auto const ix2 = ix / x1_values.size();
                            assert( table->GetXIndex(obs, ix1, ix2) == ix );

                            slice.at(x2_values.size() * ix1 + ix2) = value / factor;
                        }
                    }

                    if (non_zero)
                    {
                        non_zero_subgrid = true;
                        pineappl_subgrid_import_q2_slice(subgrid, k, slice.data());
                    }
                }

                if (non_zero_subgrid)
                {
                    pineappl_grid_replace_and_delete(pgrid, subgrid, 0, obs, subproc);
                }
                else
                {
                    pineappl_subgrid_delete(subgrid);
                }
            }
        }
    }

    return pgrid;
}

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        return EXIT_FAILURE;
    }

    std::string in(argv[1]);
    std::string out(argv[2]);

    // TODO: read this from an argument
    uint32_t alpha = 0;

    LHAPDF::setVerbosity(0);
    pdf.reset(LHAPDF::mkPDF("NNPDF31_nlo_as_0118_luxqed", 1));

    auto file = fastNLOLHAPDF(in, "NNPDF31_nlo_as_0118_luxqed");

    // TODO: for the time being only one-dimensional distributions are supported
    auto dim = file.GetNumDiffBin();
    assert( dim == 1 );

    auto const& bin_limits = file.GetDim0BinBounds();
    std::vector<double> pine_bin_limits(bin_limits.size() + 1);
    pine_bin_limits.at(0) = bin_limits.at(0).first;

    for (std::size_t i = 0; i != bin_limits.size(); ++i)
    {
        pine_bin_limits.at(i + 1) = bin_limits.at(i).second;

        // TODO: for the time being we only support consecutive bin limits
        assert( bin_limits.at(i).first == pine_bin_limits.at(i) );
    }

    auto const id_lo = file.ContrId(fastNLO::kFixedOrder, fastNLO::kLeading);
    auto const id_nlo = file.ContrId(fastNLO::kFixedOrder, fastNLO::kNextToLeading);
    auto const id_nnlo = file.ContrId(fastNLO::kFixedOrder, fastNLO::kNextToNextToLeading);

    std::vector<int> ids;

    if (id_lo >= 0)
    {
        ids.push_back(id_lo);
    }

    if (id_nlo >= 0)
    {
        ids.push_back(id_nlo);
    }

    if (id_nnlo >= 0)
    {
        ids.push_back(id_nnlo);
    }

    std::vector<pineappl_grid*> grids;

    for (auto const id : ids)
    {
        auto coeff_table = file.GetCoeffTable(id);
        assert( coeff_table != nullptr );

        auto converted = dynamic_cast <fastNLOCoeffAddFix*> (coeff_table);

        if (converted != nullptr)
        {
            grids.push_back(convert_coeff_add_fix(converted, pine_bin_limits, alpha));
        }
        else
        {
            auto converted = dynamic_cast <fastNLOCoeffAddFlex*> (coeff_table);

            if (converted != nullptr)
            {
                // TODO: NYI
                assert( false );
            }
            else
            {
                // TODO: NYI
                assert( false );
            }
        }
    }

    for (std::size_t i = 1; i < grids.size(); ++i)
    {
        pineappl_grid_merge_and_delete(grids.at(0), grids.at(i));
    }

    pineappl_grid_scale_by_order(grids.at(0), 0.5 / std::acos(-1.0), 1.0, 1.0, 1.0, 1.0);
    //pineappl_grid_optimize(grids.at(0));

    auto const& results = file.GetCrossSection();
    std::vector<double> other_results(results.size());

    pineappl_grid_convolute(
        grids.at(0),
        xfx,
        xfx,
        alphas,
        pdf.get(),
        nullptr,
        nullptr,
        1.0,
        1.0,
        other_results.data()
    );

    std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
    std::cout.precision(16);

    bool different = false;

    for (std::size_t i = 0; i != results.size(); ++i)
    {
        auto const one = results.at(i);
        auto const two = other_results.at(i);

        // catches the case where both results are zero
        if (one == two)
        {
            continue;
        }

        if (std::fabs(two / one - 1.0) > 1e-10)
        {
            std::cout << ">>> fastNLO: " << one << " PineAPPL: " << two << " fN/P: " << (one / two)
                << " P/fN: " << (two / one) << '\n';
            different = true;
        }
        else
        {
            std::cout << ">>> Success!\n";
        }
    }

    pineappl_grid_write(grids.at(0), out.c_str());
    pineappl_grid_delete(grids.at(0));
}
