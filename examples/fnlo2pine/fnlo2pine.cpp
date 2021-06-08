#include <cmath>
#include <cstdlib>
#include <numeric>
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

int32_t convert_to_pdg_id(int id)
{
    if ((id >= -6) && (id <= 6))
    {
        return (id == 0) ? 21 : id;
    }
    else
    {
        assert( false );
    }
}

pineappl_grid* convert_coeff_add_fix(
    fastNLOCoeffAddFix* table,
    fastNLOPDFLinearCombinations const& comb,
    std::size_t bins,
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
            combinations.push_back(convert_to_pdg_id(entry.first));
            combinations.push_back(convert_to_pdg_id(entry.second));
            factors.push_back(1.0);
        }

        pineappl_lumi_add(lumi, factors.size(), combinations.data(), factors.data());
    }

    // if there is no luminosity definition, we have to become creative
    if (pdf.empty())
    {
        std::vector<double> xfx1(13);
        std::vector<double> xfx2(13);

        std::vector<std::vector<int32_t>> combinations(table->GetNSubproc());
        std::vector<std::vector<double>> factors(table->GetNSubproc());

        for (int a = 0; a != 13; ++a)
        {
            xfx1.at(a) = 1.0;

            for (int b = 0; b != 13; ++b)
            {
                xfx2.at(b) = 1.0;

                auto const& lumi = comb.CalcPDFLinearCombination(table, xfx1, xfx2, false);

                assert( lumi.size() == table->GetNSubproc() );

                for (std::size_t i = 0; i != lumi.size(); ++i)
                {
                    if (lumi.at(i) == 0.0)
                    {
                        continue;
                    }

                    auto const ap = convert_to_pdg_id(a - 6);
                    auto const bp = convert_to_pdg_id(b - 6);

                    combinations.at(i).push_back(ap);
                    combinations.at(i).push_back(bp);
                    factors.at(i).push_back(lumi.at(i));
                }

                xfx2.at(b) = 0.0;
            }

            xfx1.at(a) = 0.0;
        }

        for (std::size_t i = 0; i != combinations.size(); ++i)
        {
            pineappl_lumi_add(lumi, factors.at(i).size(), combinations.at(i).data(),
                factors.at(i).data());
        }
    }

    std::vector<double> bin_limits(bins + 1);
    std::iota(bin_limits.begin(), bin_limits.end(), 0.0);
    auto* key_vals = pineappl_keyval_new();
    auto* pgrid = pineappl_grid_new(lumi, 1, order_params.data(), bins, bin_limits.data(),
        key_vals);
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

        // TODO: is this the correct assumption?
        auto const x2_values = (table->GetNxtot2(0) == -1) ? x1_values : table->GetXNodes2(obs);

        for (std::size_t subproc = 0; subproc != n_subproc; ++subproc)
        {
            //if ((pdf[subproc][0].first != 0) || (pdf[subproc][0].second != 0))
            //{
            //    continue;
            //}

            double const factor = table->GetNevt(obs, subproc);

            for (std::size_t j = 0; j != total_scalevars; ++j)
            {
                // TODO: for the time being we only extract the central scale result
                if (table->GetScaleFactor(j) != 1.0)
                {
                    continue;
                }

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

                    std::size_t ix1 = 0;
                    std::size_t ix2 = 0;

                    for (std::size_t ix = 0; ix != static_cast <std::size_t> (table->GetNxmax(obs)); ++ix)
                    {
                        assert( table->GetXIndex(obs, ix1, ix2) == ix );

                        auto const value = table->GetSigmaTilde(obs, j, k, ix, subproc);

                        if (value != 0.0)
                        {
                            non_zero = true;
                            slice.at(x2_values.size() * ix2 + ix1) = value / factor
                                * x1_values.at(ix1) * x2_values.at(ix2);
                        }

                        ++ix1;

                        switch (table->GetNPDFDim())
                        {
                        case 2:
                            if (ix1 == x1_values.size())
                            {
                                ix1 = 0;
                                ++ix2;
                            }
                            break;

                        case 1:
                            if (ix1 > ix2)
                            {
                                ix1 = 0;
                                ++ix2;
                            }
                            break;

                        default:
                            // TODO: NYI
                            assert( false );
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
    pdf.reset(LHAPDF::mkPDF("NNPDF31_nlo_as_0118_luxqed", 0));

    auto file = fastNLOLHAPDF(in, "NNPDF31_nlo_as_0118_luxqed");

    //file.SelectProcesses({ std::make_pair(0, 0) });
    //file.ActivateContribution(fastNLO::kFixedOrder, fastNLO::kNextToLeading, false);
    //file.ActivateContribution(fastNLO::kFixedOrder, fastNLO::kNextToNextToLeading, false);

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

    auto const& normalizations = file.GetBinSize();
    std::size_t const bins = normalizations.size();

    std::vector<pineappl_grid*> grids;

    for (auto const id : ids)
    {
        auto coeff_table = file.GetCoeffTable(id);
        assert( coeff_table != nullptr );

        auto converted = dynamic_cast <fastNLOCoeffAddFix*> (coeff_table);

        if (converted != nullptr)
        {
            grids.push_back(convert_coeff_add_fix(converted, file, bins, alpha));
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
    pineappl_grid_optimize(grids.at(0));

    auto const dimensions = file.GetNumDiffBin();
    std::vector<double> limits(2 * dimensions * bins);

    for (std::size_t i = 0; i != bins; ++i)
    {
        for (std::size_t j = 0; j != dimensions; ++j)
        {
            auto const& pair = file.GetObsBinDimBounds(i, j);

            limits.at(2 * (i * dimensions + j) + 0) = pair.first;
            limits.at(2 * (i * dimensions + j) + 1) = pair.second;
        }
    }

    pineappl_grid_set_remapper(grids.at(0), dimensions, normalizations.data(), limits.data());

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
        auto const two = other_results.at(i) * normalizations.at(i);

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

    if (different)
    {
        return 1;
    }
}
