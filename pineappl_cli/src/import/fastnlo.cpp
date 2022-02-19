#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <numeric>
#include <string>

#include <LHAPDF/PDF.h>
#include <pineappl_cli/src/import/fastnlo.hpp>
#include <fastnlotk/fastNLOLHAPDF.h>
#include <pineappl_capi.h>

std::unique_ptr<LHAPDF::PDF> pdf;

extern "C" double xfx(int32_t pdg_id, double x, double q2, void* state)
{
    return static_cast <LHAPDF::PDF*> (state)->xfxQ2(pdg_id, x, q2);
}

extern "C" double xfx2(int32_t pdg_id, double x, double, void*)
{
    assert( pdg_id == 11 );
    return x;
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

void create_lumi(
    fastNLOCoeffAddBase* table,
    fastNLOPDFLinearCombinations const& comb,
    pineappl_lumi* lumi
) {
    auto const& pdf = table->GetPDFCoeff();

    // TODO: set this to the right value if there's only one PDF
    int32_t const lepton_id = (table->GetNPDF() == 2) ? 0 : 11;

    // if there's a (non-empty) PDF coefficient vector reconstruct the luminosity function; the
    // advantage is that we preserve the order of the lumi entries in the PineAPPL grid
    for (auto const& pdf_entries : pdf)
    {
        std::vector<int32_t> combinations;
        std::vector<double> factors;

        for (auto const& entry : pdf_entries)
        {
            combinations.push_back(convert_to_pdg_id(entry.first));

            if (lepton_id == 0)
            {
                combinations.push_back(convert_to_pdg_id(entry.second));
            }
            else
            {
                combinations.push_back(11);
            }

            factors.push_back(1.0);
        }

        pineappl_lumi_add(lumi, factors.size(), combinations.data(), factors.data());
    }

    // if there's an empty PDF coefficient vector, we reconstruct the lumi function
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

                assert( static_cast <int> (lumi.size()) == table->GetNSubproc() );

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
}

pineappl_grid* convert_coeff_add_fix(
    fastNLOCoeffAddFix* table,
    fastNLOPDFLinearCombinations const& comb,
    std::size_t bins,
    uint32_t alpha
) {
    std::vector<uint32_t> order_params = { static_cast <uint32_t> (table->GetNpow()), alpha, 0, 0 };

    auto* lumi = pineappl_lumi_new();
    create_lumi(table, comb, lumi);

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
        auto const& x1_values = table->GetXNodes1(obs);

        // TODO: is this the correct assumption?
        auto const x2_values = (table->GetNxtot2(0) == -1) ? x1_values : table->GetXNodes2(obs);

        for (std::size_t subproc = 0; subproc != n_subproc; ++subproc)
        {
            double const factor = table->GetNevt(obs, subproc);

            for (std::size_t j = 0; j != total_scalevars; ++j)
            {
                // TODO: for the time being we only extract the central scale result
                if (table->GetScaleFactor(j) != 1.0)
                {
                    continue;
                }

                auto q2_values = table->GetScaleNodes(obs, j);
                std::vector<double> mu2_values;

                // the values are the unsquared q values, correct that
                for (auto& value : q2_values)
                {
                    value *= value;
                    mu2_values.push_back(value);
                    mu2_values.push_back(value);
                }

                auto* subgrid = pineappl_subgrid_new2(mu2_values.size() / 2, mu2_values.data(),
                    x1_values.size(), x1_values.data(), x2_values.size(), x2_values.data());

                // TODO: figure out what the general case is supposed to be
                assert( j == 0 );

                //std::cout << table->GetScaleFactor(j) << '\n';

                bool non_zero_subgrid = false;

                for (std::size_t mu2_slice = 0; mu2_slice != total_scalenodes; ++mu2_slice)
                {
                    std::vector<double> slice(x1_values.size() * x2_values.size());
                    bool non_zero = false;

                    std::size_t ix1 = 0;
                    std::size_t ix2 = 0;

                    for (int ix = 0; ix != table->GetNxmax(obs); ++ix)
                    {
                        assert( table->GetXIndex(obs, ix1, ix2) == ix );

                        auto const value = table->GetSigmaTilde(obs, j, mu2_slice, ix, subproc);

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
                        pineappl_subgrid_import_mu2_slice(subgrid, mu2_slice, slice.data());
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

pineappl_grid* convert_coeff_add_flex(
    fastNLOCoeffAddFlex* table,
    fastNLOPDFLinearCombinations const& comb,
    fastNLO::EScaleFunctionalForm mur_ff,
    fastNLO::EScaleFunctionalForm muf_ff,
    std::size_t bins,
    uint32_t alpha,
    int ipub_units
) {
    std::vector<uint32_t> order_params = { static_cast <uint32_t> (table->GetNpow()), alpha, 0, 0 };

    auto* lumi = pineappl_lumi_new();
    create_lumi(table, comb, lumi);

    std::vector<double> bin_limits(bins + 1);
    std::iota(bin_limits.begin(), bin_limits.end(), 0.0);
    auto* key_vals = pineappl_keyval_new();

    // flexible grids always have a hadron in initial state 1 ...
    pineappl_keyval_set_string(key_vals, "initial_state_1",
        std::to_string(table->GetPDFPDG(0)).c_str());
    // and something else at 2
    pineappl_keyval_set_string(key_vals, "initial_state_2", "11");

    auto* pgrid = pineappl_grid_new(lumi, 1, order_params.data(), bins, bin_limits.data(),
        key_vals);
    pineappl_keyval_delete(key_vals);
    pineappl_lumi_delete(lumi);

    std::size_t n_obs_bin = table->GetNObsBin();
    std::size_t n_subproc = table->GetNSubproc();

    auto const& sigma_tildes = table->GetSigmaTildes();

    auto const rescale = std::pow(0.1, table->GetIXsectUnits() - ipub_units);

    for (std::size_t obs = 0; obs != n_obs_bin; ++obs)
    {
        auto const& scale_nodes1 = table->GetScaleNodes1(obs);
        auto const& scale_nodes2 = table->GetScaleNodes2(obs);
        auto const& x1_values = table->GetXNodes1(obs);

        std::vector<double> mur2_values;

        switch (mur_ff)
        {
        case fastNLO::kScale1:
            for (auto const s1 : scale_nodes1)
            {
                for (std::size_t i = 0; i != scale_nodes2.size(); ++i)
                {
                    mur2_values.push_back(s1 * s1);
                }
            }
            break;

        case fastNLO::kScale2:
            for (std::size_t i = 0; i != scale_nodes1.size(); ++i)
            {
                for (auto const s2 : scale_nodes2)
                {
                    mur2_values.push_back(s2 * s2);
                }
            }
            break;

        case fastNLO::kQuadraticSum:
            for (auto const s1 : scale_nodes1)
            {
                for (auto const s2 : scale_nodes2)
                {
                    mur2_values.push_back(s1 * s1 + s2 * s2);
                }
            }
            break;

        case fastNLO::kQuadraticMean:
            for (auto const s1 : scale_nodes1)
            {
                for (auto const s2 : scale_nodes2)
                {
                    mur2_values.push_back(0.5 * (s1 * s1 + s2 * s2));
                }
            }
            break;

        default:
            // TODO: NYI
            assert( false );
        }

        std::vector<double> muf2_values;

        switch (muf_ff)
        {
        case fastNLO::kScale1:
            for (auto const s1 : scale_nodes1)
            {
                for (std::size_t i = 0; i != scale_nodes2.size(); ++i)
                {
                    muf2_values.push_back(s1 * s1);
                }
            }
            break;

        case fastNLO::kScale2:
            for (std::size_t i = 0; i != scale_nodes1.size(); ++i)
            {
                for (auto const s2 : scale_nodes2)
                {
                    muf2_values.push_back(s2 * s2);
                }
            }
            break;

        case fastNLO::kQuadraticSum:
            for (auto const s1 : scale_nodes1)
            {
                for (auto const s2 : scale_nodes2)
                {
                    muf2_values.push_back(s1 * s1 + s2 * s2);
                }
            }
            break;

        case fastNLO::kQuadraticMean:
            for (auto const s1 : scale_nodes1)
            {
                for (auto const s2 : scale_nodes2)
                {
                    muf2_values.push_back(0.5 * (s1 * s1 + s2 * s2));
                }
            }
            break;

        default:
            // TODO: NYI
            assert( false );
        }

        std::vector<double> mu2_values;

        for (std::size_t i = 0; i != scale_nodes1.size() * scale_nodes2.size(); ++i)
        {
            mu2_values.push_back(mur2_values.at(i));
            mu2_values.push_back(muf2_values.at(i));
        }

        std::vector<double> x2_values{1.0};

        for (std::size_t subproc = 0; subproc != n_subproc; ++subproc)
        {
            auto* subgrid = pineappl_subgrid_new2(mu2_values.size() / 2, mu2_values.data(),
                x1_values.size(), x1_values.data(), x2_values.size(), x2_values.data());

            auto const factor = rescale / table->GetNevt(obs, subproc);
            bool non_zero_subgrid = false;

            std::size_t mu2_slice = 0;

            for (std::size_t is1 = 0; is1 != scale_nodes1.size(); ++is1)
            {
                for (std::size_t is2 = 0; is2 != scale_nodes2.size(); ++is2)
                {
                    std::vector<double> slice(x1_values.size());
                    bool non_zero = false;
                    auto const logmur2 = std::log(mu2_values.at(2 * mu2_slice + 0));
                    auto const logmuf2 = std::log(mu2_values.at(2 * mu2_slice + 1));

                    // flexible scale grids only allow one initial-state hadron
                    for (std::size_t ix = 0; ix != sigma_tildes.at(0)->at(obs).size(); ++ix)
                    {
                        double value = sigma_tildes.at(0)->at(obs).at(ix).at(is1).at(is2)
                            .at(subproc);

                        if (table->GetNScaleDep() >= 5)
                        {
                            // mur
                            value += logmur2 *
                                sigma_tildes.at(1)->at(obs).at(ix).at(is1).at(is2).at(subproc);
                            // muf
                            value += logmuf2 *
                                sigma_tildes.at(2)->at(obs).at(ix).at(is1).at(is2).at(subproc);

                            if (table->GetNScaleDep() >= 6)
                            {
                                // mur mur
                                value += logmur2 * logmur2 *
                                    sigma_tildes.at(3)->at(obs).at(ix).at(is1).at(is2).at(subproc);
                            }

                            if (table->GetNScaleDep() >= 7)
                            {
                                // muf muf
                                value += logmuf2 * logmuf2 *
                                    sigma_tildes.at(4)->at(obs).at(ix).at(is1).at(is2).at(subproc);
                                // mur muf
                                value += logmur2 * logmuf2 *
                                    sigma_tildes.at(5)->at(obs).at(ix).at(is1).at(is2).at(subproc);
                            }
                        }

                        if (value != 0.0)
                        {
                            non_zero = true;
                            slice.at(ix) = value * factor * x1_values.at(ix);
                        }
                    }

                    if (non_zero)
                    {
                        non_zero_subgrid = true;
                        pineappl_subgrid_import_mu2_slice(subgrid, mu2_slice, slice.data());
                    }

                    ++mu2_slice;
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

    return pgrid;
}

int this_would_be_main(char const* in, char const* out)
{
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
                auto const mur_ff = file.GetMuRFunctionalForm();
                auto const muf_ff = file.GetMuFFunctionalForm();

                grids.push_back(convert_coeff_add_flex(converted, file, mur_ff, muf_ff, bins, alpha,
                    file.GetIpublunits()));
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

    pineappl_grid_convolute_with_one(
        grids.at(0),
        2212,
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

    pineappl_grid_write(grids.at(0), out);
    pineappl_grid_delete(grids.at(0));

    if (different)
    {
        return 1;
    }

    return 0;
}

template <typename T>
rust::Vec<T> std_vector_to_rust_vec(std::vector<T> vector)
{
    rust::Vec<T> result;
    result.reserve(vector.size());
    std::move(vector.begin(), vector.end(), std::back_inserter(result));
    return result;
}

rust::Vec<double> CalcPDFLinearCombination(
    fastNLOPDFLinearCombinations const& lc,
    fastNLOCoeffAddBase const& base,
    rust::Slice<double const> pdfx1,
    rust::Slice<double const> pdfx2,
    bool pdf2IsAntiParticle
) {
    std::vector<double> fx1(pdfx1.begin(), pdfx1.end());
    std::vector<double> fx2(pdfx2.begin(), pdfx2.end());

    return std_vector_to_rust_vec(lc.CalcPDFLinearCombination(&base, fx1, fx2, pdf2IsAntiParticle));
}

rust::Vec<double> GetScaleNodes(fastNLOCoeffAddFix const& coeffs, int iObs, int iSvar)
{
    return std_vector_to_rust_vec(coeffs.GetScaleNodes(iObs, iSvar));
}

rust::Vec<double> GetXNodes1(fastNLOCoeffAddBase const& coeffs, int iObsBin)
{
    return std_vector_to_rust_vec(coeffs.GetXNodes1(iObsBin));
}

rust::Vec<double> GetXNodes2(fastNLOCoeffAddBase const& coeffs, int iObsBin)
{
    return std_vector_to_rust_vec(coeffs.GetXNodes2(iObsBin));
}

std::unique_ptr<fastNLOLHAPDF> make_fastnlo_lhapdf_with_name_file_set(rust::Str name, rust::Str LHAPDFfile, int PDFSet)
{
    std::string arg0(name.begin(), name.end());
    std::string arg1(LHAPDFfile.begin(), LHAPDFfile.end());

    return std::unique_ptr<fastNLOLHAPDF>(new fastNLOLHAPDF(arg0, arg1, PDFSet));
}

rust::Vec<double> GetCrossSection(fastNLOReader& reader, bool lNorm)
{
    return std_vector_to_rust_vec(reader.GetCrossSection(lNorm));
}

rust::Vec<double> GetBinSize(fastNLOTable const& table)
{
    return std_vector_to_rust_vec(table.GetBinSize());
}

rust::Vec<double> GetScaleNodes1(fastNLOCoeffAddFlex const& coeffs, int iObsBin)
{
    return std_vector_to_rust_vec(coeffs.GetScaleNodes1(iObsBin));
}

rust::Vec<double> GetScaleNodes2(fastNLOCoeffAddFlex const& coeffs, int iObsBin)
{
    return std_vector_to_rust_vec(coeffs.GetScaleNodes2(iObsBin));
}
