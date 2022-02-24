#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>

#include <appl_grid/appl_grid.h>
#include <appl_grid/appl_igrid.h>
#include <appl_grid/lumi_pdf.h>
#include <LHAPDF/LHAPDF.h>
#include <pineappl_capi.h>

std::unique_ptr<LHAPDF::PDF> pdf;

enum flavour_map_index : std::size_t
{
    anti_top,      // -6: anti-top
    anti_bottom,   // -5: anti-bottom
    anti_charm,    // -4: anti-charm
    anti_strange,  // -3: anti-strange
    anti_up,       // -2: anti-up
    anti_down,     // -1: anti-down
    gluon,         // 21: gluon
    down,          //  1: down
    up,            //  2: up
    strange,       //  3: strange
    charm,         //  4: charm
    bottom,        //  5: bottom
    top,           //  6: top
//    photon,        // 22: photon
};

//std::array<bool, 14> flavour_map = {
std::array<bool, 13> flavour_map = {
    true,  // -6: anti-top
    true,  // -5: anti-bottom
    true,  // -4: anti-charm
    true,  // -3: anti-strange
    true,  // -2: anti-up
    true,  // -1: anti-down
    true,  // 21: gluon
    true,  //  1: down
    true,  //  2: up
    true,  //  3: strange
    true,  //  4: charm
    true,  //  5: bottom
    true,  //  6: top
//    true,  // 22: photon
};

constexpr int index_to_pdg_id(std::size_t index)
{
//    return (index == gluon) ? 21 : ((index == photon) ? 22 : (static_cast <int> (index) - 6));
    return (index == gluon) ? 21 : static_cast <int> (index) - 6;
}

extern "C" void evolvepdf(double const& x, double const& q, double* xfx)
{
    for (std::size_t i = 0; i != flavour_map.size(); ++i)
    {
        if (flavour_map.at(i))
        {
            xfx[i] = pdf->xfxQ(index_to_pdg_id(i), std::clamp(x, 0.0, 1.0), q);
        }
        else
        {
            xfx[i] = 0.0;
        }
    }
}

extern "C" double alphaspdf(double const& q)
{
    return pdf->alphasQ(q);
}

extern "C" double xfx(int32_t pdg_id, double x, double q2, void* state)
{
    return static_cast <LHAPDF::PDF*> (state)->xfxQ2(pdg_id, x, q2);
}

extern "C" double alphas(double q2, void* state)
{
    return static_cast <LHAPDF::PDF*> (state)->alphasQ2(q2);
}

void error_exit(std::string const& message)
{
    std::cerr << "Error: " << message << '\n';
    std::exit(EXIT_FAILURE);
}

int32_t convert_to_pdg_id(int id)
{
    if ((id >= -6) && (id <= 6))
    {
        return (id == 0) ? 21 : id;
    }
    else if (id == 7)
    {
        // applgridphoton extension
        return 22;
    }
    else
    {
        assert( false );
    }
}

double ckm_factors(
    int a,
    int b,
    std::vector<std::vector<double>> const& ckm2,
    std::vector<double> const& ckm_sum
) {
    if (ckm_sum.empty())
    {
        return 1.0;
    }

    if (a == 0)
    {
        if (b == 0)
        {
            return 1.0;
        }
        else if (b == 7)
        {
            return 1.0;
        }
        else
        {
            return ckm_sum.at(b + 6);
        }
    }
    else if (a == 7)
    {
        return 1.0;
    }
    else
    {
        if (b == 0)
        {
            return ckm_sum.at(a + 6);
        }
        else if (b == 7)
        {
            return 1.0;
        }
        else
        {
            return ckm2.at(a + 6).at(b + 6);
        }
    }
}

void reconstruct_luminosity_function(appl::grid& grid, int order, pineappl_lumi* lumi)
{
    auto* pdf = const_cast <appl::appl_pdf*> (static_cast <appl::grid const&> (grid).genpdf(order));

    std::vector<std::vector<int32_t>> combinations(pdf->Nproc());
    std::vector<std::vector<double>> factors(pdf->Nproc());
    std::vector<double> xfx1(14);
    std::vector<double> xfx2(14);
    std::vector<double> results(pdf->Nproc());

    for (int a = -6; a != 8; ++a)
    {
        xfx1.at(a + 6) = 1.0;

        for (int b = -6; b != 8; ++b)
        {
            xfx2.at(b + 6) = 1.0;

            pdf->evaluate(xfx1.data(), xfx2.data(), results.data());

            for (std::size_t i = 0; i != results.size(); ++i)
            {
                if (results.at(i) != 0.0)
                {
                    combinations.at(i).push_back(convert_to_pdg_id(a));
                    combinations.at(i).push_back(convert_to_pdg_id(b));
                    factors.at(i).push_back(results.at(i));
                }
            }

            xfx2.at(b + 6) = 0.0;
        }

        xfx1.at(a + 6) = 0.0;
    }

    for (std::size_t i = 0; i != combinations.size(); ++i)
    {
        assert( !factors.at(i).empty() );
        assert( !combinations.at(i).empty() );

        pineappl_lumi_add(lumi, factors.at(i).size(), combinations.at(i).data(),
            factors.at(i).data());
    }
}

std::tuple<pineappl_grid*, bool> convert_grid(appl::grid& grid, bool reweight, uint32_t alpha)
{
    std::vector<double> bin_limits(grid.Nobs_internal() + 1);
    for (std::size_t i = 0; i != bin_limits.size(); ++i)
    {
        bin_limits.at(i) = grid.obslow_internal(i);
    }

    std::vector<uint32_t> order_params;
    uint32_t const leading_order = grid.leadingOrder();
    int orders;
    double alphas_factor;

    if (grid.calculation() == appl::grid::AMCATNLO)
    {
        if (grid.nloops() == 0)
        {
            order_params = {
                leading_order, 0, alpha, 0, // LO
            };
            orders = 1;
        }
        else if (grid.nloops() == 1)
        {
            order_params = {
                leading_order + 1, alpha, 0, 0, // NLO
                leading_order + 1, alpha, 1, 0, // NLO mur
                leading_order + 1, alpha, 0, 1, // NLO muf
                leading_order    , alpha, 0, 0, // LO
            };
            orders = 4;
        }
        else
        {
            error_exit("`grid.nloops` not supported");
        }

        alphas_factor = 4.0 * std::acos(-1.0);
    }
    else if (grid.calculation() == appl::grid::STANDARD)
    {
        orders = grid.nloops() + 1;

        for (int i = 0; i != orders; ++i)
        {
            order_params.push_back(leading_order + i);
            order_params.push_back(alpha);
            order_params.push_back(0);
            order_params.push_back(0);
        }

        alphas_factor = 0.5 / std::acos(-1.0);
    }
    else
    {
        error_exit("`grid.calculation = " +
            appl::grid::_calculation(grid.calculation()) + "` not supported");
    }

    if (grid.getApplyCorrections())
    {
        error_exit("`grid.getApplyCorrections() = true` not supported");
    }

    if (grid.getDynamicScale() != 0.0)
    {
        error_exit("`grid.getDynamicScale() != 1.0` not supported");
    }

    std::vector<pineappl_grid*> grids;
    grids.reserve(orders);

    for (int i = 0; i != orders; ++i)
    {
        lumi_pdf const* lumi_ptr =
            dynamic_cast <lumi_pdf const*> (static_cast <appl::grid const&> (grid).genpdf(i));

        auto* lumi = pineappl_lumi_new();

        if (lumi_ptr == nullptr)
        {
            reconstruct_luminosity_function(grid, i, lumi);
        }
        else
        {
            for (unsigned j = 0; j != (*lumi_ptr).size(); ++j)
            {
                auto const& combination = (*lumi_ptr)[j];

                std::vector<int32_t> combinations;
                std::vector<double> factors;

                for (unsigned k = 0; k != combination.size(); ++k)
                {
                    auto const a = combination[k].first;
                    auto const b = combination[k].second;

                    combinations.push_back(convert_to_pdg_id(a));
                    combinations.push_back(convert_to_pdg_id(b));

                    auto const factor = ckm_factors(a, b, lumi_ptr->getckm2(), lumi_ptr->getckmsum());

                    factors.push_back(factor);
                }

                pineappl_lumi_add(lumi, factors.size(), combinations.data(), factors.data());
            }
        }

        auto* key_vals = pineappl_keyval_new();
        auto* pgrid = pineappl_grid_new(lumi, 1, order_params.data() + 4 * i,
            bin_limits.size() - 1, bin_limits.data(), key_vals);
        auto const lumi_size = pineappl_lumi_count(lumi);
        pineappl_keyval_delete(key_vals);
        pineappl_lumi_delete(lumi);

        grids.push_back(pgrid);

        std::vector<double> slice;

        for (int bin = 0; bin != grid.Nobs_internal(); ++bin)
        {
            auto const* igrid = grid.weightgrid(i, bin);

            std::vector<double> mu2_values(2 * igrid->Ntau());
            std::vector<double> x1_values(igrid->Ny1());
            std::vector<double> x1_weights(igrid->Ny1());
            std::vector<double> x2_values(igrid->Ny2());
            std::vector<double> x2_weights(igrid->Ny2());

            for (std::size_t i = 0; i != mu2_values.size() / 2; ++i)
            {
                mu2_values[2 * i + 0] = appl::igrid::fQ2(igrid->gettau(i));
                mu2_values[2 * i + 1] = appl::igrid::fQ2(igrid->gettau(i));
            }

            bool different_x_grids = false;

            for (std::size_t i = 0; i != x1_values.size(); ++i)
            {
                x1_values[i] = std::clamp(igrid->fx(igrid->gety1(i)), 0.0, 1.0);
                x1_weights[i] = reweight ? appl::igrid::weightfun(x1_values[i]) : 1.0;
            }

            for (std::size_t i = 0; i != x2_values.size(); ++i)
            {
                x2_values[i] = std::clamp(igrid->fx(igrid->gety2(i)), 0.0, 1.0);
                x2_weights[i] = reweight ? appl::igrid::weightfun(x2_values[i]) : 1.0;

                if ((x1_values[i] / x2_values[i]) - 1.0 > 1e-10)
                {
                    different_x_grids = true;
                }
            }

            if (different_x_grids)
            {
                std::cout << ">>> Different x1 and x2 grids!\n";
            }

            slice.resize(x1_values.size() * x2_values.size());

            for (std::size_t lumi = 0; lumi != lumi_size; ++lumi)
            {
                auto const* matrix = const_cast <appl::igrid*> (igrid)->weightgrid(lumi);

                if (matrix == nullptr)
                {
                    continue;
                }

                auto* subgrid = pineappl_subgrid_new2(mu2_values.size() / 2, mu2_values.data(),
                    x1_values.size(), x1_values.data(), x2_values.size(), x2_values.data());

                bool non_zero_subgrid = false;

                for (std::size_t itau = 0; itau != mu2_values.size() / 2; ++itau)
                {
                    bool non_zero = false;

                    for (std::size_t ix1 = 0; ix1 != x1_values.size(); ++ix1)
                    {
                        for (std::size_t ix2 = 0; ix2 != x2_values.size(); ++ix2)
                        {
                            double const value = (*matrix)(itau, ix1, ix2);

                            if (value != 0.0)
                            {
                                non_zero = true;
                            }

                            slice.at(x2_values.size() * ix1 + ix2) =
                                value * x1_weights[ix1] * x2_weights[ix2];
                        }
                    }

                    if (non_zero)
                    {
                        non_zero_subgrid = true;
                        pineappl_subgrid_import_mu2_slice(subgrid, itau, slice.data());
                    }
                }

                if (non_zero_subgrid)
                {
                    pineappl_grid_replace_and_delete(pgrid, subgrid, 0, bin, lumi);
                }
                else
                {
                    pineappl_subgrid_delete(subgrid);
                }
            }
        }
    }

    for (std::size_t i = 1; i < grids.size(); ++i)
    {
        pineappl_grid_merge_and_delete(grids.at(0), grids.at(i));
    }

    double global = 1.0;

    if (!grid.getNormalised())
    {
        double const factor = grid.run();

        if (factor != 0.0)
        {
            global = 1.0 / factor;
        }
    }

    pineappl_grid_scale_by_order(grids.at(0), alphas_factor, 1.0, 1.0, 1.0, global);
    pineappl_grid_optimize(grids.at(0));

    LHAPDF::setVerbosity(0);
    pdf.reset(LHAPDF::mkPDF("NNPDF31_nlo_as_0118_luxqed", 0));

    auto const& results = grid.vconvolute(evolvepdf, evolvepdf, alphaspdf, 1);
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
        auto const two = other_results.at(i);

        // catches the case where both results are zero
        if (one == two)
        {
            continue;
        }

        if (std::fabs(two / one - 1.0) > 1e-10)
        {
            std::cout << ">>> APPLgrid: " << one << " PineAPPL: " << two << " A/P: " << (one / two)
                << " P/A: " << (two / one) << '\n';
            different = true;
        }
        else
        {
            std::cout << ">>> Success!\n";
        }
    }

    if (different)
    {
        return std::make_tuple(grids.at(0), false);
    }

    return std::make_tuple(grids.at(0), true);
}

int main(int argc, char* argv[])
{
    if (argc < 3 && argc > 4)
    {
        return EXIT_FAILURE;
    }

    std::string in(argv[1]);
    std::string out(argv[2]);

    long alpha = (argc == 4) ? std::stoul(argv[3]) : 0;

    appl::grid grid(in);

    std::cout << ">>> Trying `reweight = true`. This may fail.\n";

    pineappl_grid* pgrid = nullptr;
    pineappl_grid* pgrid_reweight_false;
    bool success;
    std::tie(pgrid_reweight_false, success) = convert_grid(grid, true, alpha);

    if (success)
    {
        pgrid = pgrid_reweight_false;
    }
    else
    {
        std::cout << ">>> `reweight = true` didn't work. Trying `reweight = false`.\n";

        pineappl_grid* pgrid_reweight_true;
        std::tie(pgrid_reweight_true, success) = convert_grid(grid, false, alpha);

        if (!success)
        {
            auto const out_false = out + "_reweight_false";
            auto const out_true = out + "_reweight_true";

            pineappl_grid_write(pgrid_reweight_false, out_false.c_str());
            pineappl_grid_write(pgrid_reweight_true, out_true.c_str());

            pineappl_grid_delete(pgrid_reweight_false);
            pineappl_grid_delete(pgrid_reweight_true);

            error_exit("grids are different");
        }

        pineappl_grid_delete(pgrid_reweight_false);
        pgrid = pgrid_reweight_true;
    }

    pineappl_grid_write(pgrid, out.c_str());
    pineappl_grid_delete(pgrid);
}
