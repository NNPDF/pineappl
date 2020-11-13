#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>

#include <appl_grid/appl_grid.h>
#include <appl_grid/appl_igrid.h>
#include <appl_grid/lumi_pdf.h>
#include <pineappl_capi.h>

int error_exit(std::string const& message)
{
    std::cerr << "Error: " << message << '\n';
    return EXIT_FAILURE;
}

bool is_reweighting_enabled(appl::igrid const& /*igrid*/)
{
    // TODO: that is just a guess
    return true;
}

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        return EXIT_FAILURE;
    }

    std::string in(argv[1]);
    std::string out(argv[2]);

    if (appl::igrid::transformvar() != 5.0)
    {
        return error_exit("`appl::igrid::transformvar != 5` not supported");
    }

    appl::grid const grid(in);

    std::vector<double> bin_limits(grid.Nobs_internal() + 1);
    for (std::size_t i = 0; i != bin_limits.size(); ++i)
    {
        bin_limits.at(i) = grid.obslow_internal(i);
    }

    std::vector<uint32_t> order_params;
    uint32_t const leading_order = grid.leadingOrder();
    int orders;

    if (grid.calculation() == appl::grid::AMCATNLO)
    {
        if (grid.nloops() == 0)
        {
            order_params = {
                leading_order, 0, 0, 0, // LO
            };
            orders = 1;
        }
        else if (grid.nloops() == 1)
        {
            order_params = {
                leading_order + 1, 0, 0, 0, // NLO
                leading_order + 1, 0, 1, 0, // NLO mur
                leading_order + 1, 0, 0, 1, // NLO muf
                leading_order    , 0, 0, 0, // LO
            };
            orders = 4;
        }
        else
        {
            return error_exit("`grid.nloops` not supported");
        }
    }
    else if (grid.calculation() == appl::grid::STANDARD)
    {
        orders = grid.nloops() + 1;

        for (int i = 0; i != orders; ++i)
        {
            order_params.push_back(leading_order + i);
            order_params.push_back(0);
            order_params.push_back(0);
            order_params.push_back(0);
        }
    }
    else
    {
        return error_exit("`grid.calculation = " +
            appl::grid::_calculation(grid.calculation()) + "` not supported");
    }

    std::vector<pineappl_grid*> grids;
    grids.reserve(orders);

    for (int i = 0; i != orders; ++i)
    {
        lumi_pdf const* lumi_ptr = nullptr;

        try
        {
            lumi_ptr = dynamic_cast <lumi_pdf const*> (grid.genpdf(i));
        }
        catch (std::bad_cast const&)
        {
            return error_exit("could not cast into `lumi_pdf`");
        }

        auto* lumi = pineappl_lumi_new();

        for (unsigned j = 0; j != (*lumi_ptr).size(); ++j)
        {
            auto const& combination = (*lumi_ptr)[j];

            std::vector<int32_t> combinations;
            std::vector<double> factors;

            for (unsigned k = 0; k != combination.size(); ++k)
            {
                combinations.push_back(combination[k].first);
                combinations.push_back(combination[k].second);
                factors.push_back(1.0);
            }

            pineappl_lumi_add(lumi, factors.size(), combinations.data(), factors.data());
        }

        auto* key_vals = pineappl_keyval_new();

        auto* pgrid = pineappl_grid_new(lumi, 1, order_params.data() + 4 * i,
            bin_limits.size() - 1, bin_limits.data(), key_vals);
        pineappl_keyval_delete(key_vals);
        pineappl_lumi_delete(lumi);

        grids.push_back(pgrid);

        std::vector<double> slice;

        for (int bin = 0; bin != grid.Nobs_internal(); ++bin)
        {
            auto const* igrid = grid.weightgrid(i, bin);

            if (igrid->transform() != "f2")
            {
                return error_exit("`igrid.transform` not supported");
            }

            if (igrid->Ny1() != igrid->Ny2())
            {
                return error_exit("`igrid.Ny1 != igrid.Ny2`");
            }

            if (igrid->isDISgrid())
            {
                return error_exit("`igrid.isDISgrid == true` not supported");
            }

            if (igrid->y2min() != igrid->y1min())
            {
                return error_exit("`igrid.y2min != igrid.y1min` not supported");
            }

            if (igrid->y2max() != igrid->y1max())
            {
                return error_exit("`igrid.y2max != igrid.y1max` not supported");
            }

            double const q2_min = appl::igrid::fQ2(igrid->taumin());
            double const q2_max = appl::igrid::fQ2(igrid->taumax());
            double const x_min = igrid->fx(igrid->y1min());
            double const x_max = igrid->fx(igrid->y1max());

            int const q2_order = igrid->tauorder();
            int const x_order = igrid->yorder();
            int const q2_bins = igrid->Ntau();
            int const x_bins = igrid->Ny1();

            bool const reweight = is_reweighting_enabled(*igrid);

            auto* keyvals = pineappl_keyval_new();
            pineappl_keyval_set_double(keyvals, "q2_min", q2_min);
            pineappl_keyval_set_double(keyvals, "q2_max", q2_max);
            pineappl_keyval_set_double(keyvals, "x_min", x_min);
            pineappl_keyval_set_double(keyvals, "x_max", x_max);
            pineappl_keyval_set_int(keyvals, "q2_order", q2_order);
            pineappl_keyval_set_int(keyvals, "x_order", x_order);
            pineappl_keyval_set_int(keyvals, "q2_bins", q2_bins);
            pineappl_keyval_set_int(keyvals, "x_bins", x_bins);
            pineappl_keyval_set_bool(keyvals, "reweight", reweight);

            for (std::size_t lumi = 0; lumi != (*lumi_ptr).size(); ++lumi)
            {
                auto const* matrix = const_cast <appl::igrid*> (igrid)->weightgrid(lumi);

                if (matrix == nullptr)
                {
                    continue;
                }

                int itau_min = 0;
                int itau_max = 0;
                int ix1_min = 0;
                int ix1_max = 0;
                int ix2_min = 0;
                int ix2_max = 0;

                for (int itau = 0; itau != igrid->Ntau(); ++itau)
                {
                    for (int ix1 = 0; ix1 != igrid->Ny1(); ++ix1)
                    {
                        for (int ix2 = 0; ix2 != igrid->Ny2(); ++ix2)
                        {
                            if ((*matrix)(itau, ix1, ix2) != 0.0)
                            {
                                itau_min = std::min(itau_min, itau);
                                itau_max = std::max(itau_max, itau + 1);
                                ix1_min = std::min(ix1_min, ix1);
                                ix1_max = std::max(ix1_max, ix1 + 1);
                                ix2_min = std::min(ix2_min, ix2);
                                ix2_max = std::max(ix2_max, ix2 + 1);
                            }
                        }
                    }
                }

                if (itau_min == itau_max)
                {
                    continue;
                }

                auto* subgrid = pineappl_subgrid_new(keyvals);

                for (int itau = itau_min; itau != itau_max; ++itau)
                {
                    slice.assign(igrid->Ny1() * igrid->Ny2(), 0.0);

                    for (int ix1 = ix1_min; ix1 != ix1_max; ++ix1)
                    {
                        for (int ix2 = ix2_min; ix2 != ix2_max; ++ix2)
                        {
                            slice.at(igrid->Ny1() * ix1 + ix2) = (*matrix)(itau, ix1, ix2);
                        }
                    }

                    pineappl_subgrid_fill_q2_slice(subgrid, itau, slice.data());
                }

                pineappl_subgrid_replace_and_delete(pgrid, subgrid, 0, bin, lumi);
            }

            pineappl_keyval_delete(keyvals);
        }
    }

    for (std::size_t i = 1; i < grids.size(); ++i)
    {
        pineappl_grid_merge_and_delete(grids.at(0), grids.at(i));
    }

    pineappl_grid_scale_by_order(grids.at(0), 4.0 * std::acos(-1.0), 1.0, 1.0, 1.0, 1.0);

    if (!grid.getNormalised())
    {
        double const factor = const_cast <appl::grid&> (grid).run();

        if (factor != 0.0)
        {
            pineappl_grid_scale(grids.at(0), 1.0 / factor);
        }
    }

    pineappl_grid_write(grids.at(0), out.c_str());
    pineappl_grid_delete(grids.at(0));
}
