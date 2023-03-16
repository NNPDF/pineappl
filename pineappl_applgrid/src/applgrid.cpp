#include "pineappl_applgrid/src/applgrid.hpp"

#include <algorithm>
#include <appl_igrid.h>
#include <array>
#include <cassert>
#include <cmath>
#include <iterator>
#include <LHAPDF/LHAPDF.h>
#include <string>
#include <vector>

std::unique_ptr<LHAPDF::PDF> pdf;

template <typename T>
rust::Vec<T> std_vector_to_rust_vec(std::vector<T> vector)
{
    rust::Vec<T> result;
    result.reserve(vector.size());
    std::move(vector.begin(), vector.end(), std::back_inserter(result));
    return result;
}

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
    photon,        // 22: photon
};

std::array<bool, 14> flavour_map = {
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
    true,  // 22: photon
};

constexpr int index_to_pdg_id(std::size_t index)
{
    return (index == gluon) ? 21 : ((index == photon) ? 22 : (static_cast <int> (index) - 6));
}

void xfx(double const& x, double const& q, double* xfx)
{
    for (std::size_t i = 0; i != flavour_map.size(); ++i)
    {
        if (flavour_map.at(i))
        {
            xfx[i] = pdf->xfxQ(index_to_pdg_id(i), std::fmin(x, 1.0), q);
        }
        else
        {
            xfx[i] = 0.0;
        }
    }
}

double as(double const& q)
{
    return pdf->alphasQ(q);
}

std::unique_ptr<appl::grid> make_grid(rust::Str filename)
{
    std::string name(filename.begin(), filename.end());
    return std::unique_ptr<appl::grid>(new appl::grid(name));
}

std::unique_ptr<appl::grid> make_new_grid(
    rust::Slice<double const> bin_limits,
    int NQ2,
    double Q2min,
    double Q2max,
    int Q2order,
    int Nx,
    double xmin,
    double xmax,
    int xorder,
    rust::Str genpdf,
    int leading_order,
    int loops,
    rust::Str transform,
    rust::Str qtransform,
    bool is_dis
) {
    return std::unique_ptr<appl::grid>(new appl::grid(
        bin_limits.size() - 1,
        bin_limits.data(),
        NQ2,
        Q2min,
        Q2max,
        Q2order,
        Nx,
        xmin,
        xmax,
        xorder,
        static_cast <std::string> (genpdf),
        leading_order,
        loops,
        static_cast <std::string> (transform),
        static_cast <std::string> (qtransform),
        is_dis
    ));
}

std::unique_ptr<lumi_pdf> make_lumi_pdf(rust::Str s, rust::Slice<int const> combinations)
{
    return std::unique_ptr<lumi_pdf>(new lumi_pdf(
        static_cast <std::string> (s),
        std::vector<int>(combinations.begin(), combinations.end())
    ));
}

rust::Vec<int> grid_combine(appl::grid const& grid)
{
    return std_vector_to_rust_vec(grid.combine());
}

rust::Vec<double> grid_convolute(
    appl::grid& grid,
    rust::Str pdfset,
    int member,
    int nloops,
    double rscale,
    double fscale,
    double escale
) {
    pdf.reset(LHAPDF::mkPDF(std::string(pdfset.begin(), pdfset.end()), member));

    auto const results = grid.vconvolute(xfx, as, nloops, rscale, fscale, escale);

    return std_vector_to_rust_vec(results);
}

double sparse_matrix_get(SparseMatrix3d const& matrix, int x, int y, int z)
{
    return matrix(x, y, z);
}

void sparse_matrix_set(SparseMatrix3d& matrix, int x, int y, int z, double value)
{
    matrix(x, y, z) = value;
}

double weightfun(double x)
{
    return appl::igrid::weightfun(x);
}

template <typename tag, typename tag::type M>
struct access_private_member_variable
{
    friend typename tag::type access(tag)
    {
        return M;
    }
};

struct appl_igrid_m_reweight
{
    typedef bool appl::igrid::*type;
    friend type access(appl_igrid_m_reweight);
};

template class access_private_member_variable<appl_igrid_m_reweight, &appl::igrid::m_reweight>;

// we need access to `m_reweight`, but it is private
bool igrid_m_reweight(appl::igrid const& igrid)
{
	// adapted from https://stackoverflow.com/a/3173080/812178
    return igrid.*access(appl_igrid_m_reweight());
}

SparseMatrix3d& igrid_weightgrid(appl::igrid& igrid, std::size_t lumi)
{
    assert(lumi < static_cast <std::size_t> (igrid.SubProcesses()));
    return *igrid.weightgrid()[lumi];
}

struct appl_grid_m_grids
{
    using type = std::vector<appl::igrid*> (appl::grid::*) [appl::MAXGRIDS];
    friend type access(appl_grid_m_grids);
};

template class access_private_member_variable<appl_grid_m_grids, &appl::grid::m_grids>;

appl::igrid& grid_get_igrid(appl::grid& grid, std::size_t order, std::size_t bin)
{
    assert(order < static_cast <std::size_t> (appl::MAXGRIDS));

    return *(grid.*access(appl_grid_m_grids()))[order].at(bin);
}
