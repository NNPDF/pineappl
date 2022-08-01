#include "pineappl_applgrid/src/applgrid.hpp"

#include <algorithm>
#include <appl_igrid.h>
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
//    photon,        // 22: photon
};

std::array<bool, 13 /*14*/> flavour_map = {
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

lumi_pdf const* dynamic_cast_lumi_pdf(appl::appl_pdf const* pdf)
{
    return dynamic_cast <lumi_pdf const*> (pdf);
}

double weightfun(double x)
{
    return appl::igrid::weightfun(x);
}
