#include "pineappl_applgrid/src/applgrid.hpp"

// TODO: is this portable enough?
#if defined (__unix__) || defined(__unix) || defined(unix) || \
    (defined (__APPLE__) && defined (__MACH__))
#define HAVE_UNISTD_H
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <algorithm>
#include <appl_igrid.h>
#include <array>
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

std::unique_ptr<appl::grid> make_grid(rust::Str filename, bool silence)
{
    std::string name(filename.begin(), filename.end());
    std::unique_ptr<appl::grid> result;

#ifdef HAVE_UNISTD_H
    int backup_fd = -1;
    FILE *dev_null = NULL;

    if (silence)
    {
        fflush(stdout);
        backup_fd = dup(STDOUT_FILENO);
        dev_null = fopen("/dev/null", "w");
        dup2(fileno(dev_null), STDOUT_FILENO);
    }
#endif

    result.reset(new appl::grid(name));

#ifdef HAVE_UNISTD_H
    if (silence)
    {
        fflush(stdout);
        fclose(dev_null);
        dup2(backup_fd, STDOUT_FILENO);
        close(backup_fd);
    }
#endif

    return result;
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

lumi_pdf const* dynamic_cast_lumi_pdf(appl::appl_pdf const* pdf)
{
    return dynamic_cast <lumi_pdf const*> (pdf);
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
