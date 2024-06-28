#include "pineappl_applgrid/src/applgrid.hpp"

#include <algorithm>
#include <appl_igrid.h>
#include <array>
#include <cassert>
#include <cmath>
#include <iterator>
#include <string>
#include <vector>

void* user_data = nullptr;
rust::Fn<void(double const&, double const&, double*, void*)>* xfx = nullptr;
rust::Fn<double(double const&, void*)>* alphas = nullptr;

template <typename T>
rust::Vec<T> std_vector_to_rust_vec(std::vector<T> vector)
{
    rust::Vec<T> result;
    result.reserve(vector.size());
    std::move(vector.begin(), vector.end(), std::back_inserter(result));
    return result;
}

std::unique_ptr<appl::grid> make_grid(rust::Str filename)
{
    std::string name(filename.begin(), filename.end());
    return std::unique_ptr<appl::grid>(new appl::grid(name));
}

std::unique_ptr<appl::igrid> make_igrid(
    int nq2,
    double q2min,
    double q2max,
    int q2order,
    int nx,
    double xmin,
    double xmax,
    int xorder,
    rust::Str transform,
    rust::Str qtransform,
    int nproc,
    bool disflag
) {
    return std::unique_ptr<appl::igrid>(new appl::igrid(
        nq2,
        q2min,
        q2max,
        q2order,
        nx,
        xmin,
        xmax,
        xorder,
        static_cast <std::string> (transform),
        static_cast <std::string> (qtransform),
        nproc,
        disflag
    ));
}

std::unique_ptr<appl::grid> make_empty_grid(
    rust::Slice<double const> obs,
    rust::Str genpdf,
    int leading_order,
    int nloops,
    rust::Str transform,
    rust::Str qtransform
) {
    std::vector<double> limits(obs.begin(), obs.end());

    return std::unique_ptr<appl::grid>(new appl::grid(
        limits,
        static_cast <std::string> (genpdf),
        leading_order,
        nloops,
        static_cast <std::string> (transform),
        static_cast <std::string> (qtransform)
    ));
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

rust::Vec<double> grid_convolve_with_one(
    appl::grid& grid,
    rust::Fn<void(double const&, double const&, double*, void*)> xfx,
    rust::Fn<double(double const&, void*)> alphas,
    void* user_data,
    int nloops,
    double rscale,
    double fscale,
    double escale
) {
    // TODO: using global variables isn't thread-safe
    ::user_data = user_data;
    ::xfx = &xfx;
    ::alphas = &alphas;

    auto const results = grid.vconvolute(
        [](double const& x, double const& q2, double* results) {
            (*::xfx)(x, q2, results, ::user_data);
        },
        [](double const& q2) {
            return (*::alphas)(q2, ::user_data);
        },
        nloops,
        rscale,
        fscale,
        escale
    );

    ::user_data = nullptr;
    ::xfx = nullptr;
    ::alphas = nullptr;

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

template struct access_private_member_variable<appl_igrid_m_reweight, &appl::igrid::m_reweight>;

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

template struct access_private_member_variable<appl_grid_m_grids, &appl::grid::m_grids>;

appl::igrid& grid_get_igrid(appl::grid& grid, std::size_t order, std::size_t bin)
{
    assert(order < static_cast <std::size_t> (appl::MAXGRIDS));

    return *(grid.*access(appl_grid_m_grids()))[order].at(bin);
}
