#ifndef APPLGRID_HPP
#define APPLGRID_HPP

#include "pineappl_applgrid/src/lib.rs.h"
#include "rust/cxx.h"

#include <cstddef>
#include <appl_grid/appl_grid.h>
#include <appl_grid/lumi_pdf.h>
#include <appl_igrid.h>
#include <memory>

std::unique_ptr<appl::grid> make_grid(rust::Str filename);

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
);

std::unique_ptr<appl::grid> make_empty_grid(
    rust::Slice<double const> obs,
    rust::Str genpdf,
    int leading_order,
    int nloops,
    rust::Str transform,
    rust::Str qtransform
);

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
);

std::unique_ptr<lumi_pdf> make_lumi_pdf(rust::Str s, rust::Slice<int const> combinations);

rust::Vec<int> grid_combine(appl::grid const& grid);

rust::Vec<double> grid_convolve_with_one(
    appl::grid& grid,
    rust::Fn<void(double const&, double const&, double*, void*)> xfx,
    rust::Fn<double(double const&, void*)> alphas,
    void* user_data,
    int nloops,
    double rscale,
    double fscale,
    double escale
);

double sparse_matrix_get(SparseMatrix3d const& matrix, int x, int y, int z);

void sparse_matrix_set(SparseMatrix3d& matrix, int x, int y, int z, double value);

double weightfun(double x);

bool igrid_m_reweight(appl::igrid const& igrid);

SparseMatrix3d& igrid_weightgrid(appl::igrid& igrid, std::size_t lumi);

appl::igrid& grid_get_igrid(appl::grid& grid, std::size_t order, std::size_t bin);

#endif
