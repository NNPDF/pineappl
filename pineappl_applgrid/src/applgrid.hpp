#ifndef APPLGRID_HPP
#define APPLGRID_HPP

#include "pineappl_applgrid/src/lib.rs.h"
#include "rust/cxx.h"

#include <appl_grid/appl_grid.h>
#include <appl_grid/lumi_pdf.h>
#include <appl_igrid.h>
#include <memory>

std::unique_ptr<appl::grid> make_grid(rust::Str filename);

std::unique_ptr<appl::grid> make_empty_grid(
    rust::Slice<double const> bin_limits,
    rust::Str genpdf,
    int leading_order,
    int loops,
    rust::Str transform,
    rust::Str qtransform
);

std::unique_ptr<appl::igrid> make_igrid(
    int NQ2,
    double Q2min,
    double Q2max,
    int Q2order,
    int Nx,
    double xmin,
    double xmax,
    int xorder,
    rust::Str transform,
    rust::Str qtransform,
    int Nproc,
    bool disflag
);

std::unique_ptr<lumi_pdf> make_lumi_pdf(rust::Str s, rust::Slice<int const> combinations);

void grid_add_igrid(appl::grid& grid, int bin, int order, std::unique_ptr<appl::igrid> igrid);

rust::Vec<int> grid_combine(appl::grid const& grid);

rust::Vec<double> grid_convolute(
    appl::grid& grid,
    rust::Str pdfset,
    int member,
    int nloops,
    double rscale,
    double fscale,
    double escale
);

double sparse_matrix_get(SparseMatrix3d const& matrix, int x, int y, int z);

double weightfun(double x);

bool igrid_m_reweight(appl::igrid const& igrid);

#endif
