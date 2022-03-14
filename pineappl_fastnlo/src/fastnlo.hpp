#ifndef FASTNLO_HPP
#define FASTNLO_HPP

#include "pineappl_fastnlo/src/lib.rs.h"
#include "rust/cxx.h"

#include <fastnlotk/fastNLOCoeffAddBase.h>
#include <fastnlotk/fastNLOCoeffAddFix.h>
#include <fastnlotk/fastNLOCoeffAddFlex.h>
#include <fastnlotk/fastNLOLHAPDF.h>
#include <fastnlotk/fastNLOPDFLinearCombinations.h>
#include <fastnlotk/fastNLOReader.h>
#include <fastnlotk/fastNLOTable.h>
#include <memory>

std::unique_ptr<fastNLOLHAPDF> make_fastnlo_lhapdf_with_name_file_set(
    rust::Str name,
    rust::Str LHAPDFfile,
    int PDFSet
);

rust::Vec<double> CalcPDFLinearCombination(
    fastNLOPDFLinearCombinations const& lc,
    fastNLOCoeffAddBase const& base,
    rust::Slice<double const> pdfx1,
    rust::Slice<double const> pdfx2,
    bool pdf2IsAntiParticle
);

rust::Vec<double> GetScaleNodes(fastNLOCoeffAddFix const& coeffs, int iObs, int iSvar);

rust::Vec<double> GetXNodes1(fastNLOCoeffAddBase const& coeffs, int iObsBin);

rust::Vec<double> GetXNodes2(fastNLOCoeffAddBase const& coeffs, int iObsBin);

rust::Vec<double> GetCrossSection(fastNLOReader& reader, bool lNorm);

rust::Vec<double> GetBinSize(fastNLOTable const& table);

rust::Vec<double> GetScaleNodes1(fastNLOCoeffAddFlex const& coeffs, int iObsBin);

rust::Vec<double> GetScaleNodes2(fastNLOCoeffAddFlex const& coeffs, int iObsBin);

std::size_t GetPDFCoeffSize(fastNLOCoeffAddBase const& coeffs);

rust::Vec<pair_int_int> GetPDFCoeff(fastNLOCoeffAddBase const& coeffs, std::size_t index);

double GetSigmaTilde(
    fastNLOCoeffAddFlex const& coeffs,
    std::size_t,
    std::size_t,
    std::size_t,
    std::size_t,
    std::size_t,
    int
);

std::size_t GetNx(fastNLOCoeffAddFlex const& coeffs, std::size_t);

fastNLOCoeffAddBase const& downcast_coeff_add_fix_to_base(fastNLOCoeffAddFix const& coeffs);

fastNLOCoeffAddBase const& downcast_coeff_add_flex_to_base(fastNLOCoeffAddFlex const& coeffs);

fastNLOReader const& downcast_lhapdf_to_reader(fastNLOLHAPDF const& lhapdf);

fastNLOReader& downcast_lhapdf_to_reader_mut(fastNLOLHAPDF& lhapdf);

fastNLOTable const& downcast_lhapdf_to_table(fastNLOLHAPDF const& lhapdf);

fastNLOCoeffAddFix* dynamic_cast_coeff_add_fix(fastNLOCoeffBase* coeffs);

fastNLOCoeffAddFlex* dynamic_cast_coeff_add_flex(fastNLOCoeffBase* coeffs);

fastNLOPDFLinearCombinations const& downcast_reader_to_pdf_linear_combinations(
    fastNLOReader const& reader
);

pair_double_double GetObsBinDimBounds(
    fastNLOTable const& table,
    unsigned int iObs,
    unsigned int iDim
);

#endif
