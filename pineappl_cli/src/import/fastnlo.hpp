#include <memory>

#include <fastnlotk/fastNLOCoeffAddBase.h>
#include <fastnlotk/fastNLOCoeffAddFix.h>
#include <fastnlotk/fastNLOCoeffAddFlex.h>
#include <fastnlotk/fastNLOLHAPDF.h>
#include <fastnlotk/fastNLOPDFLinearCombinations.h>
#include <fastnlotk/fastNLOReader.h>
#include <fastnlotk/fastNLOTable.h>
#include <rust/cxx.h>

int this_would_be_main(char const* in, char const* out);

rust::Vec<double> CalcPDFLinearCombination(fastNLOPDFLinearCombinations const& lc, fastNLOCoeffAddBase const& base, rust::Slice<double const> pdfx1, rust::Slice<double const> pdfx2, bool pdf2IsAntiParticle);

rust::Vec<double> GetScaleNodes(fastNLOCoeffAddFix const& coeffs, int iObs, int iSvar);

rust::Vec<double> GetXNodes1(fastNLOCoeffAddBase const& coeffs, int iObsBin);

rust::Vec<double> GetXNodes2(fastNLOCoeffAddBase const& coeffs, int iObsBin);

std::unique_ptr<fastNLOLHAPDF> make_fastnlo_lhapdf_with_name_file_set(rust::Str name, rust::Str LHAPDFfile, int PDFSet);

rust::Vec<double> GetCrossSection(fastNLOReader& reader, bool lNorm);

rust::Vec<double> GetBinSize(fastNLOTable const& table);

rust::Vec<double> GetScaleNodes1(fastNLOCoeffAddFlex const& coeffs, int iObsBin);

rust::Vec<double> GetScaleNodes2(fastNLOCoeffAddFlex const& coeffs, int iObsBin);
