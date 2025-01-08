#include "pineappl_fastnlo/src/fastnlo.hpp"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <string>

template <typename T>
rust::Vec<T> std_vector_to_rust_vec(std::vector<T> vector)
{
    rust::Vec<T> result;
    result.reserve(vector.size());
    std::move(vector.begin(), vector.end(), std::back_inserter(result));
    return result;
}

rust::Vec<rust::String> std_vector_string_to_rust_vec_string(std::vector<std::string> const& vector)
{
    rust::Vec<rust::String> result;
    result.reserve(vector.size());
    std::transform(
        vector.begin(),
        vector.end(),
        std::back_inserter(result),
        [](std::string const& s) { return rust::String(s); }
    );
    return result;
}

rust::Vec<double> CalcPDFLinearCombination(
    fastNLOPDFLinearCombinations const& lc,
    fastNLOCoeffAddBase const& base,
    rust::Slice<double const> pdfx1,
    rust::Slice<double const> pdfx2,
    bool pdf2IsAntiParticle
) {
    std::vector<double> fx1(pdfx1.begin(), pdfx1.end());
    std::vector<double> fx2(pdfx2.begin(), pdfx2.end());

    return std_vector_to_rust_vec(lc.CalcPDFLinearCombination(&base, fx1, fx2, pdf2IsAntiParticle));
}

rust::Vec<double> GetScaleNodes(fastNLOCoeffAddFix const& coeffs, int iObs, int iSvar)
{
    return std_vector_to_rust_vec(coeffs.GetScaleNodes(iObs, iSvar));
}

rust::Vec<double> GetXNodes1(fastNLOCoeffAddBase const& coeffs, int iObsBin)
{
    return std_vector_to_rust_vec(coeffs.GetXNodes1(iObsBin));
}

rust::Vec<double> GetXNodes2(fastNLOCoeffAddBase const& coeffs, int iObsBin)
{
    return std_vector_to_rust_vec(coeffs.GetXNodes2(iObsBin));
}

std::unique_ptr<fastNLOLHAPDF> make_fastnlo_lhapdf_with_name_file_set(
    rust::Str name,
    rust::Str LHAPDFfile,
    int PDFSet
) {
    std::string arg0(name.begin(), name.end());
    std::string arg1(LHAPDFfile.begin(), LHAPDFfile.end());

    return std::unique_ptr<fastNLOLHAPDF>(new fastNLOLHAPDF(arg0, arg1, PDFSet));
}

std::unique_ptr<fastNLOCreate> make_fastnlo_create(
    int alphas_lo,
    rust::Slice<rust::Vec<double> const> left_bin_limits,
    rust::Slice<rust::Vec<double> const> right_bin_limits,
    rust::Slice<double const> normalizations,
    int lo_channels,
    int nlo_channels,
    int nnlo_channels,
    rust::Slice<int const> convolutions,
    rust::Slice<rust::Vec<pair_int_int> const> channels
) {
    assert(left_bin_limits.size() == right_bin_limits.size());
    auto const bins = left_bin_limits.size();
    assert(bins == normalizations.size());
    assert(bins > 0);
    auto const dimensions = left_bin_limits.at(0).size();
    assert(dimensions > 0);
    assert(convolutions.size() <= 2);
    assert(convolutions.size() >= 1);

    std::vector<std::vector<double>> bin_limits(dimensions);

    // TODO: check if this is the right ordering
    for (std::size_t i = 0; i != dimensions; ++i) {
        assert(left_bin_limits.at(i).size() == dimensions);
        assert(right_bin_limits.at(i).size() == dimensions);

        //bin_limits.at(i).resize(2 * limits);

        //for (std::size_t j = 0; j != limits; ++j) {
        //    bin_limits.at(i).at(2 * j + 0) = left_bin_limits.at(j).at(i);
        //    bin_limits.at(i).at(2 * j + 1) = right_bin_limits.at(j).at(i);
        //}
        bin_limits.at(i).resize(bins + 1);
        bin_limits.at(i).at(0) = left_bin_limits.at(0).front();

        for (std::size_t j = 0; j != bins; ++j) {
            bin_limits.at(i).at(j + 1) = right_bin_limits.at(j).at(i);
        }
    }

    fastNLO::GeneratorConstants gconst;
    // TODO: add PineAPPL's version number
    gconst.Name = "PineAPPL-fastNLO interface";

    fastNLO::ProcessConstants pconst;
    pconst.LeadingOrder = alphas_lo;
    pconst.NPDF = convolutions.size();
    pconst.NSubProcessesLO = lo_channels;
    pconst.NSubProcessesNLO = nlo_channels;
    pconst.NSubProcessesNNLO = nnlo_channels;

    if (convolutions.size() == 1) {
        pconst.IPDFdef1 = 2;
    } else {
        pconst.IPDFdef1 = 3;
    }

    // TODO: is this the correct value to set the linear combinations ourselves?
    pconst.IPDFdef2 = 0;
    pconst.IPDFdef3LO = 2;
    pconst.IPDFdef3NLO = 0;
    pconst.IPDFdef3NNLO = 0;

    if (convolutions.size() == 1) {
        // TODO: not yet implemented
        assert(false);
    } else {
        pconst.NPDFDim = 2;
    }

    std::vector<std::vector<std::pair<int, int>>> linear_combinations(channels.size());
    for (std::size_t i = 0; i != channels.size(); ++i) {
        std::vector<std::pair<int, int>> entries(channels.at(i).size());
        for (std::size_t j = 0; j != channels.at(i).size(); ++j) {
            auto const first = channels.at(i).at(j).first;
            auto const second = channels.at(i).at(j).second;
            entries.at(j) = std::make_pair(first, second);
        }
        linear_combinations.at(i) = entries;
    }
    pconst.PDFCoeffLO = linear_combinations;

    fastNLO::ScenarioConstants sconst;
    sconst.DifferentialDimension = dimensions;
    sconst.DimensionIsDifferential = std::vector<int>(dimensions, 0);
    sconst.CalculateBinSize = false;
    sconst.BinSize = std::vector<double>(normalizations.begin(), normalizations.end());

    switch (sconst.DifferentialDimension) {
        case 1:
            sconst.SingleDifferentialBinning = bin_limits.at(0);
            break;

        case 2:
            sconst.DoubleDifferentialBinning = bin_limits;
            break;

        case 3:
            sconst.TripleDifferentialBinning = bin_limits;
            break;

        default:
            // ASSERT: there are no or too many dimensions, which fastNLO doesn't support
            assert(false);
    }
    sconst.FlexibleScaleTable = true;

    if (convolutions.size() == 1) {
        sconst.PDF1 = convolutions.at(0);
        // TODO: do we leave PDF2 unchanged (set to 'proton') for DIS?
    } else {
        sconst.PDF1 = convolutions.at(0);
        sconst.PDF2 = convolutions.at(1);
    }

    sconst.ReadBinningFromSteering = true;
    sconst.IgnoreWarmupBinningCheck = true;
    sconst.X_NNodeCounting = "NodesPerBin";
    sconst.Mu1_NNodeCounting = "NodesPerBin";
    sconst.Mu2_NNodeCounting = "NodesPerBin";

    fastNLO::WarmupConstants wconst(sconst);

    for (std::size_t bin = 0; bin != bins; ++bin) {
        std::vector<double> limits;

        // TODO: the following code assumes one dimension
        assert(dimensions == 1);

        // TODO: don't know the meaning of this field
        limits.push_back(-1.0);
        // left bin limit
        limits.push_back(bin_limits.at(0).at(2 * bin + 0));
        // right bin limit
        limits.push_back(bin_limits.at(0).at(2 * bin + 1));

        wconst.Binning.push_back(limits);
    }

    // these values are probably irrelevant but must nevertheless given
    wconst.Values.resize(bins, std::vector<double>{
        // bin index
        0,
        // x-min
        2e-7,
        // x-max
        1.0,
        // scale1-min
        10.0,
        // scale1-max
        100.0,
        // scale2-min
        10.0,
        // scale2-max
        100.0
    });
    for (std::size_t i = 0; i != wconst.Values.size(); ++i) {
        wconst.Values.at(i).at(0) = static_cast <double> (i);
    }
    // wconst.headerValues = ;

    return std::unique_ptr<fastNLOCreate>(new fastNLOCreate(gconst, pconst, sconst, wconst));
}

rust::Vec<double> GetCrossSection(fastNLOReader& reader, bool lNorm)
{
    return std_vector_to_rust_vec(reader.GetCrossSection(lNorm));
}

rust::Vec<rust::String> GetDimLabels(fastNLOTable const& table)
{
    return std_vector_string_to_rust_vec_string(table.GetDimLabels());
}

rust::Vec<rust::String> GetScDescr(fastNLOTable const& table)
{
    return std_vector_string_to_rust_vec_string(table.GetScDescr());
}

rust::String GetXSDescr(fastNLOTable const& table)
{
    return rust::String(table.GetXSDescr());
}

rust::Vec<double> GetScaleNodes1(fastNLOCoeffAddFlex const& coeffs, int iObsBin)
{
    return std_vector_to_rust_vec(coeffs.GetScaleNodes1(iObsBin));
}

rust::Vec<double> GetScaleNodes2(fastNLOCoeffAddFlex const& coeffs, int iObsBin)
{
    return std_vector_to_rust_vec(coeffs.GetScaleNodes2(iObsBin));
}

std::size_t GetPDFCoeffSize(fastNLOCoeffAddBase const& coeffs)
{
    return coeffs.GetPDFCoeff().size();
}

rust::Vec<pair_int_int> GetPDFCoeff(fastNLOCoeffAddBase const& coeffs, std::size_t index)
{
    auto const& cpp_result = coeffs.GetPDFCoeff().at(index);
    rust::Vec<pair_int_int> result;
    result.reserve(cpp_result.size());

    for (auto const& pair : cpp_result)
    {
        pair_int_int res;
        res.first = pair.first;
        res.second = pair.second;
        result.push_back(res);
    }

    return result;
}

double GetSigmaTilde(
    fastNLOCoeffAddFlex const& coeffs,
    std::size_t mu,
    std::size_t obs,
    std::size_t ix,
    std::size_t is1,
    std::size_t is2,
    int subproc
) {
    return coeffs.GetSigmaTildes().at(mu)->at(obs).at(ix).at(is1).at(is2).at(subproc);
}

std::size_t GetNx(fastNLOCoeffAddFlex const& coeffs, std::size_t obs)
{
    return coeffs.GetSigmaTildes().at(0)->at(obs).size();
}

fastNLOCoeffAddBase const& downcast_coeff_add_fix_to_base(fastNLOCoeffAddFix const& coeffs)
{
    return coeffs;
}

fastNLOCoeffAddBase const& downcast_coeff_add_flex_to_base(fastNLOCoeffAddFlex const& coeffs)
{
    return coeffs;
}

fastNLOReader const& downcast_lhapdf_to_reader(fastNLOLHAPDF const& lhapdf)
{
    return lhapdf;
}

fastNLOReader& downcast_lhapdf_to_reader_mut(fastNLOLHAPDF& lhapdf)
{
    return lhapdf;
}

fastNLOTable const& downcast_lhapdf_to_table(fastNLOLHAPDF const& lhapdf)
{
    return lhapdf;
}

fastNLOCoeffAddFix const* dynamic_cast_coeff_add_fix(fastNLOCoeffBase const* coeffs)
{
    return dynamic_cast <fastNLOCoeffAddFix const*> (coeffs);
}

fastNLOCoeffAddFlex const* dynamic_cast_coeff_add_flex(fastNLOCoeffBase const* coeffs)
{
    return dynamic_cast <fastNLOCoeffAddFlex const*> (coeffs);
}

fastNLOCoeffData const* dynamic_cast_coeff_data(fastNLOCoeffBase const* coeffs)
{
    return dynamic_cast <fastNLOCoeffData const*> (coeffs);
}

fastNLOCoeffMult const* dynamic_cast_coeff_mult(fastNLOCoeffBase const* coeffs)
{
    return dynamic_cast <fastNLOCoeffMult const*> (coeffs);
}

fastNLOPDFLinearCombinations const& downcast_reader_to_pdf_linear_combinations(
    fastNLOReader const& reader
) {
    return reader;
}

pair_double_double GetObsBinDimBounds(
    fastNLOTable const& table,
    unsigned int iObs,
    unsigned int iDim
) {
    pair_double_double result;
    auto const cpp_result = table.GetObsBinDimBounds(iObs, iDim);
    result.first = cpp_result.first;
    result.second = cpp_result.second;
    return result;
}
