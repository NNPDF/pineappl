//! TODO

#![allow(clippy::missing_safety_doc)]
#![allow(clippy::must_use_candidate)]
#![allow(missing_docs)]

use std::str::FromStr;
use thiserror::Error;

#[cxx::bridge]
pub mod ffi {
    #[namespace = "fastNLO"]
    #[repr(u32)]
    enum EScaleFunctionalForm {
        kScale1,
        kScale2,
        kQuadraticSum,
        kQuadraticMean,
        kQuadraticSumOver4,
        kLinearMean,
        kLinearSum,
        kScaleMax,
        kScaleMin,
        kProd,
        kS2plusS1half,
        kPow4Sum,
        kWgtAvg,
        kS2plusS1fourth,
        kExpProd2,
        kExtern,
        kConst,
    }

    struct pair_double_double {
        first: f64,
        second: f64,
    }

    struct pair_int_int {
        first: i32,
        second: i32,
    }

    extern "C++" {
        include!("fastnlotk/fastNLOConstants.h");

        #[namespace = "fastNLO"]
        type EScaleFunctionalForm;

        #[namespace = "fastNLO"]
        type ESMCalculation;

        #[namespace = "fastNLO"]
        type ESMOrder;
    }

    unsafe extern "C++" {
        include!("fastnlotk/fastNLOCoeffAddBase.h");

        type fastNLOCoeffAddBase;

        fn GetNevt(&self, _: i32, _: i32) -> f64;
        fn GetNObsBin(&self) -> i32;
        fn GetNPDF(&self) -> i32;
        fn GetNSubproc(&self) -> i32;
        fn GetNpow(&self) -> i32;
    }

    unsafe extern "C++" {
        include!("fastnlotk/fastNLOCoeffAddFix.h");

        type fastNLOCoeffAddFix;

        fn GetNPDFDim(&self) -> i32;
        fn GetNxmax(&self, _: i32) -> i32;
        fn GetNxtot2(&self, _: i32) -> i32;
        fn GetScaleFactor(&self, _: i32) -> f64;
        fn GetSigmaTilde(&self, _: i32, _: i32, _: i32, _: i32, _: i32) -> f64;
        fn GetTotalScalenodes(&self) -> i32;
        fn GetTotalScalevars(&self) -> i32;
        fn GetXIndex(&self, _: i32, _: i32, _: i32) -> i32;
    }

    unsafe extern "C++" {
        include!("fastnlotk/fastNLOPDFLinearCombinations.h");

        type fastNLOPDFLinearCombinations;
    }

    unsafe extern "C++" {
        include!("fastnlotk/fastNLOLHAPDF.h");

        type fastNLOLHAPDF;
    }

    unsafe extern "C++" {
        include!("fastnlotk/fastNLOCoeffBase.h");

        type fastNLOCoeffBase;

        fn IsEnabled(&self) -> bool;
    }

    unsafe extern "C++" {
        include!("fastnlotk/fastNLOReader.h");

        type fastNLOReader;

        fn CalcCrossSection(self: Pin<&mut fastNLOReader>);
        unsafe fn GetIsFlexibleScaleTable(&self, _: *mut fastNLOCoeffAddBase) -> bool;
        fn GetMuRFunctionalForm(&self) -> EScaleFunctionalForm;
        fn GetMuFFunctionalForm(&self) -> EScaleFunctionalForm;
        fn SetScaleFactorsMuRMuF(self: Pin<&mut Self>, _: f64, _: f64) -> bool;
        fn SetMuRFunctionalForm(self: Pin<&mut fastNLOReader>, _: EScaleFunctionalForm);
        fn SetMuFFunctionalForm(self: Pin<&mut fastNLOReader>, _: EScaleFunctionalForm);
    }

    unsafe extern "C++" {
        include!("fastnlotk/fastNLOTable.h");

        type fastNLOTable;

        fn GetCoeffTable(&self, _: i32) -> *mut fastNLOCoeffBase;
        fn GetIpublunits(&self) -> i32;
        fn GetNObsBin(&self) -> u32;
        fn GetNumDiffBin(&self) -> u32;
    }

    unsafe extern "C++" {
        include!("fastnlotk/fastNLOCoeffAddFlex.h");

        type fastNLOCoeffAddFlex;

        fn GetIXsectUnits(&self) -> i32;
        fn GetNScaleDep(&self) -> i32;
        fn GetPDFPDG(&self, _: i32) -> i32;
    }

    unsafe extern "C++" {
        include!("fastnlotk/fastNLOCoeffData.h");

        type fastNLOCoeffData;
    }

    unsafe extern "C++" {
        include!("fastnlotk/fastNLOCoeffMult.h");

        type fastNLOCoeffMult;
    }

    unsafe extern "C++" {
        include!("pineappl_fastnlo/src/fastnlo.hpp");

        fn CalcPDFLinearCombination(
            _: &fastNLOPDFLinearCombinations,
            _: &fastNLOCoeffAddBase,
            _: &[f64],
            _: &[f64],
            _: bool,
        ) -> Vec<f64>;

        fn GetCrossSection(_: Pin<&mut fastNLOReader>, _: bool) -> Vec<f64>;
        fn GetNx(_: &fastNLOCoeffAddFlex, _: usize) -> usize;
        fn GetDimLabels(_: &fastNLOTable) -> Vec<String>;
        fn GetObsBinDimBounds(_: &fastNLOTable, _: u32, _: u32) -> pair_double_double;
        fn GetScDescr(_: &fastNLOTable) -> Vec<String>;
        fn GetXSDescr(_: &fastNLOTable) -> String;
        fn GetPDFCoeff(_: &fastNLOCoeffAddBase, index: usize) -> Vec<pair_int_int>;
        fn GetPDFCoeffSize(_: &fastNLOCoeffAddBase) -> usize;
        fn GetScaleNodes(_: &fastNLOCoeffAddFix, _: i32, _: i32) -> Vec<f64>;
        fn GetScaleNodes1(_: &fastNLOCoeffAddFlex, _: i32) -> Vec<f64>;
        fn GetScaleNodes2(_: &fastNLOCoeffAddFlex, _: i32) -> Vec<f64>;
        fn GetXNodes1(_: &fastNLOCoeffAddBase, _: i32) -> Vec<f64>;
        fn GetXNodes2(_: &fastNLOCoeffAddBase, _: i32) -> Vec<f64>;
        fn GetSigmaTilde(
            _: &fastNLOCoeffAddFlex,
            _: usize,
            _: usize,
            _: usize,
            _: usize,
            _: usize,
            _: i32,
        ) -> f64;

        fn downcast_coeff_add_fix_to_base(_: &fastNLOCoeffAddFix) -> &fastNLOCoeffAddBase;
        fn downcast_coeff_add_flex_to_base(_: &fastNLOCoeffAddFlex) -> &fastNLOCoeffAddBase;
        fn downcast_lhapdf_to_reader(_: &fastNLOLHAPDF) -> &fastNLOReader;
        fn downcast_lhapdf_to_reader_mut(_: Pin<&mut fastNLOLHAPDF>) -> Pin<&mut fastNLOReader>;
        fn downcast_lhapdf_to_table(_: &fastNLOLHAPDF) -> &fastNLOTable;
        fn downcast_reader_to_pdf_linear_combinations(
            _: &fastNLOReader,
        ) -> &fastNLOPDFLinearCombinations;

        unsafe fn dynamic_cast_coeff_add_fix(
            _: *const fastNLOCoeffBase,
        ) -> *const fastNLOCoeffAddFix;
        unsafe fn dynamic_cast_coeff_add_flex(
            _: *const fastNLOCoeffBase,
        ) -> *const fastNLOCoeffAddFlex;

        fn make_fastnlo_lhapdf_with_name_file_set(
            _: &str,
            _: &str,
            _: i32,
        ) -> UniquePtr<fastNLOLHAPDF>;
    }
}

#[derive(Debug, Error)]
#[error("invalid ScaleFunctionalForm syntax: {0}")]
pub struct ScaleFunctionalFormParseError(String);

impl FromStr for ffi::EScaleFunctionalForm {
    type Err = ScaleFunctionalFormParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "kScale1" => Ok(Self::kScale1),
            "kScale2" => Ok(Self::kScale2),
            "kQuadraticSum" => Ok(Self::kQuadraticSum),
            "kQuadraticMean" => Ok(Self::kQuadraticMean),
            "kQuadraticSumOver4" => Ok(Self::kQuadraticSumOver4),
            "kLinearMean" => Ok(Self::kLinearMean),
            "kLinearSum" => Ok(Self::kLinearSum),
            "kScaleMax" => Ok(Self::kScaleMax),
            "kScaleMin" => Ok(Self::kScaleMin),
            "kProd" => Ok(Self::kProd),
            "kS2plusS1half" => Ok(Self::kS2plusS1half),
            "kPow4Sum" => Ok(Self::kPow4Sum),
            "kWgtAvg" => Ok(Self::kWgtAvg),
            "kS2plusS1fourth" => Ok(Self::kS2plusS1fourth),
            "kExpProd2" => Ok(Self::kExpProd2),
            "kExtern" => Ok(Self::kExtern),
            "kConst" => Ok(Self::kConst),
            _ => Err(ScaleFunctionalFormParseError(s.to_string())),
        }
    }
}

impl ffi::EScaleFunctionalForm {
    pub fn compute_scale(self, s1: f64, s2: f64) -> f64 {
        match self {
            Self::kScale1 => s1 * s1,
            Self::kScale2 => s2 * s2,
            Self::kQuadraticSum => s1 * s1 + s2 * s2,
            Self::kQuadraticMean => 0.5 * (s1 * s1 + s2 * s2),
            Self::kQuadraticSumOver4 => 0.25 * (s1 * s1 + s2 * s2),
            Self::kLinearMean => 0.25 * (s1 + s2).powi(2),
            Self::kLinearSum => (s1 + s2).powi(2),
            Self::kScaleMax => s1.max(s2).powi(2),
            Self::kScaleMin => s1.min(s2).powi(2),
            Self::kProd => (s1 * s2).powi(2),
            Self::kS2plusS1half => 0.5 * (s1 * s1 + 2.0 * s2 * s2),
            Self::kPow4Sum => (s1.powi(4) + s2.powi(4)).sqrt(),
            Self::kWgtAvg => (s1.powi(4) + s2.powi(4)) / (s1 * s1 + s2 * s2),
            Self::kS2plusS1fourth => 0.25 * s1 * s1 + s2 * s2,
            Self::kExpProd2 => (s1 * (0.3 * s2).exp()).powi(2),
            Self::kExtern => todo!(),
            Self::kConst => todo!(),
            _ => unreachable!(),
        }
    }
}
