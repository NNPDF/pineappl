//! TODO

#![expect(
    clippy::missing_safety_doc,
    reason = "adding a safety section does not seem to work"
)]

use std::str::FromStr;
use thiserror::Error;

/// TODO.
#[cxx::bridge]
pub mod ffi {
    /// TODO.
    #[namespace = "fastNLO"]
    #[repr(u32)]
    enum EScaleFunctionalForm {
        /// TODO.
        kScale1,
        /// TODO.
        kScale2,
        /// TODO.
        kQuadraticSum,
        /// TODO.
        kQuadraticMean,
        /// TODO.
        kQuadraticSumOver4,
        /// TODO.
        kLinearMean,
        /// TODO.
        kLinearSum,
        /// TODO.
        kScaleMax,
        /// TODO.
        kScaleMin,
        /// TODO.
        kProd,
        /// TODO.
        kS2plusS1half,
        /// TODO.
        kPow4Sum,
        /// TODO.
        kWgtAvg,
        /// TODO.
        kS2plusS1fourth,
        /// TODO.
        kExpProd2,
        /// TODO.
        kExtern,
        /// TODO.
        kConst,
    }

    /// TODO.
    struct pair_double_double {
        /// TODO.
        first: f64,
        /// TODO.
        second: f64,
    }

    /// TODO.
    struct pair_int_int {
        /// TODO.
        first: i32,
        /// TODO.
        second: i32,
    }

    extern "C++" {
        include!("fastnlotk/fastNLOConstants.h");

        #[namespace = "fastNLO"]
        type EScaleFunctionalForm;
    }

    unsafe extern "C++" {
        include!("fastnlotk/fastNLOCoeffAddBase.h");

        /// TODO.
        type fastNLOCoeffAddBase;

        /// TODO.
        fn GetNevt(&self, _: i32, _: i32) -> f64;
        /// TODO.
        fn GetNObsBin(&self) -> i32;
        /// TODO.
        fn GetNPDF(&self) -> i32;
        /// TODO.
        fn GetNSubproc(&self) -> i32;
        /// TODO.
        fn GetNpow(&self) -> i32;
    }

    unsafe extern "C++" {
        include!("fastnlotk/fastNLOCoeffAddFix.h");

        /// TODO.
        type fastNLOCoeffAddFix;

        /// TODO.
        fn GetNPDFDim(&self) -> i32;
        /// TODO.
        fn GetNxmax(&self, _: i32) -> i32;
        /// TODO.
        fn GetNxtot2(&self, _: i32) -> i32;
        /// TODO.
        fn GetScaleFactor(&self, _: i32) -> f64;
        /// TODO.
        fn GetSigmaTilde(&self, _: i32, _: i32, _: i32, _: i32, _: i32) -> f64;
        /// TODO.
        fn GetTotalScalenodes(&self) -> i32;
        /// TODO.
        fn GetTotalScalevars(&self) -> i32;
        /// TODO.
        fn GetXIndex(&self, _: i32, _: i32, _: i32) -> i32;
    }

    unsafe extern "C++" {
        include!("fastnlotk/fastNLOPDFLinearCombinations.h");

        /// TODO.
        type fastNLOPDFLinearCombinations;
    }

    unsafe extern "C++" {
        include!("fastnlotk/fastNLOLHAPDF.h");

        /// TODO.
        type fastNLOLHAPDF;
    }

    unsafe extern "C++" {
        include!("fastnlotk/fastNLOCoeffBase.h");

        /// TODO.
        type fastNLOCoeffBase;

        /// TODO.
        fn IsEnabled(&self) -> bool;
    }

    unsafe extern "C++" {
        include!("fastnlotk/fastNLOReader.h");

        /// TODO.
        type fastNLOReader;

        /// TODO.
        fn CalcCrossSection(self: Pin<&mut fastNLOReader>);

        /// TODO.
        unsafe fn GetIsFlexibleScaleTable(&self, _: *mut fastNLOCoeffAddBase) -> bool;

        /// TODO.
        fn GetMuRFunctionalForm(&self) -> EScaleFunctionalForm;

        /// TODO.
        fn GetMuFFunctionalForm(&self) -> EScaleFunctionalForm;

        /// TODO.
        #[must_use]
        fn SetScaleFactorsMuRMuF(self: Pin<&mut Self>, _: f64, _: f64) -> bool;

        /// TODO.
        fn SetMuRFunctionalForm(self: Pin<&mut fastNLOReader>, _: EScaleFunctionalForm);

        /// TODO.
        fn SetMuFFunctionalForm(self: Pin<&mut fastNLOReader>, _: EScaleFunctionalForm);
    }

    unsafe extern "C++" {
        include!("fastnlotk/fastNLOTable.h");

        /// TODO.
        type fastNLOTable;

        /// TODO.
        fn GetCoeffTable(&self, _: i32) -> *mut fastNLOCoeffBase;
        /// TODO.
        fn GetIpublunits(&self) -> i32;
        /// TODO.
        fn GetNObsBin(&self) -> u32;
        /// TODO.
        fn GetNumDiffBin(&self) -> u32;
    }

    unsafe extern "C++" {
        include!("fastnlotk/fastNLOCoeffAddFlex.h");

        /// TODO.
        type fastNLOCoeffAddFlex;

        /// TODO.
        fn GetIXsectUnits(&self) -> i32;
        /// TODO.
        fn GetNScaleDep(&self) -> i32;
        /// TODO.
        fn GetPDFPDG(&self, _: i32) -> i32;
    }

    unsafe extern "C++" {
        include!("pineappl_fastnlo/src/fastnlo.hpp");

        /// TODO.
        fn CalcPDFLinearCombination(
            _: &fastNLOPDFLinearCombinations,
            _: &fastNLOCoeffAddBase,
            _: &[f64],
            _: &[f64],
            _: bool,
        ) -> Vec<f64>;

        /// TODO.
        #[must_use]
        fn GetCrossSection(_: Pin<&mut fastNLOReader>, _: bool) -> Vec<f64>;
        /// TODO.
        fn GetNx(_: &fastNLOCoeffAddFlex, _: usize) -> usize;
        /// TODO.
        fn GetDimLabels(_: &fastNLOTable) -> Vec<String>;
        /// TODO.
        fn GetObsBinDimBounds(_: &fastNLOTable, _: u32, _: u32) -> pair_double_double;
        /// TODO.
        fn GetScDescr(_: &fastNLOTable) -> Vec<String>;
        /// TODO.
        fn GetXSDescr(_: &fastNLOTable) -> String;
        /// TODO.
        fn GetPDFCoeff(_: &fastNLOCoeffAddBase, index: usize) -> Vec<pair_int_int>;
        /// TODO.
        fn GetPDFCoeffSize(_: &fastNLOCoeffAddBase) -> usize;
        /// TODO.
        fn GetScaleNodes(_: &fastNLOCoeffAddFix, _: i32, _: i32) -> Vec<f64>;
        /// TODO.
        fn GetScaleNodes1(_: &fastNLOCoeffAddFlex, _: i32) -> Vec<f64>;
        /// TODO.
        fn GetScaleNodes2(_: &fastNLOCoeffAddFlex, _: i32) -> Vec<f64>;
        /// TODO.
        fn GetXNodes1(_: &fastNLOCoeffAddBase, _: i32) -> Vec<f64>;
        /// TODO.
        fn GetXNodes2(_: &fastNLOCoeffAddBase, _: i32) -> Vec<f64>;
        /// TODO.
        fn GetSigmaTilde(
            _: &fastNLOCoeffAddFlex,
            _: usize,
            _: usize,
            _: usize,
            _: usize,
            _: usize,
            _: i32,
        ) -> f64;

        /// TODO.
        fn downcast_coeff_add_fix_to_base(_: &fastNLOCoeffAddFix) -> &fastNLOCoeffAddBase;
        /// TODO.
        fn downcast_coeff_add_flex_to_base(_: &fastNLOCoeffAddFlex) -> &fastNLOCoeffAddBase;
        /// TODO.
        fn downcast_lhapdf_to_reader(_: &fastNLOLHAPDF) -> &fastNLOReader;
        /// TODO.
        #[must_use]
        fn downcast_lhapdf_to_reader_mut(_: Pin<&mut fastNLOLHAPDF>) -> Pin<&mut fastNLOReader>;
        /// TODO.
        fn downcast_lhapdf_to_table(_: &fastNLOLHAPDF) -> &fastNLOTable;
        /// TODO.
        fn downcast_reader_to_pdf_linear_combinations(
            _: &fastNLOReader,
        ) -> &fastNLOPDFLinearCombinations;

        /// TODO.
        unsafe fn dynamic_cast_coeff_add_fix(
            _: *const fastNLOCoeffBase,
        ) -> *const fastNLOCoeffAddFix;
        /// TODO.
        unsafe fn dynamic_cast_coeff_add_flex(
            _: *const fastNLOCoeffBase,
        ) -> *const fastNLOCoeffAddFlex;

        /// TODO.
        #[must_use]
        fn make_fastnlo_lhapdf_with_name_file_set(
            _: &str,
            _: &str,
            _: i32,
        ) -> UniquePtr<fastNLOLHAPDF>;
    }
}

/// TODO.
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
            _ => Err(ScaleFunctionalFormParseError(s.to_owned())),
        }
    }
}

impl ffi::EScaleFunctionalForm {
    /// TODO.
    #[must_use]
    pub fn compute_scale(self, s1: f64, s2: f64) -> f64 {
        match self {
            Self::kScale1 => s1 * s1,
            Self::kScale2 => s2 * s2,
            Self::kQuadraticSum => s1.mul_add(s1, s2 * s2),
            Self::kQuadraticMean => 0.5 * s1.mul_add(s1, s2 * s2),
            Self::kQuadraticSumOver4 => 0.25 * s1.mul_add(s1, s2 * s2),
            Self::kLinearMean => 0.25 * (s1 + s2).powi(2),
            Self::kLinearSum => (s1 + s2).powi(2),
            Self::kScaleMax => s1.max(s2).powi(2),
            Self::kScaleMin => s1.min(s2).powi(2),
            Self::kProd => (s1 * s2).powi(2),
            Self::kS2plusS1half => 0.5 * s1.mul_add(s1, 2.0 * s2 * s2),
            Self::kPow4Sum => (s1.powi(4) + s2.powi(4)).sqrt(),
            Self::kWgtAvg => (s1.powi(4) + s2.powi(4)) / s1.mul_add(s1, s2 * s2),
            Self::kS2plusS1fourth => (0.25 * s1).mul_add(s1, s2 * s2),
            Self::kExpProd2 => (s1 * (0.3 * s2).exp()).powi(2),
            Self::kExtern => todo!(),
            Self::kConst => todo!(),
            _ => unreachable!(),
        }
    }
}
