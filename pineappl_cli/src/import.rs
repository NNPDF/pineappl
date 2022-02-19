use super::helpers::Subcommand;
use anyhow::{anyhow, Result};
use clap::{Parser, ValueHint};
use std::path::PathBuf;

// the indirection of two modules instead of one is needed to hide the C++ bridge behind the
// feature gate
#[cfg(feature = "fastnlo")]
mod fastnlo {
    #[cxx::bridge]
    pub mod ffi {
        #[namespace = "fastNLO"]
        #[repr(u32)] // TODO: this seems to work for me, but does it work everywhere?
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

        #[namespace = "fastNLO"]
        #[repr(u32)] // TODO: this seems to work for me, but does it work everywhere?
        enum ESMCalculation {
            kFixedOrder = 0,
            kThresholdCorrection = 1,
            kElectroWeakCorrection = 2,
            kNonPerturbativeCorrection = 3,
            kContactInteraction = 10,
        }

        #[namespace = "fastNLO"]
        #[repr(u32)] // TODO: this seems to work for me, but does it work everywhere?
        enum ESMOrder {
            kLeading,
            kNextToLeading,
            kNextToNextToLeading,
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
            include!("pineappl_cli/src/import/fastnlo.hpp");

            unsafe fn this_would_be_main(input: *const c_char, output: *const c_char) -> i32;
        }

        unsafe extern "C++" {
            include!("fastnlotk/fastNLOCoeffAddBase.h");

            type fastNLOCoeffAddBase;

            fn GetNObsBin(&self) -> i32;
            fn GetNPDF(&self) -> i32;
            fn GetNSubproc(&self) -> i32;
            fn GetNpow(&self) -> i32;

            fn GetXNodes1(_: &fastNLOCoeffAddBase, _: i32) -> Vec<f64>;
            fn GetXNodes2(_: &fastNLOCoeffAddBase, _: i32) -> Vec<f64>;
        }

        unsafe extern "C++" {
            include!("fastnlotk/fastNLOCoeffAddFix.h");

            type fastNLOCoeffAddFix;

            fn GetNPDFDim(&self) -> i32;
            fn GetNevt(&self, _: i32, _: i32) -> f64;
            fn GetNxmax(&self, _: i32) -> i32;
            fn GetNxtot2(&self, _: i32) -> i32;
            fn GetScaleFactor(&self, _: i32) -> f64;
            fn GetSigmaTilde(&self, _: i32, _: i32, _: i32, _: i32, _: i32) -> f64;
            fn GetTotalScalenodes(&self) -> i32;
            fn GetTotalScalevars(&self) -> i32;
            fn GetXIndex(&self, _: i32, _: i32, _: i32) -> i32;

            fn GetScaleNodes(_: &fastNLOCoeffAddFix, _: i32, _: i32) -> Vec<f64>;
        }

        unsafe extern "C++" {
            include!("fastnlotk/fastNLOPDFLinearCombinations.h");

            type fastNLOPDFLinearCombinations;

            fn CalcPDFLinearCombination(
                _: &fastNLOPDFLinearCombinations,
                _: &fastNLOCoeffAddBase,
                _: &[f64],
                _: &[f64],
                _: bool,
            ) -> Vec<f64>;
        }

        unsafe extern "C++" {
            include!("fastnlotk/fastNLOLHAPDF.h");

            type fastNLOLHAPDF;

            fn make_fastnlo_lhapdf_with_name_file_set(
                name: &str,
                lhapdf: &str,
                set: i32,
            ) -> UniquePtr<fastNLOLHAPDF>;
        }

        unsafe extern "C++" {
            include!("fastnlotk/fastNLOCoeffBase.h");

            type fastNLOCoeffBase;
        }

        unsafe extern "C++" {
            include!("fastnlotk/fastNLOReader.h");

            type fastNLOReader;

            fn GetMuRFunctionalForm(&self) -> EScaleFunctionalForm;
            fn GetMuFFunctionalForm(&self) -> EScaleFunctionalForm;
            fn ContrId(&self, _: ESMCalculation, _: ESMOrder) -> i32;

            fn GetCrossSection(_: Pin<&mut fastNLOReader>, _: bool) -> Vec<f64>;
        }

        unsafe extern "C++" {
            include!("fastnlotk/fastNLOTable.h");

            type fastNLOTable;

            fn GetCoeffTable(&self, _: i32) -> *mut fastNLOCoeffBase;
            fn GetIpublunits(&self) -> i32;
            fn GetNumDiffBin(&self) -> u32;

            fn GetBinSize(_: &fastNLOTable) -> Vec<f64>;

            // TODO: GetObsBinDimBounds(_: &fastNLOTable, _: u32, _: u32) -> (f64, f64);
        }

        unsafe extern "C++" {
            include!("fastnlotk/fastNLOCoeffAddFlex.h");

            type fastNLOCoeffAddFlex;

            fn GetIXsectUnits(&self) -> i32;
            fn GetNScaleDep(&self) -> i32;
            fn GetPDFPDG(&self, _: i32) -> i32;

            fn GetScaleNodes1(_: &fastNLOCoeffAddFlex, _: i32) -> Vec<f64>;
            fn GetScaleNodes2(_: &fastNLOCoeffAddFlex, _: i32) -> Vec<f64>;

            // TODO: GetSigmaTildes
        }
    }
}

/// Converts fastNLO tables to PineAPPL grids.
#[derive(Parser)]
pub struct Opts {
    /// Path to the input grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    input: PathBuf,
    /// Path to the converted grid.
    #[clap(parse(from_os_str), value_hint = ValueHint::FilePath)]
    output: PathBuf,
}

impl Subcommand for Opts {
    #[cfg(feature = "fastnlo")]
    fn run(&self) -> Result<()> {
        use std::ffi::CString;
        use std::os::unix::ffi::OsStrExt;

        // TODO: this is a correct way of converting strings, but only works on UNIX
        let input = CString::new(self.input.as_path().as_os_str().as_bytes()).unwrap();
        let output = CString::new(self.output.as_path().as_os_str().as_bytes()).unwrap();

        let errcode = unsafe { fastnlo::ffi::this_would_be_main(input.as_ptr(), output.as_ptr()) };

        if errcode == 0 {
            Ok(())
        } else {
            Err(anyhow!("something went wrong"))
        }
    }

    #[cfg(not(feature = "fastnlo"))]
    fn run(&self) -> Result<()> {
        Err(anyhow!(
            "you need to install `pineappl` with feature `fastnlo`"
        ))
    }
}

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

    const HELP_STR: &str = "pineappl-import 
Converts fastNLO tables to PineAPPL grids

USAGE:
    pineappl import <INPUT> <OUTPUT>

ARGS:
    <INPUT>     Path to the input grid
    <OUTPUT>    Path to the converted grid

OPTIONS:
    -h, --help    Print help information
";

    #[test]
    fn help() {
        Command::cargo_bin("pineappl")
            .unwrap()
            .args(&["import", "--help"])
            .assert()
            .success()
            .stdout(HELP_STR);
    }
}
