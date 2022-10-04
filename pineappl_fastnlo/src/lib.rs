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

    #[namespace = "say"]
    #[repr(i32)]
    enum Verbosity {
        DEBUG = -1000,
        MANUAL = 2,
        INFO = 0,
        WARNING = 1,
        ERROR = 2,
        SILENT = 1000,
    }

    unsafe extern "C++" {
        include!("fastnlotk/speaker.h");

        #[namespace = "say"]
        type Verbosity;

        #[namespace = "say"]
        fn SetGlobalVerbosity(_: Verbosity) -> i32;
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

        fn GetMuRFunctionalForm(&self) -> EScaleFunctionalForm;
        fn GetMuFFunctionalForm(&self) -> EScaleFunctionalForm;
    }

    unsafe extern "C++" {
        include!("fastnlotk/fastNLOTable.h");

        type fastNLOTable;

        fn GetCoeffTable(&self, _: i32) -> *mut fastNLOCoeffBase;
        fn GetIpublunits(&self) -> i32;
        fn GetNObsBin(&self) -> u32;
        fn GetNumDiffBin(&self) -> u32;
        fn IsNorm(&self) -> bool;
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
        fn GetObsBinDimBounds(_: &fastNLOTable, _: u32, _: u32) -> pair_double_double;
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

        unsafe fn dynamic_cast_coeff_add_fix(
            _: *const fastNLOCoeffBase,
        ) -> *const fastNLOCoeffAddFix;
        unsafe fn dynamic_cast_coeff_add_flex(
            _: *const fastNLOCoeffBase,
        ) -> *const fastNLOCoeffAddFlex;
        unsafe fn dynamic_cast_coeff_data(_: *const fastNLOCoeffBase) -> *const fastNLOCoeffData;
        unsafe fn dynamic_cast_coeff_mult(_: *const fastNLOCoeffBase) -> *const fastNLOCoeffMult;

        fn static_cast_lhapdf_to_reader_mut(_: Pin<&mut fastNLOLHAPDF>) -> Pin<&mut fastNLOReader>;

        fn make_fastnlo_lhapdf_with_name_file_set(
            name: &str,
            lhapdf: &str,
            set: i32,
            silence: bool,
        ) -> UniquePtr<fastNLOLHAPDF>;
    }
}
