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

    #[namespace = "fastNLO"]
    #[repr(u32)]
    enum ESMCalculation {
        kFixedOrder = 0,
        kThresholdCorrection = 1,
        kElectroWeakCorrection = 2,
        kNonPerturbativeCorrection = 3,
        kContactInteraction = 10,
    }

    #[namespace = "fastNLO"]
    #[repr(u32)]
    enum ESMOrder {
        kLeading,
        kNextToLeading,
        kNextToNextToLeading,
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

        fn GetNObsBin(&self) -> i32;
        fn GetNPDF(&self) -> i32;
        fn GetNSubproc(&self) -> i32;
        fn GetNpow(&self) -> i32;
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

        // TODO: GetSigmaTildes
    }

    unsafe extern "C++" {
        include!("pineappl_fastnlo/src/fastnlo.hpp");

        fn GetBinSize(_: &fastNLOTable) -> Vec<f64>;
        fn GetCrossSection(_: Pin<&mut fastNLOReader>, _: bool) -> Vec<f64>;
        fn GetPDFCoeff(_: &fastNLOCoeffAddBase, index: usize) -> Vec<pair_int_int>;
        fn GetPDFCoeffSize(_: &fastNLOCoeffAddBase) -> usize;
        fn GetScaleNodes(_: &fastNLOCoeffAddFix, _: i32, _: i32) -> Vec<f64>;
        fn GetScaleNodes1(_: &fastNLOCoeffAddFlex, _: i32) -> Vec<f64>;
        fn GetScaleNodes2(_: &fastNLOCoeffAddFlex, _: i32) -> Vec<f64>;
        fn GetXNodes1(_: &fastNLOCoeffAddBase, _: i32) -> Vec<f64>;
        fn GetXNodes2(_: &fastNLOCoeffAddBase, _: i32) -> Vec<f64>;
        fn GetObsBinDimBounds(_: &fastNLOTable, _: u32, _: u32) -> pair_double_double;

        fn downcast_coeff_add_fix_to_base(_: &fastNLOCoeffAddFix) -> &fastNLOCoeffAddBase;
        fn downcast_lhapdf_to_reader(_: &fastNLOLHAPDF) -> &fastNLOReader;
        fn downcast_lhapdf_to_reader_mut(_: Pin<&mut fastNLOLHAPDF>) -> Pin<&mut fastNLOReader>;
        fn downcast_lhapdf_to_table(_: &fastNLOLHAPDF) -> &fastNLOTable;
        fn downcast_reader_to_pdf_linear_combinations(
            _: &fastNLOReader,
        ) -> &fastNLOPDFLinearCombinations;
        fn make_fastnlo_lhapdf_with_name_file_set(
            name: &str,
            lhapdf: &str,
            set: i32,
        ) -> UniquePtr<fastNLOLHAPDF>;

        unsafe fn dynamic_cast_coeff_add_fix(_: *mut fastNLOCoeffBase) -> *mut fastNLOCoeffAddFix;
        unsafe fn dynamic_cast_coeff_add_flex(_: *mut fastNLOCoeffBase)
            -> *mut fastNLOCoeffAddFlex;
    }
}
