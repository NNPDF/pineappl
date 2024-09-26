use anyhow::Result;
use cxx::UniquePtr;
use lhapdf::Pdf;
use pineappl::grid::Grid;
use pineappl_fastnlo::ffi::fastNLOLHAPDF;
use std::path::Path;
use std::pin::Pin;

pub fn convert_into_fastnlo(
    _grid: &Grid,
    _output: &Path,
    _discard_non_matching_scales: bool,
) -> Result<(UniquePtr<fastNLOLHAPDF>, Vec<bool>)> {
    todo!()
}

pub fn convolve_fastnlo(_grid: Pin<&mut fastNLOLHAPDF>, conv_funs: &mut [Pdf]) -> Vec<f64> {
    // TODO: add support for convolving an APPLgrid with two functions
    assert_eq!(conv_funs.len(), 1);

    todo!()
}
