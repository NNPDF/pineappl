use pineappl::bin::BinRemapper;

use numpy::{PyArrayMethods, PyReadonlyArray1};
use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl::bin::BinRemapper <bin/struct.BinRemapper.html>`
///
/// **Usage**: `yadism`
#[pyclass]
#[derive(Clone)]
#[repr(transparent)]
pub struct PyBinRemapper {
    pub(crate) bin_remapper: BinRemapper,
}

impl PyBinRemapper {
    pub(crate) fn new(bin_remapper: BinRemapper) -> Self {
        Self { bin_remapper }
    }
}

#[pymethods]
impl PyBinRemapper {
    #[new]
    pub fn new_f64(normalizations: PyReadonlyArray1<f64>, limits: Vec<(f64, f64)>) -> Self {
        Self::new(BinRemapper::new(normalizations.to_vec().unwrap(), limits).unwrap())
    }
}
