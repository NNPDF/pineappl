use pineappl::bin::BinRemapper;

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
    pub fn new_f64(normalizations: Vec<f64>, limits: Vec<(f64, f64)>) -> Self {
        Self::new(BinRemapper::new(normalizations, limits).unwrap())
    }
}
