use pineappl::bin::BinRemapper;

//use pyo3::types::{PyList, PyTuple};
//use pyo3::{exceptions::PyRuntimeError, Python};

#[repr(transparent)]
pub struct PyBinRemapper {
    pub bin_remapper: BinRemapper,
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
