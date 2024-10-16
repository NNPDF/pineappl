//! Binnning interface.

use numpy::{PyArrayMethods, PyReadonlyArray1};
use pineappl::bin::BinRemapper;
use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl::bin::BinRemapper <bin/struct.BinRemapper.html>`.
#[pyclass(name = "BinRemapper")]
#[derive(Clone)]
#[repr(transparent)]
pub struct PyBinRemapper {
    pub(crate) bin_remapper: BinRemapper,
}

#[pymethods]
impl PyBinRemapper {
    /// Constructor.
    ///
    /// # Panics
    ///
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// normalizations : list(float)
    ///     bin normalizations
    /// limits : list(tuple(float, float))
    ///     bin limits
    #[new]
    #[must_use]
    #[allow(clippy::needless_pass_by_value)]
    pub fn new(normalizations: PyReadonlyArray1<f64>, limits: Vec<(f64, f64)>) -> Self {
        Self {
            bin_remapper: BinRemapper::new(normalizations.to_vec().unwrap(), limits).unwrap(),
        }
    }
}

/// Register submodule in parent.
/// # Errors
///
/// Raises an error if (sub)module is not found.
pub fn register(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new_bound(parent_module.py(), "bin")?;
    m.setattr(pyo3::intern!(m.py(), "__doc__"), "Binning interface.")?;
    pyo3::py_run!(
        parent_module.py(),
        m,
        "import sys; sys.modules['pineappl.bin'] = m"
    );
    m.add_class::<PyBinRemapper>()?;
    parent_module.add_submodule(&m)
}
