//! Binnning interface.

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
    /// Parameters
    /// ----------
    /// normalizations : list(float)
    ///     bin normalizations
    /// limits : list(tuple(float, float))
    ///     bin limits
    #[new]
    pub fn new(normalizations: Vec<f64>, limits: Vec<(f64, f64)>) -> Self {
        Self {
            bin_remapper: BinRemapper::new(normalizations, limits).unwrap(),
        }
    }
}

/// Register submodule in parent.
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
