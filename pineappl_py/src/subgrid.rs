//! Subgrid interface.

use pineappl::subgrid::{Subgrid, SubgridEnum};
use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl::subgrid::SubgridEnum <subgrid/struct.SubgridEnum.html>`
#[pyclass(name = "SubgridEnum")]
#[derive(Clone)]
#[repr(transparent)]
pub struct PySubgridEnum {
    pub(crate) subgrid_enum: SubgridEnum,
}

#[pymethods]
impl PySubgridEnum {
    /// Scale the subgrid by `factor`.
    ///
    /// Parameters
    /// ----------
    /// factor : float
    ///     scaling factor
    pub fn scale(&mut self, factor: f64) {
        self.subgrid_enum.scale(factor);
    }

    /// Clone.
    #[must_use]
    pub fn into(&self) -> Self {
        self.clone()
    }
}

/// Register submodule in parent.
///
/// # Errors
///
/// Raises Errors if (sub-)module is not found.
pub fn register(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new_bound(parent_module.py(), "subgrid")?;
    m.setattr(pyo3::intern!(m.py(), "__doc__"), "Subgrid interface.")?;
    pyo3::py_run!(
        parent_module.py(),
        m,
        "import sys; sys.modules['pineappl.subgrid'] = m"
    );
    m.add_class::<PySubgridEnum>()?;
    parent_module.add_submodule(&m)
}
