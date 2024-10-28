//! Subgrid interface.

use pineappl::subgrid::Mu2;
use pineappl::subgrid::{Subgrid, SubgridEnum};
use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl::subgrid::Mu2 <subgrid/struct.Mu2.html>`
#[pyclass(name = "Mu2")]
#[repr(transparent)]
pub struct PyMu2 {
    pub(crate) mu2: Mu2,
}

#[pymethods]
impl PyMu2 {
    /// Constructor.
    ///
    /// Parameters
    /// ----------
    /// ren : float
    ///     renormalization scale
    /// fac : float
    ///     factorization scale
    /// frg : float
    ///     fragmentation scale
    #[new]
    #[must_use]
    pub const fn new(ren: f64, fac: f64, frg: f64) -> Self {
        Self {
            mu2: Mu2 { ren, fac, frg },
        }
    }

    #[getter]
    const fn ren(&self) -> f64 {
        self.mu2.ren
    }

    #[setter]
    fn set_ren(&mut self, value: f64) {
        self.mu2.ren = value;
    }

    #[getter]
    const fn fac(&self) -> f64 {
        self.mu2.fac
    }

    #[setter]
    fn set_fac(&mut self, value: f64) {
        self.mu2.fac = value;
    }

    #[getter]
    const fn frg(&self) -> f64 {
        self.mu2.frg
    }

    #[setter]
    fn set_frg(&mut self, value: f64) {
        self.mu2.frg = value;
    }
}

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
    m.add_class::<PyMu2>()?;
    parent_module.add_submodule(&m)
}
