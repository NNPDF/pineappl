//! Interpolation interface.

use pineappl::interpolation::{Interp, InterpMeth, Map, ReweightMeth};
use pyo3::{prelude::*, pyclass};

/// PyO3 wrapper to :rustdoc:`pineappl::`.
#[pyclass(name = "Interp")]
#[repr(transparent)]
pub struct PyInterp {
    pub(crate) interp: Interp,
}

impl PyInterp {
    pub(crate) const fn new(interp: Interp) -> Self {
        Self { interp }
    }
}

#[pymethods]
impl PyInterp {
    /// Constructor.
    #[new]
    #[must_use]
    pub fn new_interp(min: f64, max: f64, nodes: usize, order: usize) -> Self {
        Self::new(Interp::new(
            min,
            max,
            nodes,
            order,
            ReweightMeth::ApplGridX,
            Map::ApplGridF2,
            InterpMeth::Lagrange,
        ))
    }
}

/// Register submodule in parent.
/// # Errors
///
/// Raises an error if (sub)module is not found.
pub fn register(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new_bound(parent_module.py(), "interpolation")?;
    m.setattr(pyo3::intern!(m.py(), "__doc__"), "Interpolation submodule.")?;
    pyo3::py_run!(
        parent_module.py(),
        m,
        "import sys; sys.modules['vector_length_module.interpolation'] = m"
    );
    m.add_class::<PyInterp>()?;
    parent_module.add_submodule(&m)
}
