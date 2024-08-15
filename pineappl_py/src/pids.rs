//! PIDs interface.

use pineappl::pids::PidBasis;
use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl::pids::PidBasis <pids/enum.PidBasis.html>`.
#[pyclass(name = "PidBasis")]
#[derive(Clone)]
pub enum PyPidBasis {
    /// PDG Monte Carlo IDs.
    Pdg,
    /// NNPDF's evolution basis IDs.
    Evol,
}

impl From<PyPidBasis> for PidBasis {
    fn from(basis: PyPidBasis) -> Self {
        match basis {
            PyPidBasis::Pdg => Self::Pdg,
            PyPidBasis::Evol => Self::Evol,
        }
    }
}

/// Register submodule in parent.
pub fn register(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new_bound(parent_module.py(), "pids")?;
    m.setattr(pyo3::intern!(m.py(), "__doc__"), "PIDs interface.")?;
    pyo3::py_run!(
        parent_module.py(),
        m,
        "import sys; sys.modules['pineappl.pids'] = m"
    );
    m.add_class::<PyPidBasis>()?;
    parent_module.add_submodule(&m)
}
