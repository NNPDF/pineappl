use numpy::{IntoPyArray, PyArray1};
use pineappl::evolution::EvolveInfo;

use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl::evolution::EvolveInfo <evolution/struct.EvolveInfo.html>`
#[pyclass]
#[repr(transparent)]
pub struct PyEvolveInfo {
    pub(crate) evolve_info: EvolveInfo,
}

#[pymethods]
impl PyEvolveInfo {
    /// Squared factorization scales of the `Grid`.
    #[getter]
    fn fac1<'py>(&self, py: Python<'py>) -> &'py PyArray1<f64> {
        self.evolve_info.fac1.clone().into_pyarray(py)
    }

    /// Particle identifiers of the `Grid`.
    #[getter]
    fn pids1<'py>(&self, py: Python<'py>) -> &'py PyArray1<i32> {
        self.evolve_info.pids1.clone().into_pyarray(py)
    }

    /// `x`-grid coordinates of the `Grid`.
    #[getter]
    fn x1<'py>(&self, py: Python<'py>) -> &'py PyArray1<f64> {
        self.evolve_info.x1.clone().into_pyarray(py)
    }

    /// Renormalization scales of the `Grid`.
    #[getter]
    fn ren1<'py>(&self, py: Python<'py>) -> &'py PyArray1<f64> {
        self.evolve_info.ren1.clone().into_pyarray(py)
    }
}
