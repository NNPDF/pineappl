//! Evolution interface.
use numpy::{IntoPyArray, PyArray1};
use pineappl::evolution::{EvolveInfo, OperatorSliceInfo};
use pineappl::pids::PidBasis;

use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl::pids::PidBasis <pids/struct.PidBasis.html>`.
#[pyclass(name = "PidBasis")]
#[derive(Clone)]
pub enum PyPidBasis {
    Pdg,
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

/// PyO3 wrapper to :rustdoc:`pineappl::evolution::OperatorSliceInfo <evolution/struct.OperatorSliceInfo.html>`.
#[pyclass(name = "OperatorSliceInfo")]
#[derive(Clone)]
pub struct PyOperatorSliceInfo {
    pub(crate) info: OperatorSliceInfo,
}

#[pymethods]
impl PyOperatorSliceInfo {
    /// Constructor.
    ///
    /// Parameteters
    /// ------------
    /// fac0 : float
    ///     initial factorization scale
    /// pids0 : list(int)
    ///     flavors available at the initial scale
    /// x0 : list(float)
    ///     x-grid at the initial scale
    /// fac1 : float
    ///     evolved final scale
    /// pids1 : list(int)
    ///     flavors available at the final scale
    /// x1 : list(float)
    ///     x-grid at the final scale
    /// pid_basis : PyPidBasis
    ///     flavor basis reprentation at the initial scale
    #[new]
    pub fn new(
        fac0: f64,
        pids0: Vec<i32>,
        x0: Vec<f64>,
        fac1: f64,
        pids1: Vec<i32>,
        x1: Vec<f64>,
        pid_basis: PyPidBasis,
    ) -> Self {
        Self {
            info: OperatorSliceInfo {
                fac0,
                pids0,
                x0,
                fac1,
                pids1,
                x1,
                pid_basis: pid_basis.into(),
            },
        }
    }
}

/// PyO3 wrapper to :rustdoc:`pineappl::evolution::EvolveInfo <evolution/struct.EvolveInfo.html>`.
#[pyclass(name = "EvolveInfo")]
#[repr(transparent)]
pub struct PyEvolveInfo {
    pub(crate) evolve_info: EvolveInfo,
}

#[pymethods]
impl PyEvolveInfo {
    /// Squared factorization scales of the `Grid`.
    #[getter]
    fn fac1<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.evolve_info.fac1.clone().into_pyarray_bound(py)
    }

    /// Particle identifiers of the `Grid`.
    #[getter]
    fn pids1<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<i32>> {
        self.evolve_info.pids1.clone().into_pyarray_bound(py)
    }

    /// `x`-grid coordinates of the `Grid`.
    #[getter]
    fn x1<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.evolve_info.x1.clone().into_pyarray_bound(py)
    }

    /// Renormalization scales of the `Grid`.
    #[getter]
    fn ren1<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.evolve_info.ren1.clone().into_pyarray_bound(py)
    }
}

/// Register submodule in parent.
pub fn register(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new_bound(parent_module.py(), "evolution")?;
    m.add_class::<PyEvolveInfo>()?;
    m.add_class::<PyOperatorSliceInfo>()?;
    m.add_class::<PyPidBasis>()?;
    parent_module.add_submodule(&m)
}
