//! Evolution interface.

use super::convolutions::PyConvType;
use super::pids::PyPidBasis;
use numpy::{IntoPyArray, PyArray1};
use pineappl::evolution::{EvolveInfo, OperatorSliceInfo};
use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl::evolution::OperatorSliceInfo <evolution/struct.OperatorSliceInfo.html>`.
#[pyclass(name = "OperatorSliceInfo")]
#[derive(Clone)]
#[repr(transparent)]
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
    /// conv_type : PyConvType
    ///     the type of convolution required
    #[new]
    #[must_use]
    pub fn new(
        fac0: f64,
        pids0: Vec<i32>,
        x0: Vec<f64>,
        fac1: f64,
        pids1: Vec<i32>,
        x1: Vec<f64>,
        pid_basis: PyPidBasis,
        conv_type: PyRef<PyConvType>,
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
                conv_type: conv_type.convtype,
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
    /// Constructor.
    #[new]
    #[must_use]
    pub const fn new(fac1: Vec<f64>, pids1: Vec<i32>, x1: Vec<f64>, ren1: Vec<f64>) -> Self {
        Self {
            evolve_info: EvolveInfo {
                fac1,
                pids1,
                x1,
                ren1,
            },
        }
    }

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
///
/// # Errors
///
/// Raises an error if (sub)module is not found.
pub fn register(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new_bound(parent_module.py(), "evolution")?;
    m.setattr(pyo3::intern!(m.py(), "__doc__"), "Evolution interface.")?;
    pyo3::py_run!(
        parent_module.py(),
        m,
        "import sys; sys.modules['pineappl.evolution'] = m"
    );
    m.add_class::<PyEvolveInfo>()?;
    m.add_class::<PyOperatorSliceInfo>()?;
    parent_module.add_submodule(&m)
}
