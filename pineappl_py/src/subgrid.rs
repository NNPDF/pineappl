use ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pineappl::subgrid::{Subgrid, SubgridEnum, SubgridParams};
use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

/// PyO3 wrapper to [`pineappl::subgrid::SubgridParams`]
///
/// **Usage**: `yadism`
#[pyclass]
#[repr(transparent)]
pub struct PySubgridParams {
    pub(crate) subgrid_params: SubgridParams,
}

impl PySubgridParams {
    pub(crate) fn new(subgrid_params: SubgridParams) -> Self {
        Self { subgrid_params }
    }
}

impl Clone for PySubgridParams {
    fn clone(&self) -> Self {
        let mut subgrid_params = SubgridParams::default();
        subgrid_params.set_q2_bins(self.subgrid_params.q2_bins());
        subgrid_params.set_q2_max(self.subgrid_params.q2_max());
        subgrid_params.set_q2_min(self.subgrid_params.q2_min());
        subgrid_params.set_q2_order(self.subgrid_params.q2_order());
        subgrid_params.set_reweight(self.subgrid_params.reweight());
        subgrid_params.set_x_bins(self.subgrid_params.x_bins());
        subgrid_params.set_x_max(self.subgrid_params.x_max());
        subgrid_params.set_x_min(self.subgrid_params.x_min());
        subgrid_params.set_x_order(self.subgrid_params.x_order());
        Self { subgrid_params }
    }
}

#[pymethods]
impl PySubgridParams {
    #[new]
    pub fn default() -> Self {
        let subgrid_params = SubgridParams::default();

        Self::new(subgrid_params)
    }

    /// Set reweighting.
    ///
    /// **Usage:** `yadism`
    pub fn set_reweight(&mut self, reweight: bool) {
        self.subgrid_params.set_reweight(reweight);
    }

    /// Set number of x bins.
    ///
    /// **Usage:** `yadism`
    pub fn set_x_bins(&mut self, x_bins: usize) {
        self.subgrid_params.set_x_bins(x_bins);
    }

    /// Set `x_max`.
    ///
    /// **Usage:** `yadism`
    pub fn set_x_max(&mut self, x_max: f64) {
        self.subgrid_params.set_x_max(x_max);
    }

    /// Set `x_min`.
    ///
    /// **Usage:** `yadism`
    pub fn set_x_min(&mut self, x_min: f64) {
        self.subgrid_params.set_x_min(x_min);
    }

    /// Set interpolation order for `x_grid`.
    ///
    /// **Usage:** `yadism`
    pub fn set_x_order(&mut self, x_order: usize) {
        self.subgrid_params.set_x_order(x_order);
    }
}

/// PyO3 wrapper to [`pineappl::subgrid::SubgridEnum`]
///
/// **Usage**: `yadism`, FKTable interface
#[pyclass]
#[derive(Clone)]
#[repr(transparent)]
pub struct PySubgridEnum {
    pub(crate) subgrid_enum: SubgridEnum,
}

#[pymethods]
impl PySubgridEnum {

    /// Clone `x1_grid` to allow Python to access.
    ///
    /// **Usage:** FKTable interface
    pub fn x1_grid<'a>(&self, py: Python<'a>) -> &'a PyArray1<f64> {
        self.subgrid_enum.x1_grid().into_owned().into_pyarray(py)
    }

    /// Clone `x2_grid` to allow Python to access.
    ///
    /// **Usage:** FKTable interface
    pub fn x2_grid<'a>(&self, py: Python<'a>) -> &'a PyArray1<f64> {
        self.subgrid_enum.x2_grid().into_owned().into_pyarray(py)
    }

    /// Export grid as FKTable array.
    ///
    /// **Usage:** FKTable interface
    pub fn fk_subgrid_array<'a>(&self, py: Python<'a>) -> PyResult<&'a PyArray2<f64>> {
        let mut result = Array2::zeros((
            self.subgrid_enum.x1_grid().len(),
            self.subgrid_enum.x2_grid().len(),
        ));

        if self.subgrid_enum.q2_grid().len() != 1 {
            return Err(PyValueError::new_err(format!("FK subgrid: only a single Q2 value is allowed, found:\n {:?}.", self.subgrid_enum.q2_grid())))
        }
        assert_eq!(self.subgrid_enum.q2_grid().len(), 1);

        for ((_, ix1, ix2), &value) in self.subgrid_enum.iter() {
            result[[ix1, ix2]] = value;
        }

        Ok(result.into_pyarray(py))
    }
}
