use ndarray::Array3;
use numpy::{IntoPyArray, PyArray1, PyArray3};
use pineappl::subgrid::Mu2;
use pineappl::subgrid::{Subgrid, SubgridEnum, SubgridParams};
use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl::subgrid::SubgridParams <subgrid/struct.SubgridParams.html>`
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

    /// Set number of :math:`Q^2` bins.
    ///
    /// Parameters
    /// ----------
    /// q2_bins : int
    ///     number of bins
    pub fn set_q2_bins(&mut self, q2_bins: usize) {
        self.subgrid_params.set_q2_bins(q2_bins);
    }

    /// Set the upper limit for :math:`Q^2`.
    ///
    /// Parameters
    /// ----------
    /// q2_max: float
    ///     new `q2_max`
    pub fn set_q2_max(&mut self, q2_max: f64) {
        self.subgrid_params.set_q2_max(q2_max);
    }

    /// Set the lower limit for :math:`Q^2`.
    ///
    /// Parameters
    /// ----------
    /// q2_min: float
    ///     new `q2_min`
    pub fn set_q2_min(&mut self, q2_min: f64) {
        self.subgrid_params.set_q2_min(q2_min);
    }

    /// Set interpolation order for :math:`Q^2_{grid}`.
    ///
    /// Parameters
    /// ----------
    /// q2_order : float
    ///     new `q2_order`
    pub fn set_q2_order(&mut self, q2_order: usize) {
        self.subgrid_params.set_q2_order(q2_order);
    }

    /// Set reweighting.
    ///
    /// **Usage:** `yadism`
    ///
    /// Parameters
    /// ----------
    /// reweight : bool
    ///     apply reweighting?
    pub fn set_reweight(&mut self, reweight: bool) {
        self.subgrid_params.set_reweight(reweight);
    }

    /// Set number of x bins.
    ///
    /// **Usage:** `yadism`
    ///
    /// Parameters
    /// ----------
    /// x_bins : int
    ///     number of bins
    pub fn set_x_bins(&mut self, x_bins: usize) {
        self.subgrid_params.set_x_bins(x_bins);
    }

    /// Set :math:`x_{max}`.
    ///
    /// **Usage:** `yadism`
    ///
    /// Parameters
    /// ----------
    /// x_max : float
    ///     new `x_max`
    pub fn set_x_max(&mut self, x_max: f64) {
        self.subgrid_params.set_x_max(x_max);
    }

    /// Set :math:`x_{min}`.
    ///
    /// **Usage:** `yadism`
    ///
    /// Parameters
    /// ----------
    /// x_min : float
    ///     new `x_min`
    pub fn set_x_min(&mut self, x_min: f64) {
        self.subgrid_params.set_x_min(x_min);
    }

    /// Set interpolation order for :math:`x_{grid}`.
    ///
    /// **Usage:** `yadism`
    ///
    /// Parameters
    /// ----------
    /// x_order : float
    ///     new `x_order`
    pub fn set_x_order(&mut self, x_order: usize) {
        self.subgrid_params.set_x_order(x_order);
    }
}

/// PyO3 wrapper to :rustdoc:`pineappl::subgrid::Mu2 <subgrid/struct.Mu2.html>`
#[pyclass]
#[repr(transparent)]
pub struct PyMu2 {
    pub mu2: Mu2,
}

#[pymethods]
impl PyMu2 {
    #[new]
    pub fn new(ren: f64, fac: f64) -> Self {
        Self {
            mu2: Mu2 { ren, fac },
        }
    }

    #[getter]
    fn ren(&self) -> PyResult<f64> {
        Ok(self.mu2.ren)
    }

    #[setter]
    fn set_ren(&mut self, value: f64) -> PyResult<()> {
        self.mu2.ren = value;
        Ok(())
    }

    #[getter]
    fn fac(&self) -> PyResult<f64> {
        Ok(self.mu2.fac)
    }

    #[setter]
    fn set_fac(&mut self, value: f64) -> PyResult<()> {
        self.mu2.fac = value;
        Ok(())
    }
}

/// PyO3 wrapper to :rustdoc:`pineappl::subgrid::SubgridEnum <subgrid/struct.SubgridEnum.html>`
#[pyclass]
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

    /// Return the dense array of the subgrid.
    pub fn to_array3<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray3<f64>> {
        Array3::from(&self.subgrid_enum).into_pyarray_bound(py)
    }

    pub fn into(&self) -> Self {
        self.clone()
    }
    /// Return the array of mu2 objects of a subgrid
    pub fn mu2_grid(&self) -> Vec<PyMu2> {
        self.subgrid_enum
            .mu2_grid()
            .iter()
            .cloned()
            .map(|mu2| PyMu2 { mu2 })
            .collect()
    }
    /// Return the array of x1 of a subgrid
    pub fn x1_grid<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_slice_bound(py, &self.subgrid_enum.x1_grid())
    }
    /// Return the array of x2 of a subgrid
    pub fn x2_grid<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_slice_bound(py, &self.subgrid_enum.x2_grid())
    }
}
