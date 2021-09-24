use pineappl::subgrid::{SubgridEnum, SubgridParams};
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

    /// Set reweighting.
    ///
    /// **Usage:** `yadism`
    ///
    /// Parameters
    /// ----------
    ///     reweight : bool
    ///         apply reweighting?
    pub fn set_reweight(&mut self, reweight: bool) {
        self.subgrid_params.set_reweight(reweight);
    }

    /// Set number of x bins.
    ///
    /// **Usage:** `yadism`
    ///
    /// Parameters
    /// ----------
    ///     x_bins : int
    ///         number of bins
    pub fn set_x_bins(&mut self, x_bins: usize) {
        self.subgrid_params.set_x_bins(x_bins);
    }

    /// Set :math:`x_{max}`.
    ///
    /// **Usage:** `yadism`
    ///
    /// Parameters
    /// ----------
    ///     x_max : float
    ///         new `x_max`
    pub fn set_x_max(&mut self, x_max: f64) {
        self.subgrid_params.set_x_max(x_max);
    }

    /// Set :math:`x_{min}`.
    ///
    /// **Usage:** `yadism`
    ///
    /// Parameters
    /// ----------
    ///     x_min : float
    ///         new `x_min`
    pub fn set_x_min(&mut self, x_min: f64) {
        self.subgrid_params.set_x_min(x_min);
    }

    /// Set interpolation order for `x_grid`.
    ///
    /// **Usage:** `yadism`
    ///
    /// Parameters
    /// ----------
    ///     x_order : float
    ///         new `x_order`
    pub fn set_x_order(&mut self, x_order: usize) {
        self.subgrid_params.set_x_order(x_order);
    }
}

/// PyO3 wrapper to :rustdoc:`pineappl::subgrid::SubgridEnum <subgrid/struct.SubgridEnum.html>`
///
/// **Usage**: `yadism`, FKTable interface
#[pyclass]
#[derive(Clone)]
#[repr(transparent)]
pub struct PySubgridEnum {
    pub(crate) subgrid_enum: SubgridEnum,
}
