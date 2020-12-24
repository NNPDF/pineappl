use pineappl::subgrid::{ExtraSubgridParams, SubgridParams};

use pyo3::prelude::*;

#[pyclass]
#[repr(transparent)]
pub struct PySubgridParams {
    pub subgrid_params: SubgridParams,
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

        Self { subgrid_params }
    }
}

#[pyclass]
#[repr(transparent)]
pub struct PyExtraSubgridParams {
    pub extra_subgrid_params: ExtraSubgridParams,
}

impl PyExtraSubgridParams {
    pub(crate) fn new(extra_subgrid_params: ExtraSubgridParams) -> Self {
        Self {
            extra_subgrid_params,
        }
    }
}

#[pymethods]
impl PyExtraSubgridParams {
    #[new]
    pub fn default() -> Self {
        let extra_subgrid_params = ExtraSubgridParams::default();

        Self {
            extra_subgrid_params,
        }
    }
}
