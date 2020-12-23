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
