use pineappl::lagrange_subgrid::LagrangeSubgridV2;
use pineappl::subgrid::Subgrid;

use super::subgrid::{PyExtraSubgridParams, PySubgridParams};

use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
#[repr(transparent)]
pub struct PyLagrangeSubgridV2 {
    pub lagrange_subgrid: LagrangeSubgridV2,
}

impl PyLagrangeSubgridV2 {
    pub(crate) fn new(lagrange_subgrid: LagrangeSubgridV2) -> Self {
        Self { lagrange_subgrid }
    }
}

#[pymethods]
impl PyLagrangeSubgridV2 {
    #[new]
    pub fn new_lagrange_subgrid(
        subgrid_params: &PySubgridParams,
        extra_params: &PyExtraSubgridParams,
    ) -> Self {
        Self::new(LagrangeSubgridV2::new(
            &subgrid_params.subgrid_params,
            &extra_params.extra_subgrid_params,
        ))
    }

    fn write_q2_slice(&mut self, q2_slice: usize, grid: Vec<f64>) {
        self.lagrange_subgrid
            .write_q2_slice(q2_slice, grid.as_slice());
    }
}
