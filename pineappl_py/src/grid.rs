use pineappl::grid::{Grid, Order};

use super::lumi::PyLumiEntry;
use super::subgrid::PySubgridParams;

use pyo3::prelude::*;

#[pyclass]
#[repr(transparent)]
pub struct PyOrder {
    pub order: Order,
}

impl PyOrder {
    pub(crate) fn new(order: Order) -> Self {
        Self { order }
    }
}

#[pymethods]
impl PyOrder {
    #[new]
    pub fn new_order(alphas: u32, alpha: u32, logxir: u32, logxif: u32) -> Self {
        Self::new(Order::new(alphas, alpha, logxir, logxif))
    }
}

#[pyclass]
#[repr(transparent)]
pub struct PyGrid {
    pub grid: Grid,
}

impl PyGrid {
    pub(crate) fn new(grid: Grid) -> Self {
        Self { grid }
    }
}

#[pymethods]
impl PyGrid {
    #[new]
    pub fn new_grid(
        lumi: Vec<&PyLumiEntry>,
        orders: Vec<&PyOrder>,
        bin_limits: Vec<f64>,
        subgrid_params: &PySubgridParams,
    ) -> Self {
        Self::new(Grid::new(
            lumi.iter().map(|pyl| pyl.lumi_entry).collect(),
            orders.iter().map(|pyo| pyo.order).collect(),
            bin_limits,
            subgrid_params.subgrid_params,
        ))
    }
}
