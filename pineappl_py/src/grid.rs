use pineappl::grid::{Grid, Order};
use pineappl::lumi::LumiEntry;

use pyo3::prelude::*;

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
        lumi: Vec<LumiEntry>,
        orders: Vec<Orders>,
        bin_limits: Vec<f64>,
        subgrid_params: SubgridParams,
        more_members: MoreMembers,
    ) -> Self {
        Self::new(Grid::new(lumi))
    }
}
