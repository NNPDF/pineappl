use pineappl::grid::{Grid, Order};

use super::bin::PyBinRemapper;
use super::lagrange_subgrid::PyLagrangeSubgridV2;
use super::lumi::PyLumiEntry;
use super::subgrid::PySubgridParams;

use std::fs::File;
use std::io::BufReader;

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
        lumi: Vec<PyRef<PyLumiEntry>>,
        orders: Vec<PyRef<PyOrder>>,
        bin_limits: Vec<f64>,
        subgrid_params: PySubgridParams,
    ) -> Self {
        Self::new(Grid::new(
            lumi.iter().map(|pyl| pyl.lumi_entry.clone()).collect(),
            orders.iter().map(|pyo| pyo.order.clone()).collect(),
            bin_limits,
            subgrid_params.subgrid_params,
        ))
    }

    pub fn set_key_value(&mut self, key: &str, value: &str) {
        self.grid.set_key_value(key, value);
    }

    pub fn set_subgrid(
        &mut self,
        order: usize,
        bin: usize,
        lumi: usize,
        subgrid: PyLagrangeSubgridV2,
    ) {
        self.grid
            .set_subgrid(order, bin, lumi, subgrid.lagrange_subgrid.into());
    }

    pub fn set_remapper(&mut self, remapper: PyBinRemapper) {
        self.grid.set_remapper(remapper.bin_remapper).unwrap();
    }

    #[staticmethod]
    pub fn read(path: &str) -> Self {
        Self::new(Grid::read(BufReader::new(File::open(path).unwrap())).unwrap())
    }

    pub fn write(&self, path: &str) {
        self.grid.write(File::create(path).unwrap()).unwrap();
    }
}
