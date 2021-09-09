//use super::grid::PyGrid;
use numpy::{IntoPyArray, PyArray4};
use pineappl::fk_table::FkTable;
use pyo3::prelude::*;

#[pyclass]
#[repr(transparent)]
pub struct PyFkTable {
    fk_table: FkTable,
}

#[pymethods]
impl PyFkTable {
    //#[new]
    //pub fn try_from(grid: PyGrid) -> PyFkTable {
    //    Self {
    //        fk_table: FkTable::try_from(grid.grid).unwrap(),
    //    }
    //}

    pub fn table<'a>(&self, py: Python<'a>) -> PyResult<&'a PyArray4<f64>> {
        Ok(self.fk_table.table().into_pyarray(py))
    }

    pub fn bins(&self) -> usize {
        self.fk_table.bins()
    }

    pub fn lumi(&self) -> Vec<(i32, i32)> {
        self.fk_table.lumi()
    }

    //pub fn q2(&self) -> f64 {
    //    self.fk_table.q2()
    //}

    //pub fn x_grid(&self) -> &[f64] {
    //    self.fk_table.x_grid()
    //}
}
