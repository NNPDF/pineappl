//use super::grid::PyGrid;
use numpy::{IntoPyArray, PyArray4};
use pineappl::fk_table::FkTable;
use pyo3::prelude::*;

use std::fs::File;

#[pyclass]
#[repr(transparent)]
pub struct PyFkTable {
    pub(crate) fk_table: FkTable,
}

#[pymethods]
impl PyFkTable {

    pub fn table<'a>(&self, py: Python<'a>) -> PyResult<&'a PyArray4<f64>> {
        Ok(self.fk_table.table().into_pyarray(py))
    }

    pub fn bins(&self) -> usize {
        self.fk_table.bins()
    }

    pub fn lumi(&self) -> Vec<(i32, i32)> {
        self.fk_table.lumi()
    }

    pub fn muf2(&self) -> f64 {
        self.fk_table.muf2()
    }

    pub fn x_grid(&self) -> Vec<f64> {
        self.fk_table.x_grid()
    }

    /// Write grid to file.
    ///
    /// **Usage:** `pineko`
    pub fn write(&self, path: &str) {
        self.fk_table.write(File::create(path).unwrap()).unwrap();
    }

    /// Convolute grid with pdf.
    ///
    /// **Usage:** `pineko`
    pub fn convolute(&self, xfx1: &PyAny, xfx2: &PyAny, alphas: &PyAny) -> Vec<f64> {
        self.fk_table.convolute(
            &|id, x, q2| f64::extract(xfx1.call1((id, x, q2)).unwrap()).unwrap(),
            &|id, x, q2| f64::extract(xfx2.call1((id, x, q2)).unwrap()).unwrap(),
            &|q2| f64::extract(alphas.call1((q2,)).unwrap()).unwrap(),
            &[],
            &[],
            &[],
            &[(1.0, 1.0)],
        )
    }
}
