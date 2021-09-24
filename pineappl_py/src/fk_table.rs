use numpy::{IntoPyArray, PyArray4};
use pineappl::fk_table::FkTable;
use pyo3::prelude::*;

use std::fs::File;

/// PyO3 wrapper to :rustdoc:`pineappl::fk_table::FkTable <fk_table/struct.FkTable.html>`
///
/// *Usage*: `pineko`, `yadism`
#[pyclass]
#[repr(transparent)]
pub struct PyFkTable {
    pub(crate) fk_table: FkTable,
}

#[pymethods]
impl PyFkTable {
    /// Get cross section tensor.
    ///
    /// Returns
    /// -------
    ///     numpy.ndarray :
    ///         4-dimensional tensor with indixes: bin, lumi, x1, x2
    pub fn table<'a>(&self, py: Python<'a>) -> PyResult<&'a PyArray4<f64>> {
        Ok(self.fk_table.table().into_pyarray(py))
    }

    /// Get number of bins.
    ///
    /// Returns
    /// -------
    ///     int :
    ///         number of bins
    pub fn bins(&self) -> usize {
        self.fk_table.bins()
    }

    /// Get luminsosity functions.
    ///
    /// Returns
    /// -------
    ///     list(tuple(float,float)) :
    ///         luminosity functions as pid tuples
    pub fn lumi(&self) -> Vec<(i32, i32)> {
        self.fk_table.lumi()
    }

    /// Get reference (fitting) scale.
    ///
    /// Returns
    /// -------
    ///     float :
    ///         reference scale
    pub fn muf2(&self) -> f64 {
        self.fk_table.muf2()
    }

    /// Get (unique) interpolation grid.
    ///
    /// Returns
    /// -------
    ///     x_grid : list(float)
    ///         interpolation grid
    pub fn x_grid(&self) -> Vec<f64> {
        self.fk_table.x_grid()
    }

    /// Write grid to file.
    ///
    /// **Usage:** `pineko`
    ///
    /// Parameters
    /// ----------
    ///     path : str
    ///         file path
    pub fn write(&self, path: &str) {
        self.fk_table.write(File::create(path).unwrap()).unwrap();
    }

    /// Convolute grid with pdf.
    ///
    /// **Usage:** `pineko`
    ///
    /// Parameters
    /// ----------
    ///     xfx1 : callable
    ///         lhapdf like callable with arguments `pid, x, Q2` returning x*pdf for :math:`x_1`-grid
    ///     xfx2 : callable
    ///         lhapdf like callable with arguments `pid, x, Q2` returning x*pdf for :math:`x_2`-grid
    ///     alphas : callable
    ///         lhapdf like callable with arguments `Q2` returning :math:`\alpha_s`
    ///
    /// Returns
    /// -------
    ///     list(float) :
    ///         cross sections for all bins
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
