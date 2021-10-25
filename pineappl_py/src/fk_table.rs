use numpy::{IntoPyArray, PyArray4};
use pineappl::fk_table::FkTable;
use pineappl::grid::Grid;
use pineappl::lumi::LumiCache;
use pyo3::prelude::*;
use std::convert::TryFrom;
use std::fs::File;
use std::io::BufReader;

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
    #[staticmethod]
    pub fn read(path: &str) -> Self {
        Self {
            fk_table: FkTable::try_from(
                Grid::read(BufReader::new(File::open(path).unwrap())).unwrap(),
            )
            .unwrap(),
        }
    }

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
    ///     pdg_id : integer
    ///         PDG Monte Carlo ID of the hadronic particle `xfx` is the PDF for
    ///     xfx : callable
    ///         lhapdf like callable with arguments `pid, x, Q2` returning x*pdf for :math:`x`-grid
    ///
    /// Returns
    /// -------
    ///     list(float) :
    ///         cross sections for all bins
    pub fn convolute_with_one(&self, pdg_id: i32, xfx: &PyAny) -> Vec<f64> {
        let mut xfx = |id, x, q2| f64::extract(xfx.call1((id, x, q2)).unwrap()).unwrap();
        let mut alphas = |_| 1.0;
        let mut lumi_cache = LumiCache::with_one(pdg_id, &mut xfx, &mut alphas);
        self.fk_table.convolute(&mut lumi_cache, &[], &[], &[])
    }
}
