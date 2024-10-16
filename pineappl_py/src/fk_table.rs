//! FK table interface.

use super::convolutions::PyConv;
use super::grid::PyGrid;
use numpy::{IntoPyArray, PyArray1, PyArrayDyn, PyArrayMethods, PyReadonlyArray1};
use pineappl::convolutions::ConvolutionCache;
use pineappl::fk_table::{FkAssumptions, FkTable};
use pineappl::grid::Grid;
use pyo3::prelude::*;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::str::FromStr;

/// PyO3 wrapper to :rustdoc:`pineappl::fk_table::FkTable <fk_table/struct.FkTable.html>`.
#[pyclass(name = "FkTable")]
#[repr(transparent)]
pub struct PyFkTable {
    pub(crate) fk_table: FkTable,
}

/// PyO3 wrapper to :rustdoc:`pineappl::fk_table::FkAssumptions <fk_table/enum.FkAssumptions.html>`.
#[pyclass(name = "FkAssumptions")]
#[repr(transparent)]
pub struct PyFkAssumptions {
    pub(crate) fk_assumptions: FkAssumptions,
}

#[pymethods]
impl PyFkAssumptions {
    /// Constructor.
    /// # Panics
    /// TODO
    #[new]
    #[must_use]
    pub fn new(assumption: &str) -> Self {
        Self {
            fk_assumptions: FkAssumptions::from_str(assumption).unwrap(),
        }
    }
}

#[pymethods]
impl PyFkTable {
    /// Constructor from an existing grid.
    /// # Panics
    /// TODO
    #[new]
    #[must_use]
    pub fn new(grid: PyGrid) -> Self {
        Self {
            fk_table: FkTable::try_from(grid.grid).unwrap(),
        }
    }

    /// Read from given path.
    /// # Panics
    /// TODO
    #[must_use]
    #[staticmethod]
    pub fn read(path: PathBuf) -> Self {
        Self {
            fk_table: FkTable::try_from(
                Grid::read(BufReader::new(File::open(path).unwrap())).unwrap(),
            )
            .unwrap(),
        }
    }

    /// Get cross section tensor.
    ///
    /// # Errors
    /// TODO
    ///
    /// Returns
    /// -------
    /// numpy.ndarray :
    ///     4-dimensional tensor with indixes: bin, channel, x1, x2
    pub fn table<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        Ok(self.fk_table.table().into_pyarray_bound(py))
    }

    /// Get number of bins.
    ///
    /// Returns
    /// -------
    /// int :
    ///     number of bins
    #[must_use]
    pub fn bins(&self) -> usize {
        self.fk_table.grid().bin_info().bins()
    }

    /// Extract the normalizations for each bin.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray
    ///     bin normalizations
    #[must_use]
    pub fn bin_normalizations<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.fk_table
            .grid()
            .bin_info()
            .normalizations()
            .into_pyarray_bound(py)
    }

    /// Extract the number of dimensions for bins.
    ///
    /// E.g.: two differential cross-sections will return 2.
    ///
    /// Returns
    /// -------
    /// int :
    ///     bin dimension
    #[must_use]
    pub fn bin_dimensions(&self) -> usize {
        self.fk_table.grid().bin_info().dimensions()
    }

    /// Extract the left edges of a specific bin dimension.
    ///
    /// Parameters
    /// ----------
    /// dimension : int
    ///     bin dimension
    ///
    /// Returns
    /// -------
    /// numpy.ndarray(float) :
    ///     left edges of bins
    #[must_use]
    pub fn bin_left<'py>(&self, dimension: usize, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.fk_table
            .grid()
            .bin_info()
            .left(dimension)
            .into_pyarray_bound(py)
    }

    /// Extract the right edges of a specific bin dimension.
    ///
    /// Parameters
    /// ----------
    /// dimension : int
    ///     bin dimension
    ///
    /// Returns
    /// -------
    /// numpy.ndarray(float) :
    ///     right edges of bins
    #[must_use]
    pub fn bin_right<'py>(&self, dimension: usize, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.fk_table
            .grid()
            .bin_info()
            .right(dimension)
            .into_pyarray_bound(py)
    }

    /// Get channels.
    ///
    /// Returns
    /// -------
    /// list(tuple(float,float)) :
    ///     channel functions as pid tuples
    #[must_use]
    pub fn channels(&self) -> Vec<Vec<i32>> {
        self.fk_table.channels()
    }

    /// Get reference (fitting) scale.
    ///
    /// Returns
    /// -------
    /// float :
    ///     reference scale
    #[must_use]
    pub fn muf2(&self) -> f64 {
        self.fk_table.muf2()
    }

    /// Get (unique) interpolation grid.
    ///
    /// Returns
    /// -------
    /// x_grid : numpy.ndarray(float)
    ///     interpolation grid
    #[must_use]
    pub fn x_grid<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.fk_table.x_grid().into_pyarray_bound(py)
    }

    /// Write to file.
    ///
    /// # Panics
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// path : str
    ///     file path
    pub fn write(&self, path: PathBuf) {
        self.fk_table
            .grid()
            .write(File::create(path).unwrap())
            .unwrap();
    }

    /// Write to file using lz4.
    ///
    /// # Panics
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// path : str
    ///     file path
    pub fn write_lz4(&self, path: PathBuf) {
        self.fk_table
            .grid()
            .write_lz4(File::create(path).unwrap())
            .unwrap();
    }

    /// Convolve with a single distribution.
    ///
    /// # Panics
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// pdg_id : integer
    ///     PDG Monte Carlo ID of the hadronic particle
    /// xfx : callable
    ///     lhapdf like callable with arguments `pid, x, Q2` returning x*pdf for :math:`x`-grid
    ///
    /// Returns
    /// -------
    /// numpy.ndarray(float) :
    ///     cross sections for all bins
    #[must_use]
    #[allow(clippy::needless_pass_by_value)]
    #[pyo3(signature = (pdg_conv, xfx, bin_indices = None, channel_mask= None))]
    pub fn convolve_with_one<'py>(
        &self,
        pdg_conv: PyRef<PyConv>,
        xfx: &Bound<'py, PyAny>,
        bin_indices: Option<PyReadonlyArray1<usize>>,
        channel_mask: Option<PyReadonlyArray1<bool>>,
        py: Python<'py>,
    ) -> Bound<'py, PyArray1<f64>> {
        let mut xfx = |id, x, q2| xfx.call1((id, x, q2)).unwrap().extract().unwrap();
        let mut alphas = |_| 1.0;
        let mut lumi_cache =
            ConvolutionCache::new(vec![pdg_conv.conv.clone()], vec![&mut xfx], &mut alphas);
        self.fk_table
            .convolve(
                &mut lumi_cache,
                &bin_indices.map_or(vec![], |b| b.to_vec().unwrap()),
                &channel_mask.map_or(vec![], |l| l.to_vec().unwrap()),
            )
            .into_pyarray_bound(py)
    }

    /// Convoluve grid with two different distribution.
    ///
    /// # Panics
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// pdg_id1 : integer
    ///     PDG Monte Carlo ID of the first hadronic particle
    /// xfx1 : callable
    ///     lhapdf like callable with arguments `pid, x, Q2` returning x*pdf for :math:`x`-grid
    /// pdg_id2 : integer
    ///     PDG Monte Carlo ID of the second hadronic particle
    /// xfx2 : callable
    ///     lhapdf like callable with arguments `pid, x, Q2` returning x*pdf for :math:`x`-grid
    ///
    /// Returns
    /// -------
    /// numpy.ndarray(float) :
    ///     cross sections for all bins
    #[must_use]
    #[allow(clippy::needless_pass_by_value)]
    #[pyo3(signature = (pdg_conv1, xfx1, pdg_conv2, xfx2, bin_indices = None, channel_mask= None))]
    pub fn convolve_with_two<'py>(
        &self,
        pdg_conv1: PyRef<PyConv>,
        xfx1: &Bound<'py, PyAny>,
        pdg_conv2: PyRef<PyConv>,
        xfx2: &Bound<'py, PyAny>,
        bin_indices: Option<PyReadonlyArray1<usize>>,
        channel_mask: Option<PyReadonlyArray1<bool>>,
        py: Python<'py>,
    ) -> Bound<'py, PyArray1<f64>> {
        let mut xfx1 = |id, x, q2| xfx1.call1((id, x, q2)).unwrap().extract().unwrap();
        let mut xfx2 = |id, x, q2| xfx2.call1((id, x, q2)).unwrap().extract().unwrap();
        let mut alphas = |_| 1.0;
        let mut lumi_cache = ConvolutionCache::new(
            vec![pdg_conv1.conv.clone(), pdg_conv2.conv.clone()],
            vec![&mut xfx1, &mut xfx2],
            &mut alphas,
        );
        self.fk_table
            .convolve(
                &mut lumi_cache,
                &bin_indices.map_or(vec![], |b| b.to_vec().unwrap()),
                &channel_mask.map_or(vec![], |l| l.to_vec().unwrap()),
            )
            .into_pyarray_bound(py)
    }

    /// Optimize storage.
    ///
    /// In order to perform any relevant optimization, assumptions are needed, and they are passed
    /// as parameters to the function.
    ///
    /// Parameters
    /// ----------
    /// assumptions : PyFkAssumptions
    ///     assumptions about the FkTable properties, declared by the user, deciding which
    ///     optimizations are possible
    #[allow(clippy::needless_pass_by_value)]
    pub fn optimize(&mut self, assumptions: PyRef<PyFkAssumptions>) {
        self.fk_table.optimize(assumptions.fk_assumptions);
    }
}

/// Register submodule in parent.
/// # Errors
///
/// Raises an error if (sub)module is not found.
pub fn register(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new_bound(parent_module.py(), "fk_table")?;
    m.setattr(pyo3::intern!(m.py(), "__doc__"), "FK table interface.")?;
    pyo3::py_run!(
        parent_module.py(),
        m,
        "import sys; sys.modules['pineappl.fk_table'] = m"
    );
    m.add_class::<PyFkTable>()?;
    m.add_class::<PyFkAssumptions>()?;
    parent_module.add_submodule(&m)
}
