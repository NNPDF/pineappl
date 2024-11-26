//! FK table interface.

use super::convolutions::PyConv;
use super::grid::PyGrid;
use numpy::{IntoPyArray, PyArray1, PyArrayDyn};
use pineappl::convolutions::ConvolutionCache;
use pineappl::fk_table::{FkAssumptions, FkTable};
use pineappl::grid::Grid;
use pyo3::prelude::*;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::str::FromStr;

/// PyO3 wrapper to :rustdoc:`pineappl::fk_table::FkAssumptions <fk_table/enum.FkAssumptions.html>`.
#[pyclass(name = "FkAssumptions")]
#[repr(transparent)]
pub struct PyFkAssumptions {
    pub(crate) fk_assumptions: FkAssumptions,
}

#[pymethods]
impl PyFkAssumptions {
    /// Constructor.
    ///
    /// # Panics
    ///
    /// Panics if the `assumption` is not one of the possibilities.
    #[new]
    #[must_use]
    pub fn new(assumption: &str) -> Self {
        Self {
            fk_assumptions: FkAssumptions::from_str(assumption).unwrap(),
        }
    }
}

/// PyO3 wrapper to :rustdoc:`pineappl::fk_table::FkTable <fk_table/struct.FkTable.html>`.
#[pyclass(name = "FkTable")]
#[repr(transparent)]
pub struct PyFkTable {
    pub(crate) fk_table: FkTable,
}

#[pymethods]
impl PyFkTable {
    /// Constructor from an existing grid.
    ///
    /// # Panics
    /// TODO
    #[new]
    #[must_use]
    pub fn new(grid: PyGrid) -> Self {
        Self {
            fk_table: FkTable::try_from(grid.grid).unwrap(),
        }
    }

    /// Read an FK Table from given path.
    ///
    /// # Panics
    /// TODO
    ///
    /// Parameteters
    /// ------------
    /// path : str
    ///     path to the FK table
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

    /// Get the type(s) of convolution(s) for the current FK table.
    ///
    /// Returns
    /// list(PyConv):
    ///     list of convolution type with the corresponding PIDs
    #[getter]
    #[must_use]
    pub fn convolutions(&self) -> Vec<PyConv> {
        self.fk_table
            .grid()
            .convolutions()
            .iter()
            .map(|conv| PyConv { conv: conv.clone() })
            .collect()
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
    ///
    /// Panics if the specified path is non-writeable (non-existent or missing permissions).
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
    ///
    /// Panics if the specified path is non-writeable (non-existent or missing permissions).
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

    /// Set a metadata key-value pair in the FK Table.
    ///
    /// Parameters
    /// ----------
    /// key : str
    ///     key
    /// value : str
    ///     value
    pub fn set_key_value(&mut self, key: &str, value: &str) {
        self.fk_table.set_key_value(key, value);
    }

    /// Get metadata values stored in the grid.
    ///
    ///
    /// Returns
    /// -------
    /// dict :
    ///     key, value map
    #[getter]
    #[must_use]
    pub fn key_values(&self) -> BTreeMap<String, String> {
        self.fk_table.grid().metadata().clone()
    }

    /// Convolve the FK table with as many distributions.
    ///
    /// # Panics
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// pdg_convs : list(PyConv)
    ///     list containing the types of convolutions and PID
    /// xfxs : list(callable)
    ///     list of lhapdf-like callable with arguments `pid, x, Q2` returning x*pdf
    /// bin_indices : numpy.ndarray(int)
    ///     A list with the indices of the corresponding bins that should be calculated. An
    ///     empty list means that all bins should be calculated.
    /// channel_mask : numpy.ndarray(bool)
    ///     Mask for selecting specific channels. The value `True` means the
    ///     corresponding channel is included. An empty list corresponds to all channels being
    ///     enabled.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray(float) :
    ///     cross sections for all bins
    #[must_use]
    #[pyo3(signature = (pdg_convs, xfxs, bin_indices = None, channel_mask= None))]
    pub fn convolve<'py>(
        &self,
        pdg_convs: Vec<PyRef<PyConv>>,
        xfxs: Vec<PyObject>,
        bin_indices: Option<Vec<usize>>,
        channel_mask: Option<Vec<bool>>,
        py: Python<'py>,
    ) -> Bound<'py, PyArray1<f64>> {
        let mut xfx_funcs: Vec<_> = xfxs
            .iter()
            .map(|xfx| {
                move |id: i32, x: f64, q2: f64| {
                    xfx.call1(py, (id, x, q2)).unwrap().extract(py).unwrap()
                }
            })
            .collect();

        let mut alphas = |_| 1.0;
        let mut lumi_cache = ConvolutionCache::new(
            pdg_convs.into_iter().map(|pdg| pdg.conv.clone()).collect(),
            xfx_funcs
                .iter_mut()
                .map(|fx| fx as &mut dyn FnMut(i32, f64, f64) -> f64)
                .collect(),
            &mut alphas,
        );

        self.fk_table
            .convolve(
                &mut lumi_cache,
                &bin_indices.unwrap_or_default(),
                &channel_mask.unwrap_or_default(),
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
