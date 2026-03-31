//! Subgrid interface.
//!
//! [`ImportSubgridV1`] is mainly for **injecting coefficient tables from other programs** (dumped
//! grids, external MC output, etc.) that are **meant to be convolved with parton `f(x)`**, not for the
//! usual event-by-event filling of a grid. Those coefficients use the **`f(x)`** convention (and the
//! analogous fragmentation variable), **not** **`x * f(x)`**. When you convolve with LHAPDF-style
//! callables that return **`x * f`**, the Rust library divides by `x` so the result matches. Do not
//! store **`x * f`** in an imported subgrid or you double-count `x`. See the Rust `pineappl::convolutions`
//! and `pineappl::subgrid::ImportSubgridV1` docs and <https://github.com/NNPDF/pineappl/issues/388>.

use ndarray::{ArrayD, Dimension};
use numpy::{IntoPyArray, PyArrayDyn, PyReadonlyArrayDyn};
use pineappl::packed_array::PackedArray;
use pineappl::subgrid::{ImportSubgridV1, Subgrid, SubgridEnum};
use pyo3::prelude::*;

/// Sparse subgrid for **coefficient functions dumped from outside PineAPPL**, to be convolved with **f**.
///
/// Typical use: load tabulated coefficients from another code; they must be defined to fold with
/// parton **f(x)**, not **x * f(x)**. Non-zero entries in ``array`` follow that **f** convention.
/// :class:`~pineappl.grid.Grid` ``.convolve`` still expects LHAPDF-style callables that return
/// ``x * f``; the core library divides by ``x`` when building the luminosity. See the submodule
/// docstring and `issue 388 <https://github.com/NNPDF/pineappl/issues/388>`__.
#[pyclass(name = "ImportSubgridV1")]
#[derive(Clone)]
#[repr(transparent)]
pub struct PyImportSubgridV1 {
    pub(crate) import_subgrid: ImportSubgridV1,
}

#[pymethods]
impl PyImportSubgridV1 {
    /// Constructor for **imported coefficient tables** (external dumps), convolved with ``f`` later.
    ///
    /// Stored weights must be ``f(x)``, not ``x * f(x)``, so they match :meth:`pineappl.grid.Grid.convolve`
    /// with standard LHAPDF ``xfx`` callbacks (returning ``x * f``).
    ///
    /// # Panics
    ///
    /// Panics if the array shape does not match `node_values` or indexing the packed array fails.
    ///
    /// Parameters
    /// ----------
    /// array : numpy.ndarray(float)
    ///     N-dimensional array of non-zero coefficients (``f`` convention: fold with ``f``, not ``x * f``).
    /// node_values : list(list(float))
    ///     Per-dimension node coordinates (scales, x values, etc.).
    #[new]
    #[must_use]
    #[allow(clippy::needless_pass_by_value)]
    pub fn new(array: PyReadonlyArrayDyn<f64>, node_values: Vec<Vec<f64>>) -> Self {
        let mut sparse_array: PackedArray<f64> =
            PackedArray::new(node_values.iter().map(Vec::len).collect());

        for (index, value) in array
            .as_array()
            .indexed_iter()
            .filter(|(_, value)| **value != 0.0)
        {
            sparse_array[index.as_array_view().to_slice().unwrap()] = *value;
        }

        Self {
            import_subgrid: ImportSubgridV1::new(sparse_array, node_values),
        }
    }

    /// Ensures that the subgrid has type `PySubgridEnum`.
    #[must_use]
    pub fn into(&self) -> PySubgridEnum {
        PySubgridEnum {
            subgrid_enum: self.import_subgrid.clone().into(),
        }
    }
}

/// `PyO3` wrapper to :rustdoc:`pineappl::subgrid::SubgridEnum <subgrid/struct.SubgridEnum.html>`
#[pyclass(name = "SubgridEnum")]
#[derive(Clone)]
#[repr(transparent)]
pub struct PySubgridEnum {
    pub(crate) subgrid_enum: SubgridEnum,
}

#[pymethods]
impl PySubgridEnum {
    /// Scale the subgrid by `factor`.
    ///
    /// Parameters
    /// ----------
    /// factor : float
    ///     scaling factor
    pub fn scale(&mut self, factor: f64) {
        self.subgrid_enum.scale(factor);
    }

    /// Get the values of nodes used for the subgrids
    #[getter]
    pub fn node_values(&mut self) -> Vec<Vec<f64>> {
        self.subgrid_enum.node_values()
    }

    /// Get the shape of the subgrids
    #[getter]
    pub fn shape(&mut self) -> Vec<usize> {
        self.subgrid_enum.shape().to_vec()
    }

    /// Return the dense array of the subgrid.
    #[must_use]
    pub fn to_array<'py>(
        &mut self,
        py: Python<'py>,
        shape: Vec<usize>,
    ) -> Bound<'py, PyArrayDyn<f64>> {
        let mut array_subgrid = ArrayD::<f64>::zeros(shape);

        for (index, value) in self.subgrid_enum.indexed_iter() {
            array_subgrid[index.as_slice()] = value;
        }
        array_subgrid.into_pyarray(py)
    }

    /// Clone.
    #[must_use]
    pub fn into(&self) -> Self {
        self.clone()
    }
}

/// Register submodule in parent.
///
/// # Errors
///
/// Raises Errors if (sub-)module is not found.
pub fn register(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "subgrid")?;
    m.setattr(
        pyo3::intern!(m.py(), "__doc__"),
        "Subgrid interface.\n\n\
         ImportSubgridV1 is for coefficient tables dumped from outside PineAPPL, to be convolved with f; \
         array values use the f(x) convention, not x*f(x). Grid.convolve uses LHAPDF x*f callbacks \
         (see pineappl issue 388).",
    )?;
    pyo3::py_run!(
        parent_module.py(),
        m,
        "import sys; sys.modules['pineappl.subgrid'] = m"
    );
    m.add_class::<PyImportSubgridV1>()?;
    m.add_class::<PySubgridEnum>()?;
    parent_module.add_submodule(&m)
}
