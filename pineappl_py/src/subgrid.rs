//! Subgrid interface.

use ndarray::{ArrayD, Dimension};
use numpy::{IntoPyArray, PyArrayDyn, PyReadonlyArrayDyn};
use pineappl::packed_array::PackedArray;
use pineappl::subgrid::{ImportSubgridV1, Subgrid, SubgridEnum};
use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl::subgrid::ImportSubgridV1 <subgrid/struct.ImportSubgridV1.html>`.
#[pyclass(name = "ImportSubgridV1")]
#[derive(Clone)]
#[repr(transparent)]
pub struct PyImportSubgridV1 {
    pub(crate) import_subgrid: ImportSubgridV1,
}

#[pymethods]
impl PyImportSubgridV1 {
    /// Constructor.
    ///
    /// # Panics
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// array : numpy.ndarray(float)
    ///     `N`-dimensional array with all weights
    /// node_values: list(list(float))
    ///     list containing the arrays of energy scales {q1, ..., qn} and momentum fractions
    ///     {x1, ..., xn}.
    #[new]
    #[must_use]
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

/// PyO3 wrapper to :rustdoc:`pineappl::subgrid::SubgridEnum <subgrid/struct.SubgridEnum.html>`
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
    m.setattr(pyo3::intern!(m.py(), "__doc__"), "Subgrid interface.")?;
    pyo3::py_run!(
        parent_module.py(),
        m,
        "import sys; sys.modules['pineappl.subgrid'] = m"
    );
    m.add_class::<PyImportSubgridV1>()?;
    m.add_class::<PySubgridEnum>()?;
    parent_module.add_submodule(&m)
}
